#! python

"""
Samuel Joseph Bryson
Copyright 2026

- 1 - Run fastp QC on a set of samples (r1.fq.gz & r2.fq.gz).
- 2 - Parse json output to get filtering stats.
"""
import argparse
import multiprocessing as mp
from pathlib import Path
import json
from sb_projects.wrappers import FastpQC
from sb_projects.config_manager import ConfigManager
from sb_projects.subprocess_utilities import run_check_call


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to run fastq QC on a set of paired reads.")

    parser.add_argument("--config_file", type=str, required=True, 
                        help="Path to the configuration file (str).")
    
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Directory containing input files (str).")
    
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Directory where output files will be saved (str).")
    
    parser.add_argument("--threads", type=int, default=4, 
                        help="Number of threads to use (int, default: 4).")
    
    parser.add_argument("--processes", type=int, default=2, 
                        help="Number of parallel samples to process at once (int, default: 2).")
    
    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")
    
    ## --- ToDo: --- ##
    # add sample column name
    # add r1 column name
    # add r2 column name

    return parser.parse_args()

# Run fastp on a set of r1 and r2 <sample>_<r#>.fq.gz
def run_fastp_qc(
        threads:     int,
        output_json: Path,
        output_html: Path,
        in_r1:       Path,
        in_r2:       Path,
        out_r1:      Path,
        out_r2:      Path,
        dry_run:     bool,
    ) -> None:
    # Init a FastpQC Wrapper obj.
    process = FastpQC(
        deduplicate   = True,
        gz            = 6, # default
        cut_win_size  = 5, # default
        cut_mean_qual = 25, # default
        min_length    = 100, # default
        threads       = threads, 
        output_json   = output_json,
        output_html   = output_html,
        in_r1         = in_r1,
        in_r2         = in_r2,
        out_r1        = out_r1,
        out_r2        = out_r2,
        dry_run       = dry_run,
    )
    cmd = process.build()
    run_check_call(
        formatted_command = cmd, 
        devnull = True,
        dry_run = dry_run,
    )

def parse_fastp_json(json_file: Path, dry_run: bool) -> dict:
    """Parse the fastp JSON output file and return a dict with specific QC metrics."""
    if dry_run == True:
        # fake data for dry_run tests
        metrics = {
            "duplication_rate"     : 0.1,
            "adapter_trimmed_reads": 1,
            "adapter_trimmed_bases": 1,
            "r1_preQC_total_reads" : 1,
            "r2_preQC_total_reads" : 1,
            "r1_preQC_total_bases" : 1,
            "r2_preQC_total_bases" : 1,
            "r1_postQC_total_reads": 1,
            "r2_postQC_total_reads": 1,
            "r1_postQC_total_bases": 1,
            "r2_postQC_total_bases": 1,
        }
        return metrics
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Specific json keys to extract
    metrics = {
        "duplication_rate"     : data.get("duplication", {}).get("rate"),
        "adapter_trimmed_reads": data.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
        "adapter_trimmed_bases": data.get("adapter_cutting", {}).get("adapter_trimmed_bases"),
    }

    # Dict of json keys with same internal keys and string to use in metrics keys
    groups = {
        "read1_before_filtering": "r1_preQC", 
        "read2_before_filtering": "r2_preQC",
        "read1_after_filtering" : "r1_postQC", 
        "read2_after_filtering" : "r2_postQC",
    }
    
    # Could add: "q20_bases", "q30_bases", "q40_bases", "total_cycles"
    group_keys = ["total_reads", "total_bases"] 
    
    # Add read_groups data to metrics dict
    for k,v in groups.items():
        group_data = data.get(k, {})
        for key in group_keys:
            # Creates keys like 'r1_preQC_total_reads'
            metrics[f"{v}_{key}"] = group_data.get(key)

    return metrics

# Multiprocessing worker function
def sample_worker(task_args: dict) -> dict:
    """Worker function to run Fastp and parse the resulting JSON."""
    index      = task_args["index"]
    row        = task_args["row"]
    input_dir  = task_args["input_dir"]
    output_dir = task_args["output_dir"]
    threads    = task_args["threads"]
    dry_run    = task_args["dry_run"]

    if dry_run: 
        print(f"\nWorker processing: {index}\t{row}") 

    # Setup inputs and outputs
    SAMPLE      = row["sample"]
    fastp_json  = f"{SAMPLE}.json"
    output_json = output_dir / fastp_json
    fastp_html  = f"{SAMPLE}.html"
    output_html = output_dir / fastp_html
    raw_r1      = row["r1"]
    raw_r2      = row["r2"]
    in_r1       = input_dir / raw_r1
    in_r2       = input_dir / raw_r2
    qc_r1       = f"{SAMPLE}.qc.r1.fq.gz"
    qc_r2       = f"{SAMPLE}.qc.r2.fq.gz"
    out_r1      = output_dir / qc_r1
    out_r2      = output_dir / qc_r2

    try:
        # Step 1: Run Fastp
        run_fastp_qc(
            threads     = threads,
            output_json = output_json,
            output_html = output_html,
            in_r1       = in_r1,
            in_r2       = in_r2,
            out_r1      = out_r1,
            out_r2      = out_r2,
            dry_run     = dry_run,
        )

        # Step 2: Parse JSON
        qc_stats = parse_fastp_json(json_file=output_json, dry_run=dry_run)

        # Return all results back to main process
        return {
            "index": index,
            "success": True,
            "qc_r1": qc_r1,
            "qc_r2": qc_r2,
            "fastp_json": fastp_json,
            "stats": qc_stats,
            "error": None
        }

    except Exception as e:
        return {"index": index, "success": False, "error": str(e)}

def main():
    
    # Parse args
    args       = parse_args()
    config     = Path(args.config_file)
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    threads    = args.threads
    processes  = args.processes
    dry_run    = True if args.dry_run else False

    # Load config
    proj_config = ConfigManager(config_file = config)
    
    ## FastpQC ##

    # Define output columns and fill with default value: "incomplete".
    output_columns = [
        "qc_r1", "qc_r2", "fastp_json", "duplication_rate", 
        "adapter_trimmed_reads", "adapter_trimmed_bases", 
        "r1_preQC_total_reads", "r2_preQC_total_reads", 
        "r1_preQC_total_bases", "r2_preQC_total_bases", 
        "r1_postQC_total_reads", "r2_postQC_total_reads", 
        "r1_postQC_total_bases", "r2_postQC_total_bases"
    ]
    
    for col in output_columns:
        proj_config.add_column(name = col, default_value = "incomplete")
    
    # Build worker pool task list 
    tasks = []
    for index, row in proj_config:
        tasks.append({
            "index": index,
            "row": row,
            "input_dir": input_dir,
            "output_dir": output_dir,
            "threads": threads,
            "dry_run": dry_run
        })

    # Start the multiprocessing pool
    print(f"Starting multiprocessing pool with {processes} workers...")
    with mp.Pool(processes=processes) as pool:
        for result in pool.imap_unordered(sample_worker, tasks):
            
            index = result["index"]

            if result.get("success"):
                # Update file paths
                proj_config.update_row(index, "qc_r1", result["qc_r1"])
                proj_config.update_row(index, "qc_r2", result["qc_r2"])
                proj_config.update_row(index, "fastp_json", result["fastp_json"])
                
                # Update all the stats
                for k, v in result["stats"].items():
                    proj_config.update_row(index, k, v)
            else:
                print(f"Error processing row {index}: {result['error']}")

            # Save incrementally
            proj_config.save()

if __name__ == "__main__":
    main()