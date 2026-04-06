#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026

- 1 - Run fastp QC on a set of samples (r1.fq.gz & r2.fq.gz).
- 2 - Parse json output to get filtering stats.
"""
import argparse
from pathlib import Path
from typing import Tuple
import json
from src.config_manager import ConfigManager
from src.subprocess_utilities import run_check_call
from src.file_utilities import delete_file
from wrappers import FastpQC

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
    
    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")
    
    ## add sample column name
    ## add r1 column name
    ## add r2 column name

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

def parse_fastp_json(json_file: Path) -> dict:
    """Parse the fastp JSON output file and return a dict with specific QC metrics."""
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Specific json keys to extract
    metrics = {
        "duplication_rate": data.get("duplication", {}).get("rate"),
        "adapter_trimmed_reads": data.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
        "adapter_trimmed_bases": data.get("adapter_cutting", {}).get("adapter_trimmed_bases"),
    }

    # Dict of json keys with same internal keys and string to use in metrics keys
    groups = {
        "read1_before_filtering": "r1_preQC", 
        "read2_before_filtering": "r2_preQC",
        "read1_after_filtering": "r1_postQC", 
        "read2_after_filtering": "r2_postQC"
    }
    
    group_keys = [  #"q20_bases", "q30_bases", "q40_bases", "total_cycles"
        "total_reads", 
        "total_bases", 
    ]
    # Add read_groups data to metrics dict
    for k,v in groups.items():
        group_data = data.get(k, {})
        for key in group_keys:
            # Creates keys like 'read1_before_filtering_total_reads'
            metrics[f"{v}_{key}"] = group_data.get(key)

    return metrics

def main():
    
    # Parse args
    args       = parse_args()
    config     = Path(args.config_file)
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    threads    = args.threads
    dry_run    = True if args.dry_run else False

    # Load config
    proj_config = ConfigManager(config_file = config)
    
    # Define output columns
    workflow_outputs = [
        "qc_r1",
        "qc_r2",
        "fastp_json",
        #"",
        #"",
    ]
    
    # Update config with output fields
    for output in workflow_outputs:
        proj_config.add_column(name = output, default_value = "incomplete")
    
    # Iterate over rows in the config and run fastp
    for index, row in proj_config:

        try:
            if dry_run == True: 
                print(f"\n{index}\t{row}") 
            
            # Setup inputs & outputs
            SAMPLE      = row["sample"]
            fastp_json  = f"{SAMPLE}.json",
            output_json = output_dir / fastp_json,
            fastp_html  = f"{SAMPLE}.html",
            output_html = output_dir / fastp_html,
            raw_r1      = row["r1"],
            raw_r2      = row["r2"],
            in_r1       = input_dir / raw_r1,
            in_r2       = input_dir / raw_r2,
            qc_r1       = f"{SAMPLE}.qc.r1.fq.gz",
            qc_r2       = f"{SAMPLE}.qc.r2.fq.gz",
            out_r1      = output_dir / qc_r1,
            out_r2      = output_dir / qc_r2,
            
            # Run fastp for sample
            run_fastp_qc(
                threads     =  threads,
                output_json = output_json,
                output_html = output_html,
                in_r1       = in_r1,
                in_r2       = in_r2,
                out_r1      = out_r1,
                out_r2      = out_r2,
                dry_run     = dry_run,
            )
        except Exception as e:
            print(f"Fastp Error: {e}")
            continue

        else:
            # Update config with output files
            proj_config.update_row(index, "qc_r1", qc_r1)
            proj_config.update_row(index, "qc_r2", qc_r2)
            proj_config.update_row(index, "fastp_json", fastp_json)

            # delete html file?
            #delete_file(output_html, dry_run = dry_run)

            # save updated config after processing each row
            proj_config.save()
        

if __name__ == "__main__":
    main()