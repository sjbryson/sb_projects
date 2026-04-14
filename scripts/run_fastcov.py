#! python

"""
- 1 - Run minimap2 | fastcov pipeline on a set of samples (r1.fq.gz & r2.fq.gz).
- 2 - Output coverage and alignment stats to {sample}.json.
- 3 - Don't save the SAM output from minimap2.
"""

import argparse
import multiprocessing as mp
from pathlib import Path
import json
from typing import Optional
from dataclasses import dataclass, field
from sb_projects.wrapper import Wrapper
from sb_projects.config import ConfigDf
from sb_projects.subprocesses import run_check_call

# Filtering Defaults >>> These could be converted to args.
#MIN_AP = 0.5   # float:   Alignment_Percentage = Alignment_Length / Read_Length
#MIN_PI = 0.5   # float:   percent identity = Num_Identities / Alignment_Length
MIN_AS = 150   # integer: SAM tag -> Alignment_Score
#MIN_AL = 75    # integer: CIGAR -> Alignment_Length
#MIN_SL = 1.0   # float:   SL = Alignment_Score / Alignment_Length
#MIN_MQ = 10  # integer: MAPQ score

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to run alignment based detection on a set of paired reads.")

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
    
    parser.add_argument("--db", type=str, required=True,
                        help="Pre-computed Minimap database path (str).")
    
    parser.add_argument("--sample_col", type=str, default="sample", 
                        help="Column name for samples in config file (str).")
    
    parser.add_argument("--r1_col", type=str, default="r1", 
                        help="Column name for input R1 in config file (str).")
    
    parser.add_argument("--r2_col", type=str, default="r2", 
                        help="Column name for input R2 in config file (str).")
    
    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")

    return parser.parse_args()

# Fastp Wrapper subclass
@dataclass(kw_only=True)
class FastCov(Wrapper):
    input_mmi:    Path | str      = field(metadata={'type': 'input_file'})
    r1:           Path | str      = field(metadata={'type': 'input_file'})
    r2:           Path | str      = field(metadata={'type': 'input_file'})
    sample:       str             
    map_threads:  Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cov_threads:  Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    sort_threads: Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-@ {value}'})
    min_ap:       Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--min-ap {value}'})
    min_pi:       Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--min-pi {value}'})
    min_as:       Optional[int]   = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--min-as {value}'})
    min_al:       Optional[int]   = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--min-al {value}'})
    min_sl:       Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--min-sl {value}'})
    cmd:          str             = "minimap2 -ax sr --eqx {map_threads} {input_mmi} {r1} {r2} | \
                                     fastcov {cov_threads} -r {sample} {min_as} | \
                                     samtools sort {sort_threads} - -o {sample}.sorted.bam"

def _parse_fastp_json(json_file: Path, db_name: str, dry_run: bool) -> dict:
    """Parse the fastcov JSON output file and return a dict with specific stats."""

    if dry_run:
        # fake data for dry_run tests
        metrics = {
            f"{db_name}_run_time":           100,
            f"{db_name}_total_alignments":   100,
            f"{db_name}_primary_alignments": 100,
            f"{db_name}_reference_hits":     100,
        }
        return metrics
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Specific json keys to extract
    metrics = {
        f"{db_name}_run_time":           data.get("1-run_stats", {}).get("total_run_time_seconds"),
        f"{db_name}_total_alignments":   data.get("1-run_stats", {}).get("total_alignments"),
        f"{db_name}_primary_alignments": data.get("1-run_stats", {}).get("passed_primary_alignments"),
        f"{db_name}_reference_hits":     data.get("1-run_stats", {}).get("num_refs_primary"),
    }

    return metrics

def _sample_worker(task_args: dict) -> dict:
    """Worker function to process a single sample. Returns a dictionary with the results."""
    
    index       = task_args["index"]
    row         = task_args["row"]
    sample_col  = task_args["sample_col"]
    db_name     = task_args["db_name"]
    threads     = task_args["threads"]
    dry_run     = task_args["dry_run"]

    if dry_run: 
        print(f"\nWorker processing: {index}\t{row}") 
    
    SAMPLE    = row[sample_col]
    JSON_OUT  = f"{SAMPLE}.{db_name}.json"
    BAM_OUT   = f"{SAMPLE}.{db_name}.sorted.bam"

    try:
       # 1 - Run alignment and get alignment stats.
        process = FastCov(
            input_mmi    = task_args["db"],
            r1           = task_args["input_dir"] / row[task_args["r1_col"]],
            r2           = task_args["input_dir"] / row[task_args["r2_col"]],
            sample       = task_args["output_dir"] / f"{SAMPLE}.{db_name}",
            map_threads  = threads - 4 if threads >=12 else threads - 2,
            cov_threads  = 4 if threads >=12 else 2,
            sort_threads = 4,
            min_as       = MIN_AS,
            dry_run      = dry_run,
        )

        cmd = process.build()
        run_check_call(
            formatted_command = cmd, 
            devnull = True,
            dry_run = dry_run,
        )
        
        # 2 - Parse the json output
        metrics = _parse_fastp_json(
            task_args["output_dir"] / JSON_OUT,
            db_name,
            dry_run
        )

        # Return all results back to main process
        return {
            "index":                         index,
            "success":                       True,
            f"{db_name}_stats_json":         JSON_OUT,
            f"{db_name}_bam":                BAM_OUT,
            f"{db_name}_run_time":           metrics[f"{db_name}_run_time"],
            f"{db_name}_total_alignments":   metrics[f"{db_name}_total_alignments"],
            f"{db_name}_primary_alignments": metrics[f"{db_name}_primary_alignments"],
            f"{db_name}_reference_hits":     metrics[f"{db_name}_reference_hits"],
            "error":                         None,
        }

    except Exception as e:
        return { "index": index, "success": False, "error": str(e) }

def main():
    
    # Parse args
    args       = parse_args()
    dry_run    = True if args.dry_run else False
    config     = Path(args.config_file)
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    threads    = args.threads
    processes  = args.processes
    sample_col = args.sample_col
    r1_col     = args.r1_col
    r2_col     = args.r2_col
    db         = Path(args.db)
    db_name    = db.stem
    #min_ap     = MIN_AP
    #min_pi     = MIN_PI
    #min_as     = MIN_AS
    #min_al     = MIN_AL
    #min_sl     = MIN_SL
    #min_mq     = MIN_MQ

    # Load config
    proj_config = ConfigDf(config_file = config)

    # Define output columns & add to config
    output_columns = [
        f"{db_name}_stats_json", f"{db_name}_bam",
        f"{db_name}_run_time", f"{db_name}_total_alignments",
        f"{db_name}_primary_alignments", f"{db_name}_reference_hits"
    ]

    for col in output_columns:
        proj_config.add_column(name = col, default_value = None)
    
    # Build worker pool task list 
    tasks = []
    for index, row in proj_config:
        tasks.append({
            "index":      index,
            "row":        row,
            "input_dir":  input_dir,
            "output_dir": output_dir,
            "threads":    threads,
            "dry_run":    dry_run,
            "sample_col": sample_col,
            "r1_col":     r1_col,
            "r2_col":     r2_col,
            "db":         db,
            "db_name":    db_name
        })

    # Start the multiprocessing pool
    print(f"Starting multiprocessing pool with {processes} workers...")
    with mp.Pool(processes=processes) as pool:
        for result in pool.imap_unordered(_sample_worker, tasks):
            
            index = result["index"]

            if result.get("success"):
                # Update all the columns in config_df
                for c in output_columns:
                    proj_config.update_row(index, c, result.get(c))
            else:
                print(f"Error processing row {index}: {result['error']}")

            # Save incrementally
            proj_config.save()

if __name__ == "__main__":
    main()