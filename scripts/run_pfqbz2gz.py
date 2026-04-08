#! python

"""
Samuel Joseph Bryson
Copyright 2026

- Run pfqbz2gz on a set of samples; convert bz2 compressed paired fastq to fq.gz.
"""

import argparse
from pathlib import Path
import multiprocessing as mp
from typing import Optional
from dataclasses import dataclass, field
from sb_projects.wrapper import Wrapper
from sb_projects.config import ConfigDf
from sb_projects.subprocess_utilities import run_check_call

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to convert paired fastq.bz2 to fq.gz.")

    parser.add_argument("--config_file", type=str, required=True, 
                        help="Path to the configuration file (str).")
    
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Directory containing input files (str).")
    
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Directory where output files will be saved (str).")
    
    parser.add_argument("--threads", type=int, default=4, 
                        help="Number of threads to use (int, default: 4).")
    
    parser.add_argument("--processes", type=int, default=2, 
                        help="Number of parallel subprocesses (int, default: 2).")
    
    parser.add_argument("--sample_col", type=str, default="sample", 
                        help="Column name for samples in config file (str).")
    
    parser.add_argument("--r1_col", type=str, default="r1_bz2", 
                        help="Column name for input R1 in config file (str).")
    
    parser.add_argument("--r2_col", type=str, default="r2_bz2", 
                        help="Column name for input R2 in config file (str).")
    
    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")

    ## --- ToDo: --- ##
    # add output r1 column name --> currently "r1_gz"
    # add output r2 column name --> currently "r2_gz"

    return parser.parse_args()

# Wrapper subclass for bz2 to gz conversion tool
@dataclass(kw_only=True)
class ConvertPairedFastqBZ2toGz(Wrapper):
    """Convert paired fastq with Bz2 compression to Gz compressed files."""
    r1_fastq:      Path | str    = field(metadata={'type': 'input_file'})
    r2_fastq:      Path | str    = field(metadata={'type': 'input_file'})
    output_prefix: Path | str
    threads:       Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:           str           = "pfqbz2gz {threads} --r1 {r1_fastq} --r2 {r2_fastq} -o {output_prefix}"

# Run fastp on a set of r1 and r2 <sample>_<r#>.fq.gz
def run_pfqbz2gz(
        threads:     int,
        in_r1:       Path,
        in_r2:       Path,
        out_prefix:  Path,
        dry_run:     bool,
    ) -> None:
    # Init a ConvertPairedFastqBZ2toGz Wrapper obj.
    process = ConvertPairedFastqBZ2toGz(
        r1_fastq      = in_r1,
        r2_fastq      = in_r2,
        output_prefix = out_prefix,
        threads       = threads,
        dry_run       = dry_run,
    )
    cmd = process.build()
    run_check_call(
        formatted_command = cmd,
        devnull           = True,
        dry_run           = dry_run,
    )

def sample_worker(task_args: dict) -> dict:
    """Worker function to process a single sample. Returns a dictionary with the results."""
    index      = task_args["index"]
    row        = task_args["row"]
    input_dir  = task_args["input_dir"]
    output_dir = task_args["output_dir"]
    threads    = task_args["threads"]
    dry_run    = task_args["dry_run"]
    sample_col = task_args["sample_col"]
    r1_col     = task_args["r1_col"]
    r2_col     = task_args["r2_col"]

    if dry_run: 
        print(f"\nWorker processing: {index}\t{row}") 

    SAMPLE     = row[sample_col]
    r1_bz2     = row[r1_col]
    r2_bz2     = row[r2_col]
    in_r1      = input_dir / r1_bz2
    in_r2      = input_dir / r2_bz2
    out_prefix = output_dir / SAMPLE
    r1_gz      = f"{SAMPLE}_R1.fq.gz"
    r2_gz      = f"{SAMPLE}_R2.fq.gz"
    
    try:
        run_pfqbz2gz(
            threads    = threads,
            in_r1      = in_r1,
            in_r2      = in_r2,
            out_prefix = out_prefix,
            dry_run    = dry_run,
        )
        # Return success state
        return {"index": index, "success": True, "r1_gz": r1_gz, "r2_gz": r2_gz, "error": None}
    
    except Exception as e:
        # Return failure state
        return {"index": index, "success": False, "r1_gz": None, "r2_gz": None, "error": str(e)}



def main():
    
    # Parse args
    args       = parse_args()
    config     = Path(args.config_file)
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    threads    = args.threads
    processes  = args.processes
    dry_run    = True if args.dry_run else False
    sample_col = args.sample_col
    r1_col     = args.r1_col
    r2_col     = args.r2_col

    # Load config
    proj_config = ConfigDf(config_file = config)

    # Define output columns & add to config
    fastp_outputs = ["r1_gz", "r2_gz"]
    for output in fastp_outputs:
        proj_config.add_column(name = output, default_value = None)
    
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
        })

    # Start the multiprocessing pool
    print(f"Starting multiprocessing pool with {processes} workers...")
    with mp.Pool(processes=processes) as pool:
        for result in pool.imap_unordered(sample_worker, tasks):
            
            index = result["index"]
            
            if result["success"]:
                # Update config with output files
                proj_config.update_row(index, "r1_gz", result["r1_gz"])
                proj_config.update_row(index, "r2_gz", result["r2_gz"])
            else:
                print(f"Error processing row {index}: {result['error']}")
            
            # Save updated config after processing each row
            proj_config.save()

if __name__ == "__main__":
    main()