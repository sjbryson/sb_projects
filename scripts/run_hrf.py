#! python

"""
Samuel Joseph Bryson
Copyright 2026

- Run pfqbz2gz on a set of samples; convert bz2 compressed paired fastq to fq.gz.
"""

import argparse
from pathlib import Path
import multiprocessing as mp
import json
from typing import Optional, Tuple
from dataclasses import dataclass, field
from sb_projects.wrapper import Wrapper
from sb_projects.config_manager import ConfigManager
from sb_projects.subprocess_utilities import run_check_output_to_str

# Database paths
DB_DIR      = Path.home() / "bio_db"
HG13        = DB_DIR      / "T2T_CHM13/srHG13.mmi"
HG38        = DB_DIR      / "GRCh38.p14/srHG38.mmi" 
# Filtering Defaults
MAX_AP = None # float: Alignment_Percentage = Alignment_Length / Read_Length
MAX_PI = None # float: percent identity = Num_Identities / Alignment_Length
MAX_AS = None # integer: SAM tag -> Alignment_Score
MAX_AL = None # integer: CIGAR -> Alignment_Length
MAX_SL = None # float: SL = Alignment_Score / Alignment_Length

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to run fastq QC on a set of paired reads.")

    parser.add_argument("--config_file", type=str, required=True, 
                        help="Path to the configuration file (str).")
    
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Directory containing input files (str).")
    
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Directory where output files will be saved (str).")
    
    parser.add_argument("--temp_dir", type=str, required=True, 
                        help="Temporary directory where intermediate files will be saved (str).")

    parser.add_argument("--threads", type=int, default=4, 
                        help="Number of threads to use (int, default: 4).")
    
    parser.add_argument("--processes", type=int, default=2, 
                        help="Number of parallel samples to process at once (int, default: 2).")
    
    parser.add_argument("--mmiA", type=str, required=True,
                        help="1st pre-computed Minimap database path (str).")
    
    parser.add_argument("--mmiB", type=str, required=False,
                        help="Optional: 2nd pre-computed Minimap database path (str).")
    
    parser.add_argument("--sample_col", type=str, default="sample", 
                        help="Column name for samples in config file (str).")
    
    parser.add_argument("--r1_col", type=str, default="r1", 
                        help="Column name for input R1 in config file (str).")
    
    parser.add_argument("--r2_col", type=str, default="r2", 
                        help="Column name for input R2 in config file (str).")
    
    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")

    return parser.parse_args()

    ## --- ToDo: --- ##
    # add output r1 column name --> currently "hrf_r1_gz"
    # add output r2 column name --> currently "hrf_r2_gz"

# Wrapper subclass for running HRF pipeline: minimap2 | fss2pfq 
@dataclass(kw_only=True)
class HRFsr(Wrapper):
    input_mmi:      Path | str      = field(metadata={'type': 'input_file'})
    r1:             Path | str      = field(metadata={'type': 'input_file'})
    r2:             Path | str      = field(metadata={'type': 'input_file'})
    fq_prefix:      Path | str      = field(metadata={'type': 'output_file'})
    map_threads:    Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    filter_threads: Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    max_ap:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max_ap {value}'})
    max_pi:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max_pi {value}'})
    max_as:         Optional[int]   = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max_as {value}'})
    max_al:         Optional[int]   = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max_al {value}'})
    max_sl:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max_sl {value}'})
    cmd:            str             = "minimap2 -ax sr --secondary=no {map_threads} {input_mmi} {r1} {r2} | \
                                      fss2pfq {filter_threads} {max_ap} {max_pi} {max_as} {max_al} {max_sl} -p {fq_prefix}"



# Run fastp on a set of r1 and r2 <sample>_<r#>.fq.gz
def run_hrf(
        input_mmi:   Path,
        in_r1:       Path,
        in_r2:       Path,
        out_p:       Path,
        threads:     int,



        dry_run:     bool,
    ) -> Tuple[int, int]:
    # Init a HRFsr Wrapper obj.
    process = HRFsr(
        threads       = threads, 
        in_r1         = in_r1,
        in_r2         = in_r2,
        
        dry_run       = dry_run,
    )
    
    counts_json = run_check_output_to_str(
        formatted_command = process.build(),
        dry_run           = dry_run,
    )
    if dry_run == True:
        counts_json = json.dumps(
            {'written_pairs':11, 'total_pairs':42}
        )
    data = json.loads(counts_json)
    filtered_pairs = int(data['written_pairs'])
    original_pairs = int(data['total_pairs'])
    return original_pairs, filtered_pairs


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
        original_pairs, filtered_pairs = run_hrf(
            threads    = threads,
            in_r1      = in_r1,
            in_r2      = in_r2,
            out_prefix = out_prefix,
            dry_run    = dry_run,
        )

        #######
        ## optional second db scan ##


        # Return success state
        return {"index": index, "success": True, "r1_gz": r1_gz, "r2_gz": r2_gz, "error": None}
    
    except Exception as e:
        # Return failure state
        return {"index": index, "success": False, "r1_gz": None, "r2_gz": None, "error": str(e)}





def main():
    pass

'''
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
    proj_config = ConfigManager(config_file = config)

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

'''

if __name__ == "__main__":
    main()