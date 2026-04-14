#! python

"""
- Run host read filtering on a set of samples.
- Supply one or two (for sequential) host databases in mmi format for read filtering.
- Runs minimap2 -ax sr --secondary=no {threads} {db.mmi} {r1} {r2} |
- Piped sam formatted alignment records are parsed and filtered using fastfilter.
- Reads that pass (unaligned or low quality alignments) are written to paired fastq.gz files.
"""

import argparse
from pathlib import Path
import multiprocessing as mp
import json
from typing import Optional, Tuple
from dataclasses import dataclass, field
from sb_projects.wrapper import Wrapper
from sb_projects.config import ConfigDf
from sb_projects.subprocesses import run_check_output_to_str
from sb_projects.file_utilities import delete_file


# Filtering Defaults >>> These could be converted to args.
MAX_AP = 0.5   # float:   Alignment_Percentage = Alignment_Length / Read_Length
MAX_PI = 0.5   # float:   percent identity = Num_Identities / Alignment_Length
MAX_AS = 150   # integer: SAM tag -> Alignment_Score
MAX_AL = 75    # integer: CIGAR -> Alignment_Length
MAX_SL = 1.0   # float:   SL = Alignment_Score / Alignment_Length
#MAX_MQ = 10  # integer: MAPQ score

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to run host read filtering on a set of paired reads.")

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
    
    parser.add_argument("--db1", type=str, required=True,
                        help="1st pre-computed Minimap database path (str).")
    
    parser.add_argument("--db2", type=str, default=None,
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

# Wrapper subclass for running HRF pipeline: minimap2 | fss2pfq 
@dataclass(kw_only=True)
class HRFsr(Wrapper):
    input_mmi:      Path | str      = field(metadata={'type': 'input_file'})
    r1:             Path | str      = field(metadata={'type': 'input_file'})
    r2:             Path | str      = field(metadata={'type': 'input_file'})
    fq_prefix:      Path | str      = field(metadata={'type': 'output_file'})
    map_threads:    Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    filter_threads: Optional[int]   = field(default=4, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    max_ap:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max-ap {value}'})
    max_pi:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max-pi {value}'})
    max_as:         Optional[int]   = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max-as {value}'})
    max_al:         Optional[int]   = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max-al {value}'})
    max_sl:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max-sl {value}'})
   #max_mq:         Optional[float] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '--max-mq {value}'})
    cmd:            str             = "minimap2 -ax sr --secondary=no {map_threads} {input_mmi} {r1} {r2} | \
                                      fastfilter {filter_threads} {max_ap} {max_pi} {max_as} {max_al} {max_sl} --fq-prefix {fq_prefix}"
   #cmd:            str             = "minimap2 -ax sr --secondary=no {map_threads} {input_mmi} {r1} {r2} | \
   #                                  fastfilter {filter_threads} {max_ap} {max_pi} {max_as} {max_al} {max_sl} {max_mq} --fq-prefix {fq_prefix}"

# Run HRF on a set of r1 and r2 <sample>_<r#>.fq.gz
def run_hrf(
        input_mmi: Path,
        in_r1:     Path,
        in_r2:     Path,
        fq_prefix: Path,
        threads:   int,
        dry_run:   bool,
    ) -> Tuple[int, int]:
    
    # Init a HRFsr Wrapper obj.
    process = HRFsr(
        input_mmi      = input_mmi,
        r1             = in_r1,
        r2             = in_r2,
        fq_prefix      = fq_prefix,
        map_threads    = threads - 4 if threads >=12 else threads - 2,
        filter_threads = 4 if threads >=12 else 2,
        max_ap         = MAX_AP,
        max_pi         = MAX_PI,
        max_as         = MAX_AS,
        max_al         = MAX_AL,
        max_sl         = MAX_SL,
        #max_mq         = MAX_MQ,
        dry_run        = dry_run,
    )
    
    counts_json = run_check_output_to_str(
        formatted_command = process.build(),
        dry_run           = dry_run,
    )
    
    if dry_run:
        counts_json = json.dumps(
            {"total_pairs":42, "written_pairs":11}
        )
    
    data = json.loads(counts_json)
    original_pairs = int(data["total_pairs"])
    filtered_pairs = int(data["written_pairs"])
    
    return original_pairs, filtered_pairs


def sample_worker(task_args: dict) -> dict:
    """Worker function to process a single sample. Returns a dictionary with the results."""
    dry_run     = task_args["dry_run"]
    index       = task_args["index"]
    row         = task_args["row"]
    input_dir   = task_args["input_dir"]
    output_dir  = task_args["output_dir"]
    threads     = task_args["threads"]
    sample_col  = task_args["sample_col"]
    r1_col      = task_args["r1_col"]
    r2_col      = task_args["r2_col"]
    db_list     = task_args["db_list"]
    db_names    = task_args["db_names"]
    hrf_outputs = task_args["hrf_outputs"]

    if dry_run: 
        print(f"\nWorker processing: {index}\t{row}") 
    
    SAMPLE = row[sample_col]
    
    # Set up input dicts
    if len(db_names) == 1:
        inputs = {
            1: { 
                "input_mmi": db_list[0],
                "threads"  : threads,
                "in_r1"    : input_dir / row[r1_col],
                "in_r2"    : input_dir / row[r2_col],
                "fq_prefix": output_dir / f"{SAMPLE}_hrf",
                "dry_run"  : dry_run,
            }
        }
        intermediates = []
    else:
        inputs = {
            1: { 
                "input_mmi": db_list[0],
                "threads"  : threads,
                "in_r1"    : input_dir / row[r1_col],
                "in_r2"    : input_dir / row[r2_col],
                "fq_prefix": output_dir / f"{SAMPLE}_{db_names[0]}",
                "dry_run"  : dry_run,
            },
            2: { 
                "input_mmi": db_list[1],
                "threads"  : threads,
                "in_r1"    : output_dir / f"{SAMPLE}_{db_names[0]}.r1.fq.gz",
                "in_r2"    : output_dir / f"{SAMPLE}_{db_names[0]}.r2.fq.gz",
                "fq_prefix": output_dir / f"{SAMPLE}_hrf",
                "dry_run"  : dry_run,
            },
        }
        intermediates = [
            output_dir / f"{SAMPLE}_{db_names[0]}.r1.fq.gz",
            output_dir / f"{SAMPLE}_{db_names[0]}.r2.fq.gz",
        ]

    # Set up output dict
    results = {"index": index, "success": False, "error": None}
    for o in hrf_outputs:
        results[o] = None

    # For each input call run_hrf()
    for i in range(len(db_names)):
        db_name = db_names[i]
        try:
            hrf_inputs = inputs[i+1] 
            original_pairs, filtered_pairs = run_hrf(**hrf_inputs)

        except Exception as e:
        
            results["error"] = str(e)
            # Return failure state
            return results
        
        else:
            results[f"{db_name}_unfiltered_pairs"] = original_pairs
            results[f"{db_name}_filtered_pairs"] = filtered_pairs
        
    results["success"] = True
    results["hrf_r1"] = f"{SAMPLE}_hrf.r1.fq.gz"
    results["hrf_r2"] = f"{SAMPLE}_hrf.r2.fq.gz"
    if len(intermediates) > 0:
        for f in intermediates: delete_file(f, dry_run = dry_run)
    
    return results

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
    db1        = Path(args.db1)
    db1_name   = db1.stem
    db2        = Path(args.db2) if args.db2 is not None else None
    db2_name   = db2.stem if db2 is not None else None
    #max_ap     = MAX_AP
    #max_pi     = MAX_PI
    #max_as     = MAX_AS
    #max_al     = MAX_AL
    #max_sl     = MAX_SL
    #max_mq     = MAX_MQ

    # Load config
    proj_config = ConfigDf(config_file = config)

    # Define output columns & add to config
    db_list  = [x for x in [db1, db2] if x is not None]
    db_names = [x for x in [db1_name, db2_name] if x is not None]
    # original_pairs, filtered_pairs for each db + final r1 and r2
    hrf_outputs = []
    for db_name in db_names:
        hrf_outputs.append(f"{db_name}_unfiltered_pairs")
        hrf_outputs.append(f"{db_name}_filtered_pairs")
    hrf_outputs.append("hrf_r1")
    hrf_outputs.append("hrf_r2")
            
    for output in hrf_outputs:
        proj_config.add_column(name = output, default_value = None)
    
    # Build worker pool task list 
    tasks = []
    for index, row in proj_config:
        tasks.append({
            "index":       index,
            "row":         row,
            "input_dir":   input_dir,
            "output_dir":  output_dir,
            "threads":     threads,
            "dry_run":     dry_run,
            "sample_col":  sample_col,
            "r1_col":      r1_col,
            "r2_col":      r2_col,
            "db_list":     db_list,
            "db_names":    db_names,
            "hrf_outputs": hrf_outputs,
        })

    # Start the multiprocessing pool
    print(f"Starting multiprocessing pool with {processes} workers...")
    with mp.Pool(processes=processes) as pool:
        for result in pool.imap_unordered(sample_worker, tasks):
            
            index = result["index"]
            
            if result["success"]:
                # Update config with output files and counts
                for output in hrf_outputs:
                    proj_config.update_row(index, output, result[output])

                
            else:
                print(f"Error processing row {index}: {result['error']}")
            
            # Save updated config after processing each row
            proj_config.save()

if __name__ == "__main__":
    main()