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
from sb_projects.config_manager import ConfigManager
from sb_projects.subprocess_utilities import run_check_call

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to convert paired fastq.bz2 to fq.gz.")

    parser.add_argument("--config_file", type=str, required=True, 
                        help="Path to the configuration file (str).")
    
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Directory containing input files (str).")
    
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Directory where output files will be saved (str).")
    
    parser.add_argument("--temp_dir", type=str, default="temp", 
                        help="Temporary directory for intermediate files (str).")
    
    parser.add_argument("--threads", type=int, default=4, 
                        help="Number of threads to use (int, default: 4).")
    
    parser.add_argument("--processes", type=int, default=1, 
                        help="Number of parallel subprocesses (int, default: 1).")
    
    parser.add_argument("--sample_col", type=str, default="sample", 
                        help="Column name for samples in config file (str).")
    
    parser.add_argument("--r1_col", type=str, default="r1", 
                        help="Column name for input R1 in config file (str).")
    
    parser.add_argument("--r2_col", type=str, default="r2", 
                        help="Column name for input R2 in config file (str).")

    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")
    
    ## --- ToDo: --- ##
    # add output r1 column name --> currently "hrf_r1_gz"
    # add output r2 column name --> currently "hrf_r2_gz"

    return parser.parse_args()

@dataclass(kw_only=True)
class Minimap2SRHumanDepletion(Wrapper):
    input_mmi:  Path | str    = field(metadata={'type': 'input_file'})
    r1:         Path | str    = field(metadata={'type': 'input_file'})
    r2:         Path | str    = field(metadata={'type': 'input_file'})
    output_bam: Path | str    = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:        str           = "minimap2 -ax sr --secondary=no {threads} {input_mmi} {r1} {r2} | samtools sort -o - > {output_bam}"

    ## Fix to use stream_filter --> fss2pfq
'''## 1st v T2T-CHM13
% minimap2 -ax sr --secondary=no -t 12 ~/bio_db/T2T_CHM13/srHG13.mmi \
HELTHY_B1-A_R1.fq.gz  HELTHY_B1-A_R2.fq.gz | \
stream_filter -t 12 -p HELTHY_B1-A.stream_test.m150 -m 150

{
  "max_AS": 150,
  "orphaned_reads": 206206,
  "runtime_seconds": 82.008669833,
  "total_pairs": 4998184,
  "written_pairs": 57430
}
'''