#! python

"""
- Run Penguin assembler on a set of samples.
- github: https://github.com/soedinglab/plass
- PenguiN: Jochheim A, Jochheim FA, Kolodyazhnaya A, Morice E, Steinegger M, Soeding J. 
  Strain-resolved de-novo metagenomic assembly of viral genomes and microbial 16S rRNAs. Microbiome 12, 187, (2024)
"""

import argparse
from pathlib import Path
import multiprocessing as mp
from typing import Optional
from dataclasses import dataclass, field
from sb_projects.wrapper import Wrapper
from sb_projects.config import ConfigDf
from sb_projects.subprocesses import run_check_call

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Script to run host read filtering on a set of paired reads.")

    parser.add_argument("--config_file", type=str, required=True, 
                        help="Path to the configuration file (str).")
    
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Directory containing input files (str).")
    
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Directory where output files will be saved (str).")

    parser.add_argument("--sample_col", type=str, default="sample", 
                        help="Column name for samples in config file (str).")
    
    parser.add_argument("--r1_col", type=str, default="r1", 
                        help="Column name for input R1 in config file (str).")
    
    parser.add_argument("--r2_col", type=str, default="r2", 
                        help="Column name for input R2 in config file (str).")
    
    parser.add_argument("--min-contig-len", type=int, default=4, 
                        help="Minimum contig length to keep (int, default: 500).")
    
    parser.add_argument("--dry_run", action="store_true", default=False, 
                        help="If set, perform a trial run printing commands (default: False).")

    return parser.parse_args()

# Wrapper subclass for running HRF pipeline: minimap2 | fss2pfq 
@dataclass(kw_only=True)
class PenguinAsm(Wrapper):
    r1:         Path | str    = field(metadata={'type': 'input_file'})
    r2:         Path | str    = field(metadata={'type': 'input_file'})
    assembly:   Path | str    = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(metadata={'type': 'value_flag', 'flag_fmt': '--threads {value}'})
    min_length: Optional[int] = field(default=500)
    cmd:        str           = "penguin guided_nuclassemble --min-contig-len {min_length} {r1} {r2} {assembly} tmp"

    
# Assembler Args
#--min-length 30
#--cluster-mode 1
#--min-contig-len 500
#--contig-output-mode 0


def main():
    
    # Parse args
    args       = parse_args()
    config     = Path(args.config_file)
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    dry_run    = True if args.dry_run else False
    sample_col = args.sample_col
    r1_col     = args.r1_col
    r2_col     = args.r2_col
    min_length = args.min_length

    # Load config
    proj_config = ConfigDf(config_file = config)

    # Define output column & add to config
    proj_config.add_column(name = "assembly", default_value = None)
    
   # Run Penguin assembler on each sample
    for index, row in proj_config:
        SAMPLE = row[sample_col]
        assembly = f"{SAMPLE}.fas"

        if dry_run: 
            print(f"\nWorker processing: {index}\t{row}") 

        try:
            process = PenguinAsm(
                r1         = input_dir / row[r1_col],
                r2         = input_dir / row[r2_col],
                assembly   = output_dir / assembly,
                min_length = min_length,
                dry_run    = dry_run,
            )
            cmd = process.build()
            run_check_call(
                formatted_command = cmd,
                devnull           = False,
                dry_run           = dry_run,
            )
            # define success state
            result = {"success": True, "error": None}
        
        except Exception as e:
            # define failure state
            result = {"success": False, "error": str(e)}
        
        finally:
            
            if result["success"]:
                # Update config with output files
                proj_config.update_row(index, "assembly", assembly)
            else:
                print(f"Error processing row {index}: {result['error']}")
            
            # Save updated config after processing each row
            proj_config.save()

if __name__ == "__main__":
    main()