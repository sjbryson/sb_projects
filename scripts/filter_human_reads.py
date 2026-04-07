#! python

"""
Samuel Joseph Bryson
Copyright 2026

- 1 - Run human depletion using 2 refs (hg38 and chm13) on a set of samples (r1.fq.gz & r2.fq.gz).
- 2 - Keep unmapped pairs; this yields filtered_r1.fq and filtered.r2.fq.
"""

from pathlib import Path
from typing import Tuple
import json
from sb_projects.config_manager import ConfigManager
from sb_projects.subprocess_utilities import run_check_output_to_str
from sb_projects.fastx_utilities import fastq_count_reads
from sb_projects.file_utilities import delete_file, count_lines
from sb_projects.wrapper import MinimapFilterPairedFastq

# Database paths
DB_DIR      = Path.home() / "bio_db"
HG13        = DB_DIR      / "T2T_CHM13/srHG13.mmi"
HG38        = DB_DIR      / "GRCh38.p14/srHG38.mmi" 
# Project paths
PROJ_DIR    = Path.home() / "NIH_HVP/NovogeneComp"
CONFIG      = PROJ_DIR    / "sample_config.txt"
#TEST_CONFIG = PROJ_DIR    / "sample_config copy.txt"
READS_DIR   = PROJ_DIR    / "clean_fastq"
FILT_READS  = PROJ_DIR    / "filtered_fastq"
TEMP_DIR    = PROJ_DIR    / "temp"
# Workflow vars
DRY_RUN     = False
#DRY_RUN     = True
THREADS     = 16

# Map each sample to human genome reference write filtered r1 and r2 fastq files
def filter_human_reads(
        db:        Path, 
        p:         int,
        input_r1:  Path,
        input_r2:  Path,
        output_r1: Path,
        output_r2: Path,
    ) -> Tuple[int, int, Path, Path]:
    # map paired fastq to reference human genome
    process = MinimapFilterPairedFastq(
        input_mmi   = db,
        r1          = input_r1,
        r2          = input_r2,
        filtered_r1 = output_r1,
        filtered_r2 = output_r2,
        threads     = p,
        dry_run     = DRY_RUN
    )
    counts_json = run_check_output_to_str(
        formatted_command = process.build(),
        dry_run           = DRY_RUN
    )
    if DRY_RUN:
        counts_json = json.dumps(
            {'written_pairs':11, 'total_pairs':42}
        )
    data = json.loads(counts_json)
    filtered_pairs = int(data['written_pairs'])
    original_pairs = int(data['total_pairs'])
    return original_pairs, filtered_pairs, output_r1, output_r2


def main():
    # Load config
    proj_config = ConfigManager(config_file=CONFIG)
    #proj_config = ConfigManager(config_file=TEST_CONFIG)
    workflow_outputs = [
        "clean_read_count",
        "chm13_filtered_read_count",
        "chm13_filtered_r1",
        "chm13_filtered_r2",
        "chm13_h38_filtered_read_count",
        "chm13_h38_filtered_r1",
        "chm13_h38_filtered_r2",
    ]
    for output in workflow_outputs:
        proj_config.add_column(name=output)

    for index, row in proj_config:
        if DRY_RUN == True: 
            print(f"\n{index}\t{row}") 
        SAMPLE = row["sample"]
        # filter on CHM13
        chm13_f_r1 = TEMP_DIR / f"{SAMPLE}.chm13_filtered_r1.fq.gz"
        chm13_f_r2 = TEMP_DIR / f"{SAMPLE}.chm13_filtered_r2.fq.gz"
        clean_read_count, read_count1, temp_r1, temp_r2 = filter_human_reads(
        db        = HG13, 
        p         = THREADS, 
        input_r1  = READS_DIR / f"{row["r1_fastq"]}",
        input_r2  = READS_DIR / f"{row["r2_fastq"]}",
        output_r1 = chm13_f_r1,
        output_r2 = chm13_f_r2,
        )
        proj_config.update_row(index, "clean_read_count", clean_read_count)
        proj_config.update_row(index, "chm13_filtered_read_count", str(read_count1))
        proj_config.update_row(index, "chm13_filtered_r1", temp_r1)
        proj_config.update_row(index, "chm13_filtered_r2", temp_r2)
        # filter on h38
        _, read_count2, filtered_r1, filtered_r2 = filter_human_reads(
        db        = HG38, 
        p         = THREADS, 
        input_r1  = chm13_f_r1,
        input_r2  = chm13_f_r2,
        output_r1 = FILT_READS / f"{SAMPLE}.chm13_h38_filtered_r1.fq.gz",
        output_r2 = FILT_READS / f"{SAMPLE}.chm13_h38_filtered_r2.fq.gz",
        )
        proj_config.update_row(index, "chm13_h38_filtered_read_count", str(read_count2))
        proj_config.update_row(index, "chm13_h38_filtered_r1", filtered_r1)
        proj_config.update_row(index, "chm13_h38_filtered_r2", filtered_r2)
        # save updated config after processing each row
        proj_config.save()
        # delete temp fastq files
        delete_file(file_path=chm13_f_r1, dry_run=DRY_RUN)
        delete_file(file_path=chm13_f_r2, dry_run=DRY_RUN)

if __name__ == "__main__":
    main()