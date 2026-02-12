#!/usr/bin/env python3

"""
Samuel Joseph Bryson

- 1 - Run human depletion using 2 refs (hg38 and chm13) on a set of samples (r1.fq.gz & r2.fq.gz).
- 2 - All mapped pairs to either are removed; this yields filtered_r1.fq.gz and filtered.r2.fq.gz.
- 3 - Filtered reads are mapped to a set of reference bacterial genomes of interest + one for FDR.
        - Treponema_denticola       (GCF_000008185.1_ASM818v1)      NC_002967.9
        - Tannerella_forsythia      (GCF_000238215.1_ASM23821v1)    NC_016610.1
        - Porphyromonas_gingivalis  (GCF_000010505.1_ASM1050v1)     NC_010729.1
        - Fusobacterium_nucleatum   (GCF_003019295.1_ASM301929v1)   NZ_CP028101.1
        -- Nitrosopumilus_maritimus (GCF_000018465.1_ASM1846v1)     NC_010085.1
- 4 - Mapped read counts and coverage stats are calculated for each sample x ref genome.
"""

from pathlib import Path

from src.config_manager import ConfigManager
from src.subprocess_utilities import run_check_call
from src.file_utilities import delete_file
from src.tools import (
    Minimap2SRHumanDepletion,
    Minimap2MapSRtoBam,
    SamtoolsSort,
    SamtoolsIndex,
    GetProperMappedIDs,
    UniqueIdList,
    SamtoolsMpileup,
    GetContigPerBaseCoverageFromPileup,
    SBFilterFastq
)

##############
DRY_RUN = True
##############
THREADS = 1
##############
# Database paths
DB_DIR = Path.home() / "bio_db"
HG38   = DB_DIR      / "GRCh38.p14/srHG38.mmi" 
HG13   = DB_DIR      / "T2T_CHM13/srHG13.mmi"
BACDB  = DB_DIR      / "R01_Td/srREFS.mmi"

# Data paths
DATA_DIR       = Path.home() / "NIH_HVP/HOGBV_AD/SAMPLES"
READS_DIR      = DATA_DIR    / "qc_reads"
FILT_READS_DIR = DATA_DIR    / "filtered_reads"
BAM_DIR        = DATA_DIR    / "hg_bams"

# Project paths
PROJ_DIR = Path.home() / "UCLA/R01_prelim"
TEMP_DIR = PROJ_DIR    / "temp"
CONFIG   = PROJ_DIR    / "sample_config_R01.txt"

# Map each sample to human genome reference and get mapped mate ids
def map_to_human(record: dict, db: Path, db_name: str, p: int) -> None:
    SAMPLE = record["sample"]
    # map paired fastq to reference human genome
    step1 = Minimap2SRHumanDepletion(
        input_mmi  = db,
        r1         = READS_DIR / f"{record["r1_fastq"]}",
        r2         = READS_DIR / f"{record["r2_fastq"]}",
        output_bam = TEMP_DIR / f"{SAMPLE}_{db_name}.bam",
        threads    = p,
        dry_run    = DRY_RUN,
        
    )
    run_check_call(
        formatted_command = step1.build(),
        devnull           = True,
        dry_run           = True
    )
    # sort bam
    step2 = SamtoolsSort(
        sorted_bam = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        input_bam  = TEMP_DIR / f"{SAMPLE}_{db_name}.bam",
        threads    = p,
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step2.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # index bam
    step3 = SamtoolsIndex(
        sorted_bam = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step3.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # delete the unsorted bam file
    unsorted_bam = TEMP_DIR / f"{SAMPLE}_{db_name}.bam"
    delete_file(file_path=unsorted_bam, dry_run=DRY_RUN)
    # get the ids for mapped pairs
    step4 = GetProperMappedIDs(
        sorted_bam  = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        output_file = TEMP_DIR / f"{SAMPLE}_{db_name}.mapped_ids.txt",
        threads     = p,
        dry_run     = DRY_RUN
    )
    run_check_call(
        formatted_command = step4.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )

    
def filter_fastq(record: dict, db1_name: str, db2_name: str) -> None:
    SAMPLE = record["sample"]
    # get a list of all the paired reads that mapped to either reference human genome
    step1 = UniqueIdList(
        file_a     = FILT_READS_DIR / f"{SAMPLE}_{db1_name}.mapped_ids.txt",
        file_b     = FILT_READS_DIR / f"{SAMPLE}_{db2_name}.mapped_ids.txt",
        output_txt = FILT_READS_DIR / f"{SAMPLE}_combined_human_mapped_ids.txt",
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step1.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # filter the qc_reads files, remove reads mapped to human refs
    step2 = SBFilterFastq(
        r1_fastq      = READS_DIR / f"{record["r1_fastq"]}",
        r2_fastq      = READS_DIR / f"{record["r2_fastq"]}",
        input_file    = FILT_READS_DIR / f"{SAMPLE}_combined_human_mapped_ids.txt",
        output_prefix = FILT_READS_DIR / f"{SAMPLE}_filtered",
        gz            = True,
        dry_run       = DRY_RUN
    )
    run_check_call(
        formatted_command = step2.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )

def map_filtered_fastq(record: dict, db: Path, db_name: str, p: int) -> None:
    SAMPLE = record["sample"]
    # map the human depleted reads to the bacteria reference genomes
    step1 = Minimap2MapSRtoBam(
        input_mmi = db,
        r1 = FILT_READS_DIR / f"{SAMPLE}_filtered_R1.fastq.gz",
        r2 = FILT_READS_DIR / f"{SAMPLE}_filtered_R2.fastq.gz",
        output_bam = TEMP_DIR / f"{SAMPLE}_HGF_{db_name}.bam",
        threads    = p,
        dry_run    = DRY_RUN,
    )
    run_check_call(
        formatted_command = step1.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # sort bam
    step2 = SamtoolsSort(
        sorted_bam = BAM_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.bam",
        input_bam  = TEMP_DIR / f"{SAMPLE}_HGF_{db_name}.bam",
        threads    = p,
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step2.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # index bam
    step3 = SamtoolsIndex(
        sorted_bam = BAM_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.bam",
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step3.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # delete the unsorted bam file
    unsorted_bam = TEMP_DIR / f"{SAMPLE}_HGF_{db_name}.bam",
    delete_file(file_path=unsorted_bam, dry_run=DRY_RUN)



def get_coverage_stats() -> None:
    pass

def main():
    # Load config
    proj_config = ConfigManager(config_file=CONFIG)
    for row in proj_config:
        print(f"\n{row}")
        # h38_bam -> TEMP_DIR -> sort -> BAM_DIR -> index -> id list
        map_to_human(row, HG38, "hg38", THREADS)
        # hg13_bam -> TEMP_DIR -> sort -> BAM_DIR -> index -> id list
        map_to_human(row, HG13, "hg13", THREADS)
        # extract union of mapped ids and filter r1 & r2 fastq 
        filter_fastq(row, "hg38", "hg13", THREADS)
        # map filtered fastq's to bacterial refs
        map_filtered_fastq(row, BACDB, "R01Td", THREADS)

# get coverage for refs


# plot

if __name__ == "__main__":
    main()