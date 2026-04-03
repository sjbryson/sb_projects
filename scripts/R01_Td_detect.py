#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026

- 1 - Run human depletion using 2 refs (hg38 and chm13) on a set of samples (r1.fq.gz & r2.fq.gz).
- 2 - All mapped pairs to either are removed; this yields filtered_r1.fq.gz and filtered.r2.fq.gz.
- 3 - Filtered reads are mapped to a set of 4 reference bacterial genomes of interest + 4 controls for FDR.
        - Treponema_denticola           (GCF_000008185.1_ASM818v1)          NC_002967.9
        - Tannerella_forsythia          (GCF_000238215.1_ASM23821v1)        NC_016610.1
        - Porphyromonas_gingivalis      (GCF_000010505.1_ASM1050v1)         NC_010729.1
        - Fusobacterium_nucleatum       (GCF_003019295.1_ASM301929v1)       NZ_CP028101.1
        -- Ralstonia_pickettii          (GCF_902374465.1_MGYG-HGUT-01384)   multiple
        -- Staphylococcus_epidermidis   (GCF_006094375.1_ASM609437v1)       multiple
        -- Agrobacterium tumefaciens    (GCF_013318015.2_ASM1331801v2)      multiple
        -- Nitrosopumilus_maritimus (GCF_000018465.1_ASM1846v1)             NC_010085.1
- 4 - Mapped read counts and coverage stats are calculated for each sample x ref genome.
"""

from pathlib import Path
from typing import Tuple
from src.config_manager import ConfigManager
from src.subprocess_utilities import run_check_call
from src.fastx_utilities import fastq_count_reads
from src.file_utilities import delete_file, count_lines
from wrappers import (
    Minimap2SRHumanDepletion,
    Minimap2MapSRtoBam,
    SamtoolsSort,
    SamtoolsIndex,
    GetProperMappedIDs,
    UniqueIdList,
    SamtoolsIdxstats,
    SamtoolsMpileup,
    SBFilterFastq
)

# Database paths
DB_DIR = Path.home() / "bio_db"
HG38   = DB_DIR      / "GRCh38.p14/srHG38.mmi" 
HG13   = DB_DIR      / "T2T_CHM13/srHG13.mmi"
BACDB  = DB_DIR      / "R01_Td/srREFS.mmi"
BACFA  = DB_DIR      / "R01_Td/REFS.fasta"

# Data paths
DATA_DIR       = Path.home() / "NIH_HVP/HOGBV_AD/SAMPLES"
READS_DIR      = DATA_DIR    / "qc_reads"
FILT_READS_DIR = DATA_DIR    / "filtered_reads"
BAM_DIR        = DATA_DIR    / "hg_bams"
SEQID_DIR      = DATA_DIR    / "hg_ids"

# Project paths
PROJ_DIR    = Path.home() / "UCLA/R01_prelim"
TEMP_DIR    = PROJ_DIR    / "temp"
CONFIG      = PROJ_DIR    / "sample_config_R01.txt"
TEST_CONFIG = PROJ_DIR    / "sample_config_R01 copy.txt"
OUT_DIR     = PROJ_DIR    / "data"


DRY_RUN = True
THREADS = 16


# Map each sample to human genome reference and get mapped mate ids
def map_to_human(record: dict, db: Path, db_name: str, p: int) -> Tuple[Path, Path]:
    SAMPLE     = record["sample"]
    sorted_bam = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam"
    mapped_ids = SEQID_DIR / f"{SAMPLE}_{db_name}.mapped_ids.txt"
        # map paired fastq to reference human genome
    step1 = Minimap2SRHumanDepletion(
        input_mmi  = db,
        r1         = READS_DIR / f"{record["r1_fastq"]}",
        r2         = READS_DIR / f"{record["r2_fastq"]}",
        output_bam = TEMP_DIR / f"{SAMPLE}_{db_name}.bam",
        threads    = p,
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step1.build(),
        devnull           = True,
        dry_run           = DRY_RUN
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
        output_file = SEQID_DIR / f"{SAMPLE}_{db_name}.mapped_ids.txt",
        threads     = p,
        dry_run     = DRY_RUN
    )
    run_check_call(
        formatted_command = step4.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    return sorted_bam, mapped_ids
    
def filter_fastq(record: dict, db1_name: str, db2_name: str) -> Tuple[Path, Path, Path]:
    SAMPLE = record["sample"]
    mapped_ids = SEQID_DIR / f"{SAMPLE}_union_mapped_ids.txt"
    filtered_r1 = FILT_READS_DIR / f"{SAMPLE}_filtered_R1.fastq.gz"
    filtered_r2 = FILT_READS_DIR / f"{SAMPLE}_filtered_R2.fastq.gz"
    # get a list of all the paired reads that mapped to either reference human genome
    step1 = UniqueIdList(
        file_a     = SEQID_DIR / f"{SAMPLE}_{db1_name}.mapped_ids.txt",
        file_b     = SEQID_DIR / f"{SAMPLE}_{db2_name}.mapped_ids.txt",
        output_txt = SEQID_DIR / f"{SAMPLE}_union_mapped_ids.txt",
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
        input_file    = SEQID_DIR / f"{SAMPLE}_union_mapped_ids.txt",
        output_prefix = FILT_READS_DIR / f"{SAMPLE}_filtered",
        gz            = True,
        dry_run       = DRY_RUN
    )
    run_check_call(
        formatted_command = step2.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    return mapped_ids, filtered_r1, filtered_r2


def map_filtered_fastq(record: dict, db: Path, db_name: str, p: int) -> Path:
    SAMPLE = record["sample"]
    sorted_bam = BAM_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.bam"
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
    unsorted_bam = TEMP_DIR / f"{SAMPLE}_HGF_{db_name}.bam"
    delete_file(file_path=unsorted_bam, dry_run=DRY_RUN)
    return sorted_bam


def get_coverage_stats(record: dict, db_fasta: Path, db_name: str) -> Tuple[Path, Path]:
    SAMPLE = record["sample"]
    idxstats_tsv = OUT_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.idxstats.tsv"
    pileup_txt = OUT_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.pileup.txt"
    # run idxstats
    step1 = SamtoolsIdxstats(
    sorted_bam = BAM_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.bam",
    output_tsv = OUT_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.idxstats.tsv",
    dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step1.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # run mpileup
    step2 = SamtoolsMpileup(
        refernce_fasta = db_fasta,
        sorted_bam     = BAM_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.bam",
        output_pileup  = OUT_DIR / f"{SAMPLE}_HGF_{db_name}.sorted.pileup.txt",
        dry_run        = DRY_RUN
    )
    run_check_call(
        formatted_command = step2.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    return idxstats_tsv, pileup_txt

def main():
    # Load config
    proj_config = ConfigManager(config_file=CONFIG)
    #proj_config = ConfigManager(config_file=TEST_CONFIG)
    workflow_outputs = [
        "read_count",
        "hg38_sorted_bam",
        "hg38_mapped_ids",
        "hg38_mapped_count",
        "hg13_sorted_bam",
        "hg13_mapped_ids",
        "hg13_mapped_count",
        "union_mapped_ids",
        "union_mapped_count",
        "filtered_r1",
        "filtered_r2",
        "filtered_read_count",
        "bacdb_sorted_bam",
        "idxstats_tsv",
        "pileup_txt"
    ]
    for output in workflow_outputs:
        proj_config.add_column(name=output)

    for index, row in proj_config:
        if DRY_RUN == True: 
            print(f"\n{index}\t{row}") 
        #sample_id = row["sample"]
        fastq_r1 = READS_DIR / f"{row["r1_fastq"]}"
        read_count = fastq_count_reads(file_path = fastq_r1, dry_run = DRY_RUN)
        proj_config.update_row(index, "read_count", read_count)
        # h38_bam -> TEMP_DIR -> sort -> BAM_DIR -> index -> id list
        h38_sorted_bam, h38_mapped_ids = map_to_human(row, HG38, "hg38", THREADS)
        hg38_mapped_count = count_lines(h38_mapped_ids, dry_run = DRY_RUN)
        proj_config.update_row(index, "hg38_sorted_bam", str(h38_sorted_bam))
        proj_config.update_row(index, "hg38_mapped_ids", str(h38_mapped_ids))
        proj_config.update_row(index, "hg38_mapped_count", hg38_mapped_count)
        # hg13_bam -> TEMP_DIR -> sort -> BAM_DIR -> index -> id list
        h13_sorted_bam, h13_mapped_ids = map_to_human(row, HG13, "hg13", THREADS)
        hg13_mapped_count = count_lines(h13_mapped_ids, dry_run = DRY_RUN)
        proj_config.update_row(index, "hg13_sorted_bam", str(h13_sorted_bam))
        proj_config.update_row(index, "hg13_mapped_ids", str(h13_mapped_ids))
        proj_config.update_row(index, "hg13_mapped_count", hg13_mapped_count)
        # extract union of mapped ids and filter r1 & r2 fastq 
        mapped_ids, filtered_r1, filtered_r2 = filter_fastq(row, "hg38", "hg13")
        union_mapped_count = count_lines(mapped_ids, dry_run = DRY_RUN)
        ### get filtered count 
        proj_config.update_row(index, "union_mapped_ids", str(mapped_ids))
        proj_config.update_row(index, "union_mapped_count", union_mapped_count)
        proj_config.update_row(index, "filtered_r1", str(filtered_r1))
        proj_config.update_row(index, "filtered_r2", str(filtered_r2))
        filtered_read_count = fastq_count_reads(file_path = filtered_r1, dry_run = DRY_RUN)
        proj_config.update_row(index, "filtered_read_count", filtered_read_count)
        # map filtered fastq's to bacterial refs
        sorted_bam = map_filtered_fastq(row, BACDB, "R01Td", THREADS)
        proj_config.update_row(index, "bacdb_sorted_bam", str(sorted_bam))
        # get coverage for refs
        idxstats_tsv, pileup_txt = get_coverage_stats(row, BACFA, "R01Td")
        proj_config.update_row(index, "idxstats_tsv", str(idxstats_tsv))
        proj_config.update_row(index, "pileup_txt", str(pileup_txt))
        # save updated config after processing each row
        proj_config.save()


if __name__ == "__main__":
    main()