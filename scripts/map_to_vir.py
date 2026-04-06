

"""
Samuel Joseph Bryson
Copyright 2026

- 1 - Run virus detection using minimap2 to align human depleted (used hg38 and chm13) 
    - reads against a set of virus databases: {virus_db}.mmi.
- 2 - 
- 3 - 
- 4 - Mapped read counts and coverage stats are calculated for each sample x viral genome.
"""

from pathlib import Path
from typing import Tuple
from src.config_manager import ConfigManager
from src.subprocess_utilities import run_check_call
from src.file_utilities import count_lines
from wrappers import (
    Minimap2MapSRtoSortedBam,
    SamtoolsIndex,
    GetProperMappedIDs,
    SamtoolsIdxstats,
    SamtoolsMpileup,
)

# Database paths (all precomputed with -k 15 -w 10)
DB_DIR    = Path.home() / "bio_db"
ESV       = DB_DIR      / "esviritu_v3.2.4/srVPD1510.mmi"
ESV_FA    = DB_DIR      / "esviritu_v3.2.4/virus_pathogen_database.fna"
HVP       = DB_DIR      / "hvphcphtv_v0.2.1/srHVPDB1510.mmi"
HVP_FA    = DB_DIR      / "hvphcphtv_v0.2.1/human_virus_seqs.v0.2.1.fna"
REFSEQ    = DB_DIR      / "NCBI_REFSEQ_VIRAL/srREFSEQ1510.mmi"
REFSEQ_FA = DB_DIR      / "NCBI_REFSEQ_VIRAL/viral.1.1.genomic.fna"
VHDB      = DB_DIR      / "VHDB/srVHDB1510.mmi"
VHDB_FA   = DB_DIR      / "VHDB/virushostdb.genomic.fna"
DB_PATHS  = [ESV, HVP, REFSEQ, VHDB]
DB_FASTAS = [ESV_FA, HVP_FA, REFSEQ_FA, VHDB_FA]
DB_NAMES  = ["EsViritu", "HVP", "RefSeq", "VHDB"]

# Data paths
DATA_DIR       = Path.home() / "NIH_HVP/HOGBV_AD/SAMPLES"
FILT_READS_DIR = DATA_DIR    / "filtered_reads"
BAM_DIR        = DATA_DIR    / "vir_bams"
VIR_DIR        = DATA_DIR    / "virus_detection"

# Samples config file
CONFIG         = VIR_DIR     / "vir_cfg.txt"
TEST_CONFIG    = VIR_DIR     / "test_vir_cfg.txt"

# global vars
DRY_RUN = False
THREADS = 20

# Map each sample to a viral genome reference database
def map_to_virus(
          record: dict, 
          db_path: Path, 
          db_name: str,
          db_fasta: Path,
          p: int
) -> Tuple[Path, Path, Path, Path]:
    SAMPLE       = record["sample"]
    sorted_bam   = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam"
    mapped_ids   = VIR_DIR / f"{SAMPLE}_{db_name}.mapped_ids.txt"
    idxstats_tsv = VIR_DIR / f"{SAMPLE}_{db_name}.sorted.idxstats.tsv"
    pileup_txt   = VIR_DIR / f"{SAMPLE}_{db_name}.sorted.pileup.txt"
    
    # map paired fastq to virus database
    step1 = Minimap2MapSRtoSortedBam(
        input_mmi  = db_path,
        r1         = FILT_READS_DIR / f"{record["r1_fastq"]}",
        r2         = FILT_READS_DIR / f"{record["r2_fastq"]}",
        output_bam = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        threads    = p,
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step1.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    
    # index bam
    step2 = SamtoolsIndex(
        sorted_bam = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step2.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # get the ids for mapped pairs
    step3 = GetProperMappedIDs(
        sorted_bam  = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        output_file = VIR_DIR / f"{SAMPLE}_{db_name}.mapped_ids.txt",
        threads     = p,
        dry_run     = DRY_RUN
    )
    run_check_call(
        formatted_command = step3.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )

    # run idxstats
    step4 = SamtoolsIdxstats(
    sorted_bam = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
    output_tsv = VIR_DIR / f"{SAMPLE}_{db_name}.idxstats.tsv",
    dry_run    = DRY_RUN
    )
    run_check_call(
        formatted_command = step4.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    # run mpileup
    step5 = SamtoolsMpileup(
        refernce_fasta = db_fasta,
        sorted_bam     = BAM_DIR / f"{SAMPLE}_{db_name}.sorted.bam",
        output_pileup  = VIR_DIR / f"{SAMPLE}_{db_name}.pileup.txt",
        dry_run        = DRY_RUN
    )
    run_check_call(
        formatted_command = step5.build(),
        devnull           = True,
        dry_run           = DRY_RUN
    )
    return sorted_bam, mapped_ids, idxstats_tsv, pileup_txt

def main():
    # load config
    #proj_config = ConfigManager(config_file=CONFIG)
    proj_config = ConfigManager(config_file=TEST_CONFIG)
    
    # run the workflow on each database
    for db_name, db_path, db_fasta in zip(DB_NAMES, DB_PATHS, DB_FASTAS):
        if DRY_RUN == True: 
                print(f"\n{db_name}\t{db_path}\t{db_fasta}")
        
        proj_config.add_column(name=f"{db_name}_sorted_bam")
        proj_config.add_column(name=f"{db_name}_mapped_ids")
        proj_config.add_column(name=f"{db_name}_mapped_count")
        proj_config.add_column(name=f"{db_name}_idxstats_tsv")
        proj_config.add_column(name=f"{db_name}_pileup_txt")
        
        # process each sample (config row) for the current database
        for index, row in proj_config:
            if DRY_RUN == True: 
                print(f"\n{index}\t{row}")
            # map reads and get coverage
            sorted_bam, mapped_ids, idxstats_tsv, pileup_txt  = map_to_virus(
                 row,
                 db_path, 
                 db_name,
                 db_fasta,
                 THREADS
            )
            mapped_count = count_lines(mapped_ids, dry_run = DRY_RUN)
            proj_config.update_row(index, f"{db_name}_sorted_bam", str(sorted_bam))
            proj_config.update_row(index, f"{db_name}_mapped_ids", str(mapped_ids))
            proj_config.update_row(index, f"{db_name}_mapped_count", mapped_count)
            proj_config.update_row(index, f"{db_name}_idxstats_tsv", str(idxstats_tsv))
            proj_config.update_row(index, f"{db_name}_pileup_txt", str(pileup_txt))
            # save updated config after processing each row
            proj_config.save()


if __name__ == "__main__":
    main()