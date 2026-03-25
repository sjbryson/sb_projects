#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026
"""

from dataclasses import dataclass, asdict, field, fields, InitVar
from pathlib import Path
from typing import Optional

@dataclass(kw_only=True)
class BioTool():
    dry_run: InitVar[bool] = False
    cmd:     str = field(default="", init=False, repr=False)

    def __post_init__(self, dry_run: bool):
        for f in fields(self):
            metadata_type = f.metadata.get('type', '')
            value = getattr(self, f.name)

            if value is None:
                continue

            # --- Value Flags (e.g., threads, quality scores) --- #
            if metadata_type == 'value_flag':
                flag_fmt = f.metadata.get('flag_fmt') # (e.g., "-t {value}")
                if value is not None and value is not False:
                    setattr(self, f.name, flag_fmt.format(value=value))
                else:
                    setattr(self, f.name, "")

            # --- Bool Flags --- #
            elif metadata_type == 'flag':
                flag_string = f.metadata.get('option')
                # If the attribute is True, use the flag; otherwise, empty string
                setattr(self, f.name, flag_string if value else "")

            # --- IO Paths --- #
            elif metadata_type in ['input_file', 'output_file'] and value:
                path_obj = Path(value)
                setattr(self, f.name, path_obj)
                #if metadata_type == 'input_file' and not path_obj.exists():
                #    raise FileNotFoundError(f"Missing: {path_obj}")
                if metadata_type == 'input_file' and not dry_run:
                    if not path_obj.exists():
                        raise FileNotFoundError(f"Missing: {path_obj}")
                if metadata_type == 'output_file':
                    path_obj.parent.mkdir(parents=True, exist_ok=True)

    def build(self) -> str:
        """Format self.cmd string with attribute values."""
        return " ".join(self.cmd.format(**asdict(self)).split())

####################
## -- Minimap2 -- ##
####################
@dataclass(kw_only=True)
class PreIndexMinimap2SR(BioTool):
    input_fasta: Path | str = field(metadata={'type': 'input_file'})
    output_mmi:  Path | str = field(metadata={'type': 'output_file'})
    threads:     Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:         str = "minimap2 -k 21 -w 11 {threads} -d {output_mmi} {input_fasta}"

@dataclass(kw_only=True)
class PreIndexMinimap2LR(BioTool):
    input_fasta: Path | str = field(metadata={'type': 'input_file'})
    output_mmi:  Path | str = field(metadata={'type': 'output_file'})
    threads:     Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:         str = "minimap2 -H {threads} -d {output_mmi} {input_fasta}"

@dataclass(kw_only=True)
class Minimap2MapSRtoSortedBam(BioTool):
    input_mmi:  Path | str = field(metadata={'type': 'input_file'})
    r1:         Path | str = field(metadata={'type': 'input_file'})
    r2:         Path | str = field(metadata={'type': 'input_file'})
    output_bam: Path | str = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:        str = "minimap2 -ax sr {threads} {input_mmi} {r1} {r2} | samtools sort -o - > {output_bam}"

@dataclass(kw_only=True)
class Minimap2SRHumanDepletion(BioTool):
    input_mmi:  Path | str = field(metadata={'type': 'input_file'})
    r1:         Path | str = field(metadata={'type': 'input_file'})
    r2:         Path | str = field(metadata={'type': 'input_file'})
    output_bam: Path | str = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:        str = "minimap2 -ax sr --secondary=no {threads} {input_mmi} {r1} {r2} | samtools sort -o - > {output_bam}"

@dataclass(kw_only=True)
class Minimap2MapLR5(BioTool):
    input_mmi:  Path | str = field(metadata={'type': 'input_file'})
    reads:      Path | str = field(metadata={'type': 'input_file'})
    output_bam: Path | str = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:        str = "minimap2 -ax asm5 {threads} {input_mmi} {reads} | samtools sort -o {output_bam}"

@dataclass(kw_only=True)
class Minimap2MapLR20(BioTool):
    input_mmi:  Path | str = field(metadata={'type': 'input_file'})
    reads:      Path | str = field(metadata={'type': 'input_file'})
    output_bam: Path | str = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:        str = "minimap2 -ax asm20 {threads} {input_mmi} {reads} | samtools sort -o {output_bam}"

## -- host read depletion -- ##
@dataclass(kw_only=True)
class MinimapFilterPairedFastq(BioTool):
    """Standard Host Depletion: write unmapped pairs to paired fastq."""
    input_mmi:   Path | str = field(metadata={'type': 'input_file'})
    r1:          Path | str = field(metadata={'type': 'input_file'})
    r2:          Path | str = field(metadata={'type': 'input_file'})
    filtered_r1: Path | str = field(metadata={'type': 'output_file'})
    filtered_r2: Path | str = field(metadata={'type': 'output_file'})
    mapq:        Optional[int] = field(default=5)
    threads:     Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:         str = "minimap2 -ax sr --secondary=no {threads} {input_mmi} {r1} {r2} | \
                        samtools sort -n - | \
                        filter_bam_to_paired_fastq {threads} --r1 {filtered_r1} --r2 {filtered_r2} --max-mapq {mapq}"
    

## -- virus mapping --##
@dataclass(kw_only=True)
class MinimapAlignToSortedBam(BioTool):
    """Standard virus detection: write to position sorted bam."""
    input_mmi:  Path | str = field(metadata={'type': 'input_file'})
    r1:         Path | str = field(metadata={'type': 'input_file'})
    r2:         Path | str = field(metadata={'type': 'input_file'})
    output_bam: Path | str = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-t {value}'})
    cmd:        str = "minimap2 -ax sr -N 5 {threads} {input_mmi} {r1} {r2} | \
                       samtools sort -o {output_bam} -"

## -- deduplicate bam --##
@dataclass(kw_only=True)
class DeDuplicateToSortedBam(BioTool):
    """Standard virus detection: mark and remove PCR duplicates."""
    input_bam:  Path | str = field(metadata={'type': 'input_file'})
    output_bam: Path | str = field(metadata={'type': 'output_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-@ {value}'})
    cmd:        str = "samtools sort -n {threads} {input_bam} | \
                       samtools fixmate -m - - | \
                       samtools sort {threads} - | \
                       samtools markdup -r - {output_bam}"


####################
## -- Samtools -- ##
####################
@dataclass(kw_only=True)
class SamtoolsSort(BioTool):
    sorted_bam: Path | str = field(metadata={'type': 'output_file'})
    input_bam:  Path | str = field(metadata={'type': 'input_file'})
    threads:    Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-@ {value}'})
    cmd:        str = "samtools sort {threads} -o {sorted_bam} {input_bam}"

@dataclass(kw_only=True)
class SamtoolsIndex(BioTool):
    sorted_bam: Path | str = field(metadata={'type': 'output_file'})
    cmd:        str = "samtools index {sorted_bam}"

@dataclass(kw_only=True)
class SamtoolsMpileup(BioTool):
    refernce_fasta: Path | str = field(metadata={'type': 'input_file'})
    sorted_bam:     Path | str = field(metadata={'type': 'input_file'})
    output_pileup:  Path | str = field(metadata={'type': 'output_file'})
    cmd:            str = "samtools mpileup -f {refernce_fasta} {sorted_bam} > {output_pileup}"

@dataclass(kw_only=True)
class SamtoolsStats(BioTool):
    sorted_bam:  Path | str = field(metadata={'type': 'input_file'})
    output_file: Path | str = field(metadata={'type': 'output_file'})
    cmd:         str = "samtools stats {sorted_bam} > {output_file}"

@dataclass(kw_only=True)
class SamtoolsFlagstat(BioTool):
    sorted_bam:  Path | str = field(metadata={'type': 'input_file'})
    output_file: Path | str = field(metadata={'type': 'output_file'})
    cmd:         str = "samtools flagstat {sorted_bam} > {output_file}"

@dataclass(kw_only=True)
class SamtoolsIdxstats(BioTool):
    sorted_bam: Path | str = field(metadata={'type': 'input_file'})
    output_tsv: Path | str = field(metadata={'type': 'output_file'})
    cmd:        str = "samtools idxstats {sorted_bam} > {output_tsv}"

@dataclass(kw_only=True)
class GetProperMappedIDs(BioTool):
    sorted_bam:  Path | str = field(metadata={'type': 'input_file'})
    output_file: Path | str = field(metadata={'type': 'output_file'})
    threads:     int = 1
    cmd:         str = "samtools view -f 1 -F 12 -@ {threads} {sorted_bam} | cut -f1 | sort -u > {output_file}"

@dataclass(kw_only=True)
class GetOrphanMappedIDs(BioTool):
    sorted_bam:  Path | str = field(metadata={'type': 'input_file'})
    output_file: Path | str = field(metadata={'type': 'output_file'})
    threads:     Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-@ {value}'})
    cmd:         str = "samtools view -f 4 -F 8 {threads} {sorted_bam} | cut -f1 | sort -u > {output_file}"

#######################
## -- Blast Tools -- ##
#######################
@dataclass(kw_only=True)
class MakeBlastnDB(BioTool):
    input_fasta: Path | str = field(metadata={'type': 'input_file'})
    db_path:     Path | str = field(metadata={'type': 'output_file'})
    db_name:     str
    cmd:         str = "makeblastdb -in {input_fasta} -input_type fasta -dbtype nucl -out {db_path} -parse_seqids -title {db_name}"

@dataclass(kw_only=True)
class Blastn(BioTool):
    database:    str
    input_fasta: Path | str = field(metadata={'type': 'input_file'})
    output_tsv:  Path | str = field(metadata={'type': 'output_file'})
    threads:     Optional[int] = field(default=None, metadata={'type': 'value_flag', 'flag_fmt': '-num_threads {value}'})
    cmd:         str = "blastn -db {database} -query {input_fasta} -out {output_tsv} -outfmt 6 {threads}"


##########################
## -- SHELL Commands -- ##
##########################
@dataclass(kw_only=True)
class GetContigPerBaseCoverageFromPileup(BioTool):
    contig_id:    str
    pileup_file:  Path | str = field(metadata={'type': 'input_file'})
    output_tsv:   Path | str = field(metadata={'type': 'output_file'})
    cmd:          str = "awk '$1 == \"{reference_id}\" {{print $2 \"\\t\" $4}}' {pileup_file} > {output_tsv}"

@dataclass(kw_only=True)
class UniqueIdList(BioTool):
    file_a:     Path = field(metadata={'type': 'input_file'})
    file_b:     Path = field(metadata={'type': 'input_file'})
    output_txt: Path = field(metadata={'type': 'output_file'})
    cmd:        str  = "cat {file_a} {file_b} | sort -u > {output_txt}"

@dataclass(kw_only=True)
class GzipKeepFile(BioTool):
    """Compress a file using gzip, output to a new path."""
    input_file:  Path | str = field(metadata={'type': 'input_file'})
    output_file: Path | str = field(metadata={'type': 'output_file'})
    level:       Optional[int] = field(default=6, metadata={'type': 'value_flag', 'flag_fmt': '-{value}'})
    force:       bool = field(default=False, metadata={'type': 'flag', 'option': '-f'})
    cmd:         str = "gzip -c {level} {force} {input_file} > {output_file}"

@dataclass(kw_only=True)
class GzipReplaceFile(BioTool):
    """Compress a file in-place using gzip - original is deleted."""
    input_file:  Path | str = field(metadata={'type': 'input_file'})
    level:       Optional[int] = field(default=6, metadata={'type': 'value_flag', 'flag_fmt': '-{value}'})
    force:       bool = field(default=False, metadata={'type': 'flag', 'option': '-f'})
    cmd:         str = "gzip -c {level} {force} {input_file} > {output_file}"

####################
## -- SB_Tools -- ##
####################
@dataclass(kw_only=True)
class SBFilterFastq(BioTool):
    r1_fastq:       Path | str = field(metadata={'type': 'input_file'})
    r2_fastq:       Path | str = field(metadata={'type': 'input_file'})
    input_file:     Path | str = field(metadata={'type': 'input_file'})
    output_prefix:  str
    gz:             bool = field(default=False, metadata={'type': 'flag', 'option': '--gz'})
    keep:           bool = field(default=False, metadata={'type': 'flag', 'option': '--keep'})
    exclude:        bool = field(default=True, metadata={'type': 'flag', 'option': '--exclude'})
    cmd:            str = "sb_filter_fastq {gz} {invert} --r1 {r1_fastq} --r2 {r2_fastq} --filter {input_file} --out-prefix {output_prefix}"