#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026

- custom types and validation logic
"""

from enum import Enum

from neomodel import (
    StringProperty,
)

class StringChoices(StringProperty):
    """Subclass of neomodel StringProperty with defined options."""
    def __init__(self, choices=None, **kwargs):
        self.choices = choices
        super().__init__(**kwargs)

    def validate(self, value):
        super().validate(value)
        if value not in self.choices:
            raise ValueError(f"Value {value} is not a valid choice. Expected one of {self.choices}.")


class TaxonRank(Enum):
    NULL       = "root"         # 100
    REALM      = "realm"        # 120
    SUBREALM   = "subrealm"     # 130
    KINGDOM    = "kingdom"      # 140
    SUBKINGDOM = "subkingdom"   # 150
    PHYLUM     = "phylum"       # 160
    SUBPHYLUM  = "subphylum"    # 170
    CLASS      = "class"        # 180
    SUBCLASS   = "subclass"     # 190
    ORDER      = "order"        # 200
    SUBORDER   = "suborder"     # 250
    FAMILY     = "family"       # 300
    SUBFAMILY  = "subfamily"    # 400
    GENUS      = "genus"        # 500
    SUBGENUS   = "subgenus"     # 550
    SPECIES    = "species"      # 600

class VirusType(Enum):
    ICTV_ISOLATE = "ictv_isolate"
    ASSEMBLY     = "assembly"

class IsolateType(Enum):
    ACCESSORY = "A"
    EXEMPLAR  = "E"
    UNKNOWN   = "U"
    
class HostSource(Enum):
    AIR     = "air (S)"
    ALG     = "algae"
    ARC     = "archaea"
    BAC     = "bacteria"
    BAC_INV = "bacteria, invertebrates (S)"
    BAC_VRT = "bacteria, vertebrates (S)"
    FRW     = "freshwater (S)"
    FUN     = "fungi"
    FUNS    = "fungi (S)"
    INV     = "invertebrates"
    INVS    = "invertebrates (S)"
    INV_PLT = "invertebrates, plants"
    INV_VRT = "invertebrates, vertebrates"
    MARs    = "marine (S)"
    OOM     = "oomycetes"
    OTH     = "other (specify)"
    PHY     = "phytobiome (S)"
    PLT     = "plants"
    PLTS    = "plants (S)"
    PLT_FUN = "plants, fungi"
    PRT     = "protists"
    PRTS    = "protists (S)"
    SEW     = "sewage (S)"
    SOIL    = "soil (S)"
    UNK     = "unknown (S)"
    VRT     = "vertebrates"
    VRTS    = "vertebrates (S)"

class GenomeType(Enum):
    UN       = "Unassigned"           # (0)
    DSDNA    = "dsDNA"                # (1)  - Double-stranded DNA	1	I
    SSDNA    = "ssDNA"                # (2)  - Single-stranded DNA	2	II
    DSRNA    = "dsRNA"                # (3)  - Double-stranded RNA	3	III
    SSRNA_P  = "ssRNA(+)"             # (4)  - Single-stranded RNA - Positive-sense	4	IV
    SSRNA_N  = "ssRNA(-)"             # (5)  - Single-stranded RNA - Negative-sense	5	V
    SSRNA_RT = "ssRNA-RT"             # (6)  - Single-stranded RNA - Positive-sense - replicates through a DNA intermediate	6	VI
    DSRNA_RT = "dsRNA-RT"             # (7)  - Double-stranded DNA - replicates through a single-strand RNA intermediate	7	VII
    VIROID   = "Viroid"               # (8)  - Circular Single-stranded RNA
    SSDNA_N  = "ssDNA(-)"             # (9)  - Single-stranded DNA - Negative-sense
    SSDNA_P  = "ssDNA(+)"             # (10) - Single-stranded DNA - Positive-sense
    SSDNA_A  = "ssDNA(+/-)"           # (11) - Single-stranded DNA - Ambi-sense
    SSRNA_A  = "ssRNA(+/-)"           # (12) - Single-stranded RNA - Ambi-sense
    DSSSDNA  = "dsDNA; ssDNA"         # (13) - DNA - some taxa double stranded, some taxa single standed
    SSRNA    = "ssRNA"                # (14) - Single-stranded RNA
    SSRNA_NA = "ssRNA(-); ssRNA(+/-)" # (15) - Single-stranded RNA - some Negative-sense segments, some Ambi-sense segments

class GenomeStatus(Enum):
    NG = "No entry in Genbank"             # 100
    PG = "Partial genome (non-compliant)"  # 200
    CC = "Coding-complete genome"          # 300 - aka "CCG"
    CG = "Complete genome"                 # 400
    UK = "Unknown"

class AccessionType(Enum):
    REFSEQ  = "RefSeq"
    GENBANK = "GenBank"