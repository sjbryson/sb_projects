
"""
Samuel Joseph Bryson
Copyright 2026

- edge definitions defined with neomodel
"""

from neomodel import StructuredRel, StringProperty, DateTimeProperty

class HAS_PARENT(StructuredRel):
    """Relationship from child taxonomy node to parent taxonomy node"""
    pass

class ICTV_ASSIGNMENT(StructuredRel):
    """Relationship between an isolate or virus genome to species or other ICTV taxonomy node"""
    pass 

class GENOME_ASSIGNMENT(StructuredRel):
    """Relationship between a RefSeq or GenBank accession to species or other ICTV taxonomy node"""
    pass 


class HAS_ASSEMBLY(StructuredRel):
    """Relationship from an an isolate or viral genome to a RefSeq or GenBank accession.
    Can be one-to-many such as with segmented genomes."""
    pass