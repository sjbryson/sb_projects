#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026

- node definitions defined with neomodel
"""

from neomodel import (
    StructuredNode, 
    StringProperty, 
    IntegerProperty, 
    RelationshipTo, 
    RelationshipFrom, 
    ZeroOrOne, 
    OneOrMore
)
from src.ictv_graph.relationships import (
    HAS_PARENT, 
    ICTV_ASSIGNMENT, 
    HAS_ASSEMBLY,
    GENOME_ASSIGNMENT,
)
#from ictv_graph.types import(
#    TaxonRank,
#    VirusType
#)

class ICTVtaxon(StructuredNode):
    # properties
    ##  taxnode_id = StringProperty(unique_index=True)
    name   = StringProperty(required=True)
    rank   = StringProperty(required=True)
    # relationships
    parent = RelationshipTo('ICTVtaxon', 'HAS_PARENT', model=HAS_PARENT)
    #child = RelationshipFrom('ICTVtaxon', 'HAS_PARENT', model=HAS_PARENT)
    

class ICTVisolate(StructuredNode):
    # properties
    isolate_id    = IntegerProperty(unique_index=True, required=True)
    name          = StringProperty()     
    abbreviation  = StringProperty()
    isolate_type  = StringProperty()
    genome_type   = StringProperty()
    genome_status = StringProperty()
    host_source   = StringProperty()
    # relationships
    assignment    = RelationshipTo('ICTVtaxon', 'ICTV_ASSIGNMENT', model=ICTV_ASSIGNMENT)
    assembly      = RelationshipTo('SequenceAccession', 'HAS_ASSEMBLY', model=HAS_ASSEMBLY)
    

class SequenceAccession(StructuredNode):
    accession_id = StringProperty(unique_index=True, required=True)
    type         = StringProperty(required=True) # 'genbank' or 'refseq'
    segment      = StringProperty()
    start        = IntegerProperty(required=False)
    stop         = IntegerProperty(required=False)
    # relationships 
    assignment   = RelationshipTo('ICTVtaxon', 'GENOME_ASSIGNMENT', model=GENOME_ASSIGNMENT)
