#! python

"""
Samuel Joseph Bryson
Copyright 2026
Build a neo4j graph of the ICTV viral taxonomy.

"""

from neomodel import get_config, db
import argparse
import csv
from pathlib import Path
from typing import Tuple

from sb_projects.ictv_graph.types import (
    IsolateType,
    GenomeStatus,
    GenomeType,
    HostSource
)

from sb_projects.ictv_graph.nodes import (
    ICTVtaxon,
    ICTVisolate,
    SequenceAccession
)

class ICTVBuilder:
    def __init__(self, ictv_file: str):
        self.ictv_file = ictv_file
        self.rank_list = [
            "_realm",
            "_subrealm",
            "_kingdom",
            "_subkingdom",
            "_phylum", 
            "_subphylum",
            "_class",
            "_subclass",
            "_order",
            "_suborder",
            "_family",
            "_subfamily",
            "_genus",
            "_subgenus",
            "_species",  
            ]
        self.node_cache = {
            "names":{rank: None for rank in self.rank_list},
            "nodes":{rank: None for rank in self.rank_list},
        }

          
    def build_ictv_graph(self) -> None:
        
        root_node = self._init_root_taxon()

        with open(self.ictv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                parent_node = root_node
                current_node = None
                try:
                    isolate_dict, taxonomy_dict = self._parse_ictv_entry(row)
                    accession_list  = self._parse_genbank_accs(row["GENBANK"], acc_type="genbank") 
                    with db.transaction:
                        # check and add lineage
                        for rank_type in self.rank_list:
                            rank_name = taxonomy_dict[rank_type]
                            if rank_name:   
                                rank_node = ICTVtaxon.nodes.get_or_none(name=rank_name, rank=rank_type)
                                if not rank_node:
                                    rank_node = ICTVtaxon(name=rank_name, rank=rank_type).save()
                                current_node = rank_node
                                if not current_node.parent.is_connected(parent_node):
                                        current_node.parent.connect(parent_node)
                                parent_node = current_node
                                self.node_cache["nodes"][rank_type] = rank_node
                                self.node_cache["names"][rank_type] = rank_name
                        # use the species node for the connecting isolates and accessions
                        species_node = self.node_cache["nodes"]["_species"]
                        # add the ICTVisolate node and connect it to the ICTVtaxon species node
                        isolate_node = ICTVisolate.nodes.get_or_none(
                            isolate_id = isolate_dict["isolate_id"]
                        )
                        if not isolate_node:
                            isolate_node = ICTVisolate(**isolate_dict).save()
                        if not isolate_node.assignment.is_connected(species_node):
                                    isolate_node.assignment.connect(species_node)
                        # create each SequenceAccession node and connect the isolate and species
                        if len(accession_list) > 0:
                            for acc in accession_list:
                                acc_node = SequenceAccession.nodes.get_or_none(**acc)
                                if not acc_node:
                                    acc_node = SequenceAccession(**acc).save()
                                if not isolate_node.assembly.is_connected(acc_node):
                                        isolate_node.assembly.connect(acc_node)
                                if not acc_node.assignment.is_connected(species_node):
                                        acc_node.assignment.connect(species_node)
                except Exception as e:
                    print(f"Error processing ICTV record: {row}\n {e}")
                    continue

    ## helper functions ##
    def _init_root_taxon(self) -> ICTVtaxon:
        try:
            # add root node
            with db.transaction:
                root_node = ICTVtaxon.nodes.get_or_none(name="root", rank="root_")
                if not root_node:
                    root_node = ICTVtaxon(name="root", rank="root_").save()
            self.node_cache["nodes"]["_root"] = root_node
            self.node_cache["names"]["_root"] = "root"
        except Exception as e:
            print(f"Error initializing ICTV taxonomy root node: {e}")
        else:
            return root_node

    def _parse_ictv_entry(self, row: dict) -> Tuple[dict,dict]:

        isolate_dict = {
            "isolate_id"   : int(row["Isolate_ID"]),
            "name"         : row["Name"],
            "abbreviation" : row["Abbreviation"],
            "isolate_type" : self._parse_isolate_type(row["Exemplar"]),
            "genome_type"  : self._parse_genome_type(row["Genome_type"]),
            "genome_status": self._parse_genome_status(row["Genome_status"]),
            "host_source"  : self._parse_host_source(row["Host_source"]),
        }
        taxonomy_dict = {
            "_realm"     : row["Realm"],
            "_subrealm"  : row["Subrealm"],
            "_kingdom"   : row["Kingdom"],
            "_subkingdom": row["Subkingdom"],
            "_phylum"    : row["Phylum"],
            "_subphylum" : row["Subphylum"],
            "_class"     : row["Class"],
            "_subclass"  : row["Subclass"],
            "_order"     : row["Order"],
            "_suborder"  : row["Suborder"],
            "_family"    : row["Family"],
            "_subfamily" : row["Subfamily"],
            "_genus"     : row["Genus"],
            "_subgenus"  : row["Subgenus"],
            "_species"   : row["Species"],
        }
        #ictv_id = int(row["ICTV_ID"]) # species taxonomy ID
        return isolate_dict, taxonomy_dict


    def _parse_isolate_type(self, record: str) -> str:
        try:
            return IsolateType(record).name
        except Exception:
            return IsolateType("U").name

    def _parse_genbank_accs(self, record: str, acc_type: str) -> list[dict]:
        
        if not record or record.upper() == 'NULL':
            return []
        
        results = []
        parts = [p.strip() for p in record.split(';')]
        for part in parts:
            label = "Genome"
            start = None
            stop = None
            if ':' in part:
                label, acc = part.split(':', 1)
                label = label.strip()
                acc = acc.strip()
            else:
                acc = part.strip()
            if '(' in part:
                label = "Provirus"
                acc, coords = part.split(' ', 1)
                coords = coords.replace('(', '').replace(')', '').split('.')
                start = coords[0]
                stop = coords[1]               
            
            if start != None:
                results.append({
                    'accession_id': acc,
                    'type'        : acc_type,
                    'segment'     : label,
                    'start'       : start,
                    'stop'        : stop
                })
            else:
                results.append({
                    'accession_id': acc,
                    'type'        : acc_type,
                    'segment'     : label,
                })
            
        return results

    def _parse_genome_status(self, record: str) -> str:
        """convert to codes - CCG: Coding-complete genome,	CG	Complete genome, 
        NG	No entry in Genbank, PG	Partial genome (non-compliant)"""
        try:
            return GenomeStatus(record).name
        except Exception:
            return GenomeStatus("Unknown").name
        
    def _parse_genome_type(self, record: str) -> None:
        """e.g. dsDNA, ssDNA(+), etc. - Use GenomeType Enum"""
        try:
            return GenomeType(record).name
        except Exception:
            return GenomeType("Unassigned").name
        
    def _parse_host_source(self, record: str) -> None:
        try:
            return HostSource(record).name
        except Exception:
            return HostSource("unknown (S)").name

def main():
    
    parser = argparse.ArgumentParser(description=".")
    parser.add_argument("--uri", required=False, type=str, dest="uri", default="localhost:7687",
                        help="URL for bolt connection to neo4j database.")
    parser.add_argument("--user", required=False, type=str, dest="user", default="neo4j",
                         help="DB username.")
    parser.add_argument("--pwd", required=False, type=str, dest="pwd", default="virus_db",
                         help="DB password.")
    parser.add_argument("--ictv_file", required=True, type=str, dest="ictv_file",
                         help="Path of the ICTV taxonomy input file.")
    args = parser.parse_args()
    
    try:
        config = get_config()
        config.database_url = f"bolt://{args.user}:{args.pwd}@{args.uri}"
        input = Path(args.ictv_file)  
        builder = ICTVBuilder(input)
        builder.build_ictv_graph()

    except Exception as e:
        print(f"Exception: {e}")

    finally:
        # Closes the entire driver/pool for the script
        db.close_connection()
        print(f"Connection Pool Closed.")

if __name__ == "__main__":
    main()
    