


'''
 % cat viral.1.1.genomic.metadata.json | jq -r '.reports[30] | path(..) | map(tostring) | join(".")'

accession
bioprojects
bioprojects.0
completeness
host
host.lineage
host.lineage.0
host.lineage.0.name
host.lineage.0.tax_id
host.lineage.1
host.lineage.1.name
host.lineage.1.tax_id
host.lineage.2
host.lineage.2.name
host.lineage.2.tax_id
    - last lineage entry == {host.organism_name: host.tax_id}
...
host.organism_name
host.tax_id
is_annotated
isolate
isolate.collection_date
isolate.name
length
location
location.geographic_location
location.geographic_region
nucleotide
nucleotide.sequence_hash
protein_count
release_date
segment
source_database
submitter
submitter.affiliation
submitter.country
submitter.names
submitter.names.0
submitter.names.1
submitter.names.2
submitter.names.3
update_date
virus
virus.lineage
virus.lineage.0
virus.lineage.0.name
virus.lineage.0.tax_id
virus.lineage.1
virus.lineage.1.name
virus.lineage.1.tax_id
virus.lineage.2
virus.lineage.2.name
virus.lineage.2.tax_id
    - last lineage entry == {virus.organism_name: virus.tax_id}
...
virus.organism_name
virus.tax_id
'''