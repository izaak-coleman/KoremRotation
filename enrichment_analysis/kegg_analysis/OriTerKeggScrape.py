import sys
import os
import urllib.request
import collections

# Make data struct to store ori info

Ori = collections.namedtuple('Ori', ['seq','len','start','end'])

class Organism:
  """Container class storing data scraped from KEGG, NCBI and TUBIC relating
     to a single organism present in KEGG organisms list."""

  def __init__(self, kegg_genome_id, kegg_organism_id, strain_name, taxa_string):
    self.kegg_genome_id   = kegg_genome_id   # e.g T01001
    self.kegg_organism_id = kegg_organism_id # e.g hsa
    self.strain_name      = strain_name      # e.g Homo sapiens (human)
    self.taxa_string      = taxa_string      # e.g Eukaryotes;Animals;Vertebrates;Mammals
    self.genes = dict()                      # k='gene ids', v='kegg info, gene location info, etc'


  def in_taxonomic_group(group):
    """Returns True of group in self.taxa_string, False otherwise"""
    pass

  def get_metadata():
    """Makes request to rest.kegg.jp/get/kegg_genome_id to get various metadata:
       genome length; the sequence id of the genome used for kegg's annotations 
       (either refseq or genbank, if genbank, ncbi is searched to get refseq_id
        if available); and ncbi taxonomy"""
    # store both kegg_seq_id and refseq_id
    pass

  def get_kegg_data_for_genes():
    """Makes request to rest.kegg.jp/list/kegg_organim_id to get a list of geneids
       then for each gene_id makes request to rest.kegg.jp/get/gene_id to get
       the kegg info each gene. This is then parsed and the (k=gene_id, v=kegg_data)
       pair is stored in the self.gene dictionary."""
    pass

  def set_ori_info(tubic_entry):
    """Parses the tubic_entry and extracts relevant info to make self.Ori object."""
    pass

#      T05289	srub	Salinigranum rubrum	Prokaryotes;Archaea;Euryarchaeota;Salinigranum
