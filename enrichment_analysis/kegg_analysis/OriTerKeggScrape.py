import sys
import os
import urllib.request
import collections
import re

import KeggNames.keggnames as kn

# Make data struct to store ori info

Ori = collections.namedtuple('Ori', ['seq','len','start','end'])



class Organism:
  """Container class storing data scraped from KEGG, NCBI and TUBIC relating
     to a single organism present in KEGG organisms list."""

  def __init__(self, kegg_genome_id, kegg_organism_id = None, strain_name = None, taxa_string = None):
#      T05289	srub	Salinigranum rubrum	Prokaryotes;Archaea;Euryarchaeota;Salinigranum
    self.kegg_genome_id   = kegg_genome_id   # e.g T01001
    self.kegg_organism_id = kegg_organism_id # e.g hsa
    self.strain_name      = strain_name      # e.g Homo sapiens (human)
    self.taxa_string      = taxa_string      # e.g Eukaryotes;Animals;Vertebrates;Mammals
    self.genes = dict()                      # k='gene ids', v='kegg info, gene location info, etc'


  def in_taxonomic_group(self, group):
    """Returns True of group in self.taxa_string, False otherwise"""
    return group in self.taxa_string

  def get_genome_entry_from_kegg(self):
    url = f'http://rest.kegg.jp/get/gn:{self.kegg_genome_id}'
    # Download webpage from kegg, read, and convert to regular str type
    return urllib.request.urlopen(url).read().decode('utf-8')


  def get_metadata(self):
    """Makes request to rest.kegg.jp/get/kegg_genome_id to get various metadata:
       genome length; the sequence id of the genome used for kegg's annotations 
       (either refseq or genbank, if genbank, ncbi is searched to get refseq_id
        if available); and ncbi taxonomy"""
    entry = self.get_genome_entry_from_kegg()
    # Make regex engine
    match = lambda a, b : re.findall(rf'{a}\s+(\S.*)\n',b)[0]

    # Resolve basic metadata
    self.genome_length = int(match(kn.LENGTH, entry))
    if self.kegg_organism_id == None:
      self.kegg_organism_id = match(kn.NAME, entry).split(',')[0]
    if self.strain_name == None:
      self.strain_name = match(kn.DEFINITION, entry)
    if self.taxa_string == None:
      eslf.taxa_string = match(kn.LINEAGE, entry)

    # Attempt to determine the Refseq sequence id used for this kegg organism
    

    
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
