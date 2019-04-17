import sys
import os
import urllib.request
import collections
import re
from subprocess import Popen, PIPE

import KeggNames
kn = KeggNames.keggnames

# Use my NCBI API key to set server/second server access <= 10
os.environ["NCBI_API_KEY"] = "f4142f101db95406745385d940b13c37ab07"

# Make data struct to store ori info
Ori = collections.namedtuple('Ori', ['seq','len','start','end'])

def get_kegg_organism_list():
  url = 'http://rest.kegg.jp/list/organism'
  return urllib.request.urlopen(url).read().decode('utf-8').split('\n')

def search_ncbi_for_refseq_id(gb_id):
  cmd = f"esearch -db nuccore -query {gb_id} | elink -target nuccore -name nuccore_nuccore_gbrs | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Slen Title | sed 's/,.*//' | sort -t $'\t' -k 2,2nr"
  result,err = Popen(cmd, stdout=PIPE,stderr=PIPE,shell=True).communicate()
  result = result.decode('utf-8').strip().split('\n')
  if len(result) != 1:
    log = open('multiple_refseq_ids.log','a+')
    log.write('\n'.join(result))
    log.close()
    return result[0].split('\t')[0]
  return result[0].split('\t')[0]

class Organism:
  """Container class storing data scraped from KEGG, NCBI and TUBIC relating
     to a single organism present in KEGG organisms list."""

  def __init__(self, kegg_genome_id, kegg_organism_id = None, strain_name = None, taxa_string = None):
    self.kegg_genome_id   = kegg_genome_id   # e.g T01001
    self.kegg_organism_id = kegg_organism_id # e.g hsa
    self.strain_name      = strain_name      # e.g Homo sapiens (human)
    self.taxa_string      = taxa_string      # e.g Eukaryotes;Animals;Vertebrates;Mammals
    self.genes = dict()                      # k='gene ids', v='kegg info, gene location info, etc'
    self.gff 


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
    log = open('failed_metadata_retrieval.log','a+')
    # Resolve basic metadata
    try:
      self.genome_length = int(match(kn.LENGTH, entry))
      if self.kegg_organism_id == None:
        self.kegg_organism_id = match(kn.NAME, entry).split(',')[0]
      if self.strain_name == None:
        self.strain_name = match(kn.DEFINITION, entry)
      if self.taxa_string == None:
        self.taxa_string = match(kn.LINEAGE, entry)

      # Attempt to determine the Refseq sequence id used for this kegg organism
      seq = match(kn.SEQUENCE, entry)
      if seq[0:2] == 'RS':
        self.rs_id, self.gb_id = re.findall(r'RS:(\S+)\s+\(GB:(\S+)\)',seq)[0]
      if seq[0:2] == 'GB':
        self.gb_id = re.findall(r'GB:(\S+)',seq)[0]
        self.rs_id = search_ncbi_for_refseq_id(self.gb_id)
    except:
      log.write(entry)
    
  def get_kegg_data_for_genes(self):
    """Makes request to rest.kegg.jp/list/kegg_organim_id to get a list of geneids
       then for each gene_id makes request to rest.kegg.jp/get/gene_id to get
       the kegg info each gene. This is then parsed and the (k=gene_id, v=kegg_data)
       pair is stored in the self.gene dictionary."""
    url = f'http://rest.kegg.jp/link/pathway/{self.kegg_genome_id}'
    gene_pw_list = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    gene_pw_list = [tuple(e.split('\t') for e in gene_pw_list)]
    self.genes = {gid, self.get_gene_data(gid,pw) for gid, pw in gene_pw_list}

  def get_gene_data(self, gid, pw):
    url = f'http://rest.kegg.jp/link/pathway/{self.kegg_genome_id}'

  def set_ori_info(self, tubic_entry):
    """Parses the tubic_entry and extracts relevant info to make self.Ori object."""
    pass


def main():
  olist = get_kegg_organism_list()
  # Generator produces an Organism for each organism in kegg in relevant taxanomic rank
  orgs = (Organism(*org.split('\t')) 
            for org in olist
            if ('Bacteria' in Organism(*org.split('\t')).taxa_string))


  for org in orgs:
    org.get_metadata()


if __name__ == '__main__':
  main()

