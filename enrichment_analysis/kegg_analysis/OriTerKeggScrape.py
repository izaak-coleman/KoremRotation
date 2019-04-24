import sys
import os
import urllib.request
import collections
import re
from subprocess import Popen, PIPE
import gzip


import KeggNames
kn = KeggNames.keggnames

# Use my NCBI API key to set server/second server access <= 10
os.environ["NCBI_API_KEY"] = "f4142f101db95406745385d940b13c37ab07"

# Make data struct to store ori info
Ori = collections.namedtuple('Ori', ['seq','len','start','end'])

def get_kegg_organism_list():
  url = 'http://rest.kegg.jp/list/organism'
  data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
  data.pop()
  return data



class Organism:
  """Container class storing data scraped from KEGG, NCBI and TUBIC relating
     to a single organism present in KEGG organisms list."""

  def __init__(self, kegg_genome_id, kegg_organism_id = None, strain_name = None, taxa_string = None):
    self.kegg_genome_id   = kegg_genome_id   # e.g T01001
    self.kegg_organism_id = kegg_organism_id # e.g hsa
    self.strain_name      = strain_name      # e.g Homo sapiens (human)
    self.taxa_string      = taxa_string      # e.g Eukaryotes;Animals;Vertebrates;Mammals
    self.annotated_gff = []
    self.refseq = ''
    self.genbank = ''
    self.get_metadata()
    self.doric = ''


  def in_taxonomic_group(self, group):
    """Returns True of group in self.taxa_string, False otherwise"""
    return group in self.taxa_string

  def get_genome_entry_from_kegg(self):
    url = f'http://rest.kegg.jp/get/gn:{self.kegg_genome_id}'
    # Download webpage from kegg, read, and convert to regular str type
    return urllib.request.urlopen(url).read().decode('utf-8')

  def search_ncbi_for_refseq_id(self, gb_id):
    cmd = f"esearch -db nuccore -query {gb_id} | elink -target nuccore -name nuccore_nuccore_gbrs | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Slen Title | sed 's/,.*//' | sort -t $'\t' -k 2,2nr"
    result,err = Popen(cmd, stdout=PIPE,stderr=PIPE,shell=True).communicate()
    result = result.decode('utf-8').strip().split('\n')
    if len(result) != 1:
      log = open(self.kegg_organism_id + '_multiple_refseq_ids.log','w')
      log.write('\n'.join(result))
      log.close()
      return result[0].split('\t')[0]
    return result[0].split('\t')[0]


  def get_metadata(self):
    """Makes request to rest.kegg.jp/get/kegg_genome_id to get various metadata:
       genome length; the sequence id of the genome used for kegg's annotations 
       (either refseq or genbank, if genbank, ncbi is searched to get refseq_id
        if available); and ncbi taxonomy"""
    entry = self.get_genome_entry_from_kegg()
    # Make regex engine
    match = lambda a, b : re.findall(rf'{a}\s+(\S.*)\n',b)
    # Resolve basic metadata
    self.genome_length = match(kn.LENGTH, entry)
    if len(self.genome_length) > 0:
      self.genome_length = int(self.genome_length[0])
    if self.kegg_organism_id == None:
      if len(self.organism_id) > 0:
        self.organism_id = self.organism_id[0]
      self.kegg_organism_id = match(kn.NAME, entry).split(',')[0]

    if self.strain_name == None:
      self.strain_name = match(kn.DEFINITION, entry)
      if len(self.strain_name) > 0:
        self.strain_name = self.strain_name[0]
    if self.taxa_string == None:
      self.taxa_string = match(kn.LINEAGE, entry)
      if len(self.taxa_string) > 0:
        self.taxa_string= self.taxa_string[0]

    self.assembly = match(kn.DATA_SOURCE, entry)
    if len(self.assembly) > 0:
      self.assembly = self.assembly[0]
      self.assembly = re.findall('Assembly:([\.\w]+)', self.assembly)[0]
    # Attempt to determine the Refseq sequence id used for this kegg organism
    seq = match(kn.SEQUENCE, entry)
    if len(seq) > 0:
      seq = seq[0]
    if seq[0:2] == 'RS':
      self.refseq, self.genbank = re.findall(r'RS:(\S+)\s+\(GB:(\S+)\)',seq)[0]
    if seq[0:2] == 'GB':
      self.genbank = re.findall(r'GB:(\S+)',seq)[0]
      self.refseq = self.search_ncbi_for_refseq_id(self.genbank)
    
  def get_kegg_gene_annotations(self, annotation):
    url = f'http://rest.kegg.jp/link/{annotation}/{self.kegg_genome_id}'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    data = [tuple(e.split('\t')) for e in data]
    data.pop()
    return data


def main():
  pass

if __name__ == '__main__':
  main()


# to do

#  def generate_index(self, gff, field):
#    match = lambda a,b : re.findall(rf'[^\w]{a}=([\w]+)',b)
#    d = dict()
#    for index, line in enumerate(gff):
#      tag = match(field, line)
#      if len(tag) == 1:
#        d[tag[0]] = index
#    return d
#
#  def set_gff_filename(self, file_name):
#    self.gff_filename = file_name
#
#  def kegg_annotate_gff(self):
#    """Adds KEGG pathway, Ontology and other annotations to gff"""
#    if len(self.annotated_gff) == 0:
#      with gzip.open(self.gff_filename, 'rb') as f:
#        gff = f.read().decode('utf-8').split('\n')
#        gff.pop()
#        gff = [g for g in gff if g[0] != '#']
#        gff = [g for g in gff if (g.split('\t')[2] == 'gene' and g.split('\t')[0] == 'NC_002663.1')]
#        self.annotated_gff = gff
#
#    old_tag_idx = self.generate_index(gff, 'old_locus_tag')
#    tag_idx = self.generate_index(gff, 'locus_tag')
#    annotations = self.get_kegg_gene_annotations('pathway')
#    remove_header = lambda x : re.findall(r'.*:(.*)',x)[0]
#
#    log = open(self.kegg_genome_id + '_failed_kegg_annotation.log','w')
#    for gene_id, kegg_anno in annotations:
#      gene_id = remove_header(gene_id)
#      idx = old_tag_idx.get(gene_id, -1) 
#      if idx == -1:
#        idx = tag_idx.get(gene_id, -1)
#      if idx == -1:
#        log.write(gene_id + '\n')
#        continue
#      self.annotated_gff[idx] += '\t' + kegg_anno
#    log.close()
#
#  def write_annotated_gff(self):
#    data = [(int(e.split('\t')[3]),e) for e in self.annotated_gff]
#    data = sorted(data, key=lambda x: x[0])
#    self.annotated_gff = [e for _,e in data]
#    with gzip.open(self.gff_filename, 'wb') as f:
#      f.write('\n'.join(self.annotated_gff))
