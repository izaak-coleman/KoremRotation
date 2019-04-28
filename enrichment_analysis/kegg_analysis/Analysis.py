import os
import sys
import urllib.request
import itertools as it
import re
import gzip
import numpy as np
import pickle
import bz2
import time


PATHWAY  = 'pathway'

PICKLE_STORE = './pickles/'
ASM_PATH = './assemblies/'
def pickle_organism(org,pickle_name):
  with bz2.BZ2File(pickle_name,'wb') as f:
    pickle.dump(org,f)

def unpickle_organism(filename):
  with bz2.open(filename,'rb') as f:
    return pickle.load(f)

def get_kegg_sym_instances(kegg_sym):
    url = f'http://rest.kegg.jp/list/{kegg_sym}'
    all_kegg_syms = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    all_kegg_syms.pop()
    if kegg_sym == PATHWAY:
      remove_hdr = lambda x: re.findall(r'\S*:map([0-9]*)', x)[0]
    return sorted([remove_hdr(kegg.split('\t')[0]) for kegg in all_kegg_syms])


def linear_distance(ori, genome_length, feature_pos): 
  ter = (ori + genome_length/2.0) % genome_length
  if feature_pos > ter:
    return abs(ori - feature_pos)
  else:
    return genome_length - (ori - feature_pos)

def log_distance(ori, genome_length, feature_pos): 
  ter = (ori + genome_length/2.0) % genome_length
  if feature_pos > ter:
    return abs(ori - feature_pos)
  else:
    return genome_length - (ori - feature_pos)

class Organism:

  def __init__(self, ori_id, asm_id, refseq, kegg_gn, ori_pos, asm_file, kegg_sym,genome_length):
    self.ori_id = ori_id
    self.asm_id = asm_id
    self.refseq = refseq
    self.kegg_gn = kegg_gn
    self.ori_pos = int(ori_pos)
    self.asm_file = asm_file
    self.kegg_sym = kegg_sym
    self.genome_length = genome_length

    self.positions_vec = self.get_position_vec()
    self.kegg_index = self.get_kegg_index()


  def get_position_vec(self):
    print(f'Building position vec for {self.kegg_gn}.')
    with gzip.open(self.asm_file,'r') as f:
      asm = f.read().decode('utf-8').split('\n')
      asm.pop()
      asm = [l.split('\t') for l in asm if l[0] != '#']
    return np.array([int(elem[3]) for elem in asm if elem[2] == 'gene' and elem[0] == self.refseq])

  def get_kegg_index(self):
    print(f'Building kegg index for {self.kegg_gn}.')
    # Download all genes for kegg organism along with the relevant kegg_sym.
    url = f'http://rest.kegg.jp/link/{self.kegg_sym}/{self.kegg_gn}'
    kegg_tag_idx = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    kegg_tag_idx.pop()

    # Make tuples out of kegg_tag_idx and group by kegg_sym
    remove_hdr = lambda x: re.findall(r'\S*:(\S*)', x)[0]
    remove_sp = lambda x: re.findall(r'[a-zA-Z]{3,4}([0-9]*)', x)[0]
    kegg_tag_idx = [elem.split('\t') for elem in kegg_tag_idx]
    kegg_tag_idx = [(remove_hdr(gene), remove_sp(remove_hdr(kegg_sym))) for gene, kegg_sym in kegg_tag_idx]
    kegg_tag_idx = it.groupby(sorted(kegg_tag_idx, key = lambda x: x[1]), key = lambda x : x[1]) 
    kegg_tag_idx = {kegg:[gene_tag for gene_tag,_ in group] for kegg, group in kegg_tag_idx}

    # Make locus_tag and old_locus_tag vectors
    with gzip.open(self.asm_file,'r') as f:
      asm = f.read().decode('utf-8').split('\n')
    asm.pop()
    asm = [l.split('\t') for l in asm if l[0] != '#']
    asm = ['\t'.join(elem) for elem in asm if elem[2] == 'gene' and elem[0] == self.refseq] 
    gene_tags, old_gene_tags = dict(), dict()
    get_tag = lambda a,b : re.findall(rf'[^\w]{a}=([\w]+)',b)
    for idx in range(0, len(asm)):
      line = asm[idx]
      gene_tag = get_tag('locus_tag', line)
      old_gene_tag = get_tag('old_locus_tag', line)
      if len(gene_tag) > 0:
        gene_tags[gene_tag[0]] = idx
      if len(old_gene_tag) > 0:
        old_gene_tags[old_gene_tag[0]] = idx

    # Generate kegg index
    kegg_index = dict()
    for kegg, tags in kegg_tag_idx.items():
      indices = set()
      for tag in tags:
        if tag in gene_tags:
          indices.add(gene_tags[tag])
        if tag in old_gene_tags:
          indices.add(old_gene_tags[tag])
      kegg_index[kegg] = sorted(list(indices))
    return kegg_index
 
  def compute_analysis_vector(self, ori, kegg_syms, dist, summary):
    """For each possible kegg symbol generates a summary statistic for the
       distance from ori for each gene assigned to the symbol. """
    analysis_vec = np.full(len(kegg_syms), -1, dtype='float64')
    # vectorise the input distance function
    vec_dist = np.vectorize(lambda n : dist(self.ori_pos, self.genome_length, n))
    for i in range(0, len(kegg_syms)):
      kegg_sym = kegg_syms[i]
      if kegg_sym not in self.kegg_index:
        continue
      indices = self.kegg_index[kegg_sym]
      distances = vec_dist(self.positions_vec[indices])
      print('%s: %s' % (kegg_sym, ', '.join([str(d) for d in sorted(distances)])))
      analysis_vec[i] = summary(distances) 
    return analysis_vec


def main():
  if len(sys.argv) != 3:
    print("Usage: <exe> <organism_list> <kegg_sym_of_interest>")
    sys.exit()
  organism_list = sys.argv[1]
  kegg_sym      = sys.argv[2]

  # Parse organism list
  with open(organism_list) as f:
    organism_list = [l.strip().split('\t') for l in f]

  # Generate list of organisms. 
  # If Organism has previously been generated and pickled, load pickle
  # rather than re-downloading from kegg.
  organisms = list()
  for ori_id, asm_id, refseq, hit, kegg_gn, ori_pos, fname, genome_length in organism_list:
    pickle_name = PICKLE_STORE + refseq + f'.{sys.argv[2]}' + '.pickle.bz2'
    if os.path.isfile(pickle_name): 
      organisms.append(unpickle_organism(pickle_name))
    else:
      org = Organism(ori_id = ori_id,
        asm_id = asm_id,
        refseq = refseq,
        kegg_gn = kegg_gn,
        ori_pos = ori_pos,
        asm_file = ASM_PATH + fname,
        kegg_sym = kegg_sym,
        genome_length = genome_length
      )
      organisms.append(org)
      pickle_organism(org,pickle_name)

  for org in organisms:
    key, val = list(org.kegg_index.items())[1]
    print(f'{org.kegg_gn} {org.asm_id} {key}')
    print(f'{org.positions_vec[val]}')

#
#  # Download all the kegg symbols instances for the kegg symbol of interest.
#  # e.g if looking at kegg pathways, download all possible pathways
#  kegg_sym_instances = get_kegg_sym_instances(kegg_sym)
#
#
#  # Begin analysis. 
#
#  # dist_matrix is size |organisms| X |all_kegg_symbols|. 
#  # Entry i,j contains the distance (summary statistic) of all genes
#  # of kegg_symbol j from species i to ori in species i.
#  dims = (len(organisms), len(kegg_sym_instances))
#  dist_matrix = np.zeros(dim, dtype='float64')
#
#  # Compute the real variance of the distance of each kegg symbol from ori in 
#  # all strains
#  for idx, org in enumerate(organisms):
#    dist_matrix[idx] = org.compute_analysis_vector(ori=org.ori_pos,
#      kegg_syms=kegg_sym_instances,
#      dist=linear_distance,
#      np.median
#    )
#  print(dist_matrix)
#  real_variance = np.apply_along_axis(np.var, axis=0, dist_matrix)
#
#  # Begin permutation analysis

  

#def main():
#  borg = Organism(ori_id ='ORI10030004', gcf= 'GCF_000008745.1',refseq= 'NC_000922.1', 
#                 kegg_gn='T00021',ori_pos='841385', 
#                 asm_file='assemblies/GCF_000008745.1_ASM874v1_genomic.gff.gz', 
#                 kegg_sym='pathway', genome_length=1230230)
#  pickle_organism(borg)
#  org = unpickle_organism('NC_000922.1.pickle.bz2')
#  all_kegg_syms = get_all_kegg_syms(PATHWAY)
#  av = org.generate_analysis_vector(org.ori_pos, all_kegg_syms, linear_dist, np.mean)
#  for i in range(0, len(all_kegg_syms)):
#    print(f'{all_kegg_syms[i]}: {av[i]}')
#

if __name__ == '__main__':
  main()




