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

START = 3
END = 4

def pickle_organism(org,pickle_name):
  with bz2.BZ2File(pickle_name,'wb') as f:
    pickle.dump(org,f)

def unpickle_organism(filename):
  with bz2.open(filename,'rb') as f:
    org = pickle.load(f)
    org.genome_length = int(org.genome_length)
    return org

def get_kegg_sym_instances(kegg_sym):
    url = f'http://rest.kegg.jp/list/{kegg_sym}'
    all_kegg_syms = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    all_kegg_syms.pop()
    if kegg_sym == PATHWAY:
      remove_hdr = lambda x: re.findall(r'\S*:map([0-9]*)', x)[0]
    return sorted([remove_hdr(kegg.split('\t')[0]) for kegg in all_kegg_syms])


def linear_distance(ori, genome_length, feature_pos): 
  if abs(feature_pos - ori) < (genome_length / 2.0):
    return abs(feature_pos - ori)
  return genome_length - abs(feature_pos - ori)

def normalized_linear_distance(ori, genome_length, feature_pos): 
  """ Normalise along the Ori-Ter axis. Ori =0 , Ter = 1"""
  if abs(feature_pos - ori) < (genome_length / 2.0):
    return abs(feature_pos - ori) / (genome_length / 2.0)
  return (genome_length - abs(feature_pos - ori)) / (genome_length / 2.0)

def log_distance(ori, genome_length, feature_pos): 
  ter = (ori + genome_length/2.0) % genome_length
  if feature_pos > ter:
    return abs(ori - feature_pos)
  else:
    return genome_length - (ori - feature_pos)

class Organism:

  def __init__(self, ori_id, asm_id, refseq, kegg_gn, ori_pos, asm_file, kegg_sym, genome_length):
    self.ori_id = ori_id
    self.asm_id = asm_id
    self.refseq = refseq
    self.kegg_gn = kegg_gn
    self.ori_pos = int(ori_pos)
    self.asm_file = asm_file
    self.kegg_sym = kegg_sym
    self.genome_length = int(genome_length)

    self.position_vec_start = self.get_position_vec(START)
    self.position_vec_end   = self.get_position_vec(END)
    self.dist_from_ori = np.full(len(self.position_vec_start),0)
    for i,p in enumerate(self.position_vec_start):
      self.dist_from_ori[i] = linear_distance(self.ori_pos, self.genome_length, p)
      
    self.norm_position_vec_start = self.normalize_position_vec(START)
    self.norm_position_vec_end   = self.normalize_position_vec(END)
    self.kegg_index = self.get_kegg_index()


  def get_position_vec(self, pos):
    with gzip.open(self.asm_file,'r') as f:
      asm = f.read().decode('utf-8').split('\n')
      asm.pop()
      asm = [l.split('\t') for l in asm if l[0] != '#']
    return np.array([int(elem[pos]) for elem in asm if elem[2] == 'gene' and elem[0] == self.refseq])

  def get_kegg_index(self):
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
    vec_dist = np.vectorize(lambda n : dist(ori, self.genome_length, n))
    for i in range(0, len(kegg_syms)):
      kegg_sym = kegg_syms[i]
      if kegg_sym not in self.kegg_index:
        continue
      indices = self.kegg_index[kegg_sym]
      if len(indices) == 0:
        continue
      distances = vec_dist(self.position_vec[indices])
      analysis_vec[i] = summary(distances) 
    return analysis_vec

  def normalize_position_vec(self, pos):
    """ Normalize the entire positions vector along Ori Ter axis"""
    vec_dist = np.vectorize(lambda n : normalized_linear_distance(self.ori_pos, self.genome_length, n))
    if pos == START:
      return vec_dist(self.position_vec_start)
    return vec_dist(self.position_vec_end)

  def compute_bin_vector(self, increment, kegg_to_vec):
    """Computes a boolean vector of length |all kegg syms| which has value
      True at the correct kegg symbol position if a gene within the bin
      is annotated with the kegg symbol. """

    # Associated each kegg symbol in the organism with its list of normalised
    # start and end positions for each associated gene.
    kegg_to_starts = {kegg_sym : self.norm_position_vec_start[self.kegg_index[kegg_sym]] for kegg_sym in self.kegg_index.keys()}
    kegg_to_ends = {kegg_sym : self.norm_position_vec_end[self.kegg_index[kegg_sym]] for kegg_sym in self.kegg_index.keys()}

    # initialise matrix of bin vectors. dims = |bins| X |kegg_symbols|
    dims = (int(1/increment),len(kegg_to_vec.keys()))
    bin_vectors = np.zeros(dims,dtype=bool)

    # begin loop through bins
    bv_idx = 0
    for bin_start in np.arange(0,1,increment):
      bin_end = bin_start + increment
      for kegg_sym in self.kegg_index.keys():
        start_pos = kegg_to_starts[kegg_sym]
        end_pos   = kegg_to_ends[kegg_sym]
        if (np.any( ((bin_start <= start_pos) & (start_pos < bin_end)) |
           ((bin_start <= end_pos) & (end_pos < bin_end)) |
           ((start_pos <= bin_start) & (bin_end <= end_pos)))):
          bin_vectors[bv_idx][kegg_to_vec[int(kegg_sym)]] = True
      bv_idx += 1
    return bin_vectors

  def n_kegg_syms_in_bin(self, start, increment):
    count = 0
    for kegg_sym in self.kegg_index.keys():
      positions = self.norm_position_vec_start[self.kegg_index[kegg_sym]]
      if np.any((start <= positions) & (positions < (start + increment))):
        count += 1
    return count
      



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
  for ori_id, asm_id, refseq, hit, kegg_gn, ori_pos, fname, genome_length in organism_list[:10]:
    pickle_name = PICKLE_STORE + refseq + f'.{sys.argv[2]}' + '.pickle.bz2'
    if os.path.isfile(pickle_name): 
      organisms.append(unpickle_organism(pickle_name))
    else:
      try:
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
      except:
        print(f'Failed to scrape {kegg_gn} {refseq} {asm_id}')
        continue

  kegg_sym_instances = get_kegg_sym_instances(kegg_sym)
  kegg_to_vec = {int(kegg_sym):i for i, kegg_sym in enumerate(kegg_sym_instances)}
  inc = 0.001
  # compute np array to store each bin matrix
  dims = (len(organisms), int(1/inc), len(kegg_sym_instances))
  org_matrices = np.zeros(dims, dtype=bool)
  for i, org in enumerate(organisms):
    print(f'computing matrix for {org.kegg_gn}')
    org_matrices[i] = org.compute_bin_vector(inc, kegg_to_vec)

  for vec in org_matrices[0]:
    print(vec)

##for i in range(0,len(org.position_vec)):

#    print(f'{org.position_vec[i]} {org.dist_from_ori[i]} {org.norm_position_vec[i]}')


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
#  dist_matrix = np.zeros(dims, dtype='float64')
#
#  # Compute the real variance of the distance of each kegg symbol from ori in 
#  # all strains
#  for idx, org in enumerate(organisms):
#    print(org.kegg_gn)
#    dist_matrix[idx] = org.compute_analysis_vector(ori=org.ori_pos,
#      kegg_syms=kegg_sym_instances,
#      dist=linear_distance,
#      summary=np.median
#    )
#  real_variance = np.apply_along_axis(func1d=np.var, axis=0, arr=dist_matrix)
#  gluc = dist_matrix[:,0]
#  gluc = gluc[np.where(gluc != -1)]
#  print(np.var(gluc))
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




