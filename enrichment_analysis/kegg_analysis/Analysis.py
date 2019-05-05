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
from sklearn.metrics import jaccard_similarity_score
import scipy.stats as stats


PATHWAY  = 'pathway'
PICKLE_STORE = './pickles/'
ASM_PATH = './assemblies/'

START = 3
END = 4
MIN_LABELLING = 0.9
INC = 0.25
np.random.seed()
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

def normalized_log_distance(ori, genome_length, feature_pos):
# First, compute linear distance
  if abs(feature_pos - ori) < (genome_length / 2.0):
    return np.log2(abs(feature_pos - ori)) / np.log2(genome_length / 2.0)
  return np.log2(genome_length - abs(feature_pos - ori)) / np.log2(genome_length / 2.0)


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
    self.discard_org = False
    self.gene_tags = dict()
    self.old_gene_tags = dict()

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
    get_tag = lambda a,b : re.findall(rf'[^\w]{a}=([\w]+)',b)
    for idx in range(0, len(asm)):
      line = asm[idx]
      gene_tag = get_tag('locus_tag', line)
      old_gene_tag = get_tag('old_locus_tag', line)
      if len(gene_tag) > 0:
        self.gene_tags[gene_tag[0]] = idx
      if len(old_gene_tag) > 0:
        self.old_gene_tags[old_gene_tag[0]] = idx

    # Generate kegg index
    kegg_index = dict()
    n_added_keggs = 0
    for kegg, tags in kegg_tag_idx.items():
      indices = set()
      added = False
      for tag in tags:
        if tag in self.gene_tags:
          indices.add(self.gene_tags[tag])
          added = True
        if tag in self.old_gene_tags:
          indices.add(self.old_gene_tags[tag])
          added = True
      if added:
        n_added_keggs += 1 
      kegg_index[kegg] = sorted(list(indices))
    if (n_added_keggs / float(len(kegg_tag_idx.keys()))) < MIN_LABELLING:
        print(f'Discarded {self.kegg_gn}')
        self.discard_org = True
    return kegg_index
 
  def compute_distance_vector(self, ori, kegg_syms, dist, summary):
    """For each possible kegg symbol generates a summary statistic for the
       distance from ori for each gene assigned to the symbol. """
    distance_vec = np.full(len(kegg_syms), -1, dtype='float64')
    # vectorise the input distance function
    vec_dist = np.vectorize(lambda n : dist(ori, self.genome_length, n))
    for i in range(0, len(kegg_syms)):
      kegg_sym = kegg_syms[i]
      if kegg_sym not in self.kegg_index:
        continue
      indices = self.kegg_index[kegg_sym]
      if len(indices) == 0:
        continue
      distances = vec_dist(self.position_vec_start[indices])
      distance_vec[i] = summary(distances) 
    return distance_vec 

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
  if len(sys.argv) != 4:
    print("Usage: <exe> <organism_list> <kegg_sym_of_interest> <house_keeping_gene_list>")
    sys.exit()
  organism_list = sys.argv[1]
  kegg_sym      = sys.argv[2]
  house_keeping = sys.argv[3]

  # Parse organism list
  with open(organism_list) as f:
    organism_list = [l.strip().split('\t') for l in f]

  # Generate list of organisms. 
  # If Organism has previously been generated and pickled, load pickle
  # rather than re-downloading from kegg.
  organisms = list()
  print('Loading data')
  for ori_id, asm_id, refseq, hit, kegg_gn, ori_pos, fname, genome_length in organism_list:
    sys.stdout.flush()
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
      except Exception as e:
        print(e)
        print(f'Failed to scrape {kegg_gn} {refseq} {asm_id}')
        continue
  print('Filtering organisms')
  # Discard organisms whom enough kegg symbols could not be assigned to their genome
  organisms = [o for o in organisms if o.discard_org == False]

  # Find the set of all organism that share all housekeeping genes
  # kkg_names is a list of dictionaries where the ith element of the list is k=kegg_gn v=gene_name pairs
  # and the gene_name is the name of the ith housekeeping gene in k
  hkg_names = list()
  hkg_files = list()
  remove_sp_name = lambda x : re.findall(r'[a-z]:(.*)',x)[0]
  with open(house_keeping) as f:
    hkg_files = [l.strip() for l in f]
    for hkg_file in hkg_files:
      with open(hkg_file) as f:
        dump = [l.strip().split(',') for l in f]
        hkg_names.append({kegg_genome: remove_sp_name(gene_name) for kegg_genome, gene_name in dump})

  intersect = set(hkg_names[0].keys())
  for d in hkg_names:
    intersect = intersect.intersection(set(d.keys()))
  # Make organisms and hkg dictionaries contain only the intersection
  organisms = [o for o in organisms if o.kegg_gn in intersect]
  for i,d  in enumerate(hkg_names):
    hkg_names[i] = {k:v for k,v in d.items() if k in intersect}

  # Download all the kegg symbols instances for the kegg symbol of interest.
  # e.g if looking at kegg pathways, download all possible pathways
  kegg_sym_instances = get_kegg_sym_instances(kegg_sym)

  # Begin real variance analysis. 

  # dist_matrix is size |organisms| X |all_kegg_symbols|. 
  # Entry i,j contains the distance (summary statistic) of all genes
  # of kegg_symbol j from species i to ori in species i.
  print('computing real var')
  ori_vec  = [o.ori_pos for o in organisms]
  dist_matrix = compute_dist_matrix(organisms, kegg_sym_instances, normalized_log_distance, np.median, ori_vec)

  # Compute the variance of the median distance from Ori to each KEGG group across 
  # all strains
  real_variance = compute_variance(dist_matrix, len(kegg_sym_instances))

  # Begin variance analysis for each house keeping gene
  dims = (len(hkg_names), len(kegg_sym_instances))
  hkg_variance = np.full(dims, -1, dtype='float64')
  for i, hkg_name_dict in enumerate(hkg_names):
    print(f'computing {hkg_files[i]}')
    position_vec = list()
    for o in organisms:
      index = o.gene_tags.get(hkg_name_dict[o.kegg_gn], -1)
      if index == -1:
        index = o.old_gene_tags.get(hkg_name_dict[o.kegg_gn],-1)
      if index == -1:
        position_vec.append(-1)
      else:
        # set the "origin" to be the start of the housekeeping gene
        position_vec.append(o.position_vec_start[index])
    # Compute dist matrix using housekeeping genes as origin
    dist_matrix = compute_dist_matrix(organisms, kegg_sym_instances, normalized_log_distance, np.median, position_vec)
    # Compute variance
    hkg_variance[i] = compute_variance(dist_matrix, len(kegg_sym_instances))

  # Perform mann-whitney u tests
  # Against all housekeeping genes
  all_hkg_vars = hkg_variance.ravel()
  print('real vs all house keeping genes')
  print(f'two-sided, {stats.mannwhitneyu(real_variance, all_hkg_vars, alternative="two-sided")}')
  print(f'less, {stats.mannwhitneyu(real_variance, all_hkg_vars, alternative="less")}')
  print(f'greater, {stats.mannwhitneyu(real_variance, all_hkg_vars, alternative="greater")}')

 # Against each house keeping gene
  print('real vs each house keeping gene')
  for i, hkg_var in enumerate(hkg_variance):
    print(hkg_files[i])
    print(f'two-sided, {stats.mannwhitneyu(real_variance, hkg_variance[i], alternative="two-sided")}')
    print(f'less, {stats.mannwhitneyu(real_variance, hkg_variance[i], alternative="less")}')
    print(f'greater, {stats.mannwhitneyu(real_variance, hkg_variance[i], alternative="greater")}')


def compute_variance(dist_matrix, col_len):
  filter_nulls = lambda x : x[x != -1]
  var_vec = np.full(col_len, -1, dtype='float64')
  for i in range(0,col_len):
    values = filter_nulls(dist_matrix[:,i])
    if len(values) == 0:
    # occurs if filter_nulls feeds empty list to np.var()
      var_vec[i] = -1
    else:
      var_vec[i] = np.var(values)
  return var_vec


def compute_dist_matrix(organisms, kegg_sym_instances, dist, summary, ori_vec):
  dims = (len(organisms), len(kegg_sym_instances))
  dist_matrix = np.zeros(dims, dtype='float64')

  for idx, org in enumerate(organisms):
    if ori_vec[idx] == -1:
    # Then we failed to find the housekeeping gene for this strain, therefore
    # set its dist vect to -1 so it is skipped when computing variance
      dist_matrix[idx] = np.full(len(kegg_sym_instances), -1)
    else:
      dist_matrix[idx] = org.compute_distance_vector(ori=ori_vec[idx],
        kegg_syms=kegg_sym_instances,
        dist=normalized_linear_distance,
        summary=np.median
      )
  return dist_matrix

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


#def jaccard_computation(organisms,kegg_sym_instances):
#  # compute np array to store each bin matrix
#  kegg_to_vec = {int(kegg_sym):i for i, kegg_sym in enumerate(kegg_sym_instances)}
#  dims = (len(organisms), int(1/INC), len(kegg_sym_instances))
#  org_matrices = np.zeros(dims, dtype=bool)
#  print('Computing real jaccard_indices')
#  results = open('results.csv','w')
#  for i, org in enumerate(organisms):
#    org_matrices[i] = org.compute_bin_vector(INC, kegg_to_vec)
#  and_reduce = np.bitwise_and.reduce(org_matrices,axis=0)
#  or_reduce = np.bitwise_or.reduce(org_matrices,axis=0)
#
#  for i in range(0, int(1/INC)):
#    results.write(f'real,{i},{(float(np.sum(and_reduce[i])) / np.sum(or_reduce[i]))}\n')
#  np.set_printoptions(threshold=sys.maxsize)
#  for it in range(0, 200):
#    # Permute positions
#    print(f'iteration {it}')
#    print('Randomizing genome locations')
#    for o in organisms:
#      o.norm_position_vec_start, o.norm_position_vec_end =  shuffle_in_unison(o.norm_position_vec_start, o.norm_position_vec_end)
#    print('Computing jaccard on randomized positions')
#    org_matrices = np.zeros(dims, dtype=bool)
#    for i, org in enumerate(organisms):
#      org_matrices[i] = org.compute_bin_vector(INC, kegg_to_vec)
#    and_reduce = np.bitwise_and.reduce(org_matrices,axis=0)
#    or_reduce = np.bitwise_or.reduce(org_matrices,axis=0)
#
#    for i in range(0, int(1/INC)):
#      results.write(f'rand,{i},{(float(np.sum(and_reduce[i])) / np.sum(or_reduce[i]))}\n')
#  results.close()
#
#def shuffle_in_unison(a, b):
#  assert len(a) == len(b)
#  shuffled_a = np.empty(a.shape, dtype=a.dtype)
#  shuffled_b = np.empty(b.shape, dtype=b.dtype)
#  permutation = np.random.permutation(len(a))
#  for old_index, new_index in enumerate(permutation):
#      shuffled_a[new_index] = a[old_index]
#      shuffled_b[new_index] = b[old_index]
#  return shuffled_a, shuffled_b
