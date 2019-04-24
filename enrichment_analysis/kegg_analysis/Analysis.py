import os
import urllib.request
import itertools as it
import re
import gzip
import numpy as np
import pickle
import bz2


def pickle_organism(org):
  with bz2.BZ2File(f'{org.refseq}.pickle.bz2','wb') as f:
    pickle.dump(org,f)

def unpickle_organism(filename):
  with bz2.open(filename,'rb') as f:
    return pickle.load(f)


class Organism:

  def __init__(self, ori_id, gcf, refseq, kegg_gn, ori_pos, asm_file, kegg_sym):
    self.ori_id = ori_id
    self.gcf = gcf
    self.refseq = refseq
    self.kegg_gn = kegg_gn
    self.ori_pos = int(ori_pos)
    self.asm_file = asm_file
    self.kegg_sym = kegg_sym

    self.positions_vec = self.get_position_vec()
    self.kegg_index = self.get_kegg_index()


  def get_position_vec(self):
    with gzip.open(self.asm_file,'r') as f:
      asm = f.read().decode('utf-8').split('\n')
      asm.pop()
      asm = [l.split('\t') for l in asm if l[0] != '#']
    return np.array([int(elem[3]) for elem in asm if elem[2] == 'gene' and elem[0] == self.refseq])

  def get_kegg_index(self):
    # Download all genes for kegg organism along with the relevant kegg_sym.
    url = f'http://rest.kegg.jp/link/{self.kegg_sym}/{self.kegg_gn}'
    kegg_tag_idx = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    kegg_tag_idx.pop()

    # Make tuples out of kegg_tag_idx and group by kegg_sym
    remove_hdr = lambda x: re.findall(r'\S*:(\S*)', x)[0]
    kegg_tag_idx = [elem.split('\t') for elem in kegg_tag_idx]
    kegg_tag_idx = [(remove_hdr(gene), remove_hdr(kegg_sym)) for gene, kegg_sym in kegg_tag_idx]
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
        
def main():
  borg = Organism(ori_id ='ORI10030004', gcf= 'GCF_000008745.1',refseq= 'NC_000922.1', 
                 kegg_gn='T00021',ori_pos='-1', 
                 asm_file='assemblies/GCF_000008745.1_ASM874v1_genomic.gff.gz', kegg_sym='pathway')
  pickle_organism(borg)
  org = unpickle_organism('NC_000922.1.pickle.bz2')
  print(org.positions_vec)
  print(len(org.kegg_index.values()))
  idx = org.kegg_index['cpn04122']
  for i in idx:
    print(org.positions_vec[i])
  for kegg, indices in org.kegg_index.items():
    print(kegg)
    print('\t'.join([str(e) for e in indices]))

if __name__ == '__main__':
  main()
