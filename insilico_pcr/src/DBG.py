# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

from graphviz import Digraph 
import pickle, bz2

class DBG:
  """Stores a python based De Bruijn graph implementation that allows building
     on the fly."""

  def __init__(self, sequence, k):
    self.nodes = set()
    self.edges = dict()
    self.k = int(k)
    for i in range(0, len(sequence) - k):
      self.add_edge(sequence[i : i+k+1])

  def add_edge(self, sequence):
    """Extends the De Bruijin graph by a single edge. 

       If the suffix of the edge is present in the De Brujin graph, 
       add_edge returns false."""
    if len(sequence) != self.k + 1: raise Exception(f'edge {sequence} is not of length {self.k + 1} (k + 1)')
    edge = (sequence[:self.k], sequence[1:])

    if edge in self.edges:
      self.edges[edge] += 1
    else:
      self.edges[edge] = 1

    prefix, suffix = edge
    self.nodes.add(prefix)

    if suffix in self.nodes:
      return False
    else:
      self.nodes.add(suffix)
      return True


def pickle_vec(vec,pickle_name):
  with bz2.BZ2File(pickle_name,'wb') as f:
    pickle.dump(vec,f)

def main():
  dbg = DBG('AGTCAGATCGGATCGGTCGTATCGGATCGGGCGTGTAGTATCGGATCGGGCAG', k=10)
  pickle_vec(dbg, 'simple.dbg.pickle')

if __name__ == '__main__':
  main()
