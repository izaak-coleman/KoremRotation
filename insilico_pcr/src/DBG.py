# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

from graphviz import Digraph 

DBS = 0
FREQ = 1
class DBG:
  """Stores a python-based De Bruijn graph implementation that allows building
     on the fly."""

  def __init__(self, k):
    self.nodes = set()
    self.edges = dict() # k = edge, v = (set of dbs, frequency)
    self.k = int(k)


  def add_edge(self, edge, dbs):
    """Extends the De Bruijin graph by a single edge. 

       If the suffix of the edge (prefix, suffix) is present in the De Brujin graph, 
       add_edge returns false."""
    
    prefix, suffix = (edge[:self.k], edge[1:])
    self.nodes.add(prefix)
    self.nodes.add(suffix)
    if edge in self.edges:
      self.edges[edge][DBS].add(dbs)
      self.edges[edge][FREQ] += 1
    else:
      self.edges[edge] = (set(dbs), 1)

  def get_nodes(self, sequence)
    """Takes a sequence and splits it into the prefix and suffix node
       of an edge. """
    if len(sequence) != self.k + 1: raise Exception(f'edge {sequence} is not of length {self.k + 1} (k + 1)')
    return (sequence[:self.k], sequence[1:])

  def render(self, fname):
    """Renders the dbg into a PDF format."""
    g = Digraph()
    for edge, meta in seld.edges.items():
      dbs, freq = meta
      p, s = edge[:self.k], edge[1:]
      for i in range(0, freq):
        g.edge(p,s)
    g.render(fname) 
