# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

from graphviz import Digraph 
import sys

DBS = 0
FREQ = 1
IN, OUT = 0, 1
UNDEFINED = 1

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
    
    prefix, suffix = self.get_nodes(edge)
    self.nodes.add(prefix)
    self.nodes.add(suffix)
    if edge in self.edges:
      self.edges[edge][DBS].add(dbs)
      self.edges[edge][FREQ] += 1
    else:
      self.edges[edge] = (set(dbs), 1)

  def get_nodes(self, sequence):
    """Takes a sequence and splits it into the prefix and suffix node
       of an edge. """
    if len(sequence) != self.k + 1: 
      raise Exception(f'edge {sequence} is not of length {self.k + 1} (k + 1)')
    return (sequence[:self.k], sequence[1:])

  def render(self, fname):
    print("EDGES FROM RENDER")
    for edge in self.edges:
      print(edge)
    print("NODES FROM RENDER")
    for node in self.nodes:
      print(node)
    """Renders the dbg into a PDF format."""
    g = Digraph()
    for edge, meta in self.edges.items():
      dbs, freq = meta
      p, s = edge[:self.k], edge[1:]
      for i in range(0, freq):
        g.edge(p, s, label = s[-1])
    g.render(fname) 

  def abbreviate(self, seq):
    if len(seq) > 6:
      return seq[:3] + '...' + seq[-3:]

  def compress(self):
    """Compresses the De Bruijn graph, replacing unitig paths with single nodes. """
    # Construct all unitigs of length > 1. 
    # A unitig is path of nodes v_1, ..., v_n, where indegree(v_i) = outdegree(v_i) = 1 
    # for 1 < i < j, and indegree(v_1) != 1, outdegree(v_1) = 1, indegree(v_n) = 1, 
    # outdegree(v_n) != 1. 
    unitigs = set()
    for node in self.nodes:
      unitig_path = self.get_unitig(node)
      if len(unitig_path) > 1:
        unitigs.add(tuple(unitig_path))

    # Compress each unitig in the dbg
    for unitig_path in unitigs:
      # Reconstruct all the edges from the unitig node
      # and then remove each edge and their corresponding
      # nodes from the dbg 
      edges = [unitig_path[i] + unitig_path[i+1][-1] for i in range(0, len(unitig_path) - self.k)]
      for edge in edges:
        del self.edges[edge]
      for node in unitig_path:
        self.nodes.remove(node)
      # Add unitig node to dbg
      unitig = unitig_path[0] + ''.join([node[-1] for node in unitig_path[1:]])
      self.nodes.add(unitig)
      print(f'unitig {unitig}')
      # Remove all remaining edges between unitig_path[0], unitig_path[1] and rest of dbg
      # then connect the unitig node with its in/out-adjacent nodes.
      in_unitig_nodes, out_unitig_nodes = self.get_adjacent(unitig_path[0], IN), self.get_adjacent(unitig_path[-1], OUT)
      for node in in_unitig_nodes:
        del self.edges[node + unitig[self.k-1]]
        self.edges[node[0] + unitig] =  (UNDEFINED, UNDEFINED)
      for node in out_unitig_nodes:
        del self.edges[unitig[-self.k] + node]
        self.edges[unitig + node[-1]] = (UNDEFINED, UNDEFINED)

  def get_unitig(self, node):
    """Finds the unitig of a node."""
    in_adj, out_adj = self.get_adjacent(node, IN), self.get_adjacent(node, OUT)
    if len(in_adj) == 1 and len(out_adj) == 1:
      return self.extend_unitig(in_adj[0], IN) + [node] + self.extend_unitig(out_adj[0], OUT)
    return [node]

  def get_adjacent(self, node, direction):
    ALPHA = 'ACGT'
    if direction == IN:
      return [c + node[:-1] for c in ALPHA if (c + node[:-1])  in self.nodes]
    else: # direction == OUT
      return [node[1:] + c  for c in ALPHA if (node[1:] + c)   in self.nodes]

  def extend_unitig(self, node, direction):
    in_nodes, out_nodes = self.get_adjacent(node, IN), self.get_adjacent(node, OUT)
    # if there is more than one edge out of the node in the direction the recursion
    # came from (we have reached a node that then branches into our unitig, and
    # another, do not include this node in the unitig
    if direction == IN and len(out_nodes) > 1:
      return []
    elif direction == OUT and len(in_nodes) > 1:
      return []
    # If there is one edge out of the node in the direction the recursion came from
    # (we have hit a node with no further edges in the way we are going,
    # or that in the direction we are heading) include the node in the unitig
    if direction == IN and len(in_nodes) != 1:
      return [node]
    elif direction == OUT and len(out_nodes) != 1:
      return [node]
    # Otherwise, indegree = outdegree == 1, so move in same direction to the next node
    if direction == IN:
      return self.extend_unitig(in_nodes[0], direction) + [node]
    else: # direction == OUT
      return [node] + self.extend_unitig(out_nodes[0], direction)
