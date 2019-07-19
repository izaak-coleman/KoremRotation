from cortexpy.graph.parser.random_access import RandomAccess
import sys
import array as arr
import copy

KMER_END = 8
IN_NODE = -2
OUT_NODE = 2
LEFT_MASK  = int('00111111', 2)
RIGHT_MASK = int('11111100', 2)

EdgeMap = arr.array('B',[
int('00000000',2),
int('01000000',2),
int('10000000',2),
int('11000000',2),
int('00000011',2),
int('00000010',2),
int('00000001',2),
int('00000000',2)
])

def writeAdjacencyList(graph, fname):
  """Writes an adjacency list of the input graph in .sif format"""
  sif = open(fname, 'w')
  sif.write(f'{graph.n_records}\n')

  # Aux functions
  edgeExists = lambda x: True if x == 1 else False
  get_kmer_bytes = lambda x : bytearray(x._kmer_data._data[:KMER_END][::-1])

  for kmer in graph.values():
    byte_kmer = get_kmer_bytes(kmer)
    for (idx, edge) in enumerate(kmer.edges[0]):
      if idx < 4 and edgeExists(edge):
        sif.write(f'{get_adjacent_kmer(copy.copy(byte_kmer), IN_NODE, idx)}' + ' ' + f'{to_int(byte_kmer)}\n')
      elif idx >= 4 and edgeExists(edge):
        sif.write(f'{to_int(byte_kmer)}' + ' ' + f'{get_adjacent_kmer(copy.copy(byte_kmer), OUT_NODE, idx)}\n')
  sif.close()

def to_int(kmer):
  return int.from_bytes(kmer, byteorder='big')


def get_adjacent_kmer(kmer, bitshift, idx):
  shift = lambda x,s : int(to_int(x)*(2**s)).to_bytes(8, byteorder='big')
  # Delete unneccesary base and shift
  if bitshift == 0:
    return kmer
  elif bitshift > 0:
    kmer[0] = kmer[0] & LEFT_MASK
  else: # shift < 0
    kmer[-1] = kmer[-1] & RIGHT_MASK
  kmer = bytearray(shift(kmer, bitshift))
  # Write new base to kmer 
  if bitshift > 0:
    kmer[-1] = kmer[-1] | EdgeMap[idx]
  else: # bitshift < 0
    kmer[0] = kmer[0] | EdgeMap[idx]
  return to_int(kmer)


def main():
  if len(sys.argv) != 3: 
    print("Usage: <exe> <.ctx> <fname>")
    sys.exit()

  graph = RandomAccess(open(sys.argv[1], 'rb'))
  writeAdjacencyList(graph, sys.argv[2])


if __name__ == '__main__':
  main()
