import collections
import sys
import re

Assembly = collections.namedtuple('Assembly', ['filename', 'assembly', 'refseq', 'length','contigs'])

def get_assembly_info(assembly_data):
  with open(assembly_data) as f:
    data = [l.strip() for l in f]
  delimiters = [i for i in range(0, len(data)) if 'gff.gz' == data[i][-6:]]
  entries = [data[ delimiters[i] : delimiters[i+1] ] for i in range(0, len(delimiters)-1)]
  entries.append(data[delimiters[-1]:])
  get_assembly = lambda x: re.findall(r'genome-build-accession.*:(GCF[\w\.]+)#',x)
  get_contigs = lambda x: re.findall(r'sequence-region ([\w\.]+) [0-9]+ ([0-9]+)#',x)

  assemblies = []
  for entry in entries:
    filename = entry.pop(0)
    entry = ''.join(entry) + '#'
    assembly = get_assembly(entry)
    contigs = get_contigs(entry)
    contigs = sorted([(contig, int(length)) for contig, length in contigs], key=lambda x : x[1])
    assemblies.append(Assembly(filename, assembly, '', -1, contigs))
  return assemblies

def main():
  if len(sys.argv) != 2:
    print('Usage <exe> <assembly.dat>')
    sys.exit()
  assemblies = get_assembly_info(sys.argv[1])
  for a in assemblies[:10]:
    print(f'{a.filename} {a.assembly} {a.refseq} {a.length} {a.contigs}')

if __name__ == '__main__':
  main()
