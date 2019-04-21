import collections

Assembly = collection.namedtuple('Assembly', ['filename', 'assembly', 'refseq', 'length','n_contigs'])

def get_assembly_info(assembly_data):
  with open(assembly_data) as f:
    data = [l.strip() for l in f]
  delimiters = [i for i in range(0, len(data)) if 'gff.gz' == data[i][-6:]]
  entries = [data[ delimiters[i] : delimiters[i+1] ] for i in range(0, len(delimiters)-1)]
  entries.append(data[delimiters[-1]:])
  get_assembly = lambda x: re.findall(r'genome-build-accession.*:(GCF[\w\.]+)#',x)
  get_refseqs  = lambda x: re.findall(r'sequence-region ([\w\.]+) [0-9]+ ([0-9]+)#',x)

  assemblies = []
  for entry in entries:
    filename = entry.pop(0)
    entry = ''.join(entry) + '#'
    assembly = get_assembly(entry)
    refseqs = get_refseqs(entry)
    refseqs = sorted([(refseq, int(length)) for refseq, length in refseqs], lambda x : x[1])
    refseq, length = refseq[0]
    assemblies.append(Assembly(filename, assembly, refseq, length, len(refseqs)))
  return assemblies

def main():
  pass
#test this!
