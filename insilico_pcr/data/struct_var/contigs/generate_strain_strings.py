import sys
import csv
import gzip

def main():
  if len(sys.argv) != 3:
    print("Usage: <exe> <contig_concat> <path_to_contigs>")
    sys.exit()

  with open(sys.argv[1]) as f:
    contig_list = [line for line in csv.reader(f)]
    contig_list.pop(0)

  for strain, contig_id, pos in contig_list:
    print(contig_id)
    with gzip.open(sys.argv[2] + contig_id + '.fasta.gz', 'r') as f:
      fa = [l.strip().decode('utf-8') for l in f]
      fa.pop(0)
    with open(strain + '.concat', 'a') as g:
      g.write(''.join(fa))


if __name__ == '__main__':
  main()
      
    
