import sys

FEATURE_START = 3
FEATURE_END   = 4
def main():
  if len(sys.argv) != 3:
    print('Usage <exe> <.ptx> <.gff3.>')

  with open(sys.argv[1]) as f:
    ptx = [l.strip() for l in f]
    ptx.pop(0)
    ptx = [(int(e[3]), int(e[4])) for e in ptx]

  with open(sys.argv[2]) as f:
    gff = [l.strip().split('\t') for l in f]

  gff_tox = []
  gff_tox_idx = []
  for i in range(0, len(gff)):
    for e in ptx:
      if (e[0] == gff[i][FEATURE_START]) and (e[1] == gff[i][FEATURE_END]):
        gff_tox.append(gff[i])
        gff_tox_idx.append(i)

  # remove all toxic genes from genome
  for n in gff_tox_idx:
    gff.pop(n)

  print('PANDATOX GENES')
  print('\n'.join(['\t'.(e) for e in gff_tox]))

  print('NON_PANDATOX GENES')
  print('\n'.join(['\t'.(e) for e in gff]))

if __name__ == '__main__':
  main()
