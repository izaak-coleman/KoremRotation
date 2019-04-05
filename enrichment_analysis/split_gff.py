import sys

FEATURE_START = 3
FEATURE_END   = 4

def write_gff_list(filename, gff):
  with open(filename,'w') as f:
    f.write('\n'.join(['\t'.join(e) for e in gff]))

def get_uniquely_located_features(gff):
  """Many gff features have multiple entries each with identical locations
     but different attributes, e.g gene and protein. This function
     returns a list of features that are uniquely_located."""
  seen = set()
  for f in gff:
    if (f[FEATURE_START], f[FEATURE_END]) not in seen:
      seen.add( (f[FEATURE_START], f[FEATURE_END]) )
      yield f

def compute_offset_from_ori(ori, genome_length, feature_pos):
  ter = (ori + genome_length/2.0) % genome_length
  if feature_pos > ter:
    return abs(ori - feature_pos)
  else:
    return genome_length - (ori - feature_pos)

def main():

  if len(sys.argv) != 5:
    print('Usage <exe> <.ptx> <.gff3.> <ori> <length>')
    sys.exit()

  base = sys.argv[1][:-4]
  ori = int(sys.argv[3])
  genome_length = int(sys.argv[4])

  with open(sys.argv[1]) as f:
    ptx = [l.strip().split(' \t') for l in f]
    ptx.pop(0)

  with open(sys.argv[2]) as f:
    gff = [l.strip().split('\t') for l in f]
    # Remove first two header lines and last empty line
    gff.pop(0)
    gff.pop(0)
    if gff[-1] == ['']:
      gff.pop()
    gff = list(get_uniquely_located_features(gff))


  gff_tox = []
  gff_tox_idx = []
  gff_near_tox_idx = []
  not_present = []
  for e in ptx:
    in_gff = False
    for i in range(0, len(gff)):
      if (e[FEATURE_START] == gff[i][FEATURE_START]) and (e[FEATURE_END] == gff[i][FEATURE_END]):
        gff_tox.append(gff[i])
        gff_tox_idx.append(i)
        in_gff = True
      # If feature is near a toxic region (same start or end), remove it from 
      # non-toxic gene list, but do not add it to list of pandatox genes. 
      # Hence, it still remains in the not present list
      # Why? If not a toxic gene, since there are so many non-toxic genes it is unlikely to affect result
      # If is a toxic gene, then it could affect the results.
      elif (e[FEATURE_START] == gff[i][FEATURE_START]) or (e[FEATURE_END] == gff[i][FEATURE_END]):
        gff_near_tox_idx.append(i)
    if not in_gff:
        not_present.append(e)

  # remove all toxic/near toxic genes from genome
  gff_non_tox = gff
  for n in sorted(gff_tox_idx + gff_near_tox_idx)[::-1]:
    try:
      gff_non_tox.pop(n)
    except:
      print(n)
      sys.exit()

  # Write (.pps) 
  write_gff_list(base + '.ptx_diff_gff.pps.gff3', not_present)
  write_gff_list(base + '.ptx.pps.gff3', gff_tox)
  write_gff_list(base + '.gff_diff_ptx.pps.gff3', gff_non_tox)
  print(ori)
  print(genome_length)

  for i in range(0, len(gff)):
    f = gff[i]
    toxic = 'FALSE' 
    if i in gff_tox_idx:
      toxic = 'TRUE'
    offset = compute_offset_from_ori(ori, genome_length, int(f[FEATURE_START]))
    if offset < 0 or offset > genome_length:
      print('offset less than 0 or greater than genome length')
      sys.exit()
    print('\t'.join(f + [toxic, str(offset)]))
    

if __name__ == '__main__':
  main()
