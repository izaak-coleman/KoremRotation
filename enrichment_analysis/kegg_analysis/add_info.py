import os
import sys
import DoriC
import Assembly

def main():
  if len(sys.argv) != 4:
    print("Usage: <exe> <doric_asm_matches> <tubic_bacteria> <assembly_file_list>")
    sys.exit()

  doric = DoriC.get_entries(sys.argv[2])
  assemblies = Assembly.get_entries(sys.argv[3])
  with open(sys.argv[1]) as f:
    data = [l.strip().split('\t') for l in f]
    doric_idx = {d.id:d for d in doric}
    asm_idx   = {a.assembly:a for a in assemblies}

  new_data = list()
  for elem in data:
    ori_pos = doric_idx[elem[0]].position
    asm_filename = asm_idx[elem[1]].filename
    genome_length = -1
    if elem[3] == 'ASM':
      genome_length  =  asm_idx[elem[1]].contigs[0][1]
    if elem[3] == 'REFSEQ':
      for refseq, length in asm_idx[elem[1]].contigs:
        if refseq == elem[2]:
          genome_length = length
          break
      if genome_length == -1:
        print('ERROR')
        sys.exit()
    if genome_length == -1:
      print('ERROR')
      sys.exit()
    new_data.append(elem + [str(ori_pos), asm_filename, str(genome_length)])
  print('\t'.join(['ori_id','asm_id','refseq','hit','gn:kegg','ori_pos','fname','genome_len']))
  print('\n'.join(['\t'.join(e) for e in new_data]))

if __name__ == '__main__':
  main()
