import sys
import os
from subprocess import Popen, PIPE
import KeggNames
import re
kn = KeggNames.keggnames
os.environ["NCBI_API_KEY"] = "f4142f101db95406745385d940b13c37ab07"

def search_ncbi_for_refseq_id(gb_id):
  print('searching ncbi... %s' % gb_id)
  cmd = f"esearch -db nuccore -query {gb_id} | elink -target nuccore -name nuccore_nuccore_gbrs | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Slen Title | sed 's/,.*//' | sort -t $'\t' -k 2,2nr"
  result,err = Popen(cmd, stdout=PIPE,stderr=PIPE,shell=True).communicate()
  print(result.decode('utf-8'))
  result = result.decode('utf-8').strip().split('\n')
  if len(result) != 1:
    return result[0].split('\t')[0]
  return result[0].split('\t')[0]

def main():
  if len(sys.argv) != 3:
    print('Usage: <exe> <doric_assm_matches> <kegg_file_list>')
    sys.exit()

  with open(sys.argv[1]) as f:
    doric_asm = [l.strip().split('\t') for l in f]

  with open(sys.argv[2]) as f:
    kegg_files = [l.strip() for l in f]

  match = lambda a, b : re.findall(rf'{a}\s+(\S.*)\n',b)
  gcf_dict = {elem[1]:elem for elem in doric_asm}
  refseq_dict = {elem[2]:elem for elem in doric_asm}

  hits = []
  position = -1
  try:
    for kegg_file in kegg_files:
      position += 1
      with open(kegg_file) as f:
        kegg = f.read()
      genome_id = re.findall(rf'{kn.ENTRY}\s+(T[0-9]+)\s', kegg)
      if len(genome_id) == 0:
        print(kegg)
        continue
      genome_id = genome_id[0]

      # First match based on assembly
      kegg_asm = match(kn.DATA_SOURCE, kegg)
      if len(kegg_asm) > 0:
        kegg_asm = re.findall(r'Assembly:([\.\w]+)', kegg_asm[0])
        if len(kegg_asm) == 0:
          continue
        kegg_asm = kegg_asm[0]
        hit = gcf_dict.get('GCF_' + kegg_asm[4:], None)
        if hit != None:
          with open('matched.txt', 'a+') as f:
            f.write('\t'.join(hit) + '\t' + genome_id + '\n')
          continue

      # Try matching refseq
      seq = match(kn.SEQUENCE, kegg)
      if len(seq) == 0:
        continue
      seq = seq[0]
      if seq[0:2] == 'RS':
        seq = re.findall(r'RS:(\S+)\s+\(GB:(\S+)\)',seq)
        if len(seq) == 0:
          print(kegg)
          continue
        kegg_refseq, kegg_genbank = seq[0]
      if seq[0:2] == 'GB':
        kegg_genbank = re.findall(r'GB:(\S+)',seq)[0]
      kegg_refseq = search_ncbi_for_refseq_id(kegg_genbank)
      hit = refseq_dict.get(kegg_refseq, None)
      if hit == None:
        continue
      with open('matched.txt', 'a+') as f:
        f.write('\t'.join(hit) + '\t' + genome_id + '\n')
  except:
    with open('remaning_keggs.txt','w') as f:
      f.write('\n'.join(kegg_files[position:]))
if __name__ == '__main__':
  main()

