import sys
import os
import OriTerKeggScrape
import DoriC
import Assembly
import pickle

TAXA_OF_INTEREST = 'Bacteria'

def pair_doric_with_assembly(doric, asm):
  # Need to check whether the .[0-9] extension is going to cause an with matching up genomes
  pairs = []
  for d in doric:
    added_pair = False
    for a in asm:
      if added_pair:
        break
      for a_rs, a_len in a.contigs:
        if d.refseq == a_rs:
          if 'NC_001263.1' == d.refseq:
            print(d)
            print(a)
          a.refseq = a_rs
          a.length = a_len
          if 'NC_001263.1' == d.refseq:
            print(d)
            print(a)
          pairs.append((d, a,))
          added_pair = True
          break
  return pairs

    
def main():
  if len(sys.argv) != 3:
    print('Usage: <exe> <doric> <assembly_list>')
    sys.exit()

  doric_path = sys.argv[1]
  assembly_list = sys.argv[2]
  assembly_dir = os.path.dirname(assembly_list)

 # Get all relevant organisms from Kegg
# organisms = get_kegg_organism_list()
# organisms = (Organism(*org.split('\t')) for org in organisms[:100]
#              if (TAXA_OF_INTEREST in Organism(*org.split('\t')).taxa_string))

# Load DoriC, Assemblies and pair them up
  doric_asm_pairs = pair_doric_with_assembly(
                    DoriC.get_entries(doric_path),
                    Assembly.get_assembly_info(assembly_list))
  print(doric_asm_pairs[0])
  d, a = doric_asm_pairs[0]
  if d.refseq == a.refseq:
    print(True)
  else:
    print(False)
  

##for doric, asm in doric_asm_pairs[0]:
#    print(doric)
#    print(asm)

## Filter Organisms for those in DoriC
#  assembly_index = {a.assembly:(d,a) for d,a in doric_asm_pairs}
#  refseq_index = {a.refseq:(d,a) for (d,a) in doric_asm_pairs}
#  doric_present_organisms = []
#  for org in organisms:

#    result = assembly_index.get(org.assembly, None)
#    if result == None:
#      result = refseq_index.get(org.refseq, None) 
#    if result == None:
#      continue
#      # write some error message
#    doric_present_organisms.append(org)
#
## Now pickle dump all the relevant organisms. This will be loaded 
## by the numpy program and we can then extract the data as neccessary
#  with open('doric_present_organisms.pickle','w') as f:
#    rick.dump(doric_present_organisms, f)
#

if __name__ == '__main__':
  main()
