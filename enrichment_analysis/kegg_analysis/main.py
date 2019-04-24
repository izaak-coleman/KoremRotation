import sys
import os
import OriTerKeggScrape as ot
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
          a.refseq = a_rs
          a.length = a_len
          a.contigs = ''
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

### Get all relevant organisms from Kegg
  organisms = ot.get_kegg_organism_list()
  organisms = [org for org in organisms if TAXA_OF_INTEREST in org]
  organisms = (ot.Organism(*org.split('\t')) for org in organisms[:100])

# Load DoriC, Assemblies and pair them up
  doric_asm_pairs = pair_doric_with_assembly(
                    DoriC.get_entries(doric_path),
                    Assembly.get_assembly_info(assembly_list))
# Filter Organisms for those in DoriC
  print(len(doric_asm_pairs))
  assembly_index = {a.assembly:(d,a) for d,a in doric_asm_pairs}
  refseq_index = {a.refseq[:-2]:(d,a) for (d,a) in doric_asm_pairs}
  doric_present_organisms = []
  print(assembly_index.items())
  for k,v in list(assembly_index.items())[:10]:
    print(f'{k},{v}')
  for k,v in list(refseq_index.items())[:10]:
    print(f'{k},{v}')

  for org in organisms:
    org_asm = 'GCF_' + org.assembly[4:]
    result = assembly_index.get(org_asm, None)
    if result == None:
      result = refseq_index.get(org.refseq[:-2], None) 
    if result == None:
      continue
    org.doric, org.assembly = result
    print("org.assembly %s" % org.assembly.assembly)
    print("org.doric.ass %s" % org.doric.assembly)
    print("org.refseq.ass %s" % org.refseq.assembly)
    print("org.doric.refseq %s" % org.doric.refseq)
    doric_present_organisms.append(org)


if __name__ == '__main__':
  main()
