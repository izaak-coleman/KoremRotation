import sys
import os
import OriTerKeggScrape
import DoriC
import cPickle as rick

TAXA_OF_INTEREST = 'Bacteria'

def pair_doric_with_assembly(doric, asm):
  # TESTING. 
  # Need to check whether the .[0-9] extension is going to cause an with matching up genomes
  index = {a.refseq: a for a in asm}
  return [(d, index[d.refseq]) for d in doric]

    
def main():

 # Get all relevant organisms from Kegg
 organisms = get_kegg_organism_list()
 organisms = (Organism(*org.split('\t')) for org in organism 
              if (TAXA_OF_INTEREST in Organism(*org.split('\t')).taxa_string))

# Load DoriC, Assemblies and pair them up
 doric_asm_pairs = pair_doric_with_assembly( \
     DoriC.get_entries('doric_file'), \
     Assembly.get_assembly_info('Assembly_file') \
     )

# Filter Organisms for those in DoriC
  assembly_index = {a.assembly:(d,a) for d,a in doric_asm_pairs}
  refseq_index = {a.refseq:(d,a) for (d,a) in doric_asm_pairs}
  doric_present_organisms = []
  for org in organisms:
    result = assembly_index.get(org.assembly, None)
    if result == None:
      result = refseq_index.get(org.refseq, None) 
    if result == None:
      continue
      # write some error message
    doric_present_organisms.append(org)

# Now pickle dump all the relevant organisms. This will be loaded 
# by the numpy program and we can then extract the data as neccessary
  with open('doric_present_organisms.pickle','w') as f:
    rick.dump(doric_present_organisms, f)

