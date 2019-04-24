import sys
import DoriC
import Assembly

REFSEQ_MATCH = 0
ASSM_MATCH = 1
PREFIX_MATCH = 2

def pair_doric_with_assembly(doric, assemblies):
  # Need to check whether the .[0-9] extension is going to cause an with matching up genomes
  pairs = dict()
  doric_refseq_index = {d.refseq:d for d in doric}
  doric_refseq_prefix_index = {d.refseq[:-2]:d for d in doric}
  doric_asm_index = {d.assembly: d for d in doric}

  # Match doric and assemblies based on refseq
  for asm in assemblies:
    for asm_refseq, length in asm.contigs:
      match = doric_refseq_index.get(asm_refseq, None)
      if match == None:
        continue
      asm.refseq = asm_refseq
      asm.length = length
      if match.id not in pairs:
        pairs[match.id] = (match, asm, REFSEQ_MATCH)

  # Match doric and assemblies based on assembly
  for asm in assemblies:
    match = doric_asm_index.get(asm.assembly, None)
    if match == None:
      continue
    # Consider match if a single contig is present
    if len(asm.contigs) == 1:
      asm.refseq = asm.contigs[0][0]
      asm.length = asm.contigs[0][1]
      if match.id not in pairs:
        pairs[match.id] = (match, asm, ASSM_MATCH)
    # If more than one contig present, match based on refseq prefix
    else:
      contig = None
      for contig in asm.contigs:
        contig_match = doric_refseq_prefix_index.get(contig[0][:-2], None)
        if contig_match != None:
          if contig_match.id not in pairs:
            pairs[contig_match.id] = (contig_match, asm, PREFIX_MATCH)
          break
  return pairs


def main():

  if len(sys.argv) != 3:
    print('Usage: <exe> <doric> <assembly file list>')

  doric = DoriC.get_entries(sys.argv[1])
  assemblies = Assembly.get_entries(sys.argv[2])
  pairs = pair_doric_with_assembly(doric, assemblies)
  counts = {REFSEQ_MATCH:0, ASSM_MATCH:0, PREFIX_MATCH:0}

  get_match = lambda x: ['REFSEQ','ASM','PREFIX'][x]
  for d, asm, match_type  in pairs.values():
    print('\t'.join([d.id, d.assembly, asm.refseq, get_match(match_type)]))
    

if __name__ == '__main__':
  main()
