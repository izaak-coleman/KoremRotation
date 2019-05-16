# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

import sys
import os
import Extension
import iPCR

__version__ = 0.1

ALPHA = {'A','C','G','T'}
N_INPUT_PARAMS = 7

def check_input_validity(input_params):
  if len(input_params) != N_INPUT_PARAMS:
    msg = "Usage: <exe> <p1> <p2> <path/to/mantis/executable> " +\
          "<path/to/mantis/data/datastructure> <max p1 mismatches> "+\
          "<max p2 mismatches>"
    raise Exception(msg)
  _, p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch = input_params
  try:
    max_p1_mismatch, max_p2_mismatch = int(max_p1_mismatch), int(max_p2_mismatch)
  except ValueError as e:
    print("Both the max p1 an p2 mismatch values must be integers")
    sys.exit()
  if os.path.isfile(mantis_exec) == False or os.path.isdir(mantis_ds) == False:
    raise Exception("Either the supplied mantis executable path or mantis data structure path does not exist")
  if set(p1.upper()) != ALPHA or set(p2.upper()) != ALPHA:
    raise Exception(f"Either probe p1: {p1} or p2: {p2} contains non-ATCG character, which is invalid.")
  else:
    p1, p2  = p1.upper(), p2.upper()
  return p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch

def main():
'''Assembles the genome between a pair of probes p1, p2 from a set of sequencing databases, returning
   any non-SNP variant sequences. '''
  # Check input validity
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch = check_input_validity(sys.argv)

  # Run iPCR module, which will return a list of Extension objects.
  extensions = iPCR.run(p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch)












