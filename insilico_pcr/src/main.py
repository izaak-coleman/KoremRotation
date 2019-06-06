#!/usr/bin/env python
# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

import sys
import Extension
import iPCR
import os
import timeit

__version__ = 0.1

ALPHA = {'A','C','G','T'}
N_INPUT_PARAMS = 8
OPTIONAL = 1

def check_input_validity(input_params):
  # Ensure correct number if input args
  max_extension = -1
  if len(input_params) == N_INPUT_PARAMS:
    _, p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, k = input_params
  elif len(input_params) == N_INPUT_PARAMS + OPTIONAL: 
    _, p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, k, max_extension = input_params
  else:
    msg = "Usage: <exe> <p1> <p2> <path/to/mantis/executable> " +\
          "<path/to/mantis/data/datastructure> <max p1 mismatches> "+\
          "<max p2 mismatches> <k> <max_extension (optional)"
    raise Exception(msg)

  # Correctly set input arg types
  try:
    max_p1_mismatch, max_p2_mismatch = int(max_p1_mismatch), int(max_p2_mismatch)
  except ValueError as e:
    print("Both the max p1 an p2 mismatch values must be integers.")
    sys.exit()
  try:
    max_extension = int(max_extension)
  except ValueError as e:
    print("Max extension must be an integer.")
    sys.exit()
  if set(p1.upper()) != ALPHA or set(p2.upper()) != ALPHA:
    raise Exception(f"Either probe p1: {p1} or p2: {p2} contains non-ATCG character, which is invalid.")
  else:
    p1, p2  = p1.upper(), p2.upper()

  # Check paths to Mantis exe and Mantis db are valid
  if os.path.isfile(mantis_exec) == False or os.path.isdir(mantis_ds) == False:
    raise Exception("Either the supplied mantis executable path or mantis data structure path does not exist. Check paths are valid.")

  # Check k is an integer 
  try:
    k = int(k)
  except ValueError as e:
    print(f'{k} is not an integer.')
    sys.exit()

  if len(p1) < k + 1 or len(p2) < k + 1:
    raise Exception(f'p1 and p2 must have length >= {k + 1}. Current lengths: p1 {len(p1)}, p2 {len(p2)}')
  return p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, k, max_extension

def main():
  '''Assembles the genome between a pair of probes p1, p2 from a 
     set of sequencing databases, returning any non-SNP variant sequences. '''
  # Check input validity
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, k, max_extension = check_input_validity(sys.argv)
  
  # Initialize iPCR instance.
  start = timeit.timeit()
  ipcr = iPCR.iPCR(p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, k, max_extension)
  # Run iPCR. Constructed De Bruijn graph output. 
  dbg = ipcr.run()
  # Write De Bruijn graph.
  dbg.render('testing_dbg.gv')
  dbg.compress()
  dbg.render('testing_dbg.cmp.gv')
  end = timeit.timeit()
  print(end - start)
if __name__ == '__main__':
  main()
