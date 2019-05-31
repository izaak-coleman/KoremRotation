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
    _, p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, db_file = input_params
  elif len(input_params) == N_INPUT_PARAMS + OPTIONAL: 
    _, p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, db_file, max_extension = input_params
  else:
    msg = "Usage: <exe> <p1> <p2> <path/to/mantis/executable> " +\
          "<path/to/mantis/data/datastructure> <max p1 mismatches> "+\
          "<max p2 mismatches> <max_extension (optional)"
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

  # Check .lst file exists
  if os.path.isfile(db_file) == False:
    raise Exception(f'{db_file} does not exist. Please check filename / path is valid')

  return p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, db_file, max_extension

def main():
  '''Assembles the genome between a pair of probes p1, p2 from a 
     set of sequencing databases, returning any non-SNP variant sequences. '''
  # Check input validity
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch, db_file, max_extension = check_input_validity(sys.argv)
  
  # Run iPCR module, which will return a list of Extension objects.
  iPCR.set_max_extension(max_extension)
  iPCR.set_db_dict(db_file)
  start = timeit.timeit()
  dbg = iPCR.run(p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch)
  dbg.render('testing_dbg.pdf')
  end = timeit.timeit()
  print(end - start)
if __name__ == '__main__':
  main()
