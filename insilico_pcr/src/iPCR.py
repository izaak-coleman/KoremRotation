# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com
# -*- coding: utf-8 -*-
__version__ = 0.1

import Extension as extn
import QueryMantis
import random
import string
import os

ALPHA = 'ACGT'
DUMMY_QUERY = 'A'*32
KEY_ERR = -1

# Keys to access q_result of four bases
FWD, REV = 0, 1

# set the default length of the maximum extension to 1MBp
max_extension = 1000000

def extensions_incomplete(extensions):
 return any([e.extending for e in extensions])

def set_max_extension(m):
  if m > 0:
    max_extension = m

def chunk_list(l, n):
    """Yield chunks of size n from a list."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def initialize_extensions(q_result, probe):
  """Parses the result from the probe query into a list containing the first Extension object."""
  init_ext = extn.Extension(probe) 

  # Extract the forward and reverse databasenames from q_result
  fwd, rev = q_result.pop()
  dbs = set([db for db in fwd['res'].keys()])
  dbs.union( set([db for db in rev['res'].keys()]) )

  # Add dbs to extension
  for db in dbs: 
    # Get the number of uniq kmers in query and found in mantis
    query_fwd_kmer_count = fwd['num_kmers']
    query_rev_kmer_count = rev['num_kmers']
    mantis_fwd_kmer_count = fwd['res'].get(db, KEY_ERR)
    mantis_rev_kmer_count = rev['res'].get(db, KEY_ERR)

    # Compute U_q(Q) - U_m(Q) invariants for forward and reverse
    fwd_invariant, rev_invariant = -1, -1
    if mantis_fwd_kmer_count != -1:
      fwd_invariant = query_fwd_kmer_count - mantis_fwd_kmer_count
    if mantis_rev_kmer_count != -1:
      rev_invariant = query_rev_kmer_count - mantis_rev_kmer_count
    init_ext.databases[db] = (fwd_invariant, rev_invariant)

  return [init_ext]

def update_extensions(extensions, q_results, probe, mismatch_threshold, direction):
  """Updates the list of extensions according the current round of querying.

     Function:
       - Increments extensions by a single base according to current round of querying.
       - Adds a new Extension duplicate to list if a branching evernt occurs.
       - Updates each Extension with the databases it currently exactly matches
       - Terminates extensions if max extension reached, or a probe* reached, or no
         database hits were found for the current round of extension.  
     """
  for i, q_result in enumerate(chunk_list(q_results, 4)):
    if not extensions[i].extending: 
      continue
  # If none of the eight queries ( {A,T,C,G} x {Fwd, RevComp} ) hit a database
  # terminate the extension and move to the next extension
    if sum([len(fwd['res']) + len(rev['res']) for fwd, rev in q_result]) == 0:
      extensions[i].extending = False
      continue

  # If, for some databases in the extension's set of databases prior to the
  # query, none of the current queries hit the database, then duplicated
  # the extension, set its database list to the databases that satisfy the
  # the condition, and then terminate it. (For these databases, termination
  # by lack of extendable sequence occured). 
    unhit_dbs = set()
    for db in extensions[i].databases.keys():
      db_hit = [True for fwd, rev in q_result 
                     if db in fwd['res'].keys() or db in rev['res'].keys()]
      if not any(db_hit):
        unhit_dbs.add(db)
    if len(unhit_dbs) !=  0:
      ext_duplicate = extensions[i]
      ext_duplicated.databases = {k:v for k, v in extensions[i].items() if k in unhit_dbs}
      ext_duplicated.extending = False
      extensions.append(ext_duplicated)

  # Begin extending the extension by a single base. 

  # To avoid increasing complexity, perform the first extension inplace,
  # directly modifying extensions[i]. Make proceeding modifications 
  # from an original copy of extensions[i] (ext_old)
    mutated_inplace = False
    ext_old = extensions[i] 
    base = -1 # Keeps track of which base must be added to extension
    for fwd, rev in q_result:
      base += 1
      hit_dbs = dict() # List of databases for which extension.extension + ALPHA[base] was found
      for db, invariant in ext_old.databases.items():
        if ( (db in fwd['res'].keys() and (fwd['res'][db] - fwd['num_kmers'])  == invariant[FWD]) or
             (db in rev['res'].keys() and (rev['res'][db] - rev['num_kmers'])  == invariant[REV]) ):
        # If the query "extension.extension + ALPHA[base]" hit database db, 
        # where p1 + sigma^k = extension.extension + ALPHA[base] and 
        # sigma^k is an exact match, record the database in hit_dbs
          hit_dbs[db] = invariant

      # Then the search failed to hit a db and maintain the invariant. Hence, an 
      # exact match was not found, so continue.
      if len(hit_dbs) == 0: 
        continue
      if not mutated_inplace: 
        # Then directly extend extensions[i]
        extensions[i].databases = hit_dbs
        extensions[i].extend(ALPHA[base], direction)
        if terminate_extension(extensions[i], probe, mismatch_threshold, direction):
          extensions[i].extending = False
      else:
        # BRANCH the extension: Two different variants were found during this extension.
        # Hence, a new extension must be added. A dupilcate of the original extension
        # is therefore made, extended with the new variant and added to the list of extensions
        ext_duplicate = ext_old
        ext_duplicate.databases = hit_dbs
        ext_duplicate.extend(ALPHA[base], direction)
        if terminate_extension(ext_duplicate, probe, mismatch_threshold, direction):
          ext_duplicate.extending = False
        extensions.append(ext_duplicate)

  return extensions
     
def terminate_extension(extension, probe, mismatch_threshold, direction):
  """Check whether an extension should be terminated due to it reaching the maximum 
     allowed extension or its last |probe| character matching probe with some allowed 
     mismatch threshold. Returns True if termination required."""
  if len(extension) > max_extension:
    return True
  mismatches = 0
  suffix = extension.extension[-len(probe):]
  for a, b in zip(suffix, probe):
    if a != b:
      mismatches += 1
  if mismatches < mismatch_threshold:
  # Then we have hit a probe*, so terminate extension
    if direction == extn.Forward:
      extension.hit_p2 = True
    else:
      extension.hit_p1 = True
    return True
  return False # keep truckin'

def generate_filenames():
  """Generates a unique pair of filenames for mantis query file and mantis results file
     used by QueryMantis."""
  alphanums = string.ascii_uppercase + string.digits
  query_file = ''.join([random.choice(alphanums) for i in range(0, 15)]) + '.query_file'
  result_file = ''.join([random.choice(alphanums) for i in range(0, 15)]) + '.mantis.json'
  while (query_file[:15] == result_file[:15] or 
         os.path.isfile(query_file) or 
         os.path.isfile(result_file)):
    query_file = ''.join([random.choice(alphanums) for i in range(0, 15)]) + '.query_file'
    result_file = ''.join([random.choice(alphanums) for i in range(0, 15)]) + '.mantis.json'
  return query_file, result_file

  
def run(*args):
  """Recovers a set of exact matching sequences present in a list of sequence 
     databases between two probes, p1 and p2 and their close variants (p1*, p2*)."""
  # Unpack args, expected order is as follows.
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch = args

  # Initialize QueryMantis for querying
  query_file, result_file = generate_filenames()
  qm = QueryMantis.QueryMantis(mantis_exec, mantis_ds, query_file, result_file)

  # Run the initializing extension, which is just a query search for p1
  extensions = initialize_extensions(qm.query([p1]), p1)

  # Begin forward extension (recover p1 + sigma^k + p2*)
  rounds = 15
  while extensions_incomplete(extensions):
    if rounds == 0: break
    queries = list()
    for e in extensions:
      queries = queries + [e.extension + char if e.extending else DUMMY_QUERY for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries), p2, max_p2_mismatch, extn.Forward)
    rounds -= 1
  print(p1)
  print(extensions[0].extension)
  for k, v in extensions[0].databases.items():
    print(f'db {k}, invariants {v}')
#
#  # Begin backward extension
#  for e in extensions:
#    e.cut_initializing_probe()
#    e.extending = True
#  count_down = len(p1)
#  while count_down > 0:
#    queries = list()
#    for e in extensions:
#      queries = queries + [char + e.extension if e.extending else DUMMY_QUERY for char in ALPHA]
#    extensions = update_extensions(extensions, qm.query(queries), p1, max_p1_mismatch, extn.Backward)
#  
#  # Extension process complete, return extensions
  return extensions
  
