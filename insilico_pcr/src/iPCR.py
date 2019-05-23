# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com
# -*- coding: utf-8 -*-
__version__ = 0.1

import Extension as extn
import QueryMantis
import random
import string
import os
from copy import copy

ALPHA = 'ACGT'
DUMMY_QUERY = 'A'*32
KEY_ERR = -1


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

  # Extract databasenames from the initial query
  init_query = q_result.pop()
  dbs = [db for db in init_query['res'].keys()]

  # Add dbs to extension
  for db in dbs: 
    # Get the number of uniq kmers in query and found in mantis
    query_kmer_count = init_query['num_kmers']
    mantis_kmer_count = init_query['res'].get(db, KEY_ERR)

    # Compute U_q(Q) - U_m(Q) invariant
    invariant = -1
    if mantis_kmer_count != -1:
      invariant = query_kmer_count - mantis_kmer_count
    init_ext.databases[db] = invariant

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
    # Determine which databases listed in extensions[i] failed to 
    # exactly match this round of queries. 
    unhit_dbs = set()
    for db, invariant in extensions[i].databases.items():
      db_hit = [True  for query in q_result if 
                      db in query['res'].keys() and 
                      (query['res'][db] - query['num_kmers']) == invariant]
      if not any(db_hit):
        unhit_dbs.add(db)

    if len(unhit_dbs) == len(extensions[i].databases):
    # Then all databases did not contain an exact match for any query
    # Accordingly, terminate extensions[i]
      print(f'extensions {i} failed to match any databases')
      extensions[i].extending = False
      continue

    if len(unhit_dbs) !=  0:
    # Then one or more (but not all) database did not contain an exact 
    # match for any query. Accordingly, duplicate extensions[i], set its
    # database list to the extensions that did not exact match any query
    # and terminate it. 
      ext_duplicate = copy(extensions[i])
      ext_duplicated.databases = {k:v for k, v in extensions[i].items() if k in unhit_dbs}
      ext_duplicated.extending = False
      print(f'the following databases failed to match extension {i}')
      print(ext_duplicated.databases)
      extensions.append(ext_duplicated)

  # For the databases that had an exact match against one of the queries
  # extend extensions[i] appropriately.

  # To avoid O(n) pop(), the first extension is in place (directly
  # modifying extensions[i]). 
    mutated_inplace = False
    ext_old = copy(extensions[i])
    base = -1 # Keeps track of which base must be added to extension
    for query in q_result:
      base += 1
      hit_dbs = dict() # List of databases for which extension.extension + ALPHA[base] was found
      for db, invariant in ext_old.databases.items():
        if (db in query['res'].keys() and (query['res'][db] - query['num_kmers']) == invariant):
          # Database had an exact match against this query.
          hit_dbs[db] = invariant

      if len(hit_dbs) == 0: 
      # No databases exactly matched with the current query, so continue.
        continue
      if not mutated_inplace: 
        extensions[i].databases = hit_dbs
        extensions[i].extend(ALPHA[base], direction)
        if terminate_extension(extensions[i], probe, mismatch_threshold, direction):
          extensions[i].extending = False
        mutated_inplace = True
      else:
        # BRANCH the extension: Two different variants were found during this extension.
        # Hence, a new extension must be added. A dupilcate of the original extension
        # is therefore made, extended with the new variant and added to the list of extensions
        print(ALPHA[base])
        print(q_result)
        ext_duplicate = copy(ext_old)
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
  if mismatches <= mismatch_threshold:
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

def generate_probe_list(edit_dist, probe):
  probe_list = [probe[:i] + 'ATCG'[j] + probe[i+1:] for i in range(0, len(probe)) for j in range(0, 4)]
  all_probes = set(probe_list)
  if (edit_dist - 1) > 0:
    for p in probe_list:
      all_probes = all_probes.union(generate_probe_list(edit_dist-1, p))
  return sorted(list(all_probes))

  
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
  while extensions_incomplete(extensions):
    queries = list()
    for e in extensions:
      queries = queries + [e.extension + char if e.extending else DUMMY_QUERY for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries), p2, max_p2_mismatch, extn.Forward)
    for i, e in enumerate(extensions):
      print(e.databases)
      print(f'extension {i}: {e.extension}')
    break
  print(p1)
  print(f'Mantis query time: {qm.mantis_q_time}')
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
  
