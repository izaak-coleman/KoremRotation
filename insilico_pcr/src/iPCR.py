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

def initialize_extensions(q_results, probe_list):
  """Parses the result from the probe query into a list containing the first Extension object."""

  # List of Extensions, each of which begins with a probe p1*.
  initial_extensions = list()
  for q_result, probe in zip(q_results, probe_list):
  # Initialize a potential extension from the probe.
    ext = extn.Extension(probe)
    for db, num_kmers_hit_in_db in q_result['res'].items():
      if num_kmers_hit_in_db == q_result['num_kmers']:
      # If the current probe exactly matched database db, then
      # add the database db to the extension's set of databases.
        ext.databases.add(db)
    
    if len(ext.databases) == 0:
    # If the current probe failed to exactly match any database
    # the extension is invalid, so continue without appending
    # to initial_extensions.
      continue
    # Otherwise, append.
    initial_extensions.append(ext)
  return initial_extensions

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
    for db in extensions[i].databases:
      db_hit = [True for query in q_result if 
                db in query['res'] and 
                query['res'][db] == query['num_kmers']]
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
      ext_duplicated.databases = set([db for db in extensions[i].databases if db in unhit_dbs])
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
      hit_dbs = set() # Set of databases for which extension.extension + ALPHA[base] was found
      for db in ext_old.databases:
        if (db in query['res'] and query['res'][db] == query['num_kmers']):
          # Database had an exact match against this query.
          hit_dbs.add(db)

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

def build_probe_list(edit_dist, probe):
  """Builds a list of p1* probes, of edit distance <= edit_dist from probe."""
  all_probes = set([probe])
  if edit_dist > 0:
    probe_list = [probe[:i] + ALPHA[j] + probe[i+1:] for i in range(0, len(probe)) for j in range(0, 4)]
    for p in probe_list:
      all_probes = all_probes.union(build_probe_list(edit_dist-1, p))
  return sorted(list(all_probes))

def run(*args):
  """Recovers a set of exact matching sequences present in a list of sequence 
     databases between two probes, p1 and p2 and their close variants (p1*, p2*)."""
  # Unpack args, expected order is as follows.
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch = args

  # Initialize QueryMantis for querying
  query_file, result_file = generate_filenames()
  qm = QueryMantis.QueryMantis(mantis_exec, mantis_ds, query_file, result_file)

  # Genrate the initial extensions: List of Extension objects 
  # begining with some p1* that is an edit distance <= max_p1_mismatch
  # from p1.
  p1_probe_list = build_probe_list(max_p1_mismatch, p1)
  extensions = initialize_extensions(qm.query(p1_probe_list), p1_probe_list)

  # Begin forward extension (recover p1 + sigma^k + p2*)
  while extensions_incomplete(extensions):
    queries = list()
    for e in extensions:
      queries = queries + [e.extension + char if e.extending else DUMMY_QUERY for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries), p2, max_p2_mismatch, extn.Forward)
    for i, e in enumerate(extensions):
      print(e.databases)
      print(f'extension {i}: {e.extension}')
  print(p1)
  print(f'Mantis query time: {qm.mantis_q_time}')

  # Extension process complete, return extensions
  return extensions
  
