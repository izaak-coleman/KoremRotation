
# email: izaak.coleman1@gmail.com
# -*- coding: utf-8 -*-
__version__ = 0.1

import QueryMantis
import random
import string
import os
import sys
import DBG

ALPHA = 'ACGT'
DUMMY_QUERY = 'A'*32
KEY_ERR = -1


class iPCR:
  """Constructs a De Bruijn Graph between two probes p1 and p2."""

  def __init__(self, p1=str(), p2=str(), mantis_exec=str(), 
               mantis_ds=str(), max_p1_mismatch=int(), max_p2_mismatch=int(),
               k=int(), max_extension=int()):
    self.p1, self.p2 = p1, p2
    self.max_p1_mismatch, self.max_p2_mismatch = max_p1_mismatch, max_p2_mismatch
    self.mantis_exec = mantis_exec
    self.mantis_ds = mantis_ds
    self.k = k
    self.max_extension = max_extension
    self.db_dict = self.set_db_dict(mantis_ds)

  def set_db_dict(self, mantis_ds):
    lst_file = mantis_ds + 'sampleid.lst'
    with open(lst_file) as f:
      dbs = [l.strip().split()[1] for l in f]
    return {db:idx for idx, db in enumerate(dbs)}

  def add_probe_to_dbg(self, dbg, query_results, db_dict):
    """Constructs a De Brujin graph from the queried p1* or p2* probes that
       exactly match any of the databases."""
    for q_res in query_results:
      if not self.exact_match(q_res):
        continue
      # Otherwise, add the exact matching probe to the dbg
      dbs = set([db_dict[key] for key in q_res['res'].keys()])
      probe = q_res['query']
      edges = [probe[i:i+(dbg.k+1)] for i in range(0, len(probe) - dbg.k)] # (k+1) not k because edge not node
      for edge in edges:
        dbg.add_edge(edge, dbs)
  
  def exact_match(self, q_res):
    """Returns True if query exactly matches at least one database."""
    return any([True for num_kmers_in_db in q_res['res'].values() if num_kmers_in_db == q_res['num_kmers']])
        
  def update_dbg(self, query_results, p2, max_p2_mismatch, dbg, db_dict):
    """Updates the De Bruijn graph by adding a set of edges that exactly
       match at least one database in the supplied Mantis data structure."""
    # Each queried edge that exactly matches at least one database will
    # be stored and returned (for the next update round).
    exact_matching_edges = list()
   
    # Iterate through queried edges and add exact matching edges to dbg
    for q_res in query_results:
      if not self.exact_match(q_res):
        continue
  
      # Determine whether the edge has made a loop. 
      # Loop detection reduces to determining whether the suffix
      # node of the edge is already present in the dbg. 
      edge = q_res['query']
      suffix_previously_present = False
      if edge[1:] in dbg.nodes:
        suffix_previously_present = True
  
      # Add edge (and nodes)
      dbs = [self.db_dict[db] for db in q_res['res'].keys()]
      dbg.add_edge(q_res['query'], dbs)
  
      # If the suffix was already present, 
      # the dbg path(s) constructed from the edge constructed from
      # suffix + {A, T, C, G} will have already been or is currently
      # being added to the dbg. 
      # Hence, it is not neccessary to include this edge in the
      # list of exact matching edges for the next update round. 
      # Note: This logic avoids infinite query/update interations
      # around cycles in the dbg. 
      if suffix_previously_present:
        continue
      # Otherwise, add the edge to exact matching edges
      exact_matching_edges.append(edge)
  
    return exact_matching_edges
  
  def generate_filenames(self):
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
  
  def build_probe_list(self, edit_dist, probe):
    """Builds a list of p1* probes, of edit distance <= edit_dist from probe."""
    all_probes = set([probe])
    if edit_dist > 0:
      probe_list = [probe[:i] + ALPHA[j] + probe[i+1:] for i in range(0, len(probe)) for j in range(0, 4)]
      for p in probe_list:
        all_probes = all_probes.union(self.build_probe_list(edit_dist-1, p))
    return sorted(list(all_probes))
  
  def run(self):
    """Reconstructs a coloured De Bruijn graph from a list of sequence 
       databases between two probes, p1 and p2 and their close variants (p1*, p2*)."""
    # Initialize QueryMantis for querying
    query_file, result_file = self.generate_filenames()
    qm = QueryMantis.QueryMantis(self.mantis_exec, self.mantis_ds, query_file, result_file)
  
    # Initialize De Bruijn graph by adding all edges from p1* and p2* probes
    # that match p1 and p2 probes within an edit distance of max_p1_mismatch, 
    # and max_p2_mismatch respectively. 
    dbg = DBG.DBG(self.k)
    # Add p1* probes 
    p1_probe_list = self.build_probe_list(self.max_p1_mismatch, self.p1)
    p1_query_results = qm.query(p1_probe_list)
    self.add_probe_to_dbg(dbg, p1_query_results, self.db_dict)
    # Add p2* probes
    p2_probe_list = self.build_probe_list(self.max_p2_mismatch, self.p2)
    p2_query_results = qm.query(p2_probe_list)
    self.add_probe_to_dbg(dbg, p2_query_results, self.db_dict)
  
    # Generate the set of edges from which the Dr Bruijn graph
    # edges will begin (the edges at the end of exact matching
    # p1* probes)
    edges = set()
    for q_res in p1_query_results:
       if self.exact_match(q_res):
       # If some database exactly matches the queried probe, then add its last kmer to extensions.
         query_str = q_res['query']
         edges.add(query_str[-(self.k+1):])
  
    # Begin De Brujin graph construction.
    while len(edges) > 0:
      print(edges)
      print(dbg.nodes)
      print()
      print()
      print()
      print()
      print()
      # Generate new set of edges to query. For each exact matching edge of the
      # previous iteration, generate four new edges of form edge[1:] + {A, T, C, G}. 
      edges = [e[1:] + base for e in edges for base in ALPHA]
      edges = self.update_dbg(qm.query(edges), self.p2, self.max_p2_mismatch, dbg, self.db_dict)
  
    # Construction complete.
    return dbg

#def extensions_incomplete(extensions):
# return any([e.extending for e in extensions])
#
#
#def chunk_list(l, n):
#    """Yield chunks of size n from a list."""
#    for i in range(0, len(l), n):
#        yield l[i:i + n]
#
#def initialize_extensions(q_results, probe_list):
#  """Parses the result from the probe query into a list containing the first Extension object."""
#
#  # List of Extensions, each of which begins with a probe p1*.
#  initial_extensions = list()
#  for q_result, probe in zip(q_results, probe_list):
#  # Initialize a potential extension from the probe.
#    ext = extn.Extension(probe)
#    for db, num_kmers_hit_in_db in q_result['res'].items():
#      if num_kmers_hit_in_db == q_result['num_kmers']:
#      # If the current probe exactly matched database db, then
#      # add the database db to the extension's set of databases.
#        ext.databases.add(db)
#    
#    if len(ext.databases) == 0:
#    # If the current probe failed to exactly match any database
#    # the extension is invalid, so continue without appending
#    # to initial_extensions.
#      continue
#    # Otherwise, append.
#    initial_extensions.append(ext)
#  return initial_extensions

#def terminate_extension(extension, probe, mismatch_threshold, direction):
#  """Check whether an extension should be terminated due to it reaching the maximum 
#     allowed extension or its last |probe| character matching probe with some allowed 
#     mismatch threshold. Returns True if termination required."""
#  if len(extension) > max_extension:
#    return True
#  mismatches = 0
#  suffix = extension.extension[-len(probe):]
#  for a, b in zip(suffix, probe):
#    if a != b:
#      mismatches += 1
#  if mismatches <= mismatch_threshold:
#  # Then we have hit a probe*, so terminate extension
#    if direction == extn.Forward:
#      extension.hit_p2 = True
#    else:
#      extension.hit_p1 = True
#    return True
#  return False # keep truckin'

#def update_extensions(extensions, q_results, probe, mismatch_threshold, direction):
#  """Updates the list of extensions according the current round of querying.
#
#     Function:
#       - Increments extensions by a single base according to current round of querying.
#       - Adds a new Extension duplicate to list if a branching evernt occurs.
#       - Updates each Extension with the databases it currently exactly matches
#       - Terminates extensions if max extension reached, or a probe* reached, or no
#         database hits were found for the current round of extension.  
#     """
#  for i, q_result in enumerate(chunk_list(q_results, 4)):
#    # Determine which databases listed in extensions[i] failed to 
#    # exactly match this round of queries. 
#    unhit_dbs = set()
#    for db in extensions[i].databases:
#      db_hit = [True for query in q_result if 
#                db in query['res'] and 
#                query['res'][db] == query['num_kmers']]
#      if not any(db_hit):
#        unhit_dbs.add(db)
#
#    if len(unhit_dbs) == len(extensions[i].databases):
#    # Then all databases did not contain an exact match for any query
#    # Accordingly, terminate extensions[i]
#      print(f'extensions {i} failed to match any databases')
#      extensions[i].extending = False
#      continue
#
#    if len(unhit_dbs) !=  0:
#    # Then one or more (but not all) database did not contain an exact 
#    # match for any query. Accordingly, duplicate extensions[i], set its
#    # database list to the extensions that did not exact match any query
#    # and terminate it. 
#      ext_duplicate = copy(extensions[i])
#      ext_duplicated.databases = set([db for db in extensions[i].databases if db in unhit_dbs])
#      ext_duplicated.extending = False
#      print(f'the following databases failed to match extension {i}')
#      print(ext_duplicated.databases)
#      extensions.append(ext_duplicated)
#
#  # For the databases that had an exact match against one of the queries
#  # extend extensions[i] appropriately.
#
#  # To avoid O(n) pop(), the first extension is in place (directly
#  # modifying extensions[i]). 
#    mutated_inplace = False
#    ext_old = copy(extensions[i])
#    base = -1 # Keeps track of which base must be added to extension
#    for query in q_result:
#      base += 1
#      hit_dbs = set() # Set of databases for which extension.extension + ALPHA[base] was found
#      for db in ext_old.databases:
#        if (db in query['res'] and query['res'][db] == query['num_kmers']):
#          # Database had an exact match against this query.
#          hit_dbs.add(db)
#
#      if len(hit_dbs) == 0: 
#      # No databases exactly matched with the current query, so continue.
#        continue
#      if not mutated_inplace: 
#        extensions[i].databases = hit_dbs
#        extensions[i].extend(ALPHA[base], direction)
#        if terminate_extension(extensions[i], probe, mismatch_threshold, direction):
#          extensions[i].extending = False
#        mutated_inplace = True
#      else:
#        # BRANCH the extension: Two different variants were found during this extension.
#        # Hence, a new extension must be added. A dupilcate of the original extension
#        # is therefore made, extended with the new variant and added to the list of extensions
#        print(ALPHA[base])
#        print(q_result)
#        ext_duplicate = copy(ext_old)
#        ext_duplicate.databases = hit_dbs
#        ext_duplicate.extend(ALPHA[base], direction)
#        if terminate_extension(ext_duplicate, probe, mismatch_threshold, direction):
#          ext_duplicate.extending = False
#        extensions.append(ext_duplicate)
#  return extensions
