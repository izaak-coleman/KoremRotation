# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

__version__ = 0.1

import Extension
import QueryMantis as qm

ALPHA = 'ACGT'

def run(*args):
  """Recovers a set of exact matching sequences present in a list of sequence databases between two probes, p1 and p2
     and their close variants (p1*, p2*)."""
  # Unpack args, expected order is as follows.
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch = args

  # Initialize QueryMantis for querying
  qm.set_exec_path(mantis_exec)
  qm.set_ds_path(mantis_ds)

  # Run the initializing extension, which is just a query search for p1
  # HANGOVER HELP: The general query protocol will be to 1) send ordered list of queries to qm
  # HANGOVER HELP: which will write queries to query_file in given order. 
  # HANGOVER HELP: qm will then run mantis on query_file and return an easily parseable data structure. (all done by qm.query)
  # HANGOVER HELP: update extensions will then use appropriately update the extensions where neccessary.
  # HANGOVER HELP: note add_queries should add the complement of each query in queries
  # HANGOVER HELP: note further, update extensions is going to handle branching of Extensions, keeping track of databases
  # HANGOVER HELP: and when an extensions has terminated (max extension reached, or reached a p2*)
  extensions = list()
  qm.add_queries([p1])
  extensions = update_extensions(extensions, qm.query())

  # Begin forward extension (recover p1 + sigma^k + p2*)
  while extensions_incomplete(extensions):
  # HANGOVER HELP: extensions_incomplete is going to check if at least one Extension.extending remains true
    queries = list()
    for e in extensions:
      queries = queries + [e.extension + char for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries))

  # Begin reverse extension
  for e in extensions:
    e.cut_initializing_probe()
    e.extending = True
  count_down = len(p1)
  while count_count > 0
    queries = list()
    for e in extensions:
      queries = queries + [char + e.extension for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries))
  
  # Extension process complete, return extensions
  return extensions
  

   
