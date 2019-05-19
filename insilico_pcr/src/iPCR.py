# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

__version__ = 0.1

import Extension
import QueryMantis

ALPHA = 'ACGT'
DUMMY_QUERY = 'A'*32

# Keys to access q_result of four bases
A, C, G, T = 0, 1, 2, 3
F, R = 0, 1

# set the default length of the maximum extension to 1MBp
max_extension = 1000000

def extensions_incomplete(extensions):
 return any([e.extending for e in extensions])

def set_max_extension(m)
  max_extensions = m

def chunk_list(l, n):
    """Yield chunks of size n from a list."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def update_extensions(extensions, q_results):
  """Updates the list of extensions according the current round of querying.

     Function:
       - Increments extensions by a single base according to current round of querying.
       - Adds a new Extension duplicate to list if a branching evernt occurs.
       - Updates each Extension with the databases it currently exactly matches
       - Terminates extensions if max extension reached, or a p2* reached, or no
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
    # by lack of extendable sequence occured. 
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
        

def run(*args):
  """Recovers a set of exact matching sequences present in a list of sequence databases between two probes, p1 and p2
     and their close variants (p1*, p2*)."""
  # Unpack args, expected order is as follows.
  p1, p2, mantis_exec, mantis_ds, max_p1_mismatch, max_p2_mismatch = args

  # Initialize QueryMantis for querying
  qm = QueryMantis.QueryMantis(mantis_exec, mantis_ds)
  qm.set_exec_path(mantis_exec)
  qm.set_ds_path(mantis_ds)

  # Run the initializing extension, which is just a query search for p1
  extensions = list()
  extensions = update_extensions(extensions, qm.query([p1]))

  # Begin forward extension (recover p1 + sigma^k + p2*)
  while extensions_incomplete(extensions):
    queries = list()
    for e in extensions:
      queries = queries + [e.extension + char if e.extending else DUMMY_QUERY for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries))

  # Begin reverse extension
  for e in extensions:
    e.cut_initializing_probe()
    e.extending = True
  count_down = len(p1)
  while count_down > 0:
    queries = list()
    for e in extensions:
      queries = queries + [char + e.extension for char in ALPHA]
    extensions = update_extensions(extensions, qm.query(queries))
  
  # Extension process complete, return extensions
  return extensions
  

   
