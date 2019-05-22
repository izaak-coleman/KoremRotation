# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com
import subprocess
import json
import os
import timeit

__version__ = 0.1
class QueryMantis:
  """Queries a Mantis data structure given a set of queries and parses into an
     object-based format. """

  def __init__(self, mantis_exec = str(), mantis_ds = str(), query_file = str(), result_file = str()):
    """Set the path to the the mantis executable and mantis data structure."""
    self.mantis_exec = mantis_exec
    self.mantis_ds = mantis_ds
    self.query_file = query_file
    self.result_file = result_file
    self.mantis_q_time = 0.0
    if os.path.isfile(mantis_exec) == False or os.path.isdir(mantis_ds) == False:
      raise Exception("Either the supplied mantis executable path or mantis data structure path does not exist")
    self.complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
   

  def query(self, q_list):
    """Queries query list q_list against a Mantis data structure and returns
       an object containing the query results."""
    # DEBUG
    for q in q_list:
      print(q)

    # Write queries to query file. Use random string for query filename
    # in case multiple jobs run in parallel
    with open(self.query_file, 'w') as f:
      f.write('\n'.join(q_list))

    # Run mantis query
    cmd = f'{self.mantis_exec} query -1 -j -p {self.mantis_ds} -o {self.result_file} {self.query_file}'
    start = timeit.timeit()
    result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    end = timeit.timeit()
    self.mantis_q_time += end - start
#    print(result.stdout.decode('utf-8'))
    
    # Parse query file into list of json object
    with open(self.result_file, 'r') as f:
      q_results = f.read()
    return json.loads(q_results)


