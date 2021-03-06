# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com
import subprocess
import json
import os
import timeit
import re
import collections
import sys
import psutil
import time

__version__ = 0.1

ILLEGAL_REGEX_1 = r'{(\n,)+'
ILLEGAL_REGEX_2 = r'([0-9]),\n}'
QUERY = b'Q'
TERM = b'T'

def has_handle(fpath):
  for proc in psutil.process_iter():
    try:
      for item in proc.open_files():
        if fpath == item.path:
          return True
    except Exception:
      pass
  return False

class QueryResult(collections.UserDict):
  """Stores the Mantis json object as a dictionary, and adds to the
     dictionary a query key value pair. """
  def __init__(self, query_dict, query_string):
    self.data = query_dict 
    self['query'] = query_string

  
 
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

    # Create query file (empty) to pass mantis is_file checks
    open(query_file,'w').close()
    cmd = f'{self.mantis_exec} query -1 -j -p {self.mantis_ds} -o {self.result_file} {self.query_file}'
    self.mantis_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
    result = self.mantis_proc.stdout.readline()
    print(result)
    

  def parse_result_file(self):
    """Parse Mantis output json into a json object. 
    
       Mantis' json output is syntactically illegal with respect to 
       python's json library, specifically, when the query is not hit, 
       the database list under ['res'] has illegal synatix. The illegal
       synatix is replaced with appropriate synatix.  """
    q_results = str()
    with open(self.result_file, 'r') as f:
      q_results = f.read()
      q_results = re.sub(ILLEGAL_REGEX_1, '{', q_results)
      q_results = re.sub(ILLEGAL_REGEX_2, r'\1}', q_results)
    try:
      return json.loads(q_results)
    except:
      with open(self.result_file, 'r') as f:
        print('file is empty')
        print(f.read())
        self.terminate()
        sys.exit()

  def query(self, q_list):
    """Queries query list q_list against a Mantis data structure and returns
       an object containing the query results."""

    # Write queries to query file. Use random string for query filename
    # in case multiple jobs run in parallel
    with open(self.query_file, 'w') as f:
      f.write('\n'.join(q_list))

    # Run mantis query
    start = timeit.timeit()
    self.mantis_proc.stdin.write(QUERY)
    self.mantis_proc.stdin.flush()

    while os.path.isfile(self.result_file) == False or has_handle(self.result_file) == True:
      pass
    time.sleep(0.1)

    end = timeit.timeit()
    self.mantis_q_time += end - start
    
    # Parse query file into a list of QueryResult objects.
    try:
      return [QueryResult(json_obj, query) for json_obj, query in zip(self.parse_result_file(), q_list)]
    except Exception as e:
      print(e)
      self.terminate()
    
  def terminate(self):
    try:
      self.mantis_proc.stdin.write(TERM)
      self.mantis_proc.stdin.flush()
    except:
      pass


