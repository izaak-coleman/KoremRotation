"""Insilico PCR version 0.1. This will not directly handle sequencing errors. 
   Nor will it handle coverage gaps. """

__version__ = '0.1'
__author__ = 'Izaak Coleman'

QUERY_FILE = 'ipcr.query.fa'

base_switch = {0:'A', 1:'C', 2:'G', 3:'T'}

def iPCR(p1, p2, max_iter):
  """Reconstructs the sequences between probes p1, p2."""
  # recons: k=id, v = Reconstruction object
  recons= dict()
  # set the number of current number of sequences undergoing reconstructions to 1
  num_recons = 1

  # Begin the reconstruction process
  # Initialise p1 to be the first reconstruction (this will branch if variants)
  recons[num_recons] = Reconstruction(start = p1, end = p2, terminated = False)
  num_recons += 1
  it = 0
  while(it < max_iter and still_reconstructing(sequences.values())):
    query_id = 0
    query_file =  open(QUERY_FILE, 'w')
    # Keep track of which queries relate to which reconstructions with recon_query_dict
    recon_query_dict = {recon_id:list() for recon_id in recons.keys() if recons[recon_id].terminated == False}
    for recon_id, recon in recons.items():
      if recon.terminated:
        continue
      # otherwise, write this reconstruction's next set of queries to queryfile
      query = str(recons[recon_id])
      for base in base_switch.values():
        query_file.write(query + base + '\n')
        recon_query_dict[recon_id].append(query_id)
        query_id += 1
    # Run Mantis
    query_file.close()
    query_results = mantis(QUERY_FILE)
    for recon_id, query_ids in recon_query_dict.items():
      zero_indices = [query_id for query_id in query_ids if query_results[query_id] == 0]
      if len(zero_indices) == 0:
        # BETTER HANDLE FOR THIS
        continue
      # For the first query that exactly matched, extend the existing reconstruction
      # by the corresponding base
      recons[recon_id].append(base_switch[zero_indices[0] % 4])
      zero_indices.pop(0)
      # For any remaining queries that exactly matched, consider the remaining
      # matches new variants. Accordingly, 'branch' the reconstructions by making a new
      # reconstruction with a prefix identical to the current reconstruction's sequence
      # but with a the new base as the next character in the reconstuction. 
      num_recons += 1 # making a new reconstruction
      reconst[num_recons] = 
      
