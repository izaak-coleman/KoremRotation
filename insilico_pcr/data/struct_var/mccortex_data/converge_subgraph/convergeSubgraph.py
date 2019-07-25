# Author: Izaak Coleman

import sys
import os
import subprocess
from cortexpy.graph.parser.random_access import RandomAccess

# Cluster submission imports
from time import sleep
from logging_setup import sethandlers
from qp import qp
import logging
from Utils import shell_command
le, lo = str(), str()
memory, threads = 4, 1
MONO = '.mono.ctx'

def set_logfiles(pe, po):
  global le, lo
  le  = pe
  lo = po

def set_memory(m):
  global memory
  memory = m

def set_threads(t):
  global threads
  threads = t

def write_err(data):
  if type(data) == bytes:
    data = data.decode('utf-8')
  with open(le, 'a+') as f:
    f.write(data)

def write_out(data):
  if type(data) == bytes:
    data = data.decode('utf-8')
  with open(lo, 'a+') as f:
    f.write(data)

def parse_input(argv):
  if len(argv) != 10:
    print("Usage: <exe> <mccortex_exe> <seed.fa> <ctx_list> <subgraph_name> <distance> <convergence_ratio> <log_file> <mem> <threads>")
    sys.exit(1)
  mp, seed, ctx_file, subg, d, c, l, m, t = sys.argv[1:]
  print(f'Mcortex exe {mp}',
        f'Seed {seed}',
        f'Subgraph name {subg}',
        f'CTX file list {ctx_file}',
        f'Distance {d}',
        f'Convergence {c}',
        f'Logs {l}',
        f'Memory {m}',
        f'Threads {t}',)
  if not os.path.isfile(seed):
    print("Seed kmer file invalid.")
    sys.exit()

  with open(ctx_file) as f:
    ctx_list = [l.strip() for l in f]
  for ctx in ctx_list:
    if not os.path.isfile(seed):
      print(f'.ctx file {ctx} invalid')
      sys.exit()

  d = int(d)
  c = float(c)
  set_memory(m)
  set_threads(t)
  set_logfiles(l+'.err', l+'.out')
  return mp, seed, ctx_list, subg, d, c

def run_subprocess(cmd):
  r = subprocess.run(cmd.split(), capture_output=True)
  write_out(r.stdout)
  write_err(r.stderr)
  return r.returncode

def get_num_records(ctx):
  """Returns the number of records contained in a .ctx file."""
  graph = RandomAccess(open(ctx, 'rb'))
  return graph.n_records


def converged(prev_n_records, curr_n_records, convergence_ratio):
  """Check to see whether the records have converged"""
  return (float(prev_n_records) / curr_n_records) >= convergence_ratio


def init_monochrome_graph(mccortex_exe, mono_ctx):
  """Constructs an empty monochrome graph (.ctx)."""
  # Make empty fasta.
  run_subprocess("touch init.fa")

  # Make empty monochrome .ctx.
  run_subprocess(mccortex_exe + f" build -k 31 -f -m {memory}G --sample mono --seq2 init.fa:init.fa {mono_ctx}")

  # Clean up.
  run_subprocess(f'rm init.fa')

def full_join_cluster(ctx, mono_ctx, memory):
  cmd = (mccortex_exe + f" join -f -m {memory}G --out {ctx}.merged 0:{mono_ctx} 0:{ctx}")
  shell_command(cmd)


def extract_mono_subgraph(mccortex_exe, mono_ctx, ctx_list, subgraph, seed, dist, init=False, parallel_type='thread'):
  """Merges kmer data from the .ctx's in ctx_list into a single monochrome
     subgraph (defined by seed and dist)."""

  # If not the initial monochrome subgraph extraction, join mono_ctx to
  # each sample-specific .ctx by threading or qsub
  if parallel_type == 'thread' and not init:
    pass

  elif parallel_type == 'qsub' and not init:
    with qp(jobname='FullCTXJoin', qworker='/ifs/scratch/c2b2/korem/ic2465/shared/common/qworker.py', mem_def='32G',
            time_def='1:0:0', remote=True, remote_user='ic2465', remote_pass=pwd,
            remote_cwd='/ifs/scratch/c2b2/korem/ic2465',
            remote_cwd_smb_path='/home/local/ARCS/ic2465/cluster/ifs/scratch/c2b2/korem/ic2465') as q:

      waiton = [q.method(full_join_cluster, (ctx, mono_ctx, memory,)) for ctx in ctx_list]
      # Nothing of interest returned from the full_join_cluster
      for t in waiton:
        q.waitforresult(t)

    ctx_list = [f"{ctx}.merged" for ctx in ctx_list]

  # Extract sample-specific subgraphs
  for ctx in ctx_list:
    run_subprocess(mccortex_exe + f" subgraph -f -m {memory}G -t {threads}" +
                   f" -d {dist} --out {ctx}.{subgraph}.ctx --seq {seed} {ctx}")

  # Join sample specific subgraphs into new mono_ctx
  subgraphs = [f"{ctx}.{subgraph}.ctx" for ctx in ctx_list]
  run_subprocess(mccortex_exe + f" join -f -m {memory}G --out {mono_ctx} {' '.join([f'0:{s}' for s in subgraphs])}")

  # Clean up
  run_subprocess(f"rm {' '.join(subgraphs)}")
  if not init:
    run_subprocess(f"rm {' '.join(ctx_list)}")

def extract_mono_subgraph(mccortex_exe, mono_ctx, ctx_list, seed, dist):
  """Merges kmer data from the .ctx's in ctx_list into a single monochrome
     subgraph (defined by seed and dist)"""

  # Merge each ctx into the monochrome subgraph. 
  for ctx in ctx_list:
    run_subprocess(mccortex_exe + f" join -f -m {memory}G --out {mono_ctx}.merged {mono_ctx} 0:{ctx}")
    run_subprocess(mccortex_exe + f" subgraph -f -m {memory}G -t {threads} -d {dist} --out {mono_ctx} --seq {seed} {mono_ctx}.merged")

  # Clean up. 
  run_subprocess(f'rm {mono_ctx}.merged')

def extract_coloured_subgraphs(mccortex_exe, mono_ctx, subgraph_name, ctx_list):
  """Computes the intersection .ctx for each graph in ctx_list against subgraph."""

  # --intersect {mono_ctx} writes only the kmers in each .ctx that are contained in 
  # mono_ctx but keeps the coverage/edge information in the .ctx
  for ctx in ctx_list:
    run_subprocess(mccortex_exe + f" join -f -m {memory}G  --intersect {mono_ctx} --out {ctx[:-4]}.{subgraph_name} {ctx}")

  # Clean up. 
  #run_subprocess(f'rm {mono_ctx}')

def merge_coloured_subgraphs(mccortex_exe, ctx_list, subgraph_name):
  """Merges single-sample-coloured subgraphs into a single multicoloured .ctx"""
  subgraphs = ' '.join([f'{ctx[:-4]}.{subgraph_name}' for ctx in ctx_list])
  run_subprocess(mccortex_exe + f" join -f -m {memory}G --out {subgraph_name} {subgraphs}")
  # Clean up. 
  #run_subprocess(f'rm {subgraphs}')
  
def main():
  """Build a coloured De Bruijn subgraph from multiple single colour .ctx files. 

     Iteratively builds a coloured DBG subgraph from .ctx files using a single input
     string as the center of the subgraph. Will iterate mccortex31 join commands over between
     the subgraph and each .ctx until the number of kmers added between each round converges. """

  mccortex_exe, seed, ctx_list, subgraph_name, distance, convergence_ratio  = parse_input(sys.argv)

  ## Initialize the monochrome .ctx file.
  write_out("Initializing empty monochrome subgraph.\n")
  mono_ctx = subgraph_name[:-4] + MONO
  init_monochrome_graph(mccortex_exe, mono_ctx)
  prev_num_records = get_num_records(mono_ctx)

  # Perform initial round of subgraph extraction.
  print("Performing round 1 of monochrome subgraph extraction.\n")
  extract_mono_subgraph(mccortex_exe, mono_ctx, ctx_list, seed, distance)
  curr_num_records = get_num_records(mono_ctx)
  print(f'After round 1 of subgraph extraction: \n prev_kmer / curr_kmer = {float(prev_num_records) / curr_num_records}.\n')
  print(prev_num_records) 
  print(curr_num_records) 
  iter_count = 2
  while not converged(prev_num_records, curr_num_records, convergence_ratio):
    print(f"Performing round {iter_count} of subgraph extraction.\n")
    prev_num_records = curr_num_records
    extract_mono_subgraph(mccortex_exe, mono_ctx, ctx_list, seed, distance)
    curr_num_records = get_num_records(mono_ctx)
    print(f'After round {iter_count} of subgraph extraction: \n prev_kmer / curr_kmer = {float(prev_num_records) / curr_num_records}.\n')
    print(prev_num_records) 
    print(curr_num_records) 
    iter_count += 1

  # Extract single-sample-coloured subgraphs. 
  extract_coloured_subgraphs(mccortex_exe, mono_ctx, subgraph_name, ctx_list)

  # Merge single-sample-coloured subgraphs. 
  merge_coloured_subgraphs(mccortex_exe, ctx_list, subgraph_name)

if __name__ == '__main__':
  main()

#def mc_init(mccortex_exe, subgraph_name, ctx_list):
#  """Builds an empty coloured .ctx."""
#
#  # Make empty fasta.
#  run_subprocess("touch init.fa")
#
#  # Build empty single colour graphs.
#  for ctx in ctx_list:
#    run_subprocess(mccortex_exe + f" build -k 31 -f -m {memory}G --sample {os.path.basename(ctx)[:-4]} --seq2 init.fa:init.fa {ctx[:-4]}.init.ctx")
#
#  # Make empty coloured .ctx.
#  scg = ' '.join([ctx[:-4]+'.init.ctx' for ctx in ctx_list])
#  run_subprocess(mccortex_exe + f" join -f -m {memory}G --out {subgraph_name} {scg}")
#
#  # Clean up.
#  run_subprocess(f'rm {scg} init.fa')

#def extract_subgraph(mccortex_exe, subgraph, ctx_list, seed, dist):
#  """Extracts a subgraph from each .ctx in ctx_list centered
#     on seed with a distance of dist kmers from seed. Each subgraph
#     is then joined to the multi-coloured subgraph"""
#
#  # Iteratively join subgraph to each single coloured graph
#  # and extract a new subgraph
#  for idx, ctx in enumerate(ctx_list):
#    run_subprocess(mccortex_exe + f" join -f -m {memory}G --out {subgraph[:-4]}.merged.ctx {subgraph} {idx}:{ctx}")
#    sys.exit()
#    run_subprocess(mcorrted_exe + f" subgraph -f -m {memory}G -t {threads} -d {dist} --out {subgraph} --seq {seed} {subgraph[:-4]}.merged.ctx")
