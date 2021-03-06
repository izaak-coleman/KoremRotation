import csv
from recordclass import recordclass
import collections
import sys
import re

DoriC = collections.namedtuple('DoriC', 'position refseq id assembly')

def catch(func, handle=lambda e : e, *args, **kwargs):
  try:
    return func(*args, **kwargs)
  except Exception as e:
    return handle(e)
  

def get_entries(doric_file):
  entries = list()
  with open(doric_file) as f:
    data = [l.strip() for l in f]
    data.pop(0)
 
  get_ori_id = lambda x : re.findall(r'(^ORI[0-9]+),',x)
  get_ori_coords = lambda x : re.findall(r',([0-9]+)\.\.([0-9]+),',x)
  get_refseq = lambda x : re.findall(r'((?:NC|NZ)_[\w\.]+),',x)
  get_assembly = lambda x: re.findall(r'(GCF[\w\.]+)',x)
  for entry in data:
    if 'nt*' in entry:
      continue
    ori_id = get_ori_id(entry)
    ori_coords = get_ori_coords(entry)
    refseq = get_refseq(entry)
    assembly = get_assembly(entry)
    if len(ori_coords) == 0 or len(refseq) == 0 or len(ori_id) == 0 or len(assembly) == 0:
      continue
    ori_start, ori_end = ori_coords[0]
    ori_start = int(ori_start)
    ori_end = int(ori_end)
    entries.append(DoriC(ori_start, refseq[0], ori_id[0], assembly[0]))
  return entries

def main():
  if len(sys.argv) != 2:
    print('Usage: <exe> <doric_file>')
    sys.exit()

  doric = get_entries(sys.argv[1])
  print(len(doric))
  for d in doric:
    print(f'{d.position}, {d.refseq}, {d.id}, {d.assembly}')

if __name__ == '__main__':
  main()
