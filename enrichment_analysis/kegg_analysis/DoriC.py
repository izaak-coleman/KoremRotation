import csv
import collections

DoriC = collection.namedtuple('DoriC', ['ori', 'refseq'])

def get_entries(doric_file):
  entries = list()
  with open(doric_file) as f:
    reader = csv.reader(f):
    for entry in f:
      ori_start, ori_end = [int(e) for e in entry[4]]
      entries.append(DoriC(int((ori_start + ori_end)/2), entry[1]))
  return entries

def main():
  pass
# test this!
