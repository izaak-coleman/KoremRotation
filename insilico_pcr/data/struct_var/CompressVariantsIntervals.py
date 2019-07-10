import sys
import csv
import re
import json
LEFT, RIGHT = 0, 1
def load_raw_data(fname):
  """Parses the SGVFinder variant .csv into a dictionary containing
     k = ncbi_taxonomy_id, v = list of strings of variant intervals."""
  with open(fname) as f:
    raw_data = [l.strip() for l in f]
  # remove file header
  raw_data.pop(0)
  # load raw data
  raw_dict = dict()
  for line in raw_data:
    taxa = re.search(r',([0-9]+\.[A-Za-z0-9]+):', line)
    variants = re.findall(r'([1-9][0-9]*_[1-9][0-9]*)',line)
    if taxa == None or len(variants) == 0:
      print(f'Skipped {line}')
      continue
    taxa = taxa.group(1)
    if taxa in raw_dict:
      raw_dict[taxa] += variants
    else:
      raw_dict[taxa] = variants

  return raw_dict

def overlap(intervals):
  for i in range(0, len(intervals)-1):
    if intervals[i][1] > intervals[i+1][0]:
      return True
  return False

def compress_variants(raw_dict, max_dist):
  """Join variant intervals that are within max_dist distance of one another."""
  compressed_dict = dict()
  for k, v in raw_dict.items():
    # Integerize the list of variant intervals
    intervals = list()
    for ivl in v:
      hits = re.search(r'([1-9][0-9]*)_([1-9][0-9]*)',ivl)
      if hits == None:
        print(f'Skipped compression of {k}.')
        continue
      # else, integerize 
      intervals.append( (int(hits[1]), int(hits[2])) )

    # Compress intervals
    intervals = sorted(intervals, key = lambda x: x[0])
    if overlap(intervals):
      print("Intervals overlap!!!")
   
    compressed_intervals = list()
    if len(intervals) == 1:
      compressed_dict[k] = (intervals[0][LEFT], intervals[0][RIGHT])
      continue
      
    i = 0
    while i < len(intervals):
      left, right = intervals[i][LEFT], intervals[i][RIGHT]
      j = i + 1
      if j >= len(intervals):
        break
      while (intervals[j][LEFT] - right) <= max_dist:
        right = intervals[j][RIGHT]
        j += 1
        if j >= len(intervals):
          break
      compressed_intervals.append( (left, right) )
      i = j

    # Store compressed intervals
    compressed_dict[k] = compressed_intervals

  return compressed_dict

def main():
  if len(sys.argv) != 4:
    print('Usage: <exe> <struct_file> <dist> <compressed_filename>')

  data = load_raw_data(sys.argv[1])
  data = compress_variants(data, int(sys.argv[2]))
  with open(sys.argv[3], 'w') as f:
    json.dump(data, f, indent=2, sort_keys=True)

if __name__ == '__main__':
  main()
    
