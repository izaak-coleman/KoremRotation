import sys
import gzip

def main():

  if len(sys.argv) != 8:
    print('Usage: <exe> <p1_boundary> <p2_boundary> <left_offset> <right_offset> <units> <probe_len> <concat_file>')
    sys.exit()

  p1_pos = (int(sys.argv[1]) - int(sys.argv[3]))*int(sys.argv[5])
  p2_pos = (int(sys.argv[2]) + int(sys.argv[4]) + 1)*int(sys.argv[5])

  with gzip.open(sys.argv[7]) as f:
    contig = f.read().strip().decode('utf-8')

  if p1_pos < 0 or p2_pos > len(contig):
    print(f'probe positions p1 {p1_pos} or p2 {p2_pos} are not in range')
    sys.exit()
  probe_len = int(sys.argv[6]) 
  print(p1_pos)
  print(p2_pos)
  p1 = contig[p1_pos: p1_pos + probe_len]
  p2 = contig[p2_pos - probe_len: p2_pos]
  print(f'{sys.argv[7]}, {p1}, {p2}')


if __name__ == '__main__':
  main()

