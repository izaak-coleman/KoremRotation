import sys

def main():
  if len(sys.argv) != 2:
    print("Usage <exe> <.pandatox>")
    sys.exit(1)
  with open(sys.argv[1]) as f:
    data = [l.strip() for l in f if (l.strip() != '\t' and l.strip() != '' and l.strip() != '-')]
    print(len(data))
    print('\n'.join(data))

if __name__ == '__main__':
  main()
