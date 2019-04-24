import sys
import os
import urllib.request
import gzip

def get_kegg_organism_list():
  url = 'http://rest.kegg.jp/list/organism'
  data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
  data.pop()
  return data

def get_genome_entry_from_kegg(kegg_genome_id):
  url = f'http://rest.kegg.jp/get/gn:{kegg_genome_id}'
  # Download webpage from kegg, read, and convert to regular str type
  return urllib.request.urlopen(url).read().decode('utf-8')

def main():
  if len(sys.argv) != 2:
    print('Usage: <exe> <Relevant Taxa')

  organisms = get_kegg_organism_list()
  organisms = [elem for elem in organisms if sys.argv[1] in elem]
  for org in organisms:
    genome_id, org_id, strain_name, taxa_string = org.split('\t')
    filename = genome_id + '_' + org_id + '.kegg'
    if os.path.isfile(filename):
      print('didnt rewrite %s' % filename)
      continue
    with open(filename,'w') as f:
      print('writing %s' % filename)
      f.write(get_genome_entry_from_kegg(genome_id))

if __name__ == '__main__':
  main()
