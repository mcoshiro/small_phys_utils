"""Utility to convert pileup JSON to C++ evaluator

MO 01.03.2025
"""

from argparse import ArgumentParser
from typing import List, Tuple
import json

def get_pusf_info(filename: str) -> Tuple[List[float],List[float]]:
  """Gets pileup scale-factor info from POG JSON

  Args:
    filename: LUMI POG pielup scale-factor json filename

  Returns:
    Tuple of bin edges and bin contents
  """
  json_text = ''
  with open(filename,'r') as json_file:
    json_text = json_file.read()
  json_content = json.loads(json_text)
  bin_edges = (json_content['corrections'][0]['data']['content'][0]['value']
                           ['edges'])
  sfs = (json_content['corrections'][0]['data']['content'][0]['value']
                     ['content'])
  return (bin_edges, sfs)

def write_cppevaluator(bin_edges: List[float], sfs: List[float], name: str):
  """Writes C++ evaluator for pileup scale factors

  Args:
    bin_edges: bin edges, assumed to be sorted
    sfs: scale factors
    name: name used for output function and file
  """
  filename = name+'.cpp'
  with open(filename,'w') as cpp_file:
    cpp_file.write('float {}(int npv) {{\n'.format(name))
    for ibin in range(len(bin_edges)-1):
      bedge = int(bin_edges[ibin])
      sf = sfs[ibin]
      if (ibin == 0):
        cpp_file.write('  if (npv == {}) return {};\n'.format(bedge,sf))
      elif (ibin < len(bin_edges)-2):
        cpp_file.write('  else if (npv == {}) return {};\n'.format(bedge,sf))
      else:
        cpp_file.write('  else return {};\n'.format(sf))
    cpp_file.write('}\n')

if __name__=='__main__':
  parser = ArgumentParser(description='Converts PU JSON in to C++ evaluator')
  parser.add_argument('-i', '--input_file', type=str)
  parser.add_argument('-n', '--name', type=str)
  args = parser.parse_args()
  bins, sfs = get_pusf_info(args.input_file)
  write_cppevaluator(bins, sfs, args.name)

