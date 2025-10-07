"""small script to verify root files
"""

from argparse import ArgumentParser
from glob import glob
import ROOT

def print_weight_type(fname):
  infile = ROOT.TFile(fname,'READ')
  max_evts = 10 
  ievt = 0
  avr_weight = 0.0
  pico_type = 0
  for event in infile.tree:
    avr_weight += event.weight
    pico_type = event.type
    ievt += 1
    if ievt >= max_evts: 
      break
  avr_weight /= max_evts
  print(f'File {fname}')
  print(f'avr weight: {avr_weight}')
  print(f'type: {pico_type}')

if __name__=='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='check_production',
      description='Checks ROOT files')
  argument_parser.add_argument('-i','--input_dir')
  args = argument_parser.parse_args()

  fnames = glob(args.input_dir+'/*.root')
  fnames.sort()
  for fname in fnames:
    print_weight_type(fname)

