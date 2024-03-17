#!/usr/bin/env python3
'''
Script that applies background subtraction weights to photon samples
'''
import ROOT
import subprocess
from argparse import ArgumentParser
from skim_and_slim import write_ntuples

if __name__=='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='apply_photon_bkgsubweights',
      description='Applies weights (-j) to account for background contaimination in data to a file (-i)')
  argument_parser.add_argument('-i','--input_filename')
  argument_parser.add_argument('-j','--json_filename')
  args = argument_parser.parse_args()
  
  #JIT C++ to read weights
  ROOT.gInterpreter.AddIncludePath('inc/')
  ROOT.gInterpreter.ProcessLine('#include "correction_wrapper.hpp"')
  ROOT.gSystem.Load('libSmallPhysUtils.so')
  ROOT.gInterpreter.Declare("""
  
  using std::vector;
  
  const CorrectionWrapper corr_bkgsub("""+'"'+args.json_filename+'"'+""","bkg_subtraction_weight");
  
  float get_bkgweights(float ph_et, float ph_sc_abseta, float pair_mass) {
    vector<double> eval_args;
    eval_args.push_back(ph_et);
    eval_args.push_back(ph_sc_abseta);
    eval_args.push_back(pair_mass);
    if (ph_et < 25)
      return corr_bkgsub.evaluate(eval_args);
    if (pair_mass > 101 || pair_mass < 81)
      return 0.0;
    return 1.0;
  }
  
  """)

  if (args.input_filename[-5:] != '.root'):
    raise ValueError('Input file extension should be root')
  temp_filename = args.input_filename[:-5]+'_temp.root'

  #apply weights
  ROOT.EnableImplicitMT()
  defines = [('w_bkg','get_bkgweights(ph_et,ph_sc_abseta,pair_mass)')]
  write_ntuples(
      [args.input_filename],
      [],
      temp_filename,
      defines,
      'tree')

  subprocess.run(('rm '+args.input_filename).split())
  subprocess.run(('mv '+temp_filename+' '+args.input_filename).split())
