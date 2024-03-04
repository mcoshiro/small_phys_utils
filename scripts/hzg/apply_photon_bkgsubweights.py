#!/usr/bin/env python3
'''
Script that applies background subtraction weights to photon samples
'''
import ROOT
import skim_and_slim

ROOT.gInterpreter.AddIncludePath('inc/')
ROOT.gInterpreter.ProcessLine('#include "correction_wrapper.hpp"')
ROOT.gSystem.Load('libSmallPhysUtils.so')
ROOT.gInterpreter.Declare("""

using std::vector;

const CorrectionWrapper corr_bkgsub("json/bkg_weights.json","bkg_subtraction_weight");

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

def apply_weights():
  ROOT.EnableImplicitMT()
  defines = [('w_bkgsub','get_bkgweights(ph_et,ph_sc_abseta,pair_mass)')]
  skim_and_slim.write_ntuples(
      ['/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_data_2018_new.root'],
      [],
      'photonidskim_data_2018_new_weighted.root',
      defines,
      'tree')

if __name__=='__main__':
  apply_weights()
