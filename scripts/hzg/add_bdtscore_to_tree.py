#!/usr/bin/env python3
'''
Script that adds a BDT score to a TTree
'''
import subprocess
import ROOT
import skim_and_slim

bdt_varnames = []

ROOT.gInterpreter.AddIncludePath('inc/')
ROOT.gInterpreter.ProcessLine('#include "hzg_mem.hpp"')
#ROOT.gSystem.AddLinkedLibs('-Llib -lSmallPhysUtils')
ROOT.gSystem.Load('libSmallPhysUtils.so')
ROOT.gInterpreter.ProcessLine('''
TMVA::Experimental::RReader model("dataset/weights/shuffled_phtree_ph_BDT.weights.xml");
computeModel = TMVA::Experimental::Compute<7, float>(model);
''')

if __name__=='__main__':
  #no MT with MVA? Not clear with experimental ROOT
  defines = [('photon_newmva_vec',ROOT.computeModel,ROOT.model.GetVariableNames()),
             ('photon_newmva','photon_newmva_vec[0]')]
  branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta',
      'phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity','photon_ptransverse','mllg',
      'ph_irel','ph_r9','ph_sieie','ph_hoe','photon_newmva','w_lumi_year')
  file_names = ['shuffled_phtree_sig','shuffled_phtree_bak']
  for file_name in file_names:
    skim_and_slim.write_ntuples(
        ['ntuples/'+file_name+'.root'],
        [],
        file_name+'_bdt.root',
        defines,
        'tree',
        branches)

