#!/usr/bin/env python3
"""@package docstring
Print fractions of events with Pythia-induced photon conversion
"""

import ROOT
from math import sqrt

ROOT.gInterpreter.Declare("""
using ROOT::VecOps::Any;
""")

def print_fractions(filenames):
  print('Analyzing {}'.format(filenames))

  df = ROOT.RDataFrame('tree', filenames)
  df = df.Define('w_const','1.0f')
  df = df.Define('w_lumisign','w_lumi/fabs(w_lumi)')
  total_entries_ptr = df.Sum('w_const')
  total_lumi_ptr = df.Sum('w_lumisign')
  df = df.Filter('Any(mc_id==22&&mc_mom==25&&(mc_statusflag&0x100)!=0)')
  photon_entries_ptr = df.Sum('w_const')
  photon_lumi_ptr = df.Sum('w_lumisign')
  df = df.Filter('Any(mc_id==22&&mc_mom==25&&(mc_statusflag&0x100)!=0&&mc_status==1)')
  fs_photon_entries_ptr = df.Sum('w_const')
  fs_photon_lumi_ptr = df.Sum('w_lumisign')

  total_entries = total_entries_ptr.GetValue()
  total_lumi = total_lumi_ptr.GetValue()
  photon_entries = photon_entries_ptr.GetValue()
  photon_lumi = photon_lumi_ptr.GetValue()
  fs_photon_entries = fs_photon_entries_ptr.GetValue()
  fs_photon_lumi = fs_photon_lumi_ptr.GetValue()

  no_fs_photon_entries = total_entries-fs_photon_entries
  df_du = (1.0/total_entries-fs_photon_entries/total_entries**2)
  df_dc = (fs_photon_entries/total_entries**2)
  f_unc = df_du*sqrt(fs_photon_entries)+df_dc*sqrt(no_fs_photon_entries)
  entry_fraction = fs_photon_entries/total_entries
  lumi_fraction = fs_photon_lumi/total_lumi

  print('Fraction of events with prompt photon: {}/{} ({})'.format(
      photon_entries, total_entries, photon_entries/total_entries))
  print('Fraction of lumi with prompt photon: {}/{} ({})'.format(
      photon_lumi, total_lumi, photon_lumi/total_lumi))
  print('Fraction of events with stable prompt photon: {}/{} ({}+-{})'.format(
      fs_photon_entries, total_entries, entry_fraction, f_unc))
  print('Fraction of lumi with stable prompt photon: {}/{} ({}+-{})'.format(
      fs_photon_lumi, total_lumi, lumi_fraction, 
      f_unc*lumi_fraction/entry_fraction))

def analyze_pythia_conversions():
  ROOT.EnableImplicitMT()
  folder_prefix = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/'
  all_files = []
  for year in ['2016APV/','2016/','2017/','2018/']:
    for filename in ['pico_ttHToZG_ZToLL_M-125_TuneCP5*',
                     'pico_ZH_HToZG_ZToAll_M-125_TuneCP5*',
                     'pico_WminusH_HToZG_WToAll_M-125_TuneCP5*',
                     'pico_WplusH_HToZG_WToAll_M-125_TuneCP5*',
                     'pico_VBFHToZG_ZToLL_M-125_TuneCP5*',
                     'pico_GluGluHToZG_ZToLL_M-125_TuneCP5*']:
      all_files.append('{}{}mc/unskimmed/{}'.format(folder_prefix,year,filename))
      #print_fractions('{}{}mc/unskimmed/{}'.format(folder_prefix,year,filename))
  folder_prefix = '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_lassen_v0/'
  for year in ['2022/','2022EE/']:
    for filename in ['pico_ttHtoZG_Zto2L_M-125_TuneCP5*',
                     'pico_ZH_HtoZG_ZtoAll_M-125_TuneCP5*',
                     'pico_WplusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5*',
                     'pico_WplusH_HtoZG_WtoAll_M-125_TuneCP5*',
                     'pico_WminusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5*',
                     'pico_VBFHtoZG_Zto2L_M-125_TuneCP5*',
                     'pico_GluGluHtoZG_Zto2L_M-125_TuneCP5*']:
      all_files.append('{}{}mc/unskimmed/{}'.format(folder_prefix,year,filename))
      #print_fractions('{}{}mc/unskimmed/{}'.format(folder_prefix,year,filename))
  for year in ['2023/','2023BPix/']:
    for filename in ['pico_ttHtoZG_Zto2L_M-125_TuneCP5*',
                     'pico_ZH_ZtoAll_HtoZGto2LG_M-125_TuneCP5*',
                     'pico_WplusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5*',
                     'pico_WminusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5*',
                     'pico_VBFHtoZG_Zto2L_M-125_TuneCP5*',
                     'pico_GluGluHtoZG_Zto2L_M-125_TuneCP5*']:
      all_files.append('{}{}mc/unskimmed/{}'.format(folder_prefix,year,filename))
      #print_fractions('{}{}mc/unskimmed/{}'.format(folder_prefix,year,filename))
  print_fractions(all_files)

if __name__ == '__main__':
  analyze_pythia_conversions()

