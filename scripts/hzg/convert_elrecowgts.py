#!/usr/bin/env python3
"""@package docstring
Converts electron reco weights from root to correctionlib JSON
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from math import hypot
import json
import ROOT

#constants

ETA_BINS = [-2.5, -2.0, -1.566, -1.444, -1.0, 0.0, 1.0, 1.444, 1.566, 2.0, 2.5]
PT_BINS = [10.0, 20.0, 45.0, 75.0, 100.0, 500.0]

def fix_correctionlib_json(json_texts):
  '''Fixes the format of correctionlib json created using corr.json, since 
  it is not properly formatted by default
  '''
  corr = []
  for json_text in json_texts:
    corr.append(json.loads(json_text))
  json_dict = {
    'schema_version' : 2,
    'description' : '',
    'corrections' : corr
  }
  return json.dumps(json_dict,indent=2)

def get_bin_content(hist: ROOT.TH1, values: tuple) -> float:
  '''Gets bin content of histogram for bin corresponding to variable values

  Args:
    hist: histogram to get content from
    values: input (axis) variables of histogram

  Returns:
    Content of histogram bin into which values fall
  '''
  if len(values)==1:
    return hist.GetBinContent(hist.FindBin(values[0]))
  elif len(values)==2:
    return hist.GetBinContent(hist.FindBin(values[0], values[1]))
  else:
    return hist.GetBinContent(hist.FindBin(values[0], values[1], values[2]))

def calculate_sfs(pt: float, eta: float, EGamma_SF2D: ROOT.TH2D, 
                  EGamma_EffData2D: ROOT.TH2D, EGamma_EffMC2D: ROOT.TH2D, 
                  statData: ROOT.TH2D, statMC: ROOT.TH2D, 
                  altBkgModel: ROOT.TH2D, altSignalModel: ROOT.TH2D,
                  altMCEff: ROOT.TH2D, 
                  altTag: ROOT.TH2D) -> tuple[float,float,float,float]:
  '''Calculates SFs and uncertainties from TH2Ds and pt/eta

  Args:
    pt: pt to find bin
    eta: eta to find bin
    EGamma_SF2D: pass scale factors
    EGamma_EffData2D: data efficiency
    EGamma_EffMC2D: simulation efficiency
    statData: pass SF uncertainty from data stats
    statMC: pass SF uncertainty from simulation stats
    altBkgModel: pass SF uncertainty from uncertainty in background model
    altSignalModel: pass SF uncertainty from uncertainty in signal model
    altMCEff: pass SF uncertainty from uncertainty in MC sample
    altTag: pass SF uncertainty from uncertainty in tag selection
    
  Returns:
    pass SF, pass uncertainty, fail SF, fail uncertainty
  '''

  sf_pass = get_bin_content(EGamma_SF2D, (eta, pt))
  data_eff = get_bin_content(EGamma_EffData2D, (eta, pt))
  simu_eff = get_bin_content(EGamma_EffMC2D, (eta, pt))
  unc_statdata = get_bin_content(statData, (eta, pt))
  unc_statsimu = get_bin_content(statMC, (eta, pt))
  unc_altbkg = get_bin_content(altBkgModel, (eta, pt))
  unc_altsig = get_bin_content(altSignalModel, (eta, pt))
  unc_altmc = get_bin_content(altMCEff, (eta, pt))
  unc_alttag = get_bin_content(altTag, (eta, pt))
  sf_fail = 1.0
  data_factor = 1.0
  simu_factor = 1.0
  if (simu_eff < 1.0):
    sf_fail = (1.0-data_eff)/(1.0-simu_eff)
    data_factor = simu_eff/(1.0-simu_eff)
    if (data_eff > 0.0):
      simu_factor = (simu_eff/(1.0-simu_eff))**2*((1.0-data_eff)/data_eff)
  unc_pass = hypot(unc_statdata, unc_statsimu, unc_altbkg, unc_altsig, 
                   unc_altmc, unc_alttag)
  unc_fail = hypot(unc_statdata*data_factor, unc_statsimu*simu_factor, 
                   unc_altbkg*data_factor, unc_altsig*data_factor, 
                   unc_altmc*simu_factor, unc_alttag*data_factor)
  return sf_pass, unc_pass, sf_fail, unc_fail

def make_correction(name: str, desc: str, 
                    sfs: list[float]) -> schemav2.Correction:
  '''Generates correction object

  Args:
    name: correction name
    desc: correction description
    sfs: correction values

  Returns:
    schemva2 correction object
  '''
  return schemav2.Correction(
      name=name,
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='pt'),
              schemav2.Variable(name='eta', type='real', description='eta')],
      output=schemav2.Variable(name='sf', type='real', description=desc),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','eta'],
          edges=[PT_BINS, ETA_BINS],
          content=sfs,
          flow='clamp',
          ),
      )

def convert_files(lopt_name: str, hipt_name: str, out_name: str):
  '''Converts ROOT data to correctionlib json

  Args:
    lopt_name: name of low pt weight ROOT file
    hipt_name: name of high pt weight ROOT file
    out_name: name of correctionlib json file
  '''
  eta_means = [(ETA_BINS[ieta]+ETA_BINS[ieta+1])/2.0 for ieta in 
               range(len(ETA_BINS)-1)]
  pt_means = [(PT_BINS[ipt]+PT_BINS[ipt+1])/2.0 for ipt in 
              range(len(PT_BINS)-1)]
  sfs_pass = []
  uns_pass = []
  sfs_fail = []
  uns_fail = []
  lopt_file = ROOT.TFile(lopt_name, 'READ')
  for eta_mean in eta_means:
    sf_pass, unc_pass, sf_fail, unc_fail = calculate_sfs(pt_means[0], eta_mean,
        lopt_file.EGamma_SF2D, lopt_file.EGamma_EffData2D, 
        lopt_file.EGamma_EffMC2D, lopt_file.statData, lopt_file.statMC, 
        lopt_file.altBkgModel, lopt_file.altSignalModel, lopt_file.altMCEff, 
        lopt_file.altTagSelection)
    sfs_pass.append(sf_pass)
    uns_pass.append(unc_pass)
    sfs_fail.append(sf_fail)
    uns_fail.append(unc_fail)
  lopt_file.Close()
  hipt_file = ROOT.TFile(hipt_name, 'READ')
  for pt_mean in pt_means[1:]:
    for eta_mean in eta_means:
      sf_pass, unc_pass, sf_fail, unc_fail = calculate_sfs(pt_mean, eta_mean,
          hipt_file.EGamma_SF2D, hipt_file.EGamma_EffData2D, 
          hipt_file.EGamma_EffMC2D, hipt_file.statData, hipt_file.statMC, 
          hipt_file.altBkgModel, hipt_file.altSignalModel, hipt_file.altMCEff, 
          hipt_file.altTagSelection)
      sfs_pass.append(sf_pass)
      uns_pass.append(unc_pass)
      sfs_fail.append(sf_fail)
      uns_fail.append(unc_fail)
  hipt_file.Close()
  clib_sfs_pass = make_correction('sf_pass', 'data-MC pass SF', sfs_pass)
  clib_uns_pass = make_correction('unc_pass', 'data-MC pass unc', uns_pass)
  clib_sfs_fail = make_correction('sf_fail', 'data-MC fail SF', sfs_fail)
  clib_uns_fail = make_correction('unc_fail', 'data-MC fail unc', uns_fail)
  with open(out_name,'w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_sfs_pass.json(exclude_unset=True),
         clib_uns_pass.json(exclude_unset=True),
         clib_sfs_fail.json(exclude_unset=True),
         clib_uns_fail.json(exclude_unset=True)]))

def convert_files_run3(lopt_name: str, mdpt_name: str, hipt_name: str, 
                       out_name: str):
  '''Converts ROOT data to correctionlib json

  Args:
    lopt_name: name of low pt weight ROOT file
    mdpt_name: name of mid pt weight ROOT file
    hipt_name: name of high pt weight ROOT file
    out_name: name of correctionlib json file
  '''
  eta_means = [(ETA_BINS[ieta]+ETA_BINS[ieta+1])/2.0 for ieta in 
               range(len(ETA_BINS)-1)]
  pt_means = [(PT_BINS[ipt]+PT_BINS[ipt+1])/2.0 for ipt in 
              range(len(PT_BINS)-1)]
  sfs_pass = []
  uns_pass = []
  sfs_fail = []
  uns_fail = []
  lopt_file = ROOT.TFile(lopt_name, 'READ')
  for eta_mean in eta_means:
    sf_pass, unc_pass, sf_fail, unc_fail = calculate_sfs(pt_means[0], eta_mean,
        lopt_file.EGamma_SF2D, lopt_file.EGamma_EffData2D, 
        lopt_file.EGamma_EffMC2D, lopt_file.statData, lopt_file.statMC, 
        lopt_file.altBkgModel, lopt_file.altSignalModel, lopt_file.altMCEff, 
        lopt_file.altTagSelection)
    sfs_pass.append(sf_pass)
    uns_pass.append(unc_pass)
    sfs_fail.append(sf_fail)
    uns_fail.append(unc_fail)
  lopt_file.Close()
  mdpt_file = ROOT.TFile(mdpt_name, 'READ')
  for pt_mean in pt_means[1:3]:
    for eta_mean in eta_means:
      sf_pass, unc_pass, sf_fail, unc_fail = calculate_sfs(pt_mean, eta_mean,
          mdpt_file.EGamma_SF2D, mdpt_file.EGamma_EffData2D, 
          mdpt_file.EGamma_EffMC2D, mdpt_file.statData, mdpt_file.statMC, 
          mdpt_file.altBkgModel, mdpt_file.altSignalModel, mdpt_file.altMCEff, 
          mdpt_file.altTagSelection)
      sfs_pass.append(sf_pass)
      uns_pass.append(unc_pass)
      sfs_fail.append(sf_fail)
      uns_fail.append(unc_fail)
  mdpt_file.Close()
  hipt_file = ROOT.TFile(hipt_name, 'READ')
  for pt_mean in pt_means[3:]:
    for eta_mean in eta_means:
      sf_pass, unc_pass, sf_fail, unc_fail = calculate_sfs(pt_mean, eta_mean,
          hipt_file.EGamma_SF2D, hipt_file.EGamma_EffData2D, 
          hipt_file.EGamma_EffMC2D, hipt_file.statData, hipt_file.statMC, 
          hipt_file.altBkgModel, hipt_file.altSignalModel, hipt_file.altMCEff, 
          hipt_file.altTagSelection)
      sfs_pass.append(sf_pass)
      uns_pass.append(unc_pass)
      sfs_fail.append(sf_fail)
      uns_fail.append(unc_fail)
  hipt_file.Close()
  clib_sfs_pass = make_correction('sf_pass', 'data-MC pass SF', sfs_pass)
  clib_uns_pass = make_correction('unc_pass', 'data-MC pass unc', uns_pass)
  clib_sfs_fail = make_correction('sf_fail', 'data-MC fail SF', sfs_fail)
  clib_uns_fail = make_correction('unc_fail', 'data-MC fail unc', uns_fail)
  with open(out_name,'w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_sfs_pass.json(exclude_unset=True),
         clib_uns_pass.json(exclude_unset=True),
         clib_sfs_fail.json(exclude_unset=True),
         clib_uns_fail.json(exclude_unset=True)]))

if __name__ == '__main__':
  #for year in ['2016preVFP','2016postVFP','2017','2018']:
  #  convert_files(f'egammaEffi_ptBelow20.txt_EGM2D_UL{year}.root',
  #                f'egammaEffi_ptAbove20.txt_EGM2D_UL{year}.root', 
  #                f'electron_recoSF{year}.json')
  for year in ['2022','2022EE','2023','2023BPix']:
    convert_files_run3(f'egammaEffi.txt_EGM2D_lowpt_{year}.root',
                       f'egammaEffi.txt_EGM2D_midpt_{year}.root', 
                       f'egammaEffi.txt_EGM2D_highpt_{year}.root', 
                       f'electron_recoSF{year}.json')
