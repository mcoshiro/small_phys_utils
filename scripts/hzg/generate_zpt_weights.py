#!/usr/bin/env python
"""Generates 5D photon reweighting
"""

from argparse import ArgumentParser
from array import array
from collections.abc import Callable
from glob import glob
import math
import numpy as np
from rdataframeset import *
import ROOT
from root_plot_lib import RplPlot
from typing import Tuple, List, Union
from scipy.optimize import minimize

ROOT.gInterpreter.ProcessLine('.L scripts/hzg/hzg_corrections.cpp+')

LUMI = {'2016APV' : 19.51,
        '2016' : 16.80,
        '2017' : 41.48,
        '2018' : 59.83,
        '2022' : 8.17,
        '2022EE' : 27.01,
        '2023' : 17.61,
        '2023BPix' : 9.53}

Z_PT_BINS = [0.0,4.0,8.0,12.0,16.0,20.0,25.0,30.0,35.0,40.0,60.0,80.0,100.0,
             150.0,200.0]
P_PT_BINS = [15.0,17.5,20.0,25.0,30.0,45.0,200.0]
LLPH_PT_BINS = [0.0,5.0,10.0,15.0,20.0,30.0,40.0,60.0,80.0,100.0,500.0]
NJET_BINS = [-0.5,0.5,1.5,2.5,3.5,10.5]
NJ_LLPH_PT_BINS = [[0.0,5.0,10.0,15.0,20.0,30.0,40.0,60.0,500.0],
                   [0.0,5.0,10.0,15.0,20.0,30.0,40.0,60.0,80.0,500.0],
                   [0.0,15.0,20.0,30.0,40.0,60.0,80.0,100.0,500.0],
                   [0.0,30.0,40.0,60.0,80.0,100.0,500.0],
                   [0.0,40.0,80.0,500.0]]

FAKEPH_ABSETA_BINS = [0.0, 0.8, 1.5, 2.0, 2.5]
FAKEPH_PT_BINS = [15.0, 20.0, 30.0, 500.0]

def year_to_float(year: str) -> str:
  """Converts year to string of floating point representation

  Args:
    year: data-taking period

  Returns:
    string of float representing year
  """
  if year=='2016APV':
    return '2016.0'
  elif year=='2016':
    return '2016.5'
  elif year=='2017':
    return '2017.0'
  elif year=='2018':
    return '2018.0'
  elif year=='2022':
    return '2022.0'
  elif year=='2022EE':
    return '2022.5'
  elif year=='2023':
    return '2023.0'
  elif year=='2023BPix':
    return '2023.5'
  else:
    raise ValueError('Unknown year')
  return year

def setup_dataframes(common_defines: bool=True, 
                     add_signal: bool=False) -> List[RDataFrameSet]:
  """Initializes dataframes for Z pT reweighting

  Args:
    common_defines: defines common columns and restricts to sideband
    add_signal: adds signal dataframe sets

  Returns:
    data RDataframe and simulation RDataframeSet
  """
  dfs = []
  years_list_run2 = ['2016APV','2016','2017','2018']
  years_list_run3 = ['2022','2022EE','2023','2023BPix']
  procs = [(True, get_filenames('data','llg','run2'), years_list_run2),
           (False, get_filenames('mczg','llg','run2'), years_list_run2),
           (False, get_filenames('mcdy','llg','run2'), years_list_run2),
           (True, get_filenames('data','llg','run3'), years_list_run3),
           (False, get_filenames('mczg','llg','run3'), years_list_run3),
           (False, get_filenames('mcdy','llg','run3'), years_list_run3)]
  if add_signal:
    procs = [(True, get_filenames('data','llg','run2'), years_list_run2),
             (False, get_filenames('mczg','llg','run2'), years_list_run2),
             (False, get_filenames('mcdy','llg','run2'), years_list_run2),
             (False, get_filenames('hzph','llg','run2'), years_list_run2),
             (True, get_filenames('data','llg','run3'), years_list_run3),
             (False, get_filenames('mczg','llg','run3'), years_list_run3),
             (False, get_filenames('mcdy','llg','run3'), years_list_run3),
             (False, get_filenames('hzph','llg','run3'), years_list_run3)]
  for is_data, filenames, years in procs:
    sub_dfs = []
    for filename, year in zip(filenames, years):
      sub_dfs.append(ROOT.RDataFrame('tree',filename))
      sub_dfs[-1] = sub_dfs[-1].Define('year',year_to_float(year))
      sub_dfs[-1] = sub_dfs[-1].Define('photon_res',
          'static_cast<float>(photon_energyErr[0]/(photon_pt[0]*TMath::CosH('
          'photon_eta[0])))')
      if is_data:
        sub_dfs[-1] = sub_dfs[-1].Define('w_lumiyear','1')
      else:
        sub_dfs[-1] = sub_dfs[-1].Define('w_year',str(LUMI[year]))
        sub_dfs[-1] = sub_dfs[-1].Define('weight_fix',
            '(weight/w_lumi)>10.0 ? w_lumi*10 : weight')
        if year in ['2016APV','2016','2017','2018']:
          sub_dfs[-1] = sub_dfs[-1].Define('w_phshape',('get_photon_sf_run2({'+
              'photon_pt[0], fabs(photon_eta[0]), photon_idmva[0],'+
              ' static_cast<float>(photon_res)}, photon_pflavor[0])'))
        else:
          sub_dfs[-1] = sub_dfs[-1].Define('w_phshape',('get_photon_sf_run3({'+
              'photon_pt[0], fabs(photon_eta[0]), photon_idmva[0],'+
              ' static_cast<float>(photon_res)}, photon_pflavor[0])'))
        sub_dfs[-1] = sub_dfs[-1].Define('photon_isjet',
            'photon_isjet(nphoton, photon_eta, photon_phi, mc_pt, mc_eta, '
            'mc_phi, mc_id, mc_statusflag)')
        #sub_dfs[-1] = sub_dfs[-1].Define('w_jet',
        #      'get_w_jet(year, type, njet, photon_isjet)')
        sub_dfs[-1] = sub_dfs[-1].Define('w_photon_lowpt',
              'get_w_photon_lowpt(photon_sig, photon_pflavor, photon_pt, '
              +'photon_eta, year)')
        sub_dfs[-1] = sub_dfs[-1].Define('w_fakepteta','get_w_fakephoton(year'
              +', photon_isjet, photon_pt[0], fabs(photon_eta[0]), '
              +'photon_pflavor[0])')
        sub_dfs[-1] = sub_dfs[-1].Define('w_llph_pt','get_w_llph_pt(year'
              +', type, photon_isjet, llphoton_pt[0])')
        sub_dfs[-1] = sub_dfs[-1].Define('w_llpjj_pt','get_w_llpjj_pt(year'
              +', type, llphoton_dijet_balance, njet)')
        #sub_dfs[-1] = sub_dfs[-1].Define('w_lumiyear',
        #      'weight_fix*w_year*w_photon_lowpt*w_phshape*w_fakepteta'
        #      +'*w_llph_pt*w_llpjj_pt')
        #sub_dfs[-1] = sub_dfs[-1].Define('w_lumiyear',
        #      'weight_fix*w_year*w_photon_lowpt*w_phshape*w_fakepteta')
        sub_dfs[-1] = sub_dfs[-1].Define('w_lumiyear',
              'weight_fix*w_year*w_photon_lowpt*w_phshape')
        sub_dfs[-1] = sub_dfs[-1].Filter('use_event')
    dfs.append(RDataFrameSet(sub_dfs))
    #sideband selection
    dfs[-1] = dfs[-1].Filter('nll>=1&&nphoton>=1')
    dfs[-1] = dfs[-1].Filter('pass_trigs(trig_single_el,trig_single_mu,'
        +'trig_double_el,trig_double_mu,nel,nmu,el_pt,mu_pt,year)')
    dfs[-1] = dfs[-1].Filter('ll_m[0]>80&&ll_m[0]<100')
    dfs[-1] = dfs[-1].Filter('((photon_pt[0]/llphoton_m[0])>(15.0/110.0))')
    dfs[-1] = dfs[-1].Filter('(ll_m[0]+llphoton_m[0])>185.0')
    dfs[-1] = dfs[-1].Filter('(llphoton_m[0]>100)&&(llphoton_m[0]<180)')
    if common_defines:
      dfs[-1] = dfs[-1].Filter('(llphoton_m[0]<120)||(llphoton_m[0]>130)')
    #define things
    if common_defines:
      dfs[-1] = dfs[-1].Define('z_pt','ll_pt[0]')
      dfs[-1] = dfs[-1].Define('ll_pt0','ll_pt[0]')
      dfs[-1] = dfs[-1].Define('ph_pt','photon_pt[0]')
      dfs[-1] = dfs[-1].Define('zg_pt','llphoton_pt[0]')
      dfs[-1] = dfs[-1].Define('mht0',('get_mht('
          +'photon_pt, photon_phi, photon_sig, el_pt, el_phi, el_sig, mu_pt,'
          +' mu_phi, mu_sig, jet_pt, jet_eta, jet_phi, jet_isgood)'))
      dfs[-1] = dfs[-1].Define('ht0',('get_ht('
          +'photon_pt, photon_phi, photon_sig, el_pt, el_phi, el_sig, mu_pt,'
          +' mu_phi, mu_sig, jet_pt, jet_eta, jet_phi, jet_isgood)'))
      dfs[-1] = dfs[-1].Define('photon_mht_dphi',('get_photon_mht_dphi('
          +'photon_pt, photon_phi, photon_sig, el_pt, el_phi, el_sig, mu_pt,'
          +' mu_phi, mu_sig, jet_pt, jet_eta, jet_phi, jet_isgood)'))
      dfs[-1] = dfs[-1].Define('photon_pt0','photon_pt[0]')
      dfs[-1] = dfs[-1].Define('photon_idmva0','photon_idmva[0]')
      dfs[-1] = dfs[-1].Define('photon_abseta0','fabs(photon_eta[0])')
      dfs[-1] = dfs[-1].Define('llphoton_pt0','llphoton_pt[0]')
      dfs[-1] = dfs[-1].Define('llphoton_abseta0','fabs(llphoton_eta[0])')
      dfs[-1] = dfs[-1].Define('llphoton_abscosTheta0',
                               'fabs(llphoton_cosTheta[0])')
      dfs[-1] = dfs[-1].Define('photon_zeppenfeld0','photon_zeppenfeld[0]')
      dfs[-1] = dfs[-1].Define('photon_jet1_dr0','photon_jet1_dr[0]')
      dfs[-1] = dfs[-1].Define('llphoton_dijet_balance0',
                               'llphoton_dijet_balance[0]')
      dfs[-1] = dfs[-1].Define('llphoton_dijet_absdphi0',
                               'fabs(llphoton_dijet_dphi[0])')
      dfs[-1] = dfs[-1].Define('lead_jet_pt',
                               'get_lead_jet_pt(jet_isgood, jet_pt, jet_eta)')
      dfs[-1] = dfs[-1].Define('lead_jet_abseta',
          'fabs(get_lead_jet_eta(jet_isgood, jet_pt, jet_eta))')
      dfs[-1] = dfs[-1].Define('sublead_jet_pt',
          'get_sublead_jet_pt(jet_isgood, jet_pt, jet_eta)')
      dfs[-1] = dfs[-1].Define('sublead_jet_abseta',
          'fabs(get_sublead_jet_eta(jet_isgood, jet_pt, jet_eta))')
  return dfs

def get_filenames(proc_type: str, slim_type: str, era: str) -> list[list[str]]:
  """Returns ll skim filenames for given sample type

  Args:
    proc_type: "data", "mczg", "mcdy", or "hzph"
    slim_type: "ll" or "llg"
    era: "run2", "run3", or "all"

  Returns:
    list indexed by year of list of filenames
  """
  filenames = []
  if not proc_type in ['data','mczg','mcdy','hzph']:
    raise ValueError('Unsupported proc_type')
  if not slim_type in ['ll','llg']:
    raise ValueError('Unsupported slim_type')
  if not era in ['run2','run3']:
    raise ValueError('Unsupported era')
  datamc = 'data'
  proc_names = []
  if proc_type == 'data':
    proc_names = ['*Electron*','*Muon*','*EGamma*','*DoubleEG*']
  elif proc_type == 'mczg':
    proc_names = ['*ZGToLLG_01J_5f_lowMLL*',
                  '*DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5*',
                  '*DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5*',
                  '*DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5*',
                  '*DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5*',
                  '*DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5*']
    datamc = 'mc'
  elif proc_type == 'mcdy':
    proc_names = ['*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX*',
                  '*DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX*']
    datamc = 'mc'
  elif proc_type == 'hzph':
    proc_names = [
        '*GluGluHToZG_ZToLL_M-125_TuneCP5*',
        '*VBFHToZG_ZToLL_M-125_TuneCP5_withDipoleRecoil*',
        '*WplusH_HToZG_WToAll_M-125_TuneCP5*',
        '*WminusH_HToZG_WToAll_M-125_TuneCP5*',
        '*ZH_HToZG_ZToAll_M-125_TuneCP5*',
        '*ttHToZG_ZToLL_M-125_TuneCP5*',
        '*GluGluHtoZG_Zto2L_M-125_TuneCP5*',
        '*VBFHtoZG_Zto2L_M-125_TuneCP5_withDipoleRecoil*',
        '*WplusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5*',
        '*WminusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5*',
        '*ZH_HtoZG_ZtoAll_M-125_TuneCP5*',
        '*ttHtoZG_Zto2L_M-125_TuneCP5*']
    datamc = 'mc'
  run2_yeardata = [('2016APV','NanoAODv9UCSB2','htozgamma_pinnacles_v1'),
                   ('2016','NanoAODv9','htozgamma_pinnacles_v0'),
                   ('2017','NanoAODv9','htozgamma_pinnacles_v0'),
                   ('2018','NanoAODv9','htozgamma_pinnacles_v0')]
  run3_yeardata = [('2022','NanoAODv12','htozgamma_pinnacles_v0'),
                   ('2022EE','NanoAODv12','htozgamma_pinnacles_v0'),
                   ('2023','NanoAODv12','htozgamma_pinnacles_v0'),
                   ('2023BPix','NanoAODv12','htozgamma_pinnacles_v0')]
  if proc_type!='data':
    run2_yeardata[0] = ('2016APV','NanoAODv9','htozgamma_pinnacles_v0')
  yeardata = run2_yeardata
  if era=='run3':
    yeardata = run3_yeardata
  elif era=='all':
    yeardata += run3_yeardata
  for year, nano_type, production in yeardata:
    year_filenames = []
    for proc_name in proc_names:
      year_filenames += glob(f'/net/cms11/cms11r0/pico/{nano_type}/'
          +f'{production}/{year}/{datamc}/merged_zg{datamc}_{slim_type}/'
          +f'{proc_name}.root')
    filenames.append(year_filenames)
  return filenames

def setup_ll_dataframes() -> List[RDataFrameSet]:
  """Initializes dataframes for Z pT reweighting from ll slim

  Returns:
    RDataframeSets for data and simulation
  """
  dfs = []
  years_list = ['2016APV','2016','2017','2018','2022','2022EE','2023',
                '2023BPix']
  for is_data, filenames, years in [
      (True, get_filenames('data','ll','all'), years_list),
      (False, get_filenames('mczg','llg','all'), years_list),
      (False, get_filenames('mcdy','ll','all'), years_list)]:
    sub_dfs = []
    for filename, year in zip(filenames, years):
      sub_dfs.append(ROOT.RDataFrame('tree',filename))
      sub_dfs[-1] = sub_dfs[-1].Define('year',year_to_float(year))
      if is_data:
        sub_dfs[-1] = sub_dfs[-1].Define('w_total','1')
      else:
        sub_dfs[-1] = sub_dfs[-1].Define('w_year',str(LUMI[year]))
        sub_dfs[-1] = sub_dfs[-1].Define('weight_fix',
            '(weight/w_lumi)>10.0 ? w_lumi*10 : weight')
        if year in ['2016APV','2016','2017','2018']:
          sub_dfs[-1] = sub_dfs[-1].Define('w_phshape',('nphoton > 0 ? '
              +'get_photon_sf_run2({photon_pt[0], fabs(photon_eta[0]), '
              +'photon_idmva[0], static_cast<float>(photon_res)}) : 1.0'))
        else:
          sub_dfs[-1] = sub_dfs[-1].Define('w_phshape',('nphoton > 0 ? '
              +'get_photon_sf_run3({photon_pt[0], fabs(photon_eta[0]), '
              +'photon_idmva[0], static_cast<float>(photon_res)}) : 1.0'))
        sub_dfs[-1] = sub_dfs[-1].Define('w_photon_lowpt',
              'get_w_photon_lowpt(photon_sig, photon_pflavor, photon_pt, '
              +'photon_eta, year)')
        #sub_dfs[-1] = sub_dfs[-1].Define('w_lumiyear',
        #      'weight_fix*w_year*w_phshape*w_photon_lowpt')
        sub_dfs[-1] = sub_dfs[-1].Filter('use_event')
    dfs.append(RDataFrameSet(sub_dfs))
    dfs[-1] = dfs[-1].Filter('nll>=1')
    dfs[-1] = dfs[-1].Filter('pass_trigs(trig_single_el,trig_single_mu,'
        +'trig_double_el,trig_double_mu,nel,nmu,el_pt,mu_pt,year)')
  return dfs

def get_cr_z_weights():
  """Gets Z pt weights from control regions: ZG from Z->llG, 
  DY+jet from Z->ll+jet, and DY+PU from Z->ll
  """
  df_data, df_mczg, df_mcdy = setup_ll_dataframes()
  #Z -> llph
  df_mczg = df_mczg.Filter('nphoton>=1&&llphoton_m[0]>80&&llphoton_m[0]<100')
  df_mczg = df_mczg.Define('llphoton_pt0','llphoton_pt[0]')
  df_mczg = df_mczg.Define('w_total',
                           'weight_fix*w_years*w_phshape*w_photon_lowpt')
  df_datp = df_data.Filter('nphoton>=1&&llphoton_m[0]>80&&llphoton_m[0]<100')
  df_datp = df_datp.Define('llphoton_pt0','llphoton_pt[0]')
  llph_simu_hist_ptrs = df_mczg.Histo1D(('llph_simu_hist',
      'Z#rightarrow ll#gamma;p_{T}^{ll#gamma}',len(LLPH_PT_BINS)-1,
      array('d',LLPH_PT_BINS)),'llphoton_pt0','w_total')
  llph_data_hist_ptrs = df_data.Histo1D(('llph_data_hist',
      'Data;p_{T}^{ll#gamma}',len(LLPH_PT_BINS)-1,
      array('d',LLPH_PT_BINS)),'llphoton_pt0','w_total')
  #Z -> ll
  df_mcdy = df_mcdy.Filter('ll_m[0]>80&&ll_m[0]<100')
  df_mcdy = df_mcdy.Define('ll_pt0','ll_pt[0]')
  df_mcdy = df_mcdy.Define('w_total','weight_fx*w_years')
  df_data = df_data.Filter('ll_m[0]>80&&ll_m[0]<100')
  df_data = df_data.Define('ll_pt0','ll_pt[0]')
  ll_simu_hist_ptrs = df_mcdy.Histo1D(('llph_simu_hist',
      'Z#rightarrow ll;p_{T}^{ll}',len(LLPH_PT_BINS)-1,
      array('d',LLPH_PT_BINS)),'ll_pt0','w_total')
  ll_data_hist_ptrs = df_data.Histo1D(('llph_data_hist',
      'Data;p_{T}^{ll}',len(LLPH_PT_BINS)-1,
      array('d',LLPH_PT_BINS)),'ll_pt0','w_total')
  #Z->ll+jet
  df_mcdy = df_mcdy.Filter('njet>=1')
  df_data = df_data.Filter('njet>=1')
  df_mcdy = df_mcdy.Define('llj_pt0','get_llj_pt()')
  df_data = df_data.Define('llj_pt0','get_llj_pt()')

def get_jet_weights(dfs: List[RDataFrameSet]):
  """Calculates weights for reweighting Njets distribution
  """
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, '1', 'data2'),
                ('Z#gamma (run2)', 1, '1', 'mczg2'),
                ('DY+jet (run2)', 2, 'photon_isjet', 'dyjt2'),
                ('DY+PU (run2)', 2, '!photon_isjet', 'dypu2'),
                ('Data (run3)', 3, '1', 'data3'),
                ('Z#gamma (run3)', 4, '1', 'mczg3'),
                ('DY+jet (run3)', 5, 'photon_isjet', 'dyjt3'),
                ('DY+PU (run3)', 5, '!photon_isjet', 'dypu3')]
  for proc_info in procs_info:
    hist_ptrs.append([])
    for sel_name, id_sel in [
        ('loid_lodphi','photon_idmva[0]<0.8&&photon_mht_dphi<1.5758'),
        ('loid_hidphi','photon_idmva[0]<0.8&&photon_mht_dphi>1.5758'),
        ('hiid','photon_idmva[0]>=0.8')]:
      proc_name = proc_info[0]
      df = dfs[proc_info[1]].Filter(proc_info[2])
      simple_name = proc_info[3]
      hist_ptrs[-1].append(df.Filter(id_sel).Histo1D((
          f'njet_{sel_name}_{simple_name}',
          f'{proc_name}; N_{{jet}}',
          len(NJET_BINS)-1,array('d',NJET_BINS)),
          'njet','w_lumiyear'))

  for run in [0, 4]:
    data2_hists = [merge_hist_ptrs(hist_ptrs[0+run][icat]) for icat in 
                   range(3)] 
    simu2_hists = [[merge_hist_ptrs(hist_ptrs[iproc+run][icat]) for icat in 
                   range(3)] for iproc in range(1,4)]
    inclusive_sfs = solve_hist_system_manual(data2_hists, simu2_hists, 
                                             lambda x: x.Integral())
    print(f'inclusive sfs: {inclusive_sfs}')
    for iproc in range(3):
      for icat in range(3):
        simu2_hists[iproc][icat].Scale(inclusive_sfs[iproc])
    differential_sfs = []
    for ijet in range(len(NJET_BINS)-1):
      differential_sfs.append(solve_hist_system_manual(data2_hists, 
          simu2_hists, lambda x: x.GetBinContent(ijet+1)))
    print(f'differential sfs: {differential_sfs}')

def make_minimizer_hists():
  """Generates histograms for doing simultaneous fit
  """
  dfs = setup_dataframes()
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, '1', 'data2'),
                ('Z#gamma (run2)', 1, '1', 'mczg2'),
                ('DY+jet (run2)', 2, 'photon_isjet', 'dyjt2'),
                ('DY+PU (run2)', 2, '!photon_isjet', 'dypu2'),
                ('Data (run3)', 3, '1', 'data3'),
                ('Z#gamma (run3)', 4, '1', 'mczg3'),
                ('DY+jet (run3)', 5, 'photon_isjet', 'dyjt3'),
                ('DY+PU (run3)', 5, '!photon_isjet', 'dypu3')]
  hist_ptrs = []
  SIMUL_ZPH_PT_BINS = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 30.0, 40.0, 80.0, 
                       500.0]
  SIMUL_PH_PT_BINS = [15.0, 17.5, 20.0, 30.0, 50.0]
  SIMUL_PH_ETA_BINS = [0.0, 0.8, 1.5, 2.0, 2.5]
  SIMUL_PH_ID_BINS = [0.14, 0.33, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0]
  for proc_info in procs_info:
    for izph_pt in range(len(SIMUL_ZPH_PT_BINS)-1):
      for iph_pt in range(len(SIMUL_PH_PT_BINS)-1):
        bin_selection = (f'llphoton_pt0>={SIMUL_ZPH_PT_BINS[izph_pt]}&&'
                         f'llphoton_pt0<{SIMUL_ZPH_PT_BINS[izph_pt+1]}&&'
                         f'photon_pt0>={SIMUL_PH_PT_BINS[iph_pt]}&&'
                         f'photon_pt0<{SIMUL_PH_PT_BINS[iph_pt+1]}')
        df_bin = dfs[proc_info[1]].Filter(proc_info[2]).Filter(bin_selection)
        proc_name = proc_info[0]
        simple_name = proc_info[3]
        hist_ptrs.append(df_bin.Histo2D((
          f'{simple_name}{izph_pt}_{iph_pt}',
          f'{proc_name}; Photon |#eta|; Photon IDVMA',
          len(SIMUL_PH_ETA_BINS)-1,array('d',SIMUL_PH_ETA_BINS),
          len(SIMUL_PH_ID_BINS)-1,array('d',SIMUL_PH_ID_BINS)),
          'photon_abseta0', 'photon_idmva0', 'w_lumiyear'))
  hist_file = ROOT.TFile('bkg_weight_hists.root', 'RECREATE')
  hists = []
  for hist_ptr in hist_ptrs:
    hists.append(merge_hist_ptrs(hist_ptr))
    hists[-1].Write()
    #hists[-1].Scale(1.0/hists[-1].Integral())
  hist_file.Close()

def get_fakephoton_weights():
  """Calculates weights for reweighting fake photons in pt and eta
  """
  dfs = setup_dataframes()
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, '1', 'data2'),
                ('Z#gamma (run2)', 1, '1', 'mczg2'),
                ('DY+jet (run2)', 2, 'photon_isjet', 'dyjt2'),
                ('DY+PU (run2)', 2, '!photon_isjet', 'dypu2'),
                ('Data (run3)', 3, '1', 'data3'),
                ('Z#gamma (run3)', 4, '1', 'mczg3'),
                ('DY+jet (run3)', 5, 'photon_isjet', 'dyjt3'),
                ('DY+PU (run3)', 5, '!photon_isjet', 'dypu3')]
  for proc_info in procs_info:
    hist_ptrs.append([])
    #for sel_name, id_sel in [
    #    ('loid_lodphi','photon_idmva[0]<0.85&&photon_mht_dphi<1.5758'),
    #    ('loid_hidphi','photon_idmva[0]<0.85&&photon_mht_dphi>1.5758'),
    #    ('hiid','photon_idmva[0]>=0.85')]:
    for sel_name, id_sel in [
        ('lodphi','photon_mht_dphi<1.5758'),
        ('hidphi','photon_mht_dphi>1.5758')]:
      proc_name = proc_info[0]
      df = dfs[proc_info[1]].Filter(proc_info[2])
      simple_name = proc_info[3]
      hist_ptrs[-1].append(df.Filter(id_sel).Histo2D((
          f'fakeph_{sel_name}_{simple_name}',
          f'{proc_name}; p_{{T}}^{{#gamma}} [GeV]; #eta^{{#gamma}}',
          len(FAKEPH_PT_BINS)-1,array('d',FAKEPH_PT_BINS),
          len(FAKEPH_ABSETA_BINS)-1,array('d',FAKEPH_ABSETA_BINS)),
          'photon_pt0','photon_abseta0','w_lumiyear'))

  inclusive_sfs = []
  differential_sfs = []
  for run in [0, 4]:
    data_hists = [merge_hist_ptrs(hist_ptrs[0+run][icat]) for icat in 
                   range(2)] 
    simu_hists = [[merge_hist_ptrs(hist_ptrs[iproc+run][icat]) for icat in 
                   range(2)] for iproc in range(1,4)]
    ##reweight ZG only based on photon pT
    #inclusive_sfs.append([])
    #for ipt in range(len(FAKEPH_PT_BINS)-1):
    #  inclusive_sfs[-1].append(solve_hist_system_manual(
    #      data_hists, simu_hists, 
    #      lambda x: x.Integral(ipt+1,ipt+1,1,len(FAKEPH_ABSETA_BINS)-1)))

    ##inclusive_sfs.append(solve_hist_system_manual(data_hists, simu_hists, 
    ##                                              lambda x: x.Integral()))
    #print(f'inclusive sfs: {inclusive_sfs}')
    #for icat in range(3):
    #  for ipt in range(len(FAKEPH_PT_BINS)-1):
    #    for ieta in range(len(FAKEPH_ABSETA_BINS)-1):
    #      simu_hists[0][icat].SetBinContent(ipt, ieta, 
    #          inclusive_sfs[-1][ipt][0]*simu_hists[0][icat].GetBinContent(ipt, 
    #          ieta))
    for icat in range(2):
      data_hists[icat].Add(simu_hists[0][icat],-1.0)
    data_hists_reduced = [data_hists[icat] for icat in range(2)]
    simu_hists_reduced = [[simu_hists[iproc][icat] for icat in range(2)] 
                          for iproc in range(1,3)]
    differential_sfs.append([])
    for ipt in range(len(FAKEPH_PT_BINS)-1):
      differential_sfs[-1].append([])
      for ieta in range(len(FAKEPH_ABSETA_BINS)-1):
        differential_sfs[-1][-1].append(solve_hist_system_manual(
            data_hists_reduced, simu_hists_reduced, 
            lambda x: x.GetBinContent(ipt+1,ieta+1)))
    print(f'differential sfs: {differential_sfs}')

  make_fakephoton_cpp_evaluator(differential_sfs)

def make_fakephoton_cpp_evaluator(diff_sfs: list[list[list[list[float]]]]):
  """Generates C++ evaluator code for fake photon SFs

  Args:
    differential_sfs: differential SFs indexed by run, pt, eta, and process
  """
  #indexed by process, run, pt/eta
  sfs = [[],[]]
  for irun in range(2):
    sfs[0].append([])
    sfs[1].append([])
    for ipt in range(len(FAKEPH_PT_BINS)-1):
      for ieta in range(len(FAKEPH_ABSETA_BINS)-1):
        sfs[0][-1].append(diff_sfs[irun][ipt][ieta][0])
        sfs[1][-1].append(diff_sfs[irun][ipt][ieta][1])
  print('const std::vector<float> fakeph_ptbins = {',end='')
  for ipt in range(len(FAKEPH_PT_BINS)):
    if ipt != 0:
      print(',',end='')
    print(FAKEPH_PT_BINS[ipt],end='')
  print('};')
  print('const std::vector<float> fakeph_absetabins = {',end='')
  for ieta in range(len(FAKEPH_ABSETA_BINS)):
    if ieta != 0:
      print(',',end='')
    print(FAKEPH_ABSETA_BINS[ieta],end='')
  print('};')
  for irun in range(2):
    for proc_name, iproc in [('jet',0), ('pu',1)]:
      print(
          f'const std::vector<float> fakeph_{proc_name}_sfs_run{irun+2} = {{',
          end='')
      first = True
      for isf in range(len(sfs[iproc][irun])):
        if isf != 0:
          print(',',end='')
        print('{:.8f}'.format(sfs[iproc][irun][isf]),end='')
      print('};')
  print('std::vector<float> fakeph_sfs;')
  print('if (run == 2 && photon_isjet) fakeph_sfs = fakeph_jet_sfs_run2;')
  print('else if (run == 2 && !photon_isjet) fakeph_sfs = fakeph_pu_sfs_run2;')
  print('else if (run == 3 && photon_isjet) fakeph_sfs = fakeph_jet_sfs_run3;')
  print('else if (run == 3 && !photon_isjet) fakeph_sfs = fakeph_pu_sfs_run3;')
  print('for (unsigned ipt = 0; ipt < (fakeph_ptbins.size()-1); ipt++) {')
  print('  for (unsigned ieta = 0; ieta < ',end='')
  print('(fakeph_absetabins.size()-1); ieta++) {')
  print('    if (ph_pt >= fakeph_ptbins[ipt] && ph_pt < ',end='')
  print('fakeph_ptbins[ipt+1] && ph_abseta >= fakeph_absetabins[ieta]',end='')
  print(' && ph_abseta < fakeph_absetabins[ieta+1]) {')
  print('      return fakeph_sfs[ipt*fakeph_absetabins.size()+ieta];')
  print('      break;')
  print('    }')
  print('  }')
  print('}')

def get_fakephoton_quality_weights():
  """Calculates weights for reweighting fake photon quality
  """
  dfs = setup_dataframes()
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, '1', 'data2'),
                ('Z#gamma (run2)', 1, '1', 'mczg2'),
                ('DY+jet (run2)', 2, 'photon_isjet', 'dyjt2'),
                ('DY+PU (run2)', 2, '!photon_isjet', 'dypu2'),
                ('Data (run3)', 3, '1', 'data3'),
                ('Z#gamma (run3)', 4, '1', 'mczg3'),
                ('DY+jet (run3)', 5, 'photon_isjet', 'dyjt3'),
                ('DY+PU (run3)', 5, '!photon_isjet', 'dypu3')]
  FAKEPH_ID_BINS = [0.14,0.5,0.8,0.9,1.0]
  FAKEPH_RES_BINS = [0.0,0.025,0.05,1.0]
  for proc_info in procs_info:
    proc_name = proc_info[0]
    df = dfs[proc_info[1]].Filter(proc_info[2])
    simple_name = proc_info[3]
    hist_ptrs.append(df.Histo2D((
        f'fakeph_{simple_name}',
        f'{proc_name}; Photon IDMVA; Photon #sigma_{{E}}/E',
        len(FAKEPH_ID_BINS)-1,array('d',FAKEPH_ID_BINS),
        len(FAKEPH_RES_BINS)-1,array('d',FAKEPH_RES_BINS)),
        'photon_idmva0','photon_res','w_lumiyear'))

  for run in [0, 4]:
    data_hist = merge_hist_ptrs(hist_ptrs[0+run])
    mczg_hist = merge_hist_ptrs(hist_ptrs[1+run])
    dyjt_hist = merge_hist_ptrs(hist_ptrs[2+run])
    dypu_hist = merge_hist_ptrs(hist_ptrs[3+run])
    data_hist.Add(mczg_hist, -1.0)
    dyjt_hist.Add(dypu_hist)
    data_hist.Divide(dyjt_hist)
    print('Ratios: ',end='')
    for iid in range(len(FAKEPH_ID_BINS)-1):
      for ires in range(len(FAKEPH_RES_BINS)-1):
        print(data_hist.GetBinContent(iid+1, ires+1),end=',')
    print('')

def get_llphpt_weights():
  """Calculates weights for reweighting llph system pt
  """
  dfs = setup_dataframes()
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, '1', 'data2'),
                ('Z#gamma (run2)', 1, '1', 'mczg2'),
                ('DY+jet (run2)', 2, 'photon_isjet', 'dyjt2'),
                ('DY+PU (run2)', 2, '!photon_isjet', 'dypu2'),
                ('Data (run3)', 3, '1', 'data3'),
                ('Z#gamma (run3)', 4, '1', 'mczg3'),
                ('DY+jet (run3)', 5, 'photon_isjet', 'dyjt3'),
                ('DY+PU (run3)', 5, '!photon_isjet', 'dypu3')]
  for proc_info in procs_info:
    hist_ptrs.append([])
    for sel_name, id_sel in [
        ('loid_lodphi','photon_idmva[0]<0.8&&photon_mht_dphi<1.5758'),
        ('loid_hidphi','photon_idmva[0]<0.8&&photon_mht_dphi>1.5758'),
        ('hiid','photon_idmva[0]>=0.9')]:
      proc_name = proc_info[0]
      df = dfs[proc_info[1]].Filter(proc_info[2])
      simple_name = proc_info[3]
      hist_ptrs[-1].append(df.Filter(id_sel).Histo1D((
          f'fakeph_{sel_name}_{simple_name}',
          f'{proc_name}; p_{{T}}^{{ll #gamma}} [GeV]',
          len(LLPH_PT_BINS)-1,array('d',LLPH_PT_BINS)),
          'llphoton_pt0','w_lumiyear'))

  inclusive_sfs = []
  differential_sfs = []
  for run in [0, 4]:
    data_hists = [merge_hist_ptrs(hist_ptrs[0+run][icat]) for icat in 
                   range(3)] 
    simu_hists = [[merge_hist_ptrs(hist_ptrs[iproc+run][icat]) for icat in 
                   range(3)] for iproc in range(1,4)]
    inclusive_sfs.append(solve_hist_system_manual(data_hists, simu_hists, 
                                                  lambda x: x.Integral()))
    print(f'inclusive sfs: {inclusive_sfs}')
    for iproc in range(3):
      for icat in range(3):
        simu_hists[iproc][icat].Scale(inclusive_sfs[-1][iproc])
    differential_sfs.append([])
    for ipt in range(len(LLPH_PT_BINS)-1):
      differential_sfs[-1].append(solve_hist_system_manual(
          data_hists, simu_hists, 
          lambda x: x.GetBinContent(ipt+1)))
    print(f'differential sfs: {differential_sfs}')

  make_llphoton_cpp_evaluator(differential_sfs)

def make_llphoton_cpp_evaluator(diff_sfs: list[list[list[float]]]):
  """Generates C++ evaluator code for fake photon SFs

  Args:
    differential_sfs: differential SFs indexed by run, pt, and process
  """
  #indexed by process, run, pt/eta
  print('static const std::vector<float> llph_ptbins = {',end='')
  for ipt in range(len(LLPH_PT_BINS)):
    if ipt != 0:
      print(',',end='')
    print(LLPH_PT_BINS[ipt],end='')
  print('};')
  for irun in range(2):
    for proc_name, iproc in [('zg',0), ('jt',1), ('pu',2)]:
      print(
          f'const std::vector<float> llph_{proc_name}_sfs_run{irun+2} = {{',
          end='')
      first = True
      for ipt in range(len(diff_sfs[irun])):
        if ipt != 0:
          print(',',end='')
        print('{:.8f}'.format(diff_sfs[irun][ipt][iproc]),end='')
      print('};')
  print('std::vector<float> llph_sfs;')
  print('if (run == 2 && is_dyjet) llph_sfs = llph_jt_sfs_run2;')
  print('else if (run == 2 && is_dypu) llph_sfs = llph_pu_sfs_run2;')
  print('else if (run == 2 && is_zg) llph_sfs = llph_zg_sfs_run2;')
  print('else if (run == 3 && is_dyjet) llph_sfs = llph_jt_sfs_run3;')
  print('else if (run == 3 && is_dypu) llph_sfs = llph_pu_sfs_run3;')
  print('else if (run == 3 && is_zg) llph_sfs = llph_zg_sfs_run3;')
  print('else return 1.0;')
  print('for (unsigned ipt = 0; ipt < (llph_ptbins.size()-1); ipt++) {')
  print('  if (llphoton_pt >= fakeph_ptbins[ipt] && llphoton_pt < ',end='')
  print('llph_ptbins[ipt+1]) {')
  print('      return llph_sfs[ipt];')
  print('      break;')
  print('  }')
  print('}')

def weightloss(coefs: np.array, data: np.array, simu: np.array, 
               reg_coef: float=10.0)-> float:
  """Loss function for coefficient determination

  Args:
    coefs: coefficients as horizontal array
    data: data yields as vertical array
    simu: simu yields as matrix (categories vertical)
    reg_coef: coefficient to penalize (regularize) coefficients far from 1
  """
  coefs_column = coefs[:,np.newaxis]
  data_invsqrt = data**-0.5
  pulls = data_invsqrt*(data-(simu.dot(coefs_column)))
  coef_ones = np.ones(coefs_column.shape)
  coef_penalty = (1.0-np.maximum(coefs_column, coef_ones)
                 /np.minimum(coefs_column, coef_ones))
  mse_loss = np.linalg.norm(pulls)
  coef_loss = reg_coef*np.linalg.norm(coef_penalty)
  loss = mse_loss + coef_loss
  return loss

def get_njetllph_weights():
  """Calculates weights for reweighting njet+llph system pt
  """
  dfs = setup_dataframes()
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, '1', 'data2'),
                ('Z#gamma (run2)', 1, '1', 'mczg2'),
                ('DY+jet (run2)', 2, 'photon_isjet', 'dyjt2'),
                ('DY+PU (run2)', 2, '!photon_isjet', 'dypu2'),
                ('Data (run3)', 3, '1', 'data3'),
                ('Z#gamma (run3)', 4, '1', 'mczg3'),
                ('DY+jet (run3)', 5, 'photon_isjet', 'dyjt3'),
                ('DY+PU (run3)', 5, '!photon_isjet', 'dypu3')]
  for proc_info in procs_info:
    hist_ptrs.append([])
    for sel_name, id_sel in [
        ('loid_lodphi','photon_idmva[0]<0.8&&photon_mht_dphi<1.5758'),
        ('loid_hidphi','photon_idmva[0]<0.8&&photon_mht_dphi>1.5758'),
        ('hiid','photon_idmva[0]>=0.9')]:
      hist_ptrs[-1].append([])
      proc_name = proc_info[0]
      df = dfs[proc_info[1]].Filter(proc_info[2])
      simple_name = proc_info[3]
      hist_ptrs[-1][-1].append(df.Filter(id_sel).Histo1D((
          f'zpt_njet_{sel_name}_{simple_name}',
          f'{proc_name}; N_{{jet}}',
          len(NJET_BINS)-1,array('d',NJET_BINS)),
          'njet','w_lumiyear'))
      for ijet in range(len(NJET_BINS)-1):
        njet_sel = f'njet=={ijet}'
        if ijet==(len(NJET_BINS)-2):
          njet_sel = f'njet>={ijet}'
        hist_ptrs[-1][-1].append(df.Filter(id_sel).Filter(njet_sel).Histo1D((
            f'zpt_llph_nj{ijet}_{sel_name}_{simple_name}',
            f'{proc_name}; p_{{T}}^{{ll #gamma}} [GeV]',
            len(NJ_LLPH_PT_BINS[ijet])-1,array('d',NJ_LLPH_PT_BINS[ijet])),
            'llphoton_pt0','w_lumiyear'))

  simu_coefs = [[[], []], [[],[]], [[],[]]]
  for irun in range(2):
    hists = [[[merge_hist_ptrs(hist_ptrs[iproc+irun*4][icat][ijet]) for icat 
              in range(3)] for iproc in range(4)] for ijet in 
              range(1,6)]

    #loop over bins
    for ijet in range(len(NJET_BINS)-1):
      for ipt in range(len(NJ_LLPH_PT_BINS[ijet])-1):
        data_yields = np.array([[hists[ijet][0][icat].GetBinContent(ipt+1)] 
                                for icat in range(3)])
        simu_yields = np.array([[max(
            hists[ijet][iproc][icat].GetBinContent(ipt+1),0.0) for iproc in 
            range(1,4)] for icat in range(3)])
        #coefs = np.matmul(np.linalg.inv(simu_yields),data_yields)
        # if simple equation solving looks weird, do fit with loss term
        # penalizing SFs far from 1.0
        #if (coefs[0][0] <= 0.0 or coefs[1][0] <= 0.0 or coefs[2][0] <= 0.0):
        x0 = np.array([1.0, 1.0, 1.0])
        min_fn = lambda x: weightloss(x, data_yields, simu_yields, 2.0)
        res = minimize(min_fn, x0, tol=1.0e-6)
        coefs = res.x[:,np.newaxis]
        for iproc in range(3):
          simu_coefs[iproc][irun].append(coefs[iproc][0])
        print(f'bin: {ijet} {ipt}')
        print(f'data yields: {data_yields[0][0]} {data_yields[1][0]}',end='')
        print(f' {data_yields[2][0]}')
        print(f'mczg yields: {simu_yields[0][0]} {simu_yields[1][0]}',end='')
        print(f' {simu_yields[2][0]}')
        print(f'dyjt yields: {simu_yields[0][1]} {simu_yields[1][1]}',end='')
        print(f' {simu_yields[2][1]}')
        print(f'dypu yields: {simu_yields[0][2]} {simu_yields[1][2]}',end='')
        print(f' {simu_yields[2][2]}')
    print(f'mczg sfs: {simu_coefs[0][irun]}')
    print(f'dyjt sfs: {simu_coefs[1][irun]}')
    print(f'dypu sfs: {simu_coefs[2][irun]}')

  #print output as C++ evaluator
  print('static const std::vector<std::vector<float>> llph_ptbins = {')
  for ijet in range(len(NJET_BINS)-1):
    print('    {',end='')
    for ipt in range(len(NJ_LLPH_PT_BINS[ijet])):
      if ipt != 0:
        print(',',end='')
      print(NJ_LLPH_PT_BINS[ijet][ipt],end='')
    print('}',end='')
    if ijet == (len(NJET_BINS)-2):
      print('};')
    else:
      print(',')
  for irun in range(2):
    for proc_name, iproc in [('zg',0), ('jt',1), ('pu',2)]:
      print(
          f'const std::vector<float> llph_{proc_name}_sfs_run{irun+2} = {{',
          end='')
      first = True
      for ipt in range(len(simu_coefs[iproc][irun])):
        if ipt != 0:
          print(',',end='')
        print('{:.8f}'.format(simu_coefs[iproc][irun][ipt]),end='')
      print('};')
  print('std::vector<float> llph_sfs;')
  print('if (run == 2 && is_dyjet) llph_sfs = llph_jt_sfs_run2;')
  print('else if (run == 2 && is_dypu) llph_sfs = llph_pu_sfs_run2;')
  print('else if (run == 2 && is_zg) llph_sfs = llph_zg_sfs_run2;')
  print('else if (run == 3 && is_dyjet) llph_sfs = llph_jt_sfs_run3;')
  print('else if (run == 3 && is_dypu) llph_sfs = llph_pu_sfs_run3;')
  print('else if (run == 3 && is_zg) llph_sfs = llph_zg_sfs_run3;')
  print('else return 1.0;')
  print('for (unsigned ijet = 0; ijet < llph_ptbins.size(); ipt++) {')
  print('  if ((njet==ijet) || (njet>=ijet && ijet==llph_ptbins.size()-1)) {')
  print('    for (unsigned ipt = 0; ipt < ',end='')
  print('(llph_ptbins[ijet].size()-1); ipt++) {')
  print('      if (llphoton_pt >= llph_ptbins[ijet][ipt] && ',end='')
  print('llphoton_pt < llph_ptbins[ijet][ipt+1]) {')
  print('        return llph_sfs[ipt];')
  print('      }')
  print('    }')
  print('  }')
  print('}')

def generate_dnn_datasets():
  """Makes datasets that can be used to train DNN kinematic corrections
  """
  dfs_default = setup_dataframes()
  #merge MC together
  dfs = [dfs_default[0], merge_rdataframesets([dfs_default[1],dfs_default[2]]),
         dfs_default[3], merge_rdataframesets([dfs_default[4],dfs_default[5]])]
  #get parametric weights and apply
  print('Calculating parametric weights.')
  r2_data_hist_ptr = dfs[0].Histo2D(('r2_data_hist',
      'data; IDMVA; #Delta #phi(#gamma, H_{T}^{miss})', 30, 0.14, 1.0, 30, 
      0.0, 3.1416),'photon_idmva0','photon_mht_dphi','w_lumiyear')
  r2_simu_hist_ptr = dfs[1].Histo2D(('r2_simu_hist',
      'data; IDMVA; #Delta #phi(#gamma, H_{T}^{miss})', 30, 0.14, 1.0, 30, 
      0.0, 3.1416),'photon_idmva0','photon_mht_dphi','w_lumiyear')
  r3_data_hist_ptr = dfs[2].Histo2D(('r3_data_hist',
      'data; IDMVA; #Delta #phi(#gamma, H_{T}^{miss})', 30, 0.14, 1.0, 30, 
      0.0, 3.1416),'photon_idmva0','photon_mht_dphi','w_lumiyear')
  r3_simu_hist_ptr = dfs[3].Histo2D(('r3_simu_hist',
      'data; IDMVA; #Delta #phi(#gamma, H_{T}^{miss})', 30, 0.14, 1.0, 30, 
      0.0, 3.1416),'photon_idmva0','photon_mht_dphi','w_lumiyear')
  r2_data_hist = merge_hist_ptrs(r2_data_hist_ptr)
  r2_simu_hist = merge_hist_ptrs(r2_simu_hist_ptr)
  r3_data_hist = merge_hist_ptrs(r3_data_hist_ptr)
  r3_simu_hist = merge_hist_ptrs(r3_simu_hist_ptr)
  r2_data_hist.Scale(1.0/r2_data_hist.Integral())
  r2_simu_hist.Scale(1.0/r2_simu_hist.Integral())
  r3_data_hist.Scale(1.0/r3_data_hist.Integral())
  r3_simu_hist.Scale(1.0/r3_simu_hist.Integral())
  r2_data_hist.Divide(r2_simu_hist)
  r3_data_hist.Divide(r3_simu_hist)
  ROOT.gDirectory.Add(r2_data_hist)
  ROOT.gDirectory.Add(r3_data_hist)
  ROOT.gInterpreter.Declare("""
  TH2D* r2_rw_hist = static_cast<TH2D*>(gDirectory->Get("r2_data_hist"));
  TH2D* r3_rw_hist = static_cast<TH2D*>(gDirectory->Get("r3_data_hist"));
  TRandom3 rng3;
  """)
  dfs[0] = dfs[0].Define('w_param','1')
  dfs[2] = dfs[2].Define('w_param','1')
  dfs[1] = dfs[1].Define('w_param','(r2_rw_hist->GetBinContent(r2_rw_hist->FindBin(photon_idmva0,photon_mht_dphi)))*w_lumiyear')
  dfs[3] = dfs[3].Define('w_param','(r3_rw_hist->GetBinContent(r3_rw_hist->FindBin(photon_idmva0,photon_mht_dphi)))*w_lumiyear')
  dfs[0] = dfs[0].Define('njet0','static_cast<float>(njet)')
  dfs[1] = dfs[1].Define('njet0','static_cast<float>(njet)')
  dfs[2] = dfs[2].Define('njet0','static_cast<float>(njet)')
  dfs[3] = dfs[3].Define('njet0','static_cast<float>(njet)')
  for idf in range(4):
    dfs[idf] = dfs[idf].Define('rng','rng3.Uniform()')
  #define columns and save
  print('Saving output')
  dfs[0] = dfs[0].Define('is_data','1')
  dfs[1] = dfs[1].Define('is_data','0')
  dfs[2] = dfs[2].Define('is_data','1')
  dfs[3] = dfs[3].Define('is_data','0')
  save_columns = ('w_param', 'photon_idmva0', 'photon_mht_dphi', 'njet0', 
                  'rng', 'll_pt0', 'photon_pt0', 'llphoton_pt0', 'mht0', 
                  'ht0', 'is_data')
  write_dataframe_set(dfs[0], 'dataset/kindnn_in_r2_data.root', save_columns)
  write_dataframe_set(dfs[1], 'dataset/kindnn_in_r2_simu.root', save_columns)
  write_dataframe_set(dfs[2], 'dataset/kindnn_in_r3_data.root', save_columns)
  write_dataframe_set(dfs[3], 'dataset/kindnn_in_r3_simu.root', save_columns)

def get_llpjj_weights():
  """Calculates weights for reweighting llpjj system pt
  """
  dfs = setup_dataframes()
  hist_ptrs = []
  procs_info = [('Data (run2)', 0, 'data2'),
                ('Z#gamma (run2)', 1, 'mczg2'),
                ('DY (run2)', 2, 'dyjt2'),
                ('Data (run3)', 3, 'data3'),
                ('Z#gamma (run3)', 4, 'mczg3'),
                ('DY (run3)', 5, 'dyjt3')]
  LLPJJ_PT_BINS = [0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 
                   2.0]
  for proc_info in procs_info:
    hist_ptrs.append([])
    for sel_name, id_sel in [
        ('loid','photon_idmva[0]<0.8'),
        ('hiid','photon_idmva[0]>=0.8')]:
      proc_name = proc_info[0]
      df = dfs[proc_info[1]].Filter('njet>=2')
      simple_name = proc_info[2]
      hist_ptrs[-1].append(df.Filter(id_sel).Histo1D((
          f'llpjj_{simple_name}',
          f'{proc_name}; p_{{T}}^{{ll #gamma jj}} [GeV]',
          len(LLPJJ_PT_BINS)-1,array('d',LLPJJ_PT_BINS)),
          'llphoton_dijet_balance0','w_lumiyear'))

  inclusive_sfs = []
  differential_sfs = []
  for run in [0, 3]:
    data_hists = [merge_hist_ptrs(hist_ptrs[0+run][icat]) for icat in 
                   range(2)] 
    simu_hists = [[merge_hist_ptrs(hist_ptrs[iproc+run][icat]) for icat in 
                   range(2)] for iproc in range(1,3)]
    inclusive_sfs.append(solve_hist_system_manual(data_hists, simu_hists, 
                                                  lambda x: x.Integral()))
    print(f'inclusive sfs: {inclusive_sfs[-1]}')
    for iproc in range(2):
      for icat in range(2):
        simu_hists[iproc][icat].Scale(inclusive_sfs[-1][iproc])
    differential_sfs.append([])
    for ipt in range(len(LLPJJ_PT_BINS)-1):
      differential_sfs[-1].append(solve_hist_system_manual(
          data_hists, simu_hists, 
          lambda x: x.GetBinContent(ipt+1)))
    #print(f'differential sfs: {differential_sfs}')
    for iproc in range(2):
      for ipt in range(len(LLPJJ_PT_BINS)-1):
        print('{:.8f}'.format(differential_sfs[-1][ipt][iproc]),end=',')
      print('')

def solve_hist_system_manual(data_hists: list[ROOT.TH1], 
                             simu_hists: list[list[ROOT.TH1]], 
                             hist_method : Callable[[ROOT.TH1],float]):
  """Solves system of equations to derive scale factors from histograms. 

  Args:
    data_hists: categories of data histograms
    simu_hists: MC histograms indexed by processes, then category
    hist_method: method that extracts a number from a histogram

  Returns:
    scale factors that would scale results from MC histograms to data
  """
  #hist_file = ROOT.TFile('zpt_hists.root', 'READ')
  #n_cats = n_procs to have constrained system of equations
  n_cats = len(data_hists)
  data_yields = np.array([[hist_method(data_hists[icat])] for icat in 
                          range(n_cats)])
  simu_yields = np.array([[hist_method(simu_hists[iproc][icat]) for iproc in 
                           range(n_cats)] for icat in range(n_cats)])
  coefs = np.matmul(np.linalg.inv(simu_yields),data_yields)
  return [coefs[iproc][0] for iproc in range(n_cats)]

def make_hists(dfs: List[RDataFrameSet]):
  """Generates histograms to study Z pt reweighting
  """
  hist_ptrs = []
  names = ['Data (run2)','Z#gamma (run2)',('DY+jet (run2)', 'DY+PU (run2)'),
           'Data (run3)','Z#gamma (run3)',('DY+jet (run3)', 'DY+PU (run3)')]
  for idf in range(len(dfs)):
    for sel_name, id_sel in [
        ('loid_lodphi','photon_idmva[0]<0.8&&photon_mht_dphi<1.5758'),
        ('loid_hidphi','photon_idmva[0]<0.8&&photon_mht_dphi>1.5758'),
        ('hiid','photon_idmva[0]>=0.8')]:
      is_dy = False
      if ((idf == 2) or (idf == 5)):
        is_dy = True
      proc_name = names[idf]
      if is_dy:
        #hist_ptrs.append(dfs[idf].Filter('photon_isjet').Histo1D((f'zpt_{idf}_0',
        #    f'{proc_name[0]}; p_{{T}}^{{Z}} [GeV]',125,0.0,250.0),'z_pt',
        #    'w_lumiyear'))
        #hist_ptrs.append(dfs[idf].Filter('!photon_isjet').Histo1D((
        #    f'zpt_{idf}_1',f'{proc_name[1]}; p_{{T}}^{{Z}} [GeV]',125,0.0,250.0),
        #    'z_pt','w_lumiyear'))
        hist_ptrs.append(dfs[idf].Filter(f'photon_isjet&&{id_sel}').Histo1D((
            f'zgpt_{sel_name}_{idf}_0',
            f'{proc_name[0]}; p_{{T}}^{{ll#gamma}} [GeV]',
            len(LLPH_PT_BINS)-1,array('d',LLPH_PT_BINS)),
            'zg_pt','w_lumiyear'))
        hist_ptrs.append(dfs[idf].Filter(f'!photon_isjet&&{id_sel}').Histo1D((
            f'zgpt_{sel_name}_{idf}_1',
            f'{proc_name[1]}; p_{{T}}^{{ll#gamma}} [GeV]',
            len(LLPH_PT_BINS)-1,array('d',LLPH_PT_BINS)),
            'zg_pt','w_lumiyear'))
        #hist_ptrs.append(dfs[idf].Filter(f'photon_isjet&&{id_sel}').Histo2D((
        #    f'zgpt_{sel_name}_{idf}_0',
        #    f'{proc_name}; p_{{T}}^{{Z}} [GeV]; p_{{T}}^{{#gamma}}',
        #    len(Z_PT_BINS)-1,array('d',Z_PT_BINS),len(P_PT_BINS)-1,
        #    array('d',P_PT_BINS)),'z_pt','ph_pt','w_lumiyear'))
        #hist_ptrs.append(dfs[idf].Filter(f'!photon_isjet&&{id_sel}').Histo2D((
        #    f'zgpt_{sel_name}_{idf}_1',
        #    f'{proc_name}; p_{{T}}^{{Z}} [GeV]; p_{{T}}^{{#gamma}}',
        #    len(Z_PT_BINS)-1,array('d',Z_PT_BINS),len(P_PT_BINS)-1,
        #    array('d',P_PT_BINS)),'z_pt','ph_pt','w_lumiyear'))
      else:
        #hist_ptrs.append(dfs[idf].Histo1D((f'zpt_{idf}',
        #    f'{proc_name}; p_{{T}}^{{Z}} [GeV]',125,0.0,250.0),'z_pt',
        #    'w_lumiyear'))
        #hist_ptrs.append(dfs[idf].Filter(id_sel).Histo2D((
        #    f'zgpt_{sel_name}_{idf}',
        #    f'{proc_name}; p_{{T}}^{{Z}} [GeV]; p_{{T}}^{{#gamma}}',
        #    len(Z_PT_BINS)-1,array('d',Z_PT_BINS),len(P_PT_BINS)-1,
        #    array('d',P_PT_BINS)),'z_pt','ph_pt','w_lumiyear'))
        hist_ptrs.append(dfs[idf].Filter(id_sel).Histo1D((
            f'zgpt_{sel_name}_{idf}',
            f'{proc_name}; p_{{T}}^{{ll#gamma}} [GeV]',
            len(LLPH_PT_BINS)-1,array('d',LLPH_PT_BINS)),
            'zg_pt','w_lumiyear'))

  hist_file = ROOT.TFile('zpt_hists.root', 'RECREATE')
  #plot = RplPlot()
  #plot.lumi_data = [(118,13),(62,13.6)]
  hists = []
  for hist_ptr in hist_ptrs:
    hists.append(merge_hist_ptrs(hist_ptr))
    hists[-1].Write()
    hists[-1].Scale(1.0/hists[-1].Integral())
    #plot.plot_outline(hists[-1])
  #plot.draw('plots/zpt_check.pdf')
  hist_file.Close()

def solve_hist_system(hist_method : Callable[[ROOT.TH1],float]):
  """Solves system of equations to derive scale factors from histograms. 
  Assumes root file with histograms has already been generated

  Args:
    hist_method: method that extracts a number from a histogram

  Returns:
    scale factors that would scale results from MC histograms to data
  """
  hist_file = ROOT.TFile('zpt_hists.root', 'READ')
  sel_names = ['loid_lodphi','loid_hidphi','hiid']
  mc_indices = ['1','2_0','2_1']
  data_yields = np.array([[hist_method(getattr(hist_file,f'zgpt_{sel_name}_0')
      )] for sel_name in sel_names])
  simu_yields = np.array([[hist_method(getattr(hist_file,
      f'zgpt_{sel_name}_{mc_index}')) for mc_index in mc_indices] 
      for sel_name in sel_names])
  coefs = np.matmul(np.linalg.inv(simu_yields),data_yields)
  hist_file.Close()
  return (coefs[0][0], coefs[1][0], coefs[2][0])

def make_hists_ratio():
  """Generates ratio histogram. Assumes root file with histograms has already
  been generated
  """
  norm_mczg, norm_dyjt, norm_dypu = solve_hist_system(lambda x : x.Integral())
  print(f' MC scale factors: {norm_mczg} (ZG), {norm_dyjt} (DY+jet),',
        f' {norm_dypu} (DY+PU')
  hist_file = ROOT.TFile('zpt_hists.root', 'READ')
  for sel_name in ('loid_lodphi','loid_hidphi','hiid'):
    data_hist_run2 = getattr(hist_file,f'zgpt_{sel_name}_0').Clone('data_hist')
    mczg_hist_run2 = getattr(hist_file,f'zgpt_{sel_name}_1').Clone('mczg_hist')
    dyjt_hist_run2 = getattr(hist_file,f'zgpt_{sel_name}_2_0').Clone(
        'dyjt_hist')
    dypu_hist_run2 = getattr(hist_file,f'zgpt_{sel_name}_2_1').Clone(
        'dypu_hist')
    mczg_hist_run2.Scale(norm_mczg)
    dyjt_hist_run2.Scale(norm_dyjt)
    dypu_hist_run2.Scale(norm_dypu)
    dypu_hist_run2.Add(mczg_hist_run2)
    dypu_hist_run2.Add(dyjt_hist_run2)
    dyjt_hist_run2.Add(mczg_hist_run2)
    plot = RplPlot()
    plot.lumi_data = [(118,13)]
    plot.plot_filled(dypu_hist_run2,ROOT.TColor.GetColor('#b9ac70'))
    plot.plot_filled(dyjt_hist_run2,ROOT.TColor.GetColor('#ffa90e'))
    plot.plot_filled(mczg_hist_run2,ROOT.TColor.GetColor('#3f90da'))
    plot.plot_points(data_hist_run2,ROOT.kBlack)
    plot.add_ratio('data_hist','dypu_hist')
    plot.draw(f'plots/zpt_rw_check_{sel_name}.pdf')
  hist_file.Close()

def get_sfs() -> tuple[list[float]]:
  """Gets scale factors for p_{T}^{ll#gamma}. Assumes file with histograms has
  already been generated

  Returns:
    tuple of scale factor lists
  """
  sfs_mczg = []
  sfs_dyjt = []
  sfs_dypu = []
  norm_mczg, norm_dyjt, norm_dypu = solve_hist_system(lambda x : x.Integral())
  for ibin in range(len(LLPH_PT_BINS)-1):
    sf_mczg, sf_dyjt, sf_dypu = solve_hist_system(lambda x : x.GetBinContent(
        ibin+1))
    sfs_mczg.append(sf_mczg/norm_mczg)
    sfs_dyjt.append(sf_dyjt/norm_dyjt)
    sfs_dypu.append(sf_dypu/norm_dypu)
  print(sfs_mczg)
  print(sfs_dyjt)
  print(sfs_dypu)
  return (sfs_mczg, sfs_dyjt, sfs_dypu)

def perform_correlation_studies(dfs: List[RDataFrameSet]):
  """Studies correlation between mismodelled variables

  Args:
    dfs: dataframes for processes to study
  """
  procs = [('Data (run2)','data2',dfs[0]),
           ('Z#gamma (run2)','mczg2',dfs[1]),
           ('DY+jet (run2)','dyjt2',dfs[2].Filter('photon_isjet')),
           ('DY+PU (run2)','dypu2',dfs[2].Filter('!photon_isjet')),
           ('Data (run3)','data3',dfs[3]),
           ('Z#gamma (run3)','mczg3',dfs[4]),
           ('DY+jet (run3)','dyjt3',dfs[5].Filter('photon_isjet')),
           ('DY+PU (run3)','dypu3',dfs[5].Filter('!photon_isjet'))]
  corr_vars = [('photon_pt0','p_{T}^{#gamma} [GeV]',30,0.0,80.0),
               ('photon_abseta0','|#eta^{#gamma}|',30,0.0,2.5),
               ('njet','N_{jet}',5,-0.5,4.5),
               ('llphoton_pt0','p_{T}^{ll#gamma} [GeV]',30,0.0,150.0),
               ('llphoton_abseta0','|#eta^{ll#gamma}| [GeV]',30,0.0,5.0),
               ('llphoton_abscosTheta0','|cos(#Theta)|',30,0.0,1.0),
               ('met','p_{T}^{miss} [GeV]',30,0.0,100.0)]
  corr_vars_jj = corr_vars + [('lead_jet_abseta','|#eta^{j1}|',30,0.0,4.5),
      ('lead_jet_pt','p_{T}^{j1} [GeV]',30,0.0,250.0),
      ('sublead_jet_abseta','|#eta^{j2}|',30,0.0,4.5),
      ('sublead_jet_pt','p_{T}^{j2} [GeV]',30,0.0,125.0),
      ('photon_zeppenfeld0','#Delta #eta(#gamma, Avr(jj))',30,0.0,5.0),
      ('photon_jet1_dr0','#Delta R(#gamma, j1)',30,0.0,6.0),
      ('llphoton_dijet_balance0',
       '(|#sum_{ll#gammajj} #vec{p}_T|)/(#sum_{ll#gammajj} p_{T})',30,0.0,1.0),
      ('llphoton_dijet_absdphi0','#Delta #phi(ll#gamma, jj)',30,0.0,3.1416)]
  corr_var_pairs = [(var_a, var_b) for var_a in corr_vars 
                    for var_b in corr_vars]
  corr_var_pairs_jj = [(var_a, var_b) for var_a in corr_vars_jj 
                       for var_b in corr_vars_jj]
  result_hists = []
  for proc_name_clean, proc_name, proc_df in procs:
    result_hists.append([])
    for region_sel, region_name, region_vars in [
        ('1','inclusive',corr_var_pairs),
        ('njet>=2','dijet',corr_var_pairs_jj)]:
      result_hists[-1].append([])
      region_df = proc_df.Filter(region_sel)
      for var_a, var_b in region_vars:
        result_hists[-1][-1].append(
            (f'{proc_name_clean} {region_name} {var_a[1]} {var_b[1]}',
             region_df.Histo2D((
             f'hist_{proc_name}_{region_name}_{var_a[0]}_{var_b[0]}',
             f'{proc_name_clean};{var_a[1]};{var_b[1]}',var_a[2],var_a[3],
             var_a[4],var_b[2],var_b[3],var_b[4]),var_a[0],var_b[0],
             'w_lumiyear')))
  ROOT.gStyle.SetOptStat(0)
  c = ROOT.TCanvas()
  for iproc in range(len(procs)):
    for iregion, nvars in [(0,len(corr_vars)), (1,len(corr_vars_jj))]:
      result_hist_info = result_hists[iproc][iregion]
      proc_name = procs[iproc][1]
      corr_hist = ROOT.TH2D(f'corr_hist_{proc_name}_{iregion}','Correlation',
                            nvars,0,nvars,nvars,0,nvars)
      chi2_hist = ROOT.TH2D(f'chi2_hist_{proc_name}_{iregion}',
                            'Indpendence #chi^{2}/DOF',nvars,0,nvars,nvars,0,
                            nvars)
      xidx = 0
      yidx = 0
      for result_name, result_hist_ptrs in result_hist_info:
        result_hist = merge_hist_ptrs(result_hist_ptrs)
        corr = result_hist.GetCorrelationFactor()
        chi2 = test_independence(result_hist)
        print(result_name, end=' ')
        print(corr,end=' ')
        print(chi2)
        corr_hist.SetBinContent(xidx+1, yidx+1, corr)
        chi2_hist.SetBinContent(xidx+1, yidx+1, min(chi2,10.0))
        if xidx != (nvars-1):
          xidx += 1
        else:
          yidx += 1
          xidx = 0
      corr_hist.Draw('colz text')
      c.SaveAs(f'plots/corrcheck_{corr_hist.GetName()}.pdf')
      chi2_hist.Draw('colz text')
      c.SaveAs(f'plots/corrcheck_{chi2_hist.GetName()}.pdf')

def test_independence(hist: ROOT.TH2D) -> float:
  """Returns measure of independence between two variables of histogram

  Args:
    hist: histogram to check

  Returns:
    chi2/DOF between histogram and fully independent 
  """
  x_projection = hist.ProjectionX()
  y_projection = hist.ProjectionY()
  norm = hist.Integral()
  indep_hist = hist.Clone(hist.GetName()+'_indep')
  nbinsx = hist.GetNbinsX()
  nbinsy = hist.GetNbinsY()
  for ix in range(1,nbinsx+1):
    for iy in range(1,nbinsy+1):
      indep_hist.SetBinContent(ix, iy,
                               x_projection.GetBinContent(ix)
                               *y_projection.GetBinContent(iy)/norm)
  return get_avr_chi2(hist, indep_hist)

def get_avr_chi2(hist_a: ROOT.TH2D, hist_b: ROOT.TH2D) -> float:
  """Returns chi2/DOF across bins of two histograms, ignoring *flow

  Args:
    hist_a: histogram to compare whose error is considered
    hist_b: histogram to compare

  Returns:
    average of pulls across bins of histograms
  """
  nbinsx = hist_a.GetNbinsX()
  nbinsy = hist_a.GetNbinsY()
  total_pull = 0
  for ix in range(1,nbinsx+1):
    for iy in range(1,nbinsy+1):
      yield_a = hist_a.GetBinContent(ix, iy)
      yield_b = hist_b.GetBinContent(ix, iy)
      error_a = hist_a.GetBinError(ix, iy)
      if error_a > 0.0:
        total_pull += (yield_a-yield_b)**2/error_a**2
  return total_pull/(nbinsx*nbinsy+nbinsx+nbinsy)

def get_jeteta_sfs(dfs: List[RDataFrameSet]):
  """Calculates jet eta scale factors

  Args:
    dfs: dataframes for processes to study
  """
  procs = [('Data (run2)','data2',dfs[0]),
           ('Z#gamma (run2)','mczg2',dfs[1]),
           ('DY+jet (run2)','dyjt2',dfs[2].Filter('photon_isjet')),
           ('DY+PU (run2)','dypu2',dfs[2].Filter('!photon_isjet')),
           ('Data (run3)','data3',dfs[3]),
           ('Z#gamma (run3)','mczg3',dfs[4]),
           ('DY+jet (run3)','dyjt3',dfs[5].Filter('photon_isjet')),
           ('DY+PU (run3)','dypu3',dfs[5].Filter('!photon_isjet'))]
  hist_ptrs = [[],[]]
  for proc_name_clean, proc_name, proc_df in procs:
    filter_df = proc_df.Filter('njet>=2')
    hist_ptrs[0].append(filter_df.Histo1D((f'hist_leadjeteta_{proc_name}',
        f'{proc_name_clean};|#eta^{{j1}}|',15,0.0,4.7),'lead_jet_abseta',
        'w_lumiyear'))
    hist_ptrs[1].append(filter_df.Histo1D((f'hist_subljeteta_{proc_name}',
        f'{proc_name_clean};|#eta^{{j2}}|',15,0.0,4.7),'sublead_jet_abseta',
        'w_lumiyear'))
  eta_bins = [ieta*(4.7/15.0) for ieta in range(16)]
  for run in [0,4]:
    for ijet in range(2):
      data_hist = merge_hist_ptrs(hist_ptrs[ijet][run+0])
      mczg_hist = merge_hist_ptrs(hist_ptrs[ijet][run+1])
      dyjt_hist = merge_hist_ptrs(hist_ptrs[ijet][run+2])
      dypu_hist = merge_hist_ptrs(hist_ptrs[ijet][run+3])
      mczg_hist.Add(dyjt_hist)
      mczg_hist.Add(dypu_hist)
      data_hist.Scale(1.0/data_hist.Integral())
      mczg_hist.Scale(1.0/mczg_hist.Integral())
      data_hist.Divide(mczg_hist)
      sfs = print_hist_content(data_hist)
      fit_sfs(eta_bins, sfs, f'fit_sfs_run{int(run/4+2)}_jet{ijet}', 3)

def print_hist_content(hist: ROOT.TH1D) -> list[float]:
  """Prints bin content of histogram and returns as list

  Args:
    hist: histogram to print

  Returns:
    histogram bin content as a lit
  """
  bin_content = []
  for ix in range(1,hist.GetNbinsX()+1):
    bin_content.append(hist.GetBinContent(ix))
  underflow = hist.GetBinContent(0)
  overflow = hist.GetBinContent(hist.GetNbinsX()+1)
  print(f'Flow: {underflow} {overflow}')
  print(bin_content)
  return bin_content

def fit_sfs(sf_binning: list[float], sfs: list[float], output_name: str, 
            poly_order: int, log_x: bool=False):
  """Performs a fit to scale factors

  Args:
    sf_binning: binning for SFs
    sfs: list of scale factors to fit
    output_name: name of file to save output
    log_x: whether to take logarithm of X axis
  """
  bin_centers = [((sf_binning[ipt]+sf_binning[ipt+1])/2.0) for ipt in 
      range(len(sf_binning)-1)]
  lobound = sf_binning[0]
  hibound = sf_binning[-1]
  if log_x:
    bin_centers = [math.log((sf_binning[ipt]+sf_binning[ipt+1])/2.0) for ipt in 
        range(len(sf_binning)-1)]
    if (sf_binning > 0.0):
      lobound = math.log(sf_binning[0])
    else:
      lobound = 0.0
    hibound = math.log(sf_binning[-1])
  sf_graph = ROOT.TGraph(len(sf_binning)-1,array('d',bin_centers),
                         array('d',sfs))
  fn = None
  n_params = poly_order+1
  if poly_order==1:
    fn = ROOT.TF1('cheby',
        '[0]+[1]*(x)',
        lobound, hibound)
  elif poly_order==2:
    fn = ROOT.TF1('cheby',
        '[0]+[1]*(x)+[2]*(2*x*x-1)',
        lobound, hibound)
  elif poly_order==3:
    fn = ROOT.TF1('cheby',
        '[0]+[1]*(x)+[2]*(2*x*x-1)+[3]*(4*x*x*x-3*x)',
        lobound, hibound)
  elif poly_order==4:
    fn = ROOT.TF1('cheby',
        '[0]+[1]*(x)+[2]*(2*x*x-1)+[3]*(4*x*x*x-3*x)+[4]*(8*pow(x,4)-8*x*x+1)',
        lobound, hibound)
  else:
    raise ValueError('Unimplemented polynomial order')
  fn.SetParameter(0,1.0)
  for ipar in range(1,n_params):
    fn.SetParameter(ipar,0.0)
  sf_graph.Fit('cheby')
  #print(fn.GetChisquare())
  ROOT.gStyle.SetOptStat(0)
  can = ROOT.TCanvas()
  sf_graph.SetMarkerStyle(20)
  sf_graph.SetMarkerSize(2)
  fn.SetLineWidth(2)
  fn.SetLineColor(ROOT.TColor.GetColor('#5790fc'))
  sf_graph.Draw('AP')
  fn.Draw('same')
  can.SaveAs(f'plots/{output_name}.pdf')
  coefs = [fn.GetParameter(ipar) for ipar in range(n_params)]
  print(f'Coefficients: {coefs}')

def get_hists_from_file(run: str):
  """Extracts histograms from file

  Args:
    run: run2 or run3
  """
  if not run in ['run2', 'run3']:
    raise ValueError('unrecognized run')
  hist_names = ['zgpt_loid_0', 'zgpt_loid_1', 'zgpt_loid_2_0', 'zgpt_loid_2_1',
                'zgpt_hiid_0', 'zgpt_hiid_1', 'zgpt_hiid_2_0', 'zgpt_hiid_2_1']
  if run=='run3':
    hist_names = ['zgpt_loid_3', 'zgpt_loid_4', 'zgpt_loid_5_0', 
                  'zgpt_loid_5_1', 'zgpt_hiid_3', 'zgpt_hiid_4', 
                  'zgpt_hiid_5_0', 'zgpt_hiid_5_1']
  hists = []
  hist_file = ROOT.TFile('zpt_hists.root', 'READ')
  for hist_name in hist_names:
    hists.append(getattr(hist_file, hist_name))
    hists[-1].SetDirectory(ROOT.nullptr)
  hist_file.Close()
  return hists

def normalize_hists(hists: List[ROOT.TH2D]):
  """Normalizes MC histograms so yield matches data

  Args:
    hists: histograms as returned from get_hists_from_file, modified in-place
  """
  loid_data_ratio = hists[0].Integral()/(hists[1].Integral()
      +hists[2].Integral()+hists[3].Integral())
  hiid_data_ratio = hists[4].Integral()/(hists[5].Integral()
      +hists[6].Integral()+hists[7].Integral())
  hists[1].Scale(loid_data_ratio)
  hists[2].Scale(loid_data_ratio)
  hists[3].Scale(loid_data_ratio)
  hists[5].Scale(hiid_data_ratio)
  hists[6].Scale(hiid_data_ratio)
  hists[7].Scale(hiid_data_ratio)

def update_zpt_sfs(hists: List[ROOT.TH2D], zsfs_mczg: List[float], 
                   zsfs_mcdy: List[float], psfs_mcdy: List[float], 
                   psfs_mcpu: List[float]):
  """Calculates updated Z pt scale factors

  Args:
    hists: 2D histograms as returned by get_hists_from_file
    zsfs_mczg: Z pt scale factors for DYG events, modified in-place
    zsfs_mczg: Z pt scale factors for DY events, modified in-place
    psfs_mcdy: photon pt scale factors for DY+jet fake photon events
    psfs_mcpu: photon pt scale factors for DY+PU fake photon events
  """
  # (loid, hiid)
  hist_data = (hists[0], hists[4])
  hist_mczg = (hists[1], hists[5])
  hist_mcdy = (hists[2], hists[6])
  hist_mcpu = (hists[3], hists[7])
  hists = [ROOT.TH1D('hist_loid_mczg', 'Z+#gamma; p_{T}^{Z} [GeV]', 
                     len(Z_PT_BINS)-1, array('d',Z_PT_BINS)),
           ROOT.TH1D('hist_loid_mcdy', 'Z+fake #gamma; p_{T}^{Z} [GeV]', 
                     len(Z_PT_BINS)-1, array('d',Z_PT_BINS)),
           ROOT.TH1D('hist_loid_data', 'Data; p_{T}^{Z} [GeV]', 
                     len(Z_PT_BINS)-1, array('d',Z_PT_BINS)),
           ROOT.TH1D('hist_hiid_mczg', 'Z+#gamma; p_{T}^{Z} [GeV]', 
                     len(Z_PT_BINS)-1, array('d',Z_PT_BINS)),
           ROOT.TH1D('hist_hiid_mcdy', 'Z+fake #gamma; p_{T}^{Z} [GeV]', 
                     len(Z_PT_BINS)-1, array('d',Z_PT_BINS)),
           ROOT.TH1D('hist_hiid_data', 'Data; p_{T}^{Z} [GeV]', 
                     len(Z_PT_BINS)-1, array('d',Z_PT_BINS))]
  for izpt in range(len(Z_PT_BINS)-1):
    yield_data = [0.0, 0.0]
    yield_mczg = [0.0, 0.0]
    yield_mcdy = [0.0, 0.0]
    for ippt in range(len(P_PT_BINS)-1):
      for iid in range(2):
        yield_data[iid] += hist_data[iid].GetBinContent(izpt+1, ippt+1)
        yield_mczg[iid] += hist_mczg[iid].GetBinContent(izpt+1, ippt+1)
        yield_mcdy[iid] += (hist_mcdy[iid].GetBinContent(izpt+1, ippt+1)
                            *psfs_mcdy[ippt])
        yield_mcdy[iid] += (hist_mcpu[iid].GetBinContent(izpt+1, ippt+1)
                            *psfs_mcpu[ippt])
    for iid in range(2):
      yield_mczg[iid] *= zsfs_mczg[izpt]
      yield_mcdy[iid] *= zsfs_mcdy[izpt]
      hists[iid*3+0].SetBinContent(izpt+1,yield_mczg[iid])
      hists[iid*3+1].SetBinContent(izpt+1,yield_mczg[iid]+yield_mcdy[iid])
      hists[iid*3+2].SetBinContent(izpt+1,yield_data[iid])
    #print('DEBUG #')
    #print(f'data loid: {yield_data[0]}, data hiid: {yield_data[1]}')
    #print(f'mczg loid: {yield_mczg[0]}, mczg hiid: {yield_mczg[1]}')
    #print(f'mcdy loid: {yield_mcdy[0]}, mcdy hiid: {yield_mcdy[1]}')
    zsfs_mczg[izpt] *= (
        (yield_data[0]*yield_mcdy[1]-yield_data[1]*yield_mcdy[0])/
        (yield_mczg[0]*yield_mcdy[1]-yield_mczg[1]*yield_mcdy[0]))
    zsfs_mcdy[izpt] *= (
        (yield_data[0]*yield_mczg[1]-yield_data[1]*yield_mczg[0])/
        (yield_mcdy[0]*yield_mczg[1]-yield_mcdy[1]*yield_mczg[0]))

  #id_names = ['loid','hiid']
  #for iid in range(2):
  #  id_name = id_names[iid]
  #  plot = RplPlot()
  #  plot.lumi_data = [(118,13)]
  #  plot.plot_filled(hists[iid*3+1],ROOT.TColor.GetColor('#f89c20'))
  #  plot.plot_filled(hists[iid*3+0],ROOT.TColor.GetColor('#5790fc'))
  #  plot.plot_points(hists[iid*3+2],ROOT.kBlack)
  #  plot.add_ratio(f'hist_{id_name}_data',f'hist_{id_name}_mcdy')
  #  plot.draw(f'plots/zpt_rw_{id_name}.pdf')


def update_ppt_sfs(hists: List[ROOT.TH2D], zsfs_mczg: List[float], 
                   zsfs_mcdy: List[float], psfs_mcdy: List[float], 
                   psfs_mcpu: List[float]):
  """Calculates updated photon pt scale factors

  Args:
    hists: 2D histograms as returned by get_hists_from_file
    zsfs_mczg: Z pt scale factors for DYG events
    zsfs_mczg: Z pt scale factors for DY events
    psfs_mcdy: photon scale factors for DY+jet fake photon, modified in-place
    psfs_mcpu: photon scale factors for DY+PU fake photon, modified in-place
  """
  # (loid, hiid)
  hist_data = (hists[0], hists[4])
  hist_mczg = (hists[1], hists[5])
  hist_mcdy = (hists[2], hists[6])
  hist_mcpu = (hists[3], hists[7])
  for ippt in range(len(P_PT_BINS)-1):
    #fold Z+prompt photon into data via subtraction
    yield_data = [0.0, 0.0]
    yield_mcdy = [0.0, 0.0]
    yield_mcpu = [0.0, 0.0]
    for izpt in range(len(Z_PT_BINS)-1):
      for iid in range(2):
        yield_data[iid] += hist_data[iid].GetBinContent(izpt+1, ippt+1)
        yield_data[iid] -= (hist_mczg[iid].GetBinContent(izpt+1, ippt+1)
                            *zsfs_mczg[izpt])
        yield_mcdy[iid] += (hist_mcdy[iid].GetBinContent(izpt+1, ippt+1)
                            *zsfs_mcdy[izpt])
        yield_mcpu[iid] += (hist_mcpu[iid].GetBinContent(izpt+1, ippt+1)
                            *zsfs_mcdy[izpt])
    for iid in range(2):
      yield_mcdy[iid] *= psfs_mcdy[ippt]
      yield_mcpu[iid] *= psfs_mcpu[ippt]
      if (yield_data[iid] < 0.0):
        yield_data[iid] = 0.0
    psfs_mcdy[ippt] *= (
        (yield_data[0]*yield_mcpu[1]-yield_data[1]*yield_mcpu[0])/
        (yield_mcdy[0]*yield_mcpu[1]-yield_mcdy[1]*yield_mcpu[0]))
    psfs_mcpu[ippt] *= (
        (yield_data[0]*yield_mcdy[1]-yield_data[1]*yield_mcdy[0])/
        (yield_mcpu[0]*yield_mcdy[1]-yield_mcpu[1]*yield_mcdy[0]))

def do_sf_fit(run: str): 
  """Performs fit to extract scale factors 

  Args:
    run: run2 or run3
  """
  hists = get_hists_from_file(run)
  normalize_hists(hists)
  zsfs_mczg = [1.0 for iz in range(len(Z_PT_BINS)-1)]
  zsfs_mcdy = [1.0 for iz in range(len(Z_PT_BINS)-1)]
  psfs_mcdy = [1.0 for ip in range(len(P_PT_BINS)-1)]
  psfs_mcpu = [1.0 for ip in range(len(P_PT_BINS)-1)]
  update_zpt_sfs(hists, zsfs_mczg, zsfs_mcdy, psfs_mcdy, psfs_mcpu)
  print('DEBUG 1')
  print(zsfs_mczg)
  print(zsfs_mcdy)
  update_ppt_sfs(hists, zsfs_mczg, zsfs_mcdy, psfs_mcdy, psfs_mcpu)
  print('DEBUG 2')
  print(psfs_mcdy)
  print(psfs_mcpu)
  update_zpt_sfs(hists, zsfs_mczg, zsfs_mcdy, psfs_mcdy, psfs_mcpu)
  print('DEBUG 3')
  print(zsfs_mczg)
  print(zsfs_mcdy)
  update_ppt_sfs(hists, zsfs_mczg, zsfs_mcdy, psfs_mcdy, psfs_mcpu)
  print('DEBUG 4')
  print(psfs_mcdy)
  print(psfs_mcpu)
  update_zpt_sfs(hists, zsfs_mczg, zsfs_mcdy, psfs_mcdy, psfs_mcpu)
  print('DEBUG 5')
  print(zsfs_mczg)
  print(zsfs_mcdy)
  update_ppt_sfs(hists, zsfs_mczg, zsfs_mcdy, psfs_mcdy, psfs_mcpu)
  print('DEBUG 6')
  print(psfs_mcdy)
  print(psfs_mcpu)

def generate_bdt_datasets():
  """Generates datasets for BDT evaluation
  """
  dfs = setup_dataframes(False, True)
  saved_columns = ('photon_mva','photon_res','min_dR','max_dR','pt_mass',
                   'cosTheta','costheta','phi','photon_rapidity','l1_rapidity',
                   'l2_rapidity','w_default','w_model','mllg')
  print('Defining columns')
  for iproc, proc_name in [(0, 'data'), (1, 'mczg'), (2, 'mcdy'), (3, 'hzph')]:
    dfs[iproc] = dfs[iproc].Define('photon_mva','photon_idmva[0]')
    dfs[iproc] = dfs[iproc].Define('min_dR','photon_drmin[0]')
    dfs[iproc] = dfs[iproc].Define('max_dR','photon_drmax[0]')
    dfs[iproc] = dfs[iproc].Define('pt_mass','llphoton_pt[0]/llphoton_m[0]')
    dfs[iproc] = dfs[iproc].Define('cosTheta','llphoton_cosTheta[0]')
    dfs[iproc] = dfs[iproc].Define('costheta','llphoton_costheta[0]')
    dfs[iproc] = dfs[iproc].Define('phi','llphoton_psi[0]')
    #dfs[iproc] = dfs[iproc].Define('photon_res','llphoton_psi[0]')
    dfs[iproc] = dfs[iproc].Define('photon_rapidity','photon_eta[0]')
    dfs[iproc] = dfs[iproc].Define('l1_rapidity','get_l1_rapidity(el_pt,el_eta'
        ',mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)')
    dfs[iproc] = dfs[iproc].Define('l2_rapidity','get_l2_rapidity(el_pt,el_eta'
        ',mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)')
    dfs[iproc] = dfs[iproc].Define('mllg','llphoton_m[0]')
    if iproc>0:
      dfs[iproc] = dfs[iproc].Define('w_default',
         'static_cast<float>(weight_fix*w_year*w_photon_lowpt*w_phshape)')
      dfs[iproc] = dfs[iproc].Define('w_model',
         'static_cast<float>(w_default*w_fakepteta*w_llph_pt*w_llpjj_pt)')
    else:
      dfs[iproc] = dfs[iproc].Define('w_model','1.0f')
      dfs[iproc] = dfs[iproc].Define('w_default','1.0f')

  dfs[0] = dfs[0].Filter('(llphoton_m[0]<120)||(llphoton_m[0]>130)')
  bak_dfs = merge_rdataframesets([dfs[1], dfs[2]])
  print('Writing data output')
  write_dataframe_set(dfs[0], f'ntuples/bdtmodelling_dat.root', saved_columns)
  print('Writing bak output')
  write_dataframe_set(bak_dfs, f'ntuples/bdtmodelling_bak.root', saved_columns)
  print('Writing sig output')
  write_dataframe_set(dfs[3], f'ntuples/bdtmodelling_sig.root', saved_columns)
  print('Shuffling output')
  #shuffle output
  for proc_type in ['dat','sig','bak']:
    subprocess.run(('/data1/jbkim/Linux/el7_v1/bin/python3.9 '
        f'scripts/shuffle_tree.py ntuples/bdtmodelling_{proc_type}.root '
        'tree').split())
    subprocess.run((f'rm ntuples/bdtmodelling_{proc_type}.root').split())

def test_a():
  data_yields = np.array([[190.], [163.], [225.]])
  simu_yields = np.array([[36.51965968,  1.31030893, 27.06481356],
                          [38.222505,    6.33321074, 48.69267428],
                          [66.54406422,  0.0,        62.90017488]])
  x0 = np.array([1.0, 1.0, 1.0])
  min_fn = lambda x: weightloss(x, data_yields, simu_yields, 2.0)
  print(min_fn(np.array([0.99, 1.0, 1.0])))
  print(min_fn(np.array([1.00, 1.0, 1.0])))
  print(min_fn(np.array([1.01, 1.0, 1.0])))
  res = minimize(min_fn, x0, method='BFGS', options={'gtol':1.0e-6, 
                                                     'disp':True})
  coefs = res.x[:,np.newaxis]
  print(coefs)

if __name__=='__main__':

  #ROOT.EnableImplicitMT()

  #dfs = setup_dataframes()
  #perform_correlation_studies(dfs)
  #get_jeteta_sfs(dfs)
  #get_jet_weights(dfs)
  #make_hists(dfs)
  #make_hists_ratio()
  #sfs_mczg, sfs_dyjt, sfs_dypu = get_sfs()
  #fit_sfs(LLPH_PT_BINS, sfs_mczg, 'fit_sfs_mczg', 4, True)
  #fit_sfs(LLPH_PT_BINS, sfs_dyjt, 'fit_sfs_dyjt', 4, True)
  #fit_sfs(LLPH_PT_BINS, sfs_dypu, 'fit_sfs_dypu', 4, True)
  get_fakephoton_weights()
  #get_fakephoton_quality_weights()
  #get_llphpt_weights()
  #get_llpjj_weights()
  #generate_bdt_datasets()
  #get_njetllph_weights()
  #generate_dnn_datasets()

