#!/usr/bin/env python
"""Generates 5D photon reweighting
"""

from argparse import ArgumentParser
from array import array
import math
import ROOT
from rdataframeset import RDataFrameSet, merge_hist_ptrs
from root_plot_lib import RplPlot
from typing import Tuple, List, Union

LUMI = {'2016APV' : 19.51,
        '2016' : 16.80,
        '2017' : 41.48,
        '2018' : 59.83,
        '2022' : 8.17,
        '2022EE' : 27.01,
        '2023' : 17.61,
        '2023BPix' : 9.53}

RUN2_DATA_FILENAMES = [
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleEG_zgdata_llg_nfiles_138.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_58.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/data/merged_zgdata_llg/merged_raw_pico_llg_SingleElectron_zgdata_llg_nfiles_156.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/data/merged_zgdata_llg/merged_raw_pico_llg_SingleMuon_zgdata_llg_nfiles_157.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleEG_zgdata_llg_nfiles_173.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_147.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/data/merged_zgdata_llg/merged_raw_pico_llg_SingleElectron_zgdata_llg_nfiles_245.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/data/merged_zgdata_llg/merged_raw_pico_llg_SingleMuon_zgdata_llg_nfiles_436.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_206.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma_zgdata_llg_nfiles_738.root',
    '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/data/merged_zgdata_llg/merged_raw_pico_llg_SingleMuon_zgdata_llg_nfiles_393.root'
    ]

RUN3_DATA_FILENAMES = ['/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_12.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma_zgdata_llg_nfiles_422.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/data/merged_zgdata_llg/merged_raw_pico_llg_Muon_zgdata_llg_nfiles_206.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/data/merged_zgdata_llg/merged_raw_pico_llg_SingleMuon_zgdata_llg_nfiles_35.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma_zgdata_llg_nfiles_765.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/data/merged_zgdata_llg/merged_raw_pico_llg_Muon_zgdata_llg_nfiles_594.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma0_zgdata_llg_nfiles_306.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma1_zgdata_llg_nfiles_297.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/data/merged_zgdata_llg/merged_raw_pico_llg_Muon0_zgdata_llg_nfiles_328.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/data/merged_zgdata_llg/merged_raw_pico_llg_Muon1_zgdata_llg_nfiles_244.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma0_zgdata_llg_nfiles_157.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/data/merged_zgdata_llg/merged_raw_pico_llg_EGamma1_zgdata_llg_nfiles_134.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/data/merged_zgdata_llg/merged_raw_pico_llg_Muon0_zgdata_llg_nfiles_157.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/data/merged_zgdata_llg/merged_raw_pico_llg_Muon1_zgdata_llg_nfiles_130.root']

SIMU_FILENAMES_2016 = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_198.root'
SIMU_FILENAMES_2017 = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_403.root'
SIMU_FILENAMES_2018 = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_422.root'
SIMU_FILENAMES_2022 = ['/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_49.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_341.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_23.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_28.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_21.root']
SIMU_FILENAMES_2022EE = ['/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_872.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_29.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_33.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_62.root']
SIMU_FILENAMES_2023 = ['/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_33.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_237.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_5.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_21.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_26.root']
SIMU_FILENAMES_2023BPix = ['/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_14.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_119.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_11.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_9.root',
    '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/mc/merged_zgmc_llg/merged_pico_llg_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_28.root']

#bins used for reweighting pt-abseta-npv
RW_PT_BINS = [15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,35.0,40.0,45.0,50.0,70.0,
              500.0]
RW_ABSETA_BINS = [0.0,1.2,1.4,1.5,1.6,1.8,2.1,2.3,2.4,2.45,2.5]
RW_NPV_BINS = [0.0,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,
               35.0,37.5,40.0,42.5,45.0,47.5,50.0,55.0,60.0,200.0]

#bins used for generating 5D reweighting
PT_BINS = [15.0,17.5,20.0,35.0,50.0,500.0]
ABSETA_BINS = [0.0,0.8,1.5,2.0,2.25,2.5]
NPV_BINS = [0,30,50,500]

RW_IDMVA_BINS = [0.14,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.92,0.94,0.96,1.0]
RW_PHRES_BINS = [0.0,0.015,0.02,0.03,0.04,0.05,0.06,0.08,1.0]

def setup_dataframe_helper(df: Union[ROOT.RDataFrame,RDataFrameSet]) \
                           -> Union[ROOT.RDataFrame,RDataFrameSet]:
  """Performs dataframe defines and filters

  Args:
    df: dataframe or dataframe set to perform defines and filters on

  Returns:
    dataframe (set) with defines and filters
  """
  df = df.Filter('nmu==2&&nll==1&&nphoton>=1')
  df = df.Filter('trig_single_mu||trig_double_mu')
  df = df.Filter('ll_m[0]<75')
  df = df.Filter('llphoton_m[0]>84&&llphoton_m[0]<98')
  df = df.Define('ph_pt','photon_pt[0]')
  df = df.Define('ph_abseta','fabs(photon_eta[0])')
  df = df.Define('ph_idmva','photon_idmva[0]')
  df = df.Define('ph_res',
      'photon_energyErr[0]/(photon_pt[0]*TMath::CosH(photon_eta[0]))')
  return df

def setup_dataframes(run3: bool) -> Tuple[ROOT.RDataFrame,RDataFrameSet]:
  """Initializes dataframes for 5D binned photon reweighting

  Args:
    run3: if true uses run 3 samples

  Returns:
    data RDataframe and simulation RDataframeSet
  """
  df_data = None
  df_simu = None
  if not run3:
    df_data = ROOT.RDataFrame('tree',RUN2_DATA_FILENAMES)
    df_data = df_data.Define('w_lumiyear','1')
    df_mc16 = ROOT.RDataFrame('tree',SIMU_FILENAMES_2016)
    df_mc16 = df_mc16.Define('w_year',str(LUMI['2016']))
    df_mc17 = ROOT.RDataFrame('tree',SIMU_FILENAMES_2017)
    df_mc17 = df_mc17.Define('w_year',str(LUMI['2017']))
    df_mc18 = ROOT.RDataFrame('tree',SIMU_FILENAMES_2018)
    df_mc18 = df_mc18.Define('w_year',str(LUMI['2018']))
    df_simu = RDataFrameSet([df_mc16,df_mc17,df_mc18])
  else:
    df_data = ROOT.RDataFrame('tree',RUN3_DATA_FILENAMES)
    df_data = df_data.Define('w_lumiyear','1')
    df_mc22 = ROOT.RDataFrame('tree',SIMU_FILENAMES_2022)
    df_mc22 = df_mc22.Define('w_year',str(LUMI['2022']))
    df_mcee = ROOT.RDataFrame('tree',SIMU_FILENAMES_2022EE)
    df_mcee = df_mcee.Define('w_year',str(LUMI['2022EE']))
    df_mc23 = ROOT.RDataFrame('tree',SIMU_FILENAMES_2023)
    df_mc23 = df_mc23.Define('w_year',str(LUMI['2023']))
    df_mcbp = ROOT.RDataFrame('tree',SIMU_FILENAMES_2023BPix)
    df_mcbp = df_mcbp.Define('w_year',str(LUMI['2023BPix']))
    df_simu = RDataFrameSet([df_mc22,df_mcee,df_mc23,df_mcbp])
  df_simu = df_simu.Define('w_lumiyear','fabs(w_year*w_lumi)')
  df_data = setup_dataframe_helper(df_data)
  df_simu = setup_dataframe_helper(df_simu)
  return df_data, df_simu

def study_npv(df_data: ROOT.RDataFrame, df_simu: RDataFrameSet):
  """Method used to debug strange behavior with PU reweighting

  Args:
    df_data: dataframe for data
    df_simu: dataframe set for simulation
  """
  df_data = df_data.Filter('pass')
  df_simu = df_simu.Filter('pass')
  data_hist = df_data.Histo1D(('data_hist','Data;N_{PV}^{good}',80,0.0,80.0),
                               'npv_good').GetPtr()
  df_simu = df_simu.Define('w_lumiyearpu','w_lumi*w_year*w_pu')
  simu_hist = merge_hist_ptrs(df_simu.Histo1D(('simu_hist',
                                               'MC;N_{PV}^{good}',80,
                                               0.0,80.0),'npv_good',
                                               'w_lumiyear'))
  simu_hist_pu = merge_hist_ptrs(df_simu.Histo1D(('simu_hist_pu',
      'MC (PU weights);N_{PV}^{good}',80,0.0,80.0),'npv_good','w_lumiyear'))
  data_hist.Scale(1.0/data_hist.Integral())
  simu_hist.Scale(1.0/simu_hist.Integral())
  simu_hist_pu.Scale(1.0/simu_hist_pu.Integral())
  plot = RplPlot()
  plot.lumi_data = [(118,13)]
  plot.plot_outline(data_hist)
  plot.plot_outline(simu_hist)
  plot.plot_outline(simu_hist_pu)
  plot.add_ratio('data_hist','simu_hist')
  plot.add_ratio('data_hist','simu_hist_pu')
  plot.draw('plots/npv_check.pdf')

def apply_reweighting(df: ROOT.RDataFrame, 
    hist: Union[ROOT.TH1D, ROOT.TH2D, ROOT.TH3D], 
    columns: List[str], 
    rw_column_name: str) -> Union[ROOT.RDataFrame,RDataFrameSet]:
  """Applies reweighting to a RDataFrame from a THND histogram

  Args:
    df: dataframe to apply reweighting to
    hist: histogram with weights
    columns: branches in dataframe used to do reweighting
    rw_column_name: name of generated weight

  Returns:
    dataframe with reweighting column added
  """
  hist_name = hist.GetName()
  joined_columns = ''
  for column_name in columns:
    joined_columns += (column_name+',')
  joined_columns = joined_columns[:-1]
  hist_dim = len(columns)
  ROOT.gDirectory.Add(hist)
  ROOT.gInterpreter.Declare(f'TH{hist_dim}D* {hist_name} = '+
      f'static_cast<TH{hist_dim}D*>(gDirectory->Get("{hist_name}"));')
  return df.Define(rw_column_name, 
      f'{hist_name}->GetBinContent({hist_name}->FindBin({joined_columns}))')

def study_pteta_reweighting(df_data: ROOT.RDataFrame, 
    df_simu: RDataFrameSet):
  """Produces fine-grained histogram to study pt-eta reweighting

  Args:
    df_data: dataframe for data
    df_simu: dataframe set for simulation
  """
  fine_pt_bins = [15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,70.0,100.0]
  fine_eta_bins = [0.0,0.5,1.0,1.2,1.4,1.6,1.8,2.0,2.1,2.2,2.3,2.4,2.45,2.5]
  hist_pteta_data = df_data.Histo2D(('hist_pteta_rw',
      '; Photon p_{T} [GeV];Photon |#eta|',
      len(fine_pt_bins)-1,array('d',fine_pt_bins),
      len(fine_eta_bins)-1,array('d',fine_eta_bins)),
      'ph_pt','ph_abseta','w_lumiyear').GetPtr()
  hist_pteta_simu = merge_hist_ptrs(df_simu.Histo2D(('','',
      len(fine_pt_bins)-1,array('d',fine_pt_bins),
      len(fine_eta_bins)-1,array('d',fine_eta_bins)),
      'ph_pt','ph_abseta','w_lumiyear'))
  hist_pteta_data.Scale(1.0/hist_pteta_data.Integral())
  hist_pteta_simu.Scale(1.0/hist_pteta_simu.Integral())
  hist_pteta_data.Divide(hist_pteta_simu)
  plot = RplPlot()
  plot.lumi_data = [(118,13)]
  plot.z_title = 'Data/MC'
  plot.z_min = 0.75
  plot.z_max = 1.25
  plot.plot_colormap(hist_pteta_data)
  plot.draw('plots/mumuphoton_finepteta.pdf')

def perform_pteta_reweighting(df_data: ROOT.RDataFrame, 
    df_simu: RDataFrameSet) -> Tuple[ROOT.RDataFrame,RDataFrameSet,ROOT.TH2D]:
  """Generates branches than can reweight simulated dataframe to data in pT/eta

  Args:
    df_data: dataframe for data
    df_simu: dataframe set for simulation

  Returns:
    data dataframe and simulation dataframe with w_pteta column defined as well
    as histogram to prevent python garbage collecting it
  """
  hist_pteta_data = df_data.Histo2D(('hist_pteta_rw','',len(RW_PT_BINS)-1, 
      array('d',RW_PT_BINS), len(RW_ABSETA_BINS)-1, 
      array('d', RW_ABSETA_BINS)), 'ph_pt','ph_abseta','w_lumiyear').GetPtr()
  hist_pteta_simu = merge_hist_ptrs(df_simu.Histo2D(('','',len(RW_PT_BINS)-1, 
      array('d',RW_PT_BINS), len(RW_ABSETA_BINS)-1, 
      array('d', RW_ABSETA_BINS)), 'ph_pt','ph_abseta','w_lumiyear'))
  hist_pteta_data.Scale(1.0/hist_pteta_data.Integral())
  hist_pteta_simu.Scale(1.0/hist_pteta_simu.Integral())
  hist_pteta_data.Divide(hist_pteta_simu)
  df_simu = apply_reweighting(df_simu, hist_pteta_data, ['ph_pt','ph_abseta'],
                              'w_pteta')
  return df_data, df_simu, hist_pteta_data

def perform_npv_reweighting(df_data: ROOT.RDataFrame, 
    df_simu: RDataFrameSet) -> Tuple[ROOT.RDataFrame,RDataFrameSet,ROOT.TH2D]:
  """Generates branches than can reweight simulated dataframe to data in PU

  Args:
    df_data: dataframe for data
    df_simu: dataframe set for simulation

  Returns:
    data dataframe and simulation dataframe with w_pu column defined as well
    as histogram to prevent python garbage collecting it
  """
  hist_npv_data = df_data.Histo1D(('hist_npv_rw','',len(RW_NPV_BINS)-1, 
      array('d',RW_NPV_BINS)), 'npv', 'w_lumiyear').GetPtr()
  hist_npv_simu = merge_hist_ptrs(df_simu.Histo1D(('','',len(RW_NPV_BINS)-1, 
      array('d',RW_NPV_BINS)), 'npv', 'w_lumiyear'))
  hist_npv_data.Scale(1.0/hist_npv_data.Integral())
  hist_npv_simu.Scale(1.0/hist_npv_simu.Integral())
  hist_npv_data.Divide(hist_npv_simu)
  df_simu = apply_reweighting(df_simu, hist_npv_data, ['npv'] ,'w_npv')
  return df_data, df_simu, hist_npv_data

def make_rw_bins() -> List[str]:
  """Generates binning for 2D weights
  """
  bins = []
  for ipt in range(len(PT_BINS)-1):
    pt_loedge = PT_BINS[ipt]
    pt_hiedge = PT_BINS[ipt+1]
    for ieta in range(len(ABSETA_BINS)-1):
      eta_loedge = ABSETA_BINS[ieta]
      eta_hiedge = ABSETA_BINS[ieta+1]
      #for inpv in range(len(NPV_BINS)-1):
      #  npv_loedge = NPV_BINS[inpv]
      #  npv_hiedge = NPV_BINS[inpv+1]
      bins.append((f'ph_pt>{pt_loedge}&&ph_pt<{pt_hiedge}'
                  +f'&&ph_abseta>{eta_loedge}&&ph_abseta<{eta_hiedge}'))
                  #+f'&&npv>{npv_loedge}&&npv<{npv_hiedge}'))
  return bins


def generate_shape_weights(df_data: ROOT.RDataFrame, 
                           df_simu: RDataFrameSet, 
                           bins: List[str]) -> List[ROOT.TH2D]:
  """Generates histograms that can be used to reweight IDMVA and photon res.
  
  Args:
    df_data: dataframe for data
    df_simu: dataframe for simulation
    bins: selections defining categories for 2D histograms

  Returns:
    list of histograms, one per bin
  """
  data_hist_ptrs = []
  simu_hist_ptrs = []
  ibin = 0
  for rw_bin in bins:
    df_data_bin = df_data.Filter(rw_bin)
    df_simu_bin = df_simu.Filter(rw_bin)
    data_hist_ptrs.append(df_data_bin.Histo2D((f'hist_rw_bin{ibin}','',
      len(RW_IDMVA_BINS)-1,array('d',RW_IDMVA_BINS),len(RW_PHRES_BINS)-1,
      array('d',RW_PHRES_BINS)),'ph_idmva','ph_res','w_all'))
    simu_hist_ptrs.append(df_simu_bin.Histo2D((f'hist_simu_bin{ibin}','',
      len(RW_IDMVA_BINS)-1,array('d',RW_IDMVA_BINS),len(RW_PHRES_BINS)-1,
      array('d',RW_PHRES_BINS)),'ph_idmva','ph_res','w_all'))
    ibin += 1
  hists = []
  for irw in range(len(bins)):
    data_hist = data_hist_ptrs[irw].GetPtr()
    simu_hist = merge_hist_ptrs(simu_hist_ptrs[irw])
    data_hist.Scale(1.0/data_hist.Integral())
    simu_hist.Scale(1.0/simu_hist.Integral())
    data_hist.Divide(simu_hist)
    hists.append(data_hist)
  return hists

def save_hists(hists: List[ROOT.TH2D], filename: str):
  """Saves a list of histrograms (or any other ROOT object) to a TFile

  Args:
    hists: list of objects to save
    filename: name of file to save them to
  """
  hist_file = ROOT.TFile(filename, 'RECREATE')
  for hist in hists:
    hist.Write()
  hist_file.Close()

def make_cpp_evaluator(hists: List[ROOT.TH2D], filename: str):
  """Saves reweight evaluation code to a C++ file

  Args:
    hists: histograms storing reweight values
    filename: name of file to save them to
  """
  with open(filename, 'w') as cpp_file:
    cpp_file.write('float get_photon_weight(float pt, float abseta, ')
    cpp_file.write('float idmva, float res) {\n')
    ibin = 0
    for ipt in range(len(PT_BINS)-1):
      pt_loedge = PT_BINS[ipt]
      pt_hiedge = PT_BINS[ipt+1]
      cpp_file.write('  ')
      if ipt != 0:
        cpp_file.write('else ')
      cpp_file.write(f'if (pt>={pt_loedge} && pt<{pt_hiedge}) '+'{\n')
      for ieta in range(len(ABSETA_BINS)-1):
        eta_loedge = ABSETA_BINS[ieta]
        eta_hiedge = ABSETA_BINS[ieta+1]
        cpp_file.write('    ')
        if ieta != 0:
          cpp_file.write('else ')
        cpp_file.write(f'if (abseta>={eta_loedge} && abseta<{eta_hiedge}) ')
        cpp_file.write('{\n')
        #for inpv in range(len(NPV_BINS)-1):
        #  npv_loedge = NPV_BINS[inpv]
        #  npv_hiedge = NPV_BINS[inpv+1]
        #  cpp_file.write('      ')
        #  if inpv != 0:
        #    cpp_file.write('else ')
        #  cpp_file.write(f'if (npv>={npv_loedge} && npv<{npv_hiedge}) '+'{\n')
        hist = hists[ibin]
        for iidmva in range(len(RW_IDMVA_BINS)-1):
          idmva_loedge = RW_IDMVA_BINS[iidmva]
          idmva_hiedge = RW_IDMVA_BINS[iidmva+1]
          cpp_file.write('      ')
          if iidmva != 0:
            cpp_file.write('else ')
          cpp_file.write(f'if (idmva>={idmva_loedge}')
          cpp_file.write(f' && idmva<{idmva_hiedge}) '+'{\n')
          for ires in range(len(RW_PHRES_BINS)-1):
            res_loedge = RW_PHRES_BINS[ires]
            res_hiedge = RW_PHRES_BINS[ires+1]
            correction = hist.GetBinContent(iidmva+1,ires+1)
            cpp_file.write('        ')
            if ires != 0:
              cpp_file.write('else ')
            cpp_file.write(f'if (res>={res_loedge} && res<{res_hiedge})')
            cpp_file.write(f' return {correction};'+'\n')
          cpp_file.write('      }\n')
        ibin += 1
        cpp_file.write('    }\n')
        #cpp_file.write('    }\n')
      cpp_file.write('  }\n')
    cpp_file.write('  return 0.0;\n')
    cpp_file.write('}')

def save_ntuples_for_dnn(df_data: ROOT.RDataFrame, df_simu: RDataFrameSet,
                         output_filename: str, run3: bool=False):
  """Saves ntuples in a format that can be used to train DNN

  Args:
    df_data: dataframe for data
    df_simu: dataframes for simulation
    output_filename: base filename for output
  """
  if output_filename[-5:]=='.root':
    output_filename = output_filename[:-5]
  ROOT.gInterpreter.Declare('TRandom3 rng;')
  df_data = df_data.Define('rng','rng.Integer(10)')
  df_simu = df_simu.Define('rng','rng.Integer(10)')
  save_columns = ['ph_pt','ph_abseta','npv','ph_idmva','ph_res','w_all','rng']
  df_data.Snapshot('tree',output_filename+'_data.root',save_columns)
  if not run3:
    df_simu.dataframes[0].Snapshot('tree',output_filename+'_simu16.root',
                                   save_columns)
    df_simu.dataframes[1].Snapshot('tree',output_filename+'_simu17.root',
                                   save_columns)
    df_simu.dataframes[2].Snapshot('tree',output_filename+'_simu18.root',
                                   save_columns)
  else:
    df_simu.dataframes[0].Snapshot('tree',output_filename+'_simu22.root',
                                   save_columns)
    df_simu.dataframes[1].Snapshot('tree',output_filename+'_simu22EE.root',
                                   save_columns)
    df_simu.dataframes[2].Snapshot('tree',output_filename+'_simu23.root',
                                   save_columns)
    df_simu.dataframes[3].Snapshot('tree',output_filename+'_simu23BPix.root',
                                   save_columns)
  df_simu_merged = ROOT.RDataFrame('tree',output_filename+'_simu*.root')
  df_simu_merged.Snapshot('tree',output_filename+'_simu.root')

def add_binnedweight_to_ntuple(input_filename, output_filename, cpp_filename):
  """Adds column with binned reweighting to an n-tuple

  Args:
    input_filename: filename of input ROOT file with tree tree
    output_filename: filename of output ROOT file
    cpp_filename: filename of C++ evaluator
  """
  #ROOT.gInterpreter.Load('binned_ph_reweight.cpp')
  ROOT.gInterpreter.ProcessLine(f'#include "{cpp_filename}"')
  df = ROOT.RDataFrame('tree',input_filename)
  df = df.Define('w_binned',
      'get_photon_weight(ph_pt, ph_abseta, ph_idmva, ph_res)')
  #df = df.Define('w_binned',
  #    'ph_pt+ph_abseta+static_cast<float>(npv)+ph_idmva+ph_res')
  df.Snapshot('tree',output_filename)

if __name__=='__main__':

  ROOT.EnableImplicitMT()

  #df_data, df_simu = setup_dataframes(True)
  ####study_npv(df_data, df_simu)
  ####study_pteta_reweighting(df_data, df_simu)
  #df_data, df_simu, hist_pteta_rw = perform_pteta_reweighting(df_data, df_simu)
  ####df_data, df_simu, hist_npv_rw = perform_npv_reweighting(df_data, df_simu)
  #df_data = df_data.Define('w_all','w_lumiyear')
  #df_simu = df_simu.Define('w_all','w_lumiyear*w_pteta')
  #save_ntuples_for_dnn(df_data, df_simu, 'dataset/mmp_weighttrain3.root')
  #reweight_hists = generate_shape_weights(df_data, df_simu, make_rw_bins())
  #save_hists(reweight_hists, 'binned_ph_reweight4.root')
  #make_cpp_evaluator(reweight_hists, 'binned_ph_reweight4.cpp')
     
  add_binnedweight_to_ntuple('dataset/mmp_weighttrain3_simu.root',
                             'dataset/mmp_weighttrain3_simuw.root',
                             'binned_ph_reweight4.cpp')
