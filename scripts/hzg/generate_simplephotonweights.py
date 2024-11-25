#!/usr/bin/env python3
"""@package docstring
Generates photon preselection corrections using tag-and-probe
"""

from array import array
from argparse import ArgumentParser
from correctionlib import schemav2
from correctionlib.convert import ndpolyfit
from root_plot_lib import RplPlot
from math import hypot
import numpy as np
import ROOT
import json

#constants
#PT_BINS = [15.0,20.0,35.0,50.0,80.0,500.0]
PT_BINS = [15.0,20.0,35.0,55.0,500.0]
ETA_BINS = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
ABSETA_BINS = [0.0,0.8,1.5,2.0,2.5]
MVA_BINS_EB = [0.42,0.5,0.6,0.7,0.733,0.767,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0]
MVA_BINS_EE = [0.14,0.28,0.42,0.5,0.6,0.7,0.75,0.8,0.825,0.85,0.875,0.9,0.92,0.94,0.96,0.98,1.0]
NMVA_BINS = 25
MVA_HI = 1.0
MVA_LO = 0.14
#MVA_BINS = [i*((MVA_HI-MVA_LO)/NMVA_BINS)+MVA_LO for i in range(NMVA_BINS+1)]
LUMI = {'2016APV' : 19.51,
        '2016' : 16.80,
        '2017' : 41.48,
        '2018' : 59.83,
        '2022' : 8.17,
        '2022EE' : 27.01,
        '2023' : 17.61,
        '2023BPix' : 9.53}
ENERGY = {'2016APV' : 13,
          '2016' : 13,
          '2017' : 13,
          '2018' : 13,
          '2022' : 13.6,
          '2022EE' : 13.6,
          '2023' : 13.6,
          '2023BPix' : 13.6}

def bin_per_unit(hist, bins):
  '''Modifies histogram contents so its yields are per unit of measured 
  variable (divides each bin content by the bin width)

  hist  TH1 to modify
  bins  list of floats specifying histogram binning
  '''
  for ibin in range(1,len(bins)):
    bin_width = (bins[ibin]-bins[ibin-1])
    hist.SetBinContent(ibin, hist.GetBinContent(ibin)/bin_width)
    hist.SetBinError(ibin, hist.GetBinError(ibin)/bin_width)

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

def make_correction(name, eff_values, unc_values):
  '''Generates correctionlib corrections with given name and values
  '''
  effs = schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','eta','mva'],
          edges=[ETA_BINS,PT_BINS,MVA_BINS],
          content=eff_values,
          flow='clamp',
          )
  uncs = schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','eta','mva'],
          edges=[ETA_BINS,PT_BINS,MVA_BINS],
          content=unc_values,
          flow='clamp',
          )
  return schemav2.Correction(
      name=name,
      version=1,
      inputs=[schemav2.Variable(name='value', type='string', description='sf/unc'),
              schemav2.Variable(name='pt', type='real', description='photon pt'),
              schemav2.Variable(name='eta', type='real', description='photon eta'),
              schemav2.Variable(name='mva', type='real', description='photon IMDVA')],
      output=schemav2.Variable(name='sf', type='real', description=name),
      data=schemav2.Category(
          nodetype='category',
          input='value',
          content=[
              schemav2.CategoryItem(
                  key='sf',
                  value=effs,
              ),
              schemav2.CategoryItem(
                  key='unc',
                  value=uncs,
              ),
          ],
      ))

def get_filename(year, sample_type):
  '''Returns pico filename for given year and sample type

  year         string, data taking era
  sample_type  string, "data", "zg", or "dy"
  '''
  if sample_type=='data':
    if year == '2016APV':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_109.root'
    if year == '2016':
      #return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon*.root'
      return '/data2/oshiro/ntuples/2016/merged_raw_pico_llg_DoubleMuon*.root'
    if year == '2017':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2017/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon*.root'
    if year == '2018':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon*.root'

  if sample_type=='zg':
    if year == '2016APV':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_28.root'
    if year == '2016':
      #return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL*.root'
      return '/data2/oshiro/ntuples/2016/merged_pico_llg_ZGToLLG*.root'
    if year == '2017':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2017/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL*.root'
    if year == '2018':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL*.root'

  if sample_type=='dy':
    if year == '2016APV':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_62.root'
    if year == '2016':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX*.root'
    if year == '2017':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2017/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX*.root'
    if year == '2018':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX*.root'

  raise ValueError('Unknown sample')

def generate_validation_histogram1(year, sample_type):
  '''Generates validation plots

  year         string, data taking era
  sample_type  string, "data", "zg", or "dy"
  '''
  if not sample_type in ['data', 'zg', 'dy']:
    raise ValueError('Unsupported sample type.')
  filename = get_filename(year, sample_type)
  df = ROOT.RDataFrame('tree', filename)
  df = df.Filter('nmu>=2&&nphoton>=1&&trig_double_mu')
  df = df.Define('lead_photon_idmva','photon_idmva[0]')
  title = 'Data'
  if sample_type == 'dy':
    df = df.Filter('photon_pflavor[0]!=1')
  if sample_type == 'zg':
    title = 'MC'
  if sample_type != 'data':
    df = df.Filter('use_event')
    #df = df.Define('w_lumi_year','w_lumi*{}'.format(LUMI[year]))
    df = df.Define('w_lumi_year','1')
  else:
    df = df.Define('w_lumi_year','1')
  df_onz = df.Filter('llphoton_m[0]>80&&llphoton_m[0]<100')
  df_sideband = df.Filter('llphoton_m[0]>100&&llphoton_m[0]<115')

  mva_bins = [0.0,0.5,0.75,1.0]
  onz_hists = []
  sideband_hists = []
  for imva in range(len(mva_bins)-1):
    filter_string = ('photon_idmva[0]>{}&&photon_idmva[0]<{}').format(
                     mva_bins[imva],mva_bins[imva+1])
    df_onz_bin = df_onz.Filter(filter_string)
    df_sideband_bin = df_sideband.Filter(filter_string)
    onz_hists.append(df_onz_bin.Histo2D((
        '{}_onz_hist_bin{}'.format(sample_type,imva),
        ';Photon #eta;Photon p_{T} [GeV]',
        10,-2.4,2.4,10,0.0,100.0),'photon_eta','photon_pt','w_lumi_year'))
    sideband_hists.append(df_sideband_bin.Histo2D(
        ('{}_sideband_hist_bin{}'.format(sample_type,imva),
        ';Photon #eta;Photon p_{T} [GeV]',10,-2.4,2.4,10,0.0,100.0),
        'photon_eta','photon_pt','w_lumi_year'))

  for imva in range(len(mva_bins)-1):
    onz_plot = RplPlot()
    onz_plot.lumi_data = [(ENERGY[year],LUMI[year])]
    onz_plot.z_title = 'Events/bin'
    onz_plot.plot_colormap(onz_hists[imva])
    onz_plot.draw('plots/photon_shape2d_onz_{}_{}_bin{}.pdf'.format(
                  sample_type,year,imva))
    sdb_plot = RplPlot()
    sdb_plot.lumi_data = [(ENERGY[year],LUMI[year])]
    sdb_plot.z_title = 'Events/bin'
    sdb_plot.plot_colormap(sideband_hists[imva])
    sdb_plot.draw('plots/photon_shape2d_sdb_{}_{}_bin{}.pdf'.format(
                  sample_type,year,imva))

  return onz_hists, sideband_hists

def generate_validation_histograms1(year):
  generate_validation_histogram1(year, 'data')
  generate_validation_histogram1(year, 'zg')
  generate_validation_histogram1(year, 'dy')

def generate_validation_histogram2(year, sample_type):
  '''Generates validation plots

  year         string, data taking era
  sample_type  string, "data", "zg", or "dy"
  '''
  if not sample_type in ['data', 'zg', 'dy']:
    raise ValueError('Unsupported sample type.')
  filename = get_filename(year, sample_type)
  df = ROOT.RDataFrame('tree', filename)
  df = df.Filter('nmu>=2&&nphoton>=1&&trig_double_mu')
  df = df.Define('lead_photon_idmva','photon_idmva[0]')
  onz_title = 'Data 80<m_{ll#gamma}<100 GeV'
  sdb_title = 'Data 80<m_{ll#gamma}<115 GeV'
  if sample_type == 'dy':
    df = df.Filter('photon_pflavor[0]!=1')
  if sample_type == 'zg':
    onz_title = 'MC'
    sdb_title = 'MC'
  if sample_type != 'data':
    df = df.Filter('use_event')
    df = df.Define('w_lumi_year','w_lumi*{}'.format(LUMI[year]))
  else:
    df = df.Define('w_lumi_year','1')
  df_onz = df.Filter('llphoton_m[0]>80&&llphoton_m[0]<100')
  df_sideband = df.Filter('llphoton_m[0]>80&&llphoton_m[0]<115')

  onz_hist_ptrs = []
  sdb_hist_ptrs = []
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ABSETA_BINS)-1):
      mva_bins = MVA_BINS_EE
      if (ABSETA_BINS[ieta+1] <= 1.567):
        mva_bins = MVA_BINS_EB
      filter_string = ('fabs(photon_eta[0])>{}&&fabs(photon_eta[0])<{}'
                       +'&&photon_pt[0]>{}&&photon_pt[0]<{}').format(
                       ABSETA_BINS[ieta],ABSETA_BINS[ieta+1],PT_BINS[ipt],
                       PT_BINS[ipt+1])
      df_onz_bin = df_onz.Filter(filter_string)
      df_sideband_bin = df_sideband.Filter(filter_string)
      onz_hist_ptrs.append(df_onz_bin.Histo1D((
          '{}_onz_hist_{}_bin{}'.format(sample_type,year,ibin),
          '{};IDMVA'.format(onz_title),len(mva_bins)-1,array('d',mva_bins)),
          'lead_photon_idmva','w_lumi_year'))
      sdb_hist_ptrs.append(df_sideband_bin.Histo1D(
          ('{}_sdb_hist_{}_bin{}'.format(sample_type,year,ibin),
          '{};IDMVA'.format(sdb_title),len(mva_bins)-1,
          array('d',mva_bins)),'lead_photon_idmva','w_lumi_year'))
      ibin += 1

  return onz_hist_ptrs, sdb_hist_ptrs

def plot_yearcomp():
  dt_hists_onz_2016, dt_hists_sdb_2016 = generate_validation_histogram2('2016', 'data')
  zg_hists_onz_2016, zg_hists_sdb_2016 = generate_validation_histogram2('2016', 'zg')
  dt_hists_onz_2017, dt_hists_sdb_2017 = generate_validation_histogram2('2017', 'data')
  zg_hists_onz_2017, zg_hists_sdb_2017 = generate_validation_histogram2('2017', 'zg')
  dt_hists_onz_2018, dt_hists_sdb_2018 = generate_validation_histogram2('2018', 'data')
  zg_hists_onz_2018, zg_hists_sdb_2018 = generate_validation_histogram2('2018', 'zg')
  print('Performing event loop, this may take a while.')
  #cross-year comparison
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ABSETA_BINS)-1):
      name_2016 = 'data_onz_hist_2016_bin{}'.format(ibin)
      name_2017 = 'data_onz_hist_2017_bin{}'.format(ibin)
      name_2018 = 'data_onz_hist_2018_bin{}'.format(ibin)
      dt_onz_hist_2016 = dt_hists_onz_2016[ibin].GetPtr()
      zg_onz_hist_2016 = zg_hists_onz_2016[ibin].GetPtr()
      dt_onz_hist_2017 = dt_hists_onz_2017[ibin].GetPtr()
      zg_onz_hist_2017 = zg_hists_onz_2017[ibin].GetPtr()
      dt_onz_hist_2018 = dt_hists_onz_2018[ibin].GetPtr()
      zg_onz_hist_2018 = zg_hists_onz_2018[ibin].GetPtr()
      dt_onz_hist_2016.SetTitle('2016;IDMVA')
      dt_onz_hist_2017.SetTitle('2017;IDMVA')
      dt_onz_hist_2018.SetTitle('2018;IDMVA')
      dt_onz_hist_2016.Sumw2()
      dt_onz_hist_2017.Sumw2()
      dt_onz_hist_2018.Sumw2()
      dt_onz_hist_2016.Scale(1.0/dt_onz_hist_2016.Integral())
      zg_onz_hist_2016.Scale(1.0/zg_onz_hist_2016.Integral())
      dt_onz_hist_2017.Scale(1.0/dt_onz_hist_2017.Integral())
      zg_onz_hist_2017.Scale(1.0/zg_onz_hist_2017.Integral())
      dt_onz_hist_2018.Scale(1.0/dt_onz_hist_2018.Integral())
      zg_onz_hist_2018.Scale(1.0/zg_onz_hist_2018.Integral())
      bin_per_unit(dt_onz_hist_2016, MVA_BINS)
      bin_per_unit(zg_onz_hist_2016, MVA_BINS)
      bin_per_unit(dt_onz_hist_2017, MVA_BINS)
      bin_per_unit(zg_onz_hist_2017, MVA_BINS)
      bin_per_unit(dt_onz_hist_2018, MVA_BINS)
      bin_per_unit(zg_onz_hist_2018, MVA_BINS)
      dt_onz_hist_2016.Divide(zg_onz_hist_2016)
      dt_onz_hist_2017.Divide(zg_onz_hist_2017)
      dt_onz_hist_2018.Divide(zg_onz_hist_2018)
      var_plot = RplPlot()
      var_plot.lumi_data = [('118','13')]
      var_plot.y_title = 'Data/MC'
      var_plot.y_title_lower = 'Year/2016'
      var_plot.plot_outline(dt_onz_hist_2016,error=True)
      var_plot.plot_outline(dt_onz_hist_2017,error=True)
      var_plot.plot_outline(dt_onz_hist_2018,error=True)
      #var_plot.add_ratio(name_2017, name_2016, True)
      #var_plot.add_ratio(name_2018, name_2016, True)
      var_plot.draw('plots/photon_yearcomp_bin{}.pdf'.format(ibin))
      ibin += 1

def generate_validation_histograms2(year):
  dt_hists_onz, dt_hists_sdb = generate_validation_histogram2(year, 'data')
  zg_hists_onz, zg_hists_sdb = generate_validation_histogram2(year, 'zg')
  print('Performing event loop, this may take a while.')
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ABSETA_BINS)-1):
      dt_onz_name = 'data_onz_hist_bin{}'.format(ibin)
      dt_sdb_name = 'data_sdb_hist_bin{}'.format(ibin)
      zg_onz_name = 'zg_onz_hist_bin{}'.format(ibin)
      dt_onz_hist = dt_hists_onz[ibin].GetPtr()
      dt_sdb_hist = dt_hists_sdb[ibin].GetPtr()
      zg_onz_hist = zg_hists_onz[ibin].GetPtr()
      dt_onz_hist.Sumw2()
      dt_sdb_hist.Sumw2()
      dt_onz_hist.Scale(1.0/dt_onz_hist.Integral())
      dt_sdb_hist.Scale(1.0/dt_sdb_hist.Integral())
      zg_onz_hist.Scale(1.0/zg_onz_hist.Integral())
      bin_per_unit(dt_onz_hist, MVA_BINS)
      bin_per_unit(dt_sdb_hist, MVA_BINS)
      bin_per_unit(zg_onz_hist, MVA_BINS)
      var_plot = RplPlot()
      var_plot.lumi_data = [(ENERGY[year],LUMI[year])]
      var_plot.y_title = '% Events/IDMVA'
      var_plot.y_max_lower = 1.35
      var_plot.y_min_lower = 0.65
      var_plot.plot_outline(dt_onz_hist)
      var_plot.plot_outline(dt_sdb_hist)
      var_plot.plot_outline(zg_onz_hist)
      var_plot.add_ratio(dt_onz_name, zg_onz_name, True)
      var_plot.add_ratio(dt_sdb_name, zg_onz_name, True)
      var_plot.draw('plots/photon_shapecomp_{}_bin{}.pdf'.format(year,ibin))
      ibin += 1

def save_rootfile(year):
  dt_hists_onz, dt_hists_sdb = generate_validation_histogram2(year, 'data')
  zg_hists_onz, zg_hists_sdb = generate_validation_histogram2(year, 'zg')
  print('Performing event loop, this may take a while.')
  ibin = 0
  output_file = ROOT.TFile('temp_photonreweighthists.root','RECREATE')
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ABSETA_BINS)-1):
      mva_bins = MVA_BINS_EE
      if (ABSETA_BINS[ieta+1] <= 1.567):
        mva_bins = MVA_BINS_EB
      dt_onz_hist = dt_hists_onz[ibin].GetPtr()
      zg_onz_hist = zg_hists_onz[ibin].GetPtr()
      dt_onz_hist.Sumw2()
      dt_onz_hist.Scale(1.0/dt_onz_hist.Integral())
      zg_onz_hist.Scale(1.0/zg_onz_hist.Integral())
      bin_per_unit(dt_onz_hist, mva_bins)
      bin_per_unit(zg_onz_hist, mva_bins)
      dt_onz_hist.Divide(zg_onz_hist)
      dt_onz_hist.Write()
      ibin += 1
  output_file.Close()

def generate_polycorrections():
  input_file = ROOT.TFile('temp_photonreweighthists.root','READ')
  for ieta in range(len(ABSETA_BINS)-1):
    pt_vals = []
    mva_vals = []
    y_vals = []
    y_wgts = [] #1/var
    for ipt in range(len(PT_BINS)-1):
      ibin = ipt*(len(ABSETA_BINS)-1)+ieta
      pt_center = (PT_BINS[ipt]+PT_BINS[ipt+1])/2.0
      if (ipt==3):
        pt_center = 70
      eta_center = (ETA_BINS[ieta]+ETA_BINS[ieta+1])/2.0
      ratio_hist = input_file.Get('data_onz_hist_2016_bin{}'.format(ibin))
      for imva in range(len(MVA_BINS)-1):
        pt_vals.append(pt_center)
        mva_vals.append(ratio_hist.GetXaxis().GetBinCenter(imva+1))
        y_vals.append(ratio_hist.GetBinContent(imva+1))
        y_wgts.append(ratio_hist.GetBinError(imva+1))
        if y_wgts[-1]>0.0:
          y_wgts[-1] = 1.0/(y_wgts[-1]**2)
    corr, fit = ndpolyfit([np.array(pt_vals),np.array(mva_vals)], np.array(y_vals), 
                          np.array(y_wgts), ['pt','idmva'], (2,5))
    print('New bin')
    print(corr)
    print(fit)
  input_file.Close()

def generate_polycorrections2():
  input_file = ROOT.TFile('temp_photonreweighthists.root','READ')
  #ROOT.gInterpreter.Declare('''
  #''')
  ibin = 0
  for ieta in range(len(ABSETA_BINS)-1):
    for ipt in range(len(PT_BINS)-1):
      ratio_hist = input_file.Get('data_onz_hist_2016_bin{}'.format(ibin))
      fn = ROOT.TF1('gaussline','1.0+[2]/sqrt(2.0*3.1416*[1])*exp(-0.5*pow(x-[0],2)/[1]/[1])+[5]/sqrt(2.0*3.1416*[4])*exp(-0.5*pow(x-[3],2)/[4]/[4])',0.0,1.0)
      fn.SetParameter(0,0.0)
      fn.SetParameter(1,0.4)
      fn.SetParameter(2,0.3)
      fn.SetParameter(3,0.95)
      fn.SetParameter(4,0.02)
      fn.SetParameter(5,-0.05)
      fn.SetParLimits(1,0.3,2.0) #don't allow first gaussian to be too small
      fn.SetParLimits(4,0.01,2.0) #don't allow second gaussian to be too small
      fn.SetParLimits(2,-0.5,0.5) #don't allow very large gaussians
      fn.SetParLimits(5,-0.5,0.5) #don't allow very large gaussians
      ratio_hist.Fit('gaussline')
      print(fn.GetChisquare())
      ROOT.gStyle.SetOptStat(0)
      can = ROOT.TCanvas()
      ratio_hist.Draw()
      fn.Draw('same')
      can.SaveAs('plots/photonfittest_{}.pdf'.format(ibin))
      ibin += 1
  input_file.Close()

def generate_histograms(year, sample_type):
  '''Generates histograms used to derive weights, returns photon ID MVA 
  distribution histograms as two lists, one of on Z histograms binned in pT
  and eta, and one of high sideband (100-115 GeV) histograms, again binned in
  pT and eta

  year         string, data taking era
  sample_type  string, "data", "zg", or "dy"
  '''
  if not sample_type in ['data', 'zg', 'dy']:
    raise ValueError('Unsupported sample type.')
  filename = get_filename(year, sample_type)
  df = ROOT.RDataFrame('tree', filename)
  df = df.Filter('nmu>=2&&nphoton>=1&&trig_double_mu')
  df = df.Define('lead_photon_idmva','photon_idmva[0]')
  title = 'Data'
  if sample_type == 'dy':
    df = df.Filter('photon_pflavor[0]!=1')
  if sample_type == 'zg':
    title = 'MC'
  if sample_type != 'data':
    df = df.Filter('use_event')
    df = df.Define('w_lumi_year','w_lumi*{}'.format(LUMI[year]))
  else:
    df = df.Define('w_lumi_year','1')
  df_onz = df.Filter('llphoton_m[0]>80&&llphoton_m[0]<100')
  df_sideband = df.Filter('llphoton_m[0]>100&&llphoton_m[0]<115')

  onz_hists = []
  sideband_hists = []
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      filter_string = ('photon_eta[0]>{}&&photon_eta[0]<{}'
                       +'&&photon_pt[0]>{}&&photon_pt[0]<{}').format(
                       ETA_BINS[ieta],ETA_BINS[ieta+1],PT_BINS[ipt],
                       PT_BINS[ipt+1])
      df_onz_bin = df_onz.Filter(filter_string)
      df_sideband_bin = df_sideband.Filter(filter_string)
      onz_hists.append(df_onz_bin.Histo1D((
          '{}_onz_hist_bin{}'.format(sample_type,ibin),'{};IDMVA'.format(title),
          NMVA_BINS,MVA_LO,MVA_HI),'lead_photon_idmva','w_lumi_year'))
      sideband_hists.append(df_sideband_bin.Histo1D(
          ('{}_sideband_hist_bin{}'.format(sample_type,ibin),
          '{};IDMVA'.format(title),NMVA_BINS,
          MVA_LO,MVA_HI),'lead_photon_idmva','w_lumi_year'))
      ibin += 1

  return onz_hists, sideband_hists

def get_corrected_histograms(dt_z_hists, dt_sb_hists, zg_z_hists, zg_sb_hists,
                             dy_z_hists, dy_sb_hists):
  '''Estimates the true distribution of ZG events from measured data and 
  ZG/fake MC samples
  '''
  est_z_hists = []
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      est_z_hists.append(ROOT.TH1D('est_onz_hist_bin{}'.format(ibin),
          'Data-driven prompt photon estimate;IDMVA',NMVA_BINS,MVA_LO,MVA_HI))
      print('Bin: {}'.format(ibin))
      for imva in range(1,NMVA_BINS+1):
        dt_z_yield = dt_z_hists[ibin].GetBinContent(imva)
        dt_sb_yield = dt_sb_hists[ibin].GetBinContent(imva)
        zg_z_yield = zg_z_hists[ibin].GetBinContent(imva)
        zg_sb_yield = zg_sb_hists[ibin].GetBinContent(imva)
        dy_z_yield = dy_z_hists[ibin].GetBinContent(imva)
        dy_sb_yield = dy_sb_hists[ibin].GetBinContent(imva)
        yield_den = (dy_z_yield*zg_sb_yield-dy_sb_yield*zg_z_yield)
        est_yield = 0.0
        if (dy_z_yield <= 0.0):
          est_yield = dt_z_yield
        elif (yield_den <= 0.0):
          est_yield = zg_z_yield*dt_z_yield/(zg_z_yield+dy_z_yield)
        else:
          est_yield = ((dy_z_yield*dt_sb_yield-dy_sb_yield*dt_z_yield)
                       /yield_den*zg_z_yield)
        est_z_hists[-1].SetBinContent(imva,est_yield)
        print('Data yields: {}, {}'.format(dt_z_yield,dt_sb_yield))
        print('MC yields: {}, {}'.format(zg_z_yield,zg_sb_yield))
        print('DY yields: {}, {}'.format(dy_z_yield,dy_sb_yield))
        print('Estimated true data yield: {}'.format(est_yield))
      if (est_z_hists[-1].Integral() <= 0.0):
        raise RuntimeError('No integral')
      ibin += 1
  return est_z_hists

def write_sfs(dt_z_hists, zg_z_hists, dy_z_hists, sy_z_hists, year):
  '''Generates SFs, associated uncertainties, and writes them to a json
  '''
  values = []
  systs = []
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      for imva in range(1,NMVA_BINS+1):
        dt_frac = (dt_z_hists[ibin].GetBinContent(imva)
                   /dt_z_hists[ibin].Integral())
        zg_frac = (zg_z_hists[ibin].GetBinContent(imva)
                   /zg_z_hists[ibin].Integral())
        sy_frac = (sy_z_hists[ibin].GetBinContent(imva)
                   /sy_z_hists[ibin].Integral())
        dt_err = dt_z_hists[ibin].GetBinError(imva)/dt_z_hists[ibin].Integral()
        zg_err = zg_z_hists[ibin].GetBinError(imva)/zg_z_hists[ibin].Integral()
        if (zg_frac <= 0.0 or dt_frac <= 0.0):
          values.append(1.0)
          systs.append(0.0)
        else:
          sf = dt_frac/zg_frac
          sf_alt = sy_frac/zg_frac
          stat_err = sf*hypot(zg_err/zg_frac, dt_err/dt_frac)
          sys_err = abs(sf-sf_alt)
          values.append(sf)
          systs.append(hypot(stat_err,sys_err))
      ibin += 1

  json_content = make_correction('photon_shapecorr',values,systs)
  out_filename = 'dataset/photon_shapecorr_{}.json'.format(year)
  with open(out_filename,'w') as out_file:
    out_file.write(fix_correctionlib_json([
        json_content.json(exclude_unset=True)]))

def save_histograms(dt_z_hists, zg_z_hists, dy_z_hists, sy_z_hists, year):
  '''Makes nice plots
  '''
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      dt_name = 'data_onz_hist_bin{}'.format(ibin)
      zg_name = 'zg_onz_hist_bin{}'.format(ibin)
      sy_name = 'est_onz_hist_bin{}'.format(ibin)
      var_plot = RplPlot()
      var_plot.lumi_data = [(ENERGY[year],LUMI[year])]
      var_plot.y_title = 'Events/bin'
      var_plot.plot_outline(dt_z_hists[ibin].GetPtr())
      var_plot.plot_outline(zg_z_hists[ibin].GetPtr())
      var_plot.plot_outline(sy_z_hists[ibin])
      var_plot.add_ratio(dt_name, zg_name, True)
      var_plot.add_ratio(sy_name, zg_name, True)
      var_plot.draw('plots/photon_shape_{}_bin{}.pdf'.format(year,ibin))
      ibin += 1

def generate_weights(year):
  #get histograms
  dt_hists_onz, dt_hists_sideband = generate_histograms(year, 'data')
  zg_hists_onz, zg_hists_sideband = generate_histograms(year, 'zg')
  dy_hists_onz, dy_hists_sideband = generate_histograms(year, 'dy')
  print('Generating histograms, this may take a while.')

  #get corrected histograms
  data_zg_hists_onz = get_corrected_histograms(dt_hists_onz,dt_hists_sideband,
                                               zg_hists_onz,zg_hists_sideband,
                                               dy_hists_onz,dy_hists_sideband)

  #write text output and make summary plots
  write_sfs(dt_hists_onz, zg_hists_onz, dy_hists_onz, data_zg_hists_onz, 
            args.year)
  save_histograms(dt_hists_onz, zg_hists_onz, dy_hists_onz, data_zg_hists_onz,
                  args.year)

if __name__=='__main__':
  parser = ArgumentParser(prog='generate_simplephotonweights',
                          description='script that generates photon weights')
  parser.add_argument('-y','--year',choices=['2016APV','2016','2017','2018',
                                             '2022','2022EE','2023','2023BPix'])
  args = parser.parse_args()
  save_rootfile(args.year)
  generate_polycorrections2()


