#!/usr/bin/env python3
"""@package docstring
Generates photon preselection corrections using tag-and-probe
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from root_plot_lib import RplPlot
from math import hypot
import ROOT
import json

#constants
PT_BINS = [15.0,20.0,35.0,50.0,80.0,500.0]
ETA_BINS = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
NMVA_BINS = 25
MVA_HI = 1.0
MVA_LO = 0.14
MVA_BINS = [i*((MVA_HI-MVA_LO)/NMVA_BINS)+MVA_LO for i in range(NMVA_BINS+1)]
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

#ROOT JIT'ed C++ definitions

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

  if sample_type=='zg':
    if year == '2016APV':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_28.root'

  if sample_type=='dy':
    if year == '2016APV':
      return '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_62.root'

  raise ValueError('Unknown sample')

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

if __name__=='__main__':
  parser = ArgumentParser(prog='generate_simplephotonweights',
                          description='script that generates photon weights')
  parser.add_argument('-y','--year',choices=['2016APV','2016','2017','2018',
                                             '2022','2022EE','2023','2023BPix'])
  args = parser.parse_args()

  #get histograms
  dt_hists_onz, dt_hists_sideband = generate_histograms(args.year, 'data')
  zg_hists_onz, zg_hists_sideband = generate_histograms(args.year, 'zg')
  dy_hists_onz, dy_hists_sideband = generate_histograms(args.year, 'dy')
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
