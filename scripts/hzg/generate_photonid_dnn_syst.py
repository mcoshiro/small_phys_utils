#!/usr/bin/env python3
"""@package docstring
Generates photon preselection corrections using tag-and-probe
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from tnp_utils import do_tnp_fit, integrate_tgraph
from root_plot_lib import RplPlot
from math import sqrt
import ROOT
import json
import subprocess
import ctypes

#constants
idmva_bins = [float(i)*(1.0+0.58)/20.0-0.58 for i in range(21)]
#idmva_bins = [float(i)*(1.0+0.58)/40.0-0.58 for i in range(41)]

def book_hist1d_zpeak_bkgsub(dataframe, hist_desc, branch_name, weight_name):
  '''Return pointers to histograms used to extract yield after nonresonant background subtraction'''
  peak_ptr = dataframe.Filter('mmy_mass>80&&mmy_mass<98').Histo1D(hist_desc,branch_name,weight_name)
  subt_ptr = dataframe.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo1D(hist_desc,branch_name,weight_name)
  return (peak_ptr, subt_ptr)

def book_hist2d_zpeak_bkgsub(dataframe, hist_desc, branch1_name, branch2_name, weight_name):
  '''Return pointers to histograms used to extract yield after nonresonant background subtraction'''
  peak_ptr = dataframe.Filter('mmy_mass>80&&mmy_mass<98').Histo2D(hist_desc,branch1_name,branch2_name,weight_name)
  subt_ptr = dataframe.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo2D(hist_desc,branch1_name,branch2_name,weight_name)
  return (peak_ptr, subt_ptr)

def book_hist3d_zpeak_bkgsub(dataframe, hist_desc, branch1_name, branch2_name, branch3_name, weight_name):
  '''Return pointers to histograms used to extract yield after nonresonant background subtraction'''
  peak_ptr = dataframe.Filter('mmy_mass>80&&mmy_mass<98').Histo3D(hist_desc,branch1_name,branch2_name,branch3_name,weight_name)
  subt_ptr = dataframe.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo3D(hist_desc,branch1_name,branch2_name,branch3_name,weight_name)
  return (peak_ptr, subt_ptr)

def get_bkgsub_hist(bkgsub_hist_ptrs):
  '''Returns histogram with background subtracted'''
  tot_hist = bkgsub_hist_ptrs[0].GetValue()
  #sub_hist = bkgsub_hist_ptrs[1].GetValue()
  #tot_hist.Add(sub_hist, -4.0/18.0)
  return tot_hist

def debug_plots(data_input_file, simu_input_file, region): 
  '''temp debugging
  '''
  region_sel_string = 'photon_isScEtaEB[0]'
  #region_sel_string = 'photon_isScEtaEE[0]'
  weight_definition = 'get_sf({photon_r9[0],photon_s4[0],photon_etawidth[0],photon_phiwidth[0],photon_sieie[0],photon_sieip[0],photon_phiso[0],photon_chiso[0],photon_chiso_worst[0],photon_energy_raw[0],photon_origin_eta[0],photon_pt[0],fabs(photon_origin_eta[0]),rho})'

  zg_filename = '/net/cms11/cms11r0/pico/NanoAODv9UCSB1/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_77.root'
  dy_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_llg/merged_pico_llg_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_204.root'
  data_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_206.root'

  #setup dataframe and add tag/probe preselections
  df_data = ROOT.RDataFrame('tree',data_filename)
  df_data = df_data.Filter('llphoton_m[0]>50.0&&llphoton_m[0]<110.0')
  df_data = df_data.Filter('nmu==2&&nphoton==1&&trig_double_mu')
  df_data = df_data.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_data = df_data.Filter(region_sel_string)
  df_data = df_data.Define('mmy_mass','llphoton_m[0]')
  df_data = df_data.Define('lead_photon_idmva','photon_idmva[0]')

  df_zg = ROOT.RDataFrame('tree',zg_filename)
  df_zg = df_zg.Filter('llphoton_m[0]>50.0&&llphoton_m[0]<110.0')
  df_zg = df_zg.Filter('use_event&&nmu==2&&nphoton==1&&trig_double_mu')
  df_zg = df_zg.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_zg = df_zg.Filter(region_sel_string)
  df_zg = df_zg.Define('w_photon_dnn','photon_pflavor[0]==1 ? '+weight_definition+' : 1')
  df_zg = df_zg.Define('mmy_mass','llphoton_m[0]')
  df_zg = df_zg.Define('w_lumi_photon','w_lumi*w_photon_dnn')
  df_zg = df_zg.Define('w_lumi_photon_year','w_lumi*w_photon_dnn*60.0')
  df_zg = df_zg.Define('lead_photon_idmva','photon_idmva[0]')
  df_zg = df_zg.Define('w_lumi_year','w_lumi*60.0')

  df_dy = ROOT.RDataFrame('tree',dy_filename)
  df_dy = df_dy.Filter('llphoton_m[0]>50.0&&llphoton_m[0]<110.0')
  df_dy = df_dy.Filter('use_event&&nmu==2&&nphoton==1&&trig_double_mu')
  df_dy = df_dy.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_dy = df_dy.Filter(region_sel_string)
  df_dy = df_dy.Filter('photon_pflavor[0]!=1')
  df_dy = df_dy.Define('mmy_mass','llphoton_m[0]')
  df_dy = df_dy.Define('lead_photon_idmva','photon_idmva[0]')
  df_dy = df_dy.Define('w_lumi_year','w_lumi*60.0')

  dt_hist_ptrs = []
  zg_hist_ptrs = []
  dy_hist_ptrs = []
 
  for iid in range(len(idmva_bins)-1):
    df_dt_bin = df_data.Filter('lead_photon_idmva>'+str(idmva_bins[iid])+'&&lead_photon_idmva<'+str(idmva_bins[iid+1]))
    df_dy_bin = df_dy.Filter('lead_photon_idmva>'+str(idmva_bins[iid])+'&&lead_photon_idmva<'+str(idmva_bins[iid+1]))
    df_zg_bin = df_zg.Filter('lead_photon_idmva>'+str(idmva_bins[iid])+'&&lead_photon_idmva<'+str(idmva_bins[iid+1]))
    dt_hist_ptrs.append(df_dt_bin.Histo1D(('hist_mllg_data'+str(iid),'data;m_{#mu#mu#gamma} [GeV]',60,50.0,110.0),'mmy_mass'))
    dy_hist_ptrs.append(df_dy_bin.Histo1D(('hist_mllg_dy'+str(iid),'Z+fake photon;m_{#mu#mu#gamma} [GeV]',60,50.0,110.0),'mmy_mass','w_lumi_year'))
    zg_hist_ptrs.append(df_zg_bin.Histo1D(('hist_mllg_zg'+str(iid),'Z+#gamma;m_{#mu#mu#gamma} [GeV]',60,50.0,110.0),'mmy_mass','w_lumi_photon_year'))

  for iid in range(len(idmva_bins)-1):
    dt_hist = dt_hist_ptrs[iid].GetValue()
    dy_hist = dy_hist_ptrs[iid].GetValue()
    zg_hist = zg_hist_ptrs[iid].GetValue()
    dy_hist.Add(zg_hist) #stack

    var_plot = RplPlot()
    var_plot.lumi_data = [(60,13)]
    var_plot.y_title = 'Events/bin'
    var_plot.plot_filled(dy_hist,ROOT.TColor.GetColor('#f89c20'))
    var_plot.plot_filled(zg_hist,ROOT.TColor.GetColor('#5790fc'))
    var_plot.plot_points(dt_hist,ROOT.kBlack)
    var_plot.add_ratio('hist_mllg_data'+str(iid),'hist_mllg_dy'+str(iid), True)
    var_plot.draw('plots/check_mllg_idbin'+str(iid)+'.pdf')

  #var_defs = ['photon_r9[0]','photon_s4[0]','photon_etawidth[0]',
  #            'photon_phiwidth[0]','photon_sieie[0]','photon_sieip[0]',
  #            'photon_phiso[0]',
  #            'photon_chiso[0]','photon_chiso_worst[0]',
  #            'photon_energy_raw[0]','photon_origin_eta[0]',
  #            'photon_esoversc[0]','photon_essigmarr[0]']
  #var_names = ['lead_photon_r9','lead_photon_s4','lead_photon_etawidth',
  #             'lead_photon_phiwidth','lead_photon_sieie','lead_photon_sieiep',
  #             'lead_photon_phiso',
  #             'lead_photon_chiso','lead_photon_chiso_worst',
  #             'lead_photon_energy_raw','lead_photon_origin_eta',
  #             'lead_photon_esoversc','lead_photon_essigmarr']
  #var_descs = ['R_{9}','s_{4}','#eta width','#phi width','#sigma_{i#eta i#eta}',
  #            '#sigma_{i#eta i#phi}','I_{abs,photons} [GeV]',
  #            'I_{abs,charged} [GeV]',
  #            'I_{abs,charged} w.r.t. Worst Vertex [GeV]','E_{raw} [GeV]',
  #            '#eta_{origin}','E_{ES}/E_{raw,EE}','#sigma^{ES}_{RR}']
  #var_lowers = [0.4,0.4,0.004,0.005,0.006,-0.002,0.0,0.0,0.00,0.000,-2.5,0.00,0.00]
  #var_uppers = [1.0,1.0,0.023,0.100,0.030,0.0020,6.0,6.0,25.0,100.0,2.50,0.16,10.0]

  #var_hist_ptrs = []

  #for var_def, var_name in zip(var_defs, var_names):
  #  df_zg = df_zg.Define(var_name, var_def)

  #var_names.append('rho')
  #var_descs.append('#rho [GeV]')
  #var_uppers.append(0.00)
  #var_lowers.append(45.0)

  #idmva_broad_bins = [-0.6,-0.33,-0.07,0.2,0.47,0.73,1.0]

  #for var_name, var_desc, var_lower, var_upper in zip(var_names,var_descs,var_lowers,var_uppers):
  #  for iid in range(len(idmva_broad_bins)-1):
  #    df_zg_bin = df_zg.Filter('lead_photon_idmva>'+str(idmva_broad_bins[iid])+'&&lead_photon_idmva<'+str(idmva_broad_bins[iid+1]))
  #    var_hist_ptrs.append(df_zg_bin.Histo1D(('hist_'+var_name+'_'+str(iid),
  #                         'ID>'+str(idmva_broad_bins[iid])+',ID<'+str(idmva_broad_bins[iid+1])+';'+var_desc,
  #                         30,var_lower,var_upper),var_name,'w_lumi_year'))

  #ibin = 0
  #for ivar in range(len(var_names)):
  #  var_hists = []
  #  var_plot = RplPlot()
  #  var_plot.lumi_data = [(60,13)]
  #  var_plot.y_title = '% Events/bin'
  #  for iid in range(len(idmva_broad_bins)-1):
  #    var_hists.append(var_hist_ptrs[ibin].GetValue())
  #    var_hists[iid].Scale(1.0/var_hists[iid].Integral())
  #    var_plot.plot_outline(var_hists[iid])
  #    ibin += 1
  #  var_plot.draw('plots/check_'+var_names[ivar]+'_idcorr.pdf')

def make_photon_electron_plots():
  #Generate plots (in MC) comparing distributions of MVA inputs for electrons and photons
  #Only have data reference for electrons for most variables

  df_ph = ROOT.RDataFrame('tree','/net/cms11/cms11r0/pico/NanoAODv9UCSB1/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_77.root')
  df_ph = df_ph.Filter('photon_isScEtaEB[0]||photon_isScEtaEE[0]')
  df_ph = df_ph.Filter('llphoton_m[0]>80.0&&llphoton_m[0]<100.0')
  df_ph = df_ph.Filter('use_event&&nmu==2&&nphoton==1&&trig_double_mu')
  df_ph = df_ph.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_ph = df_ph.Define('w_lumi_year','w_lumi*60.0')
  df_ph = df_ph.Define('lead_photon_pt','photon_pt[0]')
  df_ph = df_ph.Define('lead_photon_r9','photon_r9[0]')
  df_ph = df_ph.Define('lead_photon_sieie','photon_sieie[0]')
  df_ph = df_ph.Define('lead_photon_hoe','photon_hoe[0]')
  df_ph = df_ph.Define('lead_photon_phiso','photon_phiso[0]')
  df_ph = df_ph.Define('lead_photon_chiso','photon_chiso[0]')
  df_ph = df_ph.Define('lead_photon_abseta','fabs(photon_origin_eta[0])')
  df_ph = df_ph.Define('mmy_mass','llphoton_m[0]')
  #df_ph = df_ph.Filter('lead_photon_pt>30')
  #df_ph = df_ph.Filter('lead_photon_r9>0.8&&lead_photon_sieie<0.015&&lead_photon_phiso<4.0&&lead_photon_chiso<6.0')
  #df_ph = df_ph.Filter('lead_photon_r9>0.9&&lead_photon_r9<0.95&&lead_photon_pt>30.0&&lead_photon_pt<32.5&&lead_photon_abseta<0.4')
  #df_ph = df_ph.Filter('lead_photon_pt>50.0&&lead_photon_pt<55.0&&lead_photon_abseta<0.4')
  df_ph = df_ph.Filter('lead_photon_pt>55.0&&lead_photon_abseta<1.5')

  #df_ph_data = ROOT.RDataFrame('tree','/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_206.root')
  #df_ph_data = df_ph_data.Filter('photon_isScEtaEE[0]')
  #df_ph_data = df_ph_data.Filter('nmu==2&&nphoton==1&&trig_double_mu')
  #df_ph_data = df_ph_data.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  #df_ph_data = df_ph_data.Define('lead_photon_pt','photon_pt[0]')
  #df_ph_data = df_ph_data.Define('lead_photon_abseta','fabs(photon_eta[0])')
  #df_ph_bkgd = df_ph_data.Filter('(llphoton_m[0]>78.0&&llphoton_m[0]<80.0)||(llphoton_m[0]>98.0&&llphoton_m[0]<100.0)')
  #df_ph_data = df_ph_data.Filter('llphoton_m[0]>80.0&&llphoton_m[0]<98.0')

  df_el = ROOT.RDataFrame('tree','/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root')
  #df_el = df_el.Filter('ph_sc_abseta<1.5')
  df_el = df_el.Filter('pair_mass>80.0&&pair_mass<100.0')
  df_el = df_el.Define('w_lumi_year','1')
  df_el = df_el.Define('mmy_mass','pair_mass')
  #df_el = df_el.Filter('ph_et>30')
  #df_el = df_el.Filter('ph_r9>0.8&&ph_sieie<0.015&&ph_phoIso<4.0&&ph_chIso<6.0')
  #df_el = df_el.Filter('ph_r9>0.9&&ph_r9<0.95&&ph_et>30.0&&ph_et<32.5&&ph_sc_abseta<0.4')
  #df_el = df_el.Filter('ph_et>50.0&&ph_et<55.0&&ph_sc_abseta<0.4')
  df_el = df_el.Filter('ph_et>55.0&&ph_sc_abseta<1.5')

  #df_el_data = ROOT.RDataFrame('tree','/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root')
  #df_el_data = df_el_data.Filter('ph_sc_abseta>1.5')
  #df_el_data = df_el_data.Filter('pair_mass>80.0&&pair_mass<100.0')

  #prew_ph_hist_ptr = df_ph.Histo2D(('prew_ph_hist','',40,15.0,75.0,10,1.5,2.5),'lead_photon_pt','lead_photon_abseta','w_lumi_year')
  #prew_el_hist_ptr = df_el.Histo2D(('prew_el_hist','',40,15.0,75.0,10,1.5,2.5),'ph_et','ph_sc_abseta')
  #prew_ph_data_hist_ptr = df_ph_data.Histo2D(('prew_ph_data_hist','',40,15.0,75.0,10,1.5,2.5),'lead_photon_pt','lead_photon_abseta')
  #prew_ph_bkgd_hist_ptr = df_ph_bkgd.Histo2D(('prew_ph_bkgd_hist','',40,15.0,75.0,10,1.5,2.5),'lead_photon_pt','lead_photon_abseta')
  ##prew_el_data_hist_ptr = df_el_data.Histo2D(('prew_el_data_hist','',40,15.0,75.0,10,1.5,2.5),'ph_et','ph_sc_abseta','w_bkg')
  #prew_rho_ph_hist_ptr = df_ph.Histo1D(('prew_rho_ph_hist','',40,0.0,45.0),'rho','w_lumi_year')
  #prew_rho_el_hist_ptr = df_el.Histo1D(('prew_rho_el_hist','',40,0.0,45.0),'event_rho')
  #prew_rho_ph_data_hist_ptr = df_ph_data.Histo1D(('prew_rho_ph_data_hist','',40,0.0,45.0),'rho')
  #prew_rho_ph_bkgd_hist_ptr = df_ph_bkgd.Histo1D(('prew_rho_ph_bkgd_hist','',40,0.0,45.0),'rho')
  #prew_rho_el_data_hist_ptr = df_el_data.Histo1D(('prew_rho_el_data_hist','',40,0.0,45.0),'event_rho','w_bkg')

  #reweight kinematics
  #prew_ph_hist_ptrs = book_hist3d_zpeak_bkgsub(df_ph, 
  #    ('prew_ph_hist','',30,15.0,75.0,15,0.0,2.5,40,0.0,1.0), 'lead_photon_pt',
  #    'lead_photon_abseta', 'lead_photon_r9', 'w_lumi_year')
  #prew_el_hist_ptrs = book_hist3d_zpeak_bkgsub(df_el, 
  #    ('prew_el_hist','',30,15.0,75.0,15,0.0,2.5,40,0.0,1.0), 'ph_et',
  #    'ph_sc_abseta', 'ph_r9', 'w_lumi_year')
  prew_ph_hist_ptrs = book_hist2d_zpeak_bkgsub(df_ph, 
      ('prew_ph_hist','',30,15.0,75.0,15,0.0,2.5), 'lead_photon_pt',
      'lead_photon_abseta', 'w_lumi_year')
  prew_el_hist_ptrs = book_hist2d_zpeak_bkgsub(df_el, 
      ('prew_el_hist','',30,15.0,75.0,15,0.0,2.5), 'ph_et',
      'ph_sc_abseta', 'w_lumi_year')
  prew_ph_rho_hist_ptrs = book_hist1d_zpeak_bkgsub(df_ph, 
      ('prew_ph_rho_hist','',45,0.0,45.0), 'rho', 'w_lumi_year')
  prew_el_rho_hist_ptrs = book_hist1d_zpeak_bkgsub(df_el, 
      ('prew_el_rho_hist','',45,0.0,45.0), 'event_rho', 'w_lumi_year')

  prew_ph_hist = get_bkgsub_hist(prew_ph_hist_ptrs)
  prew_el_hist = get_bkgsub_hist(prew_el_hist_ptrs)
  prew_ph_rho_hist = get_bkgsub_hist(prew_ph_rho_hist_ptrs)
  prew_el_rho_hist = get_bkgsub_hist(prew_el_rho_hist_ptrs)

  prew_ph_hist.Scale(1.0/prew_ph_hist.Integral())
  prew_el_hist.Scale(1.0/prew_el_hist.Integral())
  prew_ph_rho_hist.Scale(1.0/prew_ph_rho_hist.Integral())
  prew_el_rho_hist.Scale(1.0/prew_el_rho_hist.Integral())

  prew_ph_hist.Divide(prew_el_hist)
  prew_ph_rho_hist.Divide(prew_el_rho_hist)
  ROOT.gDirectory.Add(prew_ph_hist)
  ROOT.gDirectory.Add(prew_ph_rho_hist)

  ROOT.gInterpreter.Declare("""
  TH1D* prew_ph_hist = static_cast<TH1D*>(gDirectory->Get("prew_ph_hist"));
  TH1D* prew_ph_rho_hist = static_cast<TH1D*>(gDirectory->Get("prew_ph_rho_hist"));
  """)

  df_el = df_el.Define('w_pteta','(prew_ph_hist->GetBinContent(prew_ph_hist->FindBin(ph_et,ph_sc_abseta)))*(prew_ph_rho_hist->GetBinContent(prew_ph_rho_hist->FindBin(event_rho)))')
  #df_el = df_el.Define('w_pteta','(prew_ph_hist->GetBinContent(prew_ph_hist->FindBin(ph_et,ph_sc_abseta,ph_r9)))*(prew_ph_rho_hist->GetBinContent(prew_ph_rho_hist->FindBin(event_rho)))')
  df_el = df_el.Define('w_lumi_pteta','w_pteta')

  #prew_ph_hist = prew_ph_hist_ptr.GetValue()
  #prew_el_hist = prew_el_hist_ptr.GetValue()
  #prew_ph_data_hist = prew_ph_data_hist_ptr.GetValue()
  #prew_ph_bkgd_hist = prew_ph_bkgd_hist_ptr.GetValue()
  #prew_el_data_hist = prew_el_data_hist_ptr.GetValue()
  #prew_ph_hist.Scale(1.0/prew_ph_hist.Integral())
  #prew_el_hist.Scale(1.0/prew_el_hist.Integral())
  #prew_ph_data_hist.Add(prew_ph_bkgd_hist,-1.0*18.0/4.0)
  #prew_ph_data_hist.Scale(1.0/prew_ph_data_hist.Integral())
  #prew_el_data_hist.Scale(1.0/prew_el_data_hist.Integral())
  #prew_ph_hist_scale_phdata = prew_ph_hist.Clone('prew_ph_hist_phdata')
  #prew_ph_hist_scale_elsimu = prew_ph_hist.Clone('prew_ph_hist_elsimu')
  #prew_ph_hist_scale_eldata = prew_ph_hist.Clone('prew_ph_hist_eldata')
  #ROOT.gDirectory.Add(prew_ph_hist_scale_phdata)
  #ROOT.gDirectory.Add(prew_ph_hist_scale_elsimu)
  #ROOT.gDirectory.Add(prew_ph_hist_scale_eldata)

  #prew_rho_ph_hist = prew_rho_ph_hist_ptr.GetValue()
  #prew_rho_el_hist = prew_rho_el_hist_ptr.GetValue()
  #prew_rho_ph_data_hist = prew_rho_ph_data_hist_ptr.GetValue()
  #prew_rho_ph_bkgd_hist = prew_rho_ph_bkgd_hist_ptr.GetValue()
  #prew_rho_el_data_hist = prew_rho_el_data_hist_ptr.GetValue()
  #prew_rho_ph_hist.Scale(1.0/prew_rho_ph_hist.Integral())
  #prew_rho_el_hist.Scale(1.0/prew_rho_el_hist.Integral())
  #prew_rho_ph_data_hist.Add(prew_rho_ph_bkgd_hist,-1.0*18.0/4.0)
  #prew_rho_ph_data_hist.Scale(1.0/prew_rho_ph_data_hist.Integral())
  #prew_rho_el_data_hist.Scale(1.0/prew_rho_el_data_hist.Integral())
  #prew_rho_ph_hist_scale_phdata = prew_rho_ph_hist.Clone('prew_rho_ph_hist_phdata')
  #prew_rho_ph_hist_scale_elsimu = prew_rho_ph_hist.Clone('prew_rho_ph_hist_elsimu')
  #prew_rho_ph_hist_scale_eldata = prew_rho_ph_hist.Clone('prew_rho_ph_hist_eldata')
  #ROOT.gDirectory.Add(prew_rho_ph_hist_scale_phdata)
  #ROOT.gDirectory.Add(prew_rho_ph_hist_scale_elsimu)
  #ROOT.gDirectory.Add(prew_rho_ph_hist_scale_eldata)

  #ROOT.gInterpreter.Declare("""
  #TH1D* prew_ph_hist_phdata = static_cast<TH1D*>(gDirectory->Get("prew_ph_hist_phdata"));
  #TH1D* prew_ph_hist_elsimu = static_cast<TH1D*>(gDirectory->Get("prew_ph_hist_elsimu"));
  #TH1D* prew_ph_hist_eldata = static_cast<TH1D*>(gDirectory->Get("prew_ph_hist_eldata"));

  #TH1D* prew_rho_ph_hist_phdata = static_cast<TH1D*>(gDirectory->Get("prew_rho_ph_hist_phdata"));
  #TH1D* prew_rho_ph_hist_elsimu = static_cast<TH1D*>(gDirectory->Get("prew_rho_ph_hist_elsimu"));
  #TH1D* prew_rho_ph_hist_eldata = static_cast<TH1D*>(gDirectory->Get("prew_rho_ph_hist_eldata"));
  #""")

  #df_el = df_el.Define('w_pteta','prew_ph_hist_elsimu->GetBinContent(prew_ph_hist_elsimu->FindBin(ph_et,ph_sc_abseta))')
  #df_el = df_el.Define('w_rho','prew_rho_ph_hist_elsimu->GetBinContent(prew_rho_ph_hist_elsimu->FindBin(event_rho))')
  #df_el = df_el.Define('w_pteta_rho','w_pteta*w_rho')
  #df_ph_data = df_ph_data.Define('w_pteta','prew_ph_hist_phdata->GetBinContent(prew_ph_hist_phdata->FindBin(photon_pt[0],photon_eta[0]))')
  ##df_ph_data = df_ph_data.Define('w_rho','prew_rho_ph_hist_phdata->GetBinContent(prew_rho_ph_hist_phdata->FindBin(rho))')
  ##df_ph_data = df_ph_data.Define('w_pteta_rho','w_pteta*w_rho')
  #df_ph_bkgd = df_ph_bkgd.Define('w_pteta','prew_ph_hist_phdata->GetBinContent(prew_ph_hist_phdata->FindBin(photon_pt[0],photon_eta[0]))')
  ##df_ph_bkgd = df_ph_bkgd.Define('w_rho','prew_rho_ph_hist_phdata->GetBinContent(prew_rho_ph_hist_phdata->FindBin(rho))')
  ##df_ph_bkgd = df_ph_bkgd.Define('w_pteta_rho','w_pteta*w_rho')
  #df_el_data = df_el_data.Define('w_pteta','prew_ph_hist_eldata->GetBinContent(prew_ph_hist_eldata->FindBin(ph_et,ph_sc_abseta))')
  ##df_el_data = df_el_data.Define('w_rho','prew_rho_ph_hist_eldata->GetBinContent(prew_rho_ph_hist_eldata->FindBin(event_rho))')
  ##df_el_data = df_el_data.Define('w_pteta_rho_bkg','w_pteta*w_rho*w_bkg')
  #df_el_data = df_el_data.Define('w_pteta_bkg','w_pteta*w_bkg')

  pico_vars = ['lead_photon_r9','lead_photon_s4','lead_photon_etawidth',
               'lead_photon_phiwidth','lead_photon_sieie','lead_photon_sieip',
               'lead_photon_phiso','lead_photon_chiso',
               'lead_photon_chiso_worst','lead_photon_essigmarr','lead_photon_esoversc',
               'lead_photon_pt','lead_photon_eta',
               'event_rho','lead_photon_idmva']
  pico_defs = ['photon_r9[0]','photon_s4[0]','photon_etawidth[0]',
               'photon_phiwidth[0]','photon_sieie[0]','photon_sieip[0]',
               'photon_phiso[0]','photon_chiso[0]','photon_chiso_worst[0]',
               'photon_essigmarr[0]','photon_esoversc[0]','photon_pt[0]',
               'photon_origin_eta[0]','rho','photon_idmva[0]']
  tnp_vars = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie',
              'ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_ESsigma','ph_esEnergyOverRawE',
              'ph_et','ph_sc_eta','event_rho',
              'ph_mva94XV2']
  xtitle = ['Photon R_{9}','Photon s_{4}','Photon supercluster #eta width',
            'Photon supercluster #phi width','Photon #sigma_{i#eta i#eta}',
            'Photon #sigma_{i#eta i#phi}','Photon I_{abs}(photon) [GeV]',
            'Photon I_{abs}(charged, PV) [GeV]',
            'Photon I_{abs}(charged, worst vertex) [GeV]',
            'Photon preshower #sigma_{eff}',
            'Photon E_{ES}/E_{raw}',
            'Photon p_{T} [GeV]','Photon #eta_{SC}','#rho [GeV]',
            'Photon IDMVA']
  ranges = [(0.4,1.0),(0.4,0.975),(0.004,0.0225),(0.005,0.1),(0.006,0.03),
            (-0.002,0.002),(0.0,6.0),(0.0,6.0),(0.0,25.0),(0.0,10.0),(0.0,0.16),
            (0.0,100.0),(-2.5,2.5),(0.0,45.0),(-1.0,1.0)]

  ph_hist_ptrs = []
  el_hist_ptrs = []

  for ivar in range(len(pico_vars)):
    if not (pico_vars[ivar] in ['lead_photon_pt','lead_photon_r9','lead_photon_sieie','lead_photon_chiso','lead_photon_phiso']):
      df_ph = df_ph.Define(pico_vars[ivar],pico_defs[ivar])

    ph_hist_ptrs.append(book_hist1d_zpeak_bkgsub(df_ph, 
        ('ph_hist','Photon simulation;'+xtitle[ivar],60,ranges[ivar][0],ranges[ivar][1]),
        pico_vars[ivar], 'w_lumi_year'))
    el_hist_ptrs.append(book_hist1d_zpeak_bkgsub(df_el, 
        ('el_hist','Electron simulation;'+xtitle[ivar],60,ranges[ivar][0],ranges[ivar][1]), 
        tnp_vars[ivar], 'w_lumi_pteta'))

  for ivar in range(len(pico_vars)):
    ph_hist = get_bkgsub_hist(ph_hist_ptrs[ivar])
    el_hist = get_bkgsub_hist(el_hist_ptrs[ivar])
    ph_hist.Scale(1.0/ph_hist.Integral())
    el_hist.Scale(1.0/el_hist.Integral())
    phel_plot = RplPlot()
    phel_plot.title_type = 'cms simulation'
    phel_plot.y_title = '% Events/bin'
    phel_plot.plot_outline(ph_hist)
    phel_plot.plot_outline(el_hist)
    phel_plot.add_ratio('ph_hist','el_hist',True)
    phel_plot.y_title_lower = '#frac{Photon}{Electron}'
    phel_plot.draw('plots/check_phel_'+tnp_vars[ivar]+'.pdf')

  #r9_bins = [0.4+i*0.6/30.0 for i in range(31)]
  #r9_hist_ptrs = []
  #for ir9 in range(len(r9_bins)-1):
  #  df_bin = df_el.Filter('(ph_r9>'+str(r9_bins[ir9])+')&&(ph_r9<'+str(r9_bins[ir9+1])+')')
  #  r9_hist_ptrs.append(df_bin.Histo1D(('r9_hist',';m_{l#gamma} [GeV]',35,75.0,110.0),'pair_mass'))
  #for ir9 in range(len(r9_bins)-1):
  #  r9_hist = r9_hist_ptrs[ir9].GetValue()
  #  r9_plot = RplPlot()
  #  r9_plot.title_type = 'cms simulation'
  #  r9_plot.y_title = 'Events/bin'
  #  r9_plot.plot_outline(r9_hist)
  #  r9_plot.draw('plots/check_phel_r9bin_'+str(ir9)+'.pdf')


def make_plots(data_input_file, simu_input_file, region):
  '''Generate mll spectra to fit in the next step
  '''
  region_sel_string = 'photon_isScEtaEB[0]'
  weight_definition = 'get_sf({photon_r9[0],photon_s4[0],photon_etawidth[0],photon_phiwidth[0],photon_sieie[0],photon_sieip[0],photon_phiso[0],photon_chiso[0],photon_chiso_worst[0],photon_energy_raw[0],photon_origin_eta[0],photon_pt[0],fabs(photon_origin_eta[0]),rho})'
  if (region=='endcap'):
    region_sel_string = 'photon_isScEtaEE[0]'
    weight_definition = 'get_sf({photon_r9[0],photon_s4[0],photon_etawidth[0],photon_phiwidth[0],photon_sieie[0],photon_sieip[0],photon_phiso[0],photon_chiso[0],photon_chiso_worst[0],photon_origin_eta[0],photon_esoversc[0],photon_pt[0],fabs(photon_origin_eta[0]),rho,photon_essigmarr[0]})'

  #setup dataframe and add tag/probe preselections
  df_data = ROOT.RDataFrame('tree',data_input_file)
  df_data = df_data.Filter('llphoton_m[0]>50.0&&llphoton_m[0]<110.0')
  #df_data = df_data.Filter('photon_pt[0]>25.0')
  df_data = df_data.Filter('nmu==2&&nphoton==1&&trig_double_mu')
  df_data = df_data.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_data = df_data.Filter(region_sel_string)
  df_data = df_data.Define('mmy_mass','llphoton_m[0]')
  df_data = df_data.Define('lead_photon_idmva','photon_idmva[0]')
  df_data = df_data.Define('lead_photon_pt','photon_pt[0]')
  df_data = df_data.Define('lead_photon_eta','photon_eta[0]')
  df_data = df_data.Define('lead_photon_abseta','fabs(photon_eta[0])')
  df_data = df_data.Define('w_lumi','1.0')

  df_simu = ROOT.RDataFrame('tree',simu_input_file)
  df_simu = df_simu.Filter('llphoton_m[0]>50.0&&llphoton_m[0]<110.0')
  #df_simu = df_simu.Filter('photon_pt[0]>25.0')
  df_simu = df_simu.Filter('use_event&&nmu==2&&nphoton==1&&trig_double_mu')
  df_simu = df_simu.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_simu = df_simu.Filter(region_sel_string)
  df_simu = df_simu.Define('w_photon_dnn','photon_pflavor[0]==1 ? '+weight_definition+' : 1')
  df_simu = df_simu.Define('mmy_mass','llphoton_m[0]')
  df_simu = df_simu.Define('w_lumi_photon','w_lumi*w_photon_dnn*w_pu')
  df_simu = df_simu.Define('w_lumi_year','w_lumi*60.0*w_pu')
  df_simu = df_simu.Define('lead_photon_idmva','photon_idmva[0]')
  df_simu = df_simu.Define('lead_photon_pt','photon_pt[0]')
  df_simu = df_simu.Define('lead_photon_eta','photon_eta[0]')
  df_simu = df_simu.Define('lead_photon_abseta','fabs(photon_eta[0])')

  #reweight kinematics
  prew_data_hist_ptrs = book_hist2d_zpeak_bkgsub(df_data, 
      ('prew_data_hist','',45,15.0,75.0,20,0.0,2.5), 'lead_photon_pt',
      'lead_photon_abseta', 'w_lumi')
  prew_simu_hist_ptrs = book_hist2d_zpeak_bkgsub(df_simu, 
      ('prew_data_hist','',45,15.0,75.0,20,0.0,2.5), 'lead_photon_pt',
      'lead_photon_abseta', 'w_lumi_year')
  prew_data_hist = get_bkgsub_hist(prew_data_hist_ptrs)
  prew_simu_hist = get_bkgsub_hist(prew_simu_hist_ptrs)
  #prew_data_hist_ptr = df_data.Filter('mmy_mass>80&&mmy_mass<98').Histo2D(('prew_data_hist','',45,15.0,75.0,20,0.0,2.5),'lead_photon_pt','lead_photon_abseta')
  #prew_bkgd_hist_ptr = df_data.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo2D(('prew_bkgd_hist','',45,15.0,75.0,20,0.0,2.5),'lead_photon_pt','lead_photon_abseta')
  #prew_simu_hist_ptr = df_simu.Filter('mmy_mass>80&&mmy_mass<98').Histo2D(('prew_simu_hist','',45,15.0,75.0,20,0.0,2.5),'lead_photon_pt','lead_photon_abseta','w_lumi_year')
  #prew_smbg_hist_ptr = df_simu.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo2D(('prew_smbg_hist','',45,15.0,75.0,20,0.0,2.5),'lead_photon_pt','lead_photon_abseta','w_lumi_year')
  #prew_data_hist = prew_data_hist_ptr.GetValue()
  #prew_bkgd_hist = prew_bkgd_hist_ptr.GetValue()
  #prew_simu_hist = prew_simu_hist_ptr.GetValue()
  #prew_smbg_hist = prew_smbg_hist_ptr.GetValue()
  #prew_data_hist.Add(prew_bkgd_hist,-1.0*18.0/4.0)
  #prew_simu_hist.Add(prew_smbg_hist,-1.0*18.0/4.0)
  prew_data_hist.Scale(1.0/prew_data_hist.Integral())
  prew_simu_hist.Scale(1.0/prew_simu_hist.Integral())
  prew_data_hist.Divide(prew_simu_hist)
  ROOT.gDirectory.Add(prew_data_hist)

  ROOT.gInterpreter.Declare("""
  TH1D* prew_data_hist = static_cast<TH1D*>(gDirectory->Get("prew_data_hist"));
  """)

  df_simu = df_simu.Define('w_pteta','prew_data_hist->GetBinContent(prew_data_hist->FindBin(photon_pt[0],fabs(photon_eta[0])))')
  #df_simu = df_simu.Define('w_pteta','w_ptetaog*w_ptetaog')
  #df_simu = df_simu.Define('w_pteta','1')
  df_simu = df_simu.Define('w_lumi_pteta','w_lumi*w_pteta*w_pu*60.0')
  df_simu = df_simu.Define('w_lumi_pteta_photon','w_lumi_photon*w_pteta*60.0')

  data_hist_ptrs = []
  simu_hist_ptrs = []
  rwgt_hist_ptrs = []

  #df_data = df_data.Filter('photon_pt[0]>30')
  #df_simu = df_simu.Filter('photon_pt[0]>30')

  data_hist_ptrs = book_hist1d_zpeak_bkgsub(df_data, 
      ('phid_data','',20,-0.58,1.0), 'lead_photon_idmva', 'w_lumi')
  simu_hist_ptrs = book_hist1d_zpeak_bkgsub(df_simu, 
      ('phid_data','',20,-0.58,1.0), 'lead_photon_idmva', 'w_lumi_pteta')
  rwgt_hist_ptrs = book_hist1d_zpeak_bkgsub(df_simu, 
      ('phid_data','',20,-0.58,1.0), 'lead_photon_idmva', 'w_lumi_pteta_photon')

  #for iid in range(len(idmva_bins)-1):
  #  df_data_bin = df_data.Filter('photon_idmva[0]>'+str(idmva_bins[iid])
  #                               +'&&photon_idmva[0]<='+str(idmva_bins[iid+1]))
  #  df_simu_bin = df_simu.Filter('photon_idmva[0]>'+str(idmva_bins[iid])
  #                               +'&&photon_idmva[0]<='+str(idmva_bins[iid+1]))
  #  bin_string = '_idbin'+str(iid)
  #  data_hist_ptrs.append(df_data_bin.Histo1D(
  #      ('validation_mmumuph_data'+bin_string,'',35,75.0,110.0),'mmy_mass'))
  #  simu_hist_ptrs.append(df_simu_bin.Histo1D(
  #      ('validation_mmumuph_simu'+bin_string,'',35,75.0,110.0),'mmy_mass','w_lumi_pteta'))
  #  rwgt_hist_ptrs.append(df_simu_bin.Histo1D(
  #      ('validation_mmumuph_rwgt'+bin_string,'',35,75.0,110.0),'mmy_mass','w_lumi_pteta_photon'))

  ##sanity check plots to debug
  #debug_simu_hist1_ptr = df_simu.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('simu_idmva_hist','MC (unweighted);IDMVA',40,-0.6,1.0),'lead_photon_idmva','w_lumi_pteta')
  #debug_rwgt_hist1_ptr = df_simu.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('rwgt_idmva_hist','MC (weighted);IDMVA',40,-0.6,1.0),'lead_photon_idmva','w_lumi_pteta_photon')
  #debug_data_hist1_ptr = df_data.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('data_idmva_hist','Data;IDMVA',40,-0.6,1.0),'lead_photon_idmva')

  #debug_simu_hist2_ptr = df_simu.Filter('mmy_mass>80&&mmy_mass<98').Histo1D(('simu_pt_hist','MC (unweighted);p_{T} [GeV]',45,15.0,75.0),'lead_photon_pt','w_lumi_pteta')
  ##debug_rwgt_hist2_ptr = df_simu.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('rwgt_pt_hist','MC (weighted);p_{T} [GeV]',45,15.0,75.0),'lead_photon_pt','w_lumi_pteta_photon')
  #debug_data_hist2_ptr = df_data.Filter('mmy_mass>80&&mmy_mass<98').Histo1D(('data_pt_hist','Data;p_{T} [GeV]',45,15.0,75.0),'lead_photon_pt')
  #debug_simu_hist4_ptr = df_simu.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo1D(('simu_pt_hist','MC (unweighted);p_{T} [GeV]',45,15.0,75.0),'lead_photon_pt','w_lumi_pteta')
  ##debug_rwgt_hist2_ptr = df_simu.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo1D(('rwgt_pt_hist','MC (weighted);p_{T} [GeV]',45,15.0,75.0),'lead_photon_pt','w_lumi_pteta_photon')
  #debug_data_hist4_ptr = df_data.Filter('(mmy_mass>78&&mmy_mass<80)||(mmy_mass>98&&mmy_mass<100)').Histo1D(('data_pt_hist','Data;p_{T} [GeV]',45,15.0,75.0),'lead_photon_pt')

  #debug_simu_hist3_ptr = df_simu.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('simu_eta_hist','MC (unweighted);#eta',40,-2.5,2.5),'lead_photon_eta','w_lumi_pteta')
  #debug_rwgt_hist3_ptr = df_simu.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('rwgt_eta_hist','MC (weighted);#eta',40,-2.5,2.5),'lead_photon_eta','w_lumi_pteta_photon')
  #debug_data_hist3_ptr = df_data.Filter('mmy_mass>85&&mmy_mass<95').Histo1D(('data_eta_hist','Data;#eta',40,-2.5,2.5),'lead_photon_eta')

  #debug_simu_hist1 = debug_simu_hist1_ptr.GetValue()
  #debug_rwgt_hist1 = debug_rwgt_hist1_ptr.GetValue()
  #debug_data_hist1 = debug_data_hist1_ptr.GetValue()
  #debug_simu_hist1.Scale(1.0/debug_simu_hist1.Integral())
  #debug_rwgt_hist1.Scale(1.0/debug_rwgt_hist1.Integral())
  #debug_data_hist1.Scale(1.0/debug_data_hist1.Integral())

  #debug_simu_hist2 = debug_simu_hist2_ptr.GetValue()
  ##debug_rwgt_hist2 = debug_rwgt_hist2_ptr.GetValue()
  #debug_data_hist2 = debug_data_hist2_ptr.GetValue()
  #debug_simu_hist4 = debug_simu_hist4_ptr.GetValue()
  #debug_data_hist4 = debug_data_hist4_ptr.GetValue()
  #debug_simu_hist2.Add(debug_simu_hist4, -18.0/4.0)
  #debug_data_hist2.Add(debug_data_hist4, -18.0/4.0)
  #debug_simu_hist2.Scale(1.0/debug_simu_hist2.Integral())
  ##debug_rwgt_hist2.Scale(1.0/debug_rwgt_hist2.Integral())
  #debug_data_hist2.Scale(1.0/debug_data_hist2.Integral())

  #debug_simu_hist3 = debug_simu_hist3_ptr.GetValue()
  #debug_rwgt_hist3 = debug_rwgt_hist3_ptr.GetValue()
  #debug_data_hist3 = debug_data_hist3_ptr.GetValue()
  #debug_simu_hist3.Scale(1.0/debug_simu_hist3.Integral())
  #debug_rwgt_hist3.Scale(1.0/debug_rwgt_hist3.Integral())
  #debug_data_hist3.Scale(1.0/debug_data_hist3.Integral())

  #debug1_plot = RplPlot()
  #debug1_plot.lumi_data = [(60,13)]
  #debug1_plot.y_title = '% Events/bin'
  #debug1_plot.plot_outline(debug_simu_hist1)
  #debug1_plot.plot_outline(debug_rwgt_hist1)
  #debug1_plot.plot_outline(debug_data_hist1)
  #debug1_plot.add_ratio('data_idmva_hist','simu_idmva_hist')
  #debug1_plot.add_ratio('data_idmva_hist','rwgt_idmva_hist')
  #debug1_plot.draw('plots/check_syst_idmva.pdf')

  #debug2_plot = RplPlot()
  #debug2_plot.lumi_data = [(60,13)]
  #debug2_plot.y_title = '% Events/bin'
  #debug2_plot.plot_outline(debug_simu_hist2)
  ##debug2_plot.plot_outline(debug_rwgt_hist2)
  #debug2_plot.plot_outline(debug_data_hist2)
  #debug2_plot.draw('plots/check_syst_pt.pdf')

  #debug3_plot = RplPlot()
  #debug3_plot.lumi_data = [(60,13)]
  #debug3_plot.y_title = '% Events/bin'
  #debug3_plot.plot_outline(debug_simu_hist3)
  #debug3_plot.plot_outline(debug_rwgt_hist3)
  #debug3_plot.plot_outline(debug_data_hist3)
  #debug3_plot.draw('plots/check_syst_eta.pdf')

  #data_fit_yield = []
  #simu_fit_yield = []
  #rwgt_fit_yield = []

  #for iid in range(len(idmva_bins)-1):
  #  data_hist = data_hist_ptrs[iid].GetValue()
  #  simu_hist = simu_hist_ptrs[iid].GetValue()
  #  rwgt_hist = rwgt_hist_ptrs[iid].GetValue()
  #  lo_bin = data_hist.FindBin(80.0)
  #  hi_bin = data_hist.FindBin(98.0)
  #  data_yield_lo = (data_hist.GetBinContent(lo_bin)+data_hist.GetBinContent(lo_bin-1))/2.0
  #  data_yield_hi = (data_hist.GetBinContent(hi_bin)+data_hist.GetBinContent(hi_bin+1))/2.0
  #  data_yield_z = data_hist.Integral(lo_bin+1,hi_bin-1)-(data_yield_lo+data_yield_hi)/2.0*(hi_bin-lo_bin-1)
  #  simu_yield_lo = (simu_hist.GetBinContent(lo_bin)+simu_hist.GetBinContent(lo_bin-1))/2.0
  #  simu_yield_hi = (simu_hist.GetBinContent(hi_bin)+simu_hist.GetBinContent(hi_bin+1))/2.0
  #  simu_yield_z = simu_hist.Integral(lo_bin+1,hi_bin-1)-(simu_yield_lo+simu_yield_hi)/2.0*(hi_bin-lo_bin-1)
  #  #simu_yield_z = simu_hist.Integral(lo_bin+1,hi_bin-1)
  #  rwgt_yield_lo = (rwgt_hist.GetBinContent(lo_bin)+rwgt_hist.GetBinContent(lo_bin-1))/2.0
  #  rwgt_yield_hi = (rwgt_hist.GetBinContent(hi_bin)+rwgt_hist.GetBinContent(hi_bin+1))/2.0
  #  rwgt_yield_z = rwgt_hist.Integral(lo_bin+1,hi_bin-1)-(rwgt_yield_lo+rwgt_yield_hi)/2.0*(hi_bin-lo_bin-1)
  #  #rwgt_yield_z = rwgt_hist.Integral(lo_bin+1,hi_bin-1)
  #  data_fit_yield.append(data_yield_z)
  #  simu_fit_yield.append(simu_yield_z)
  #  rwgt_fit_yield.append(rwgt_yield_z)

  #  debug_plot = RplPlot()
  #  debug_plot.lumi_data = [(60,13)]
  #  debug_plot.y_title = 'Events/bin'
  #  debug_plot.plot_outline(simu_hist)
  #  debug_plot.plot_outline(rwgt_hist)
  #  debug_plot.plot_outline(data_hist)
  #  debug_plot.draw('plots/debug_mllyfit_'+str(iid)+'.pdf')

  #make_comp_plot(data_fit_yield,simu_fit_yield,rwgt_fit_yield,region)

  ##save mlly histograms
  #output_file = ROOT.TFile('temp_zpeak_hists.root','RECREATE')

  #for iid in range(len(idmva_bins)-1):
  #  data_hist = data_hist_ptrs[iid].GetValue()
  #  simu_hist = simu_hist_ptrs[iid].GetValue()
  #  rwgt_hist = rwgt_hist_ptrs[iid].GetValue()
  #  data_hist.Write()
  #  simu_hist.Write()
  #  rwgt_hist.Write()

  #output_file.Close()

def make_comp_plot(data_fit_yield, simu_fit_yield, rwgt_fit_yield, region):
  '''Generate histograms
  '''
  data_hist = ROOT.TH1D('data_hist','Data;IDMVA;Events/bin',len(idmva_bins)-1,-0.58,1.0)
  simu_hist = ROOT.TH1D('simu_hist','Simulation (no photon ID weights);IDMVA;Events/bin',len(idmva_bins)-1,-0.58,1.0)
  rwgt_hist = ROOT.TH1D('rwgt_hist','Simulation (DNN photon ID weights);IDMVA;Events/bin',len(idmva_bins)-1,-0.58,1.0)
  for iid in range(len(idmva_bins)-1):
    data_hist.SetBinContent(iid,data_fit_yield[iid])
    simu_hist.SetBinContent(iid,simu_fit_yield[iid])
    rwgt_hist.SetBinContent(iid,rwgt_fit_yield[iid])
  data_hist.Scale(1.0/data_hist.Integral())
  rwgt_hist.Scale(1.0/rwgt_hist.Integral()) #check norm for now
  simu_hist.Scale(1.0/simu_hist.Integral())
  var_plot = RplPlot()
  var_plot.lumi_data = [(60,13)]
  var_plot.y_title = '% Events/bin'
  var_plot.y_max_lower = 1.24
  var_plot.y_min_lower = 0.76
  var_plot.plot_outline(simu_hist)
  var_plot.plot_outline(rwgt_hist)
  var_plot.plot_outline(data_hist)
  var_plot.add_ratio('data_hist','simu_hist')
  var_plot.add_ratio('data_hist','rwgt_hist')
  var_plot.draw('plots/validate_mumuph_'+region+'.pdf')
  var_plot_zoom = RplPlot()
  var_plot_zoom.lumi_data = [(60,13)]
  var_plot_zoom.y_title = '% Events/bin'
  var_plot_zoom.y_max_lower = 1.125
  var_plot_zoom.y_min_lower = 0.875
  var_plot_zoom.plot_outline(simu_hist)
  var_plot_zoom.plot_outline(rwgt_hist)
  var_plot_zoom.plot_outline(data_hist)
  var_plot_zoom.add_ratio('data_hist','simu_hist')
  var_plot_zoom.add_ratio('data_hist','rwgt_hist')
  var_plot_zoom.draw('plots/validate_mumuph_'+region+'_zoom.pdf')

def fit_plots(region):
  '''Use histograms from the previous step, and perform fits to the Z-peak to
  get efficiencies for data and MC before and after reweighting
  '''

  data_fit_yield = []
  simu_fit_yield = []
  rwgt_fit_yield = []

  #Do Z-peak fits to extract data uncertainties
  input_file = ROOT.TFile('temp_zpeak_hists.root','READ')

  for iid in range(len(idmva_bins)-1):
    bin_string = '_idbin'+str(iid)
    data_hist = input_file.Get('validation_mmumuph_data'+bin_string)
    simu_hist = input_file.Get('validation_mmumuph_simu'+bin_string)
    rwgt_hist = input_file.Get('validation_mmumuph_rwgt'+bin_string)

    if (data_hist.Integral() != 0.0):
      fit_result_data = do_tnp_fit(data_hist,'Gauss+CMSmmy',
          'plots/photon_corr/validation_syst_data_mumuph_tnp'+bin_string+'.pdf',
          'curve',True)
      integral_data_sb = integrate_tgraph(fit_result_data[1],76,99)
      integral_data_bk = integrate_tgraph(fit_result_data[2],76,99)
      data_fit_yield.append(integral_data_sb-integral_data_bk)
    else:
      data_fit_yield.append(0.0)

    if (simu_hist.Integral() != 0.0):
      fit_result_simu = do_tnp_fit(simu_hist,'Gauss+CMSmmy',
          'plots/photon_corr/validation_syst_simu_mumuph_tnp'+bin_string+'.pdf',
          'curve',True)
      integral_simu_sb = integrate_tgraph(fit_result_simu[1],76,99)
      integral_simu_bk = integrate_tgraph(fit_result_simu[2],76,99)
      simu_fit_yield.append(integral_simu_sb-integral_simu_bk)
    else:
      simu_fit_yield.append(0.0)

    if (rwgt_hist.Integral() != 0.0):
      fit_result_rwgt = do_tnp_fit(rwgt_hist,'Gauss+CMSmmy',
          'plots/photon_corr/validation_syst_rwgt_mumuph_tnp'+bin_string+'.pdf',
          'curve',True)
      integral_rwgt_sb = integrate_tgraph(fit_result_rwgt[1],76,99)
      integral_rwgt_bk = integrate_tgraph(fit_result_rwgt[2],76,99)
      rwgt_fit_yield.append(integral_rwgt_sb-integral_rwgt_bk)
    else:
      rwgt_fit_yield.append(0.0)

  #left here, print stuff in case plot script fails, then call
  print(data_fit_yield)
  print(simu_fit_yield)
  print(rwgt_fit_yield)

  make_comp_plot(data_fit_yield, simu_fit_yield, rwgt_fit_yield,region)

if __name__ == '__main__':

  make_photon_electron_plots()
  exit()

  argument_parser = ArgumentParser(prog='generate_photonid_dnn_syst',
  description='Script to validate DNN performance')
  argument_parser.add_argument('-s','--steps',default='make_plots,fit_plots,clean')
  argument_parser.add_argument('-d','--data_input_file')
  argument_parser.add_argument('-m','--mc_input_file')
  argument_parser.add_argument('-n','--nn_filename')
  argument_parser.add_argument('-r','--region',choices=['barrel','endcap'])

  args = argument_parser.parse_args()
  steps = args.steps.split(',')

  dir_loc = args.nn_filename.rfind('/')
  class_name = args.nn_filename[:-4]
  if (dir_loc != -1):
    class_name = class_name[dir_loc+1:]
  ROOT.gInterpreter.ProcessLine('#include "'+args.nn_filename+'"')
  ROOT.gInterpreter.Declare("""
  const """+class_name+""" dnn;

  float get_sf(std::vector<float> vars) {
    float dnn_output = dnn.evaluate(vars);
    if (dnn_output >= 1.0)
      return 1.0;
    float sf = dnn_output/(1.0-dnn_output);
    if (sf > 5.0) return 5.0;
    return sf;
  }
  """)

  if 'make_photon_electron_plots' in steps:
    make_photon_electron_plots()

  if 'make_plots' in steps:
    make_plots(args.data_input_file, args.mc_input_file, args.region)
    #debug_plots(args.data_input_file, args.mc_input_file, args.region)

  if 'fit_plots' in steps:
    temp_debug(args.mc_input_file)
    #fit_plots(args.region)
    #data_fit_yield = [0.0, 0.0, 579.0648207271413, 784.4951034775838, 1144.2956918699285, 964.4596853130047, 1260.6594769097753, 1340.5484190846332, 1557.3681716839894, 1601.051539440746, 1996.316876420574, 2333.703534428538, 2798.162609953051, 3737.14424445255, 4695.40474590765, 5859.345806403799, 8617.844150502482, 13731.534036297991, 26138.76381445412, 34675.245456166595]
    #simu_fit_yield = [0.0, 0.0, 490.4160594099652, 817.1847566524309, 901.317393175006, 1046.9478151260214, 1110.0570512321076, 1185.9233544646024, 1411.4135592923894, 1532.235715347085, 1844.5087476760903, 2098.014880061689, 2632.080297372834, 3317.783891271697, 4271.505606971456, 5683.369276378562, 8589.814310676786, 14019.66311278803, 26943.797974815232, 38437.87217340166]
    #rwgt_fit_yield = [0.0, 0.0, 459.97550774590025, 789.3253232500574, 877.5097329358783, 1145.392581573074, 1080.0026242466174, 1164.9796515933886, 1395.0005364493204, 1506.0032584503751, 1826.56510222929, 2066.619580607604, 2607.1780322314216, 3291.510933783713, 4219.334888855836, 5612.650037610523, 8421.56161060106, 13569.111208999859, 25268.85377772824, 33430.49231438058]
    #make_comp_plot(data_fit_yield,simu_fit_yield,rwgt_fit_yield)

  if 'clean' in steps:
    subprocess.run('rm temp_zpeak_hists.root'.split())


