#!/usr/bin/env python3
"""@package docstring
Make plots related to photon weights
"""

from root_plot_lib import RplPlot
from array import array
import ROOT
import gc

MC_FILE = '/data2/oshiro/ntuples/triggervalidate/*.root'

DATA_FILES = ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_SingleMuon*.root',
              '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root']

LEAD_BINS = [0.0,20.0,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,47.5,50.0,52.5,55.0,57.5,60.0,65.0,70.0,90.0]
SUBL_BINS = [0.0,10.0,15.0,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,50.0,60.0,90.0]

if __name__ == '__main__':

  #generate histograms
  ROOT.EnableImplicitMT()
  df_sim = ROOT.RDataFrame('tree', MC_FILE)
  df_dat = ROOT.RDataFrame('tree', DATA_FILES)

  mc_norm = df_sim.Sum('w_lumi')
  mc_norm_str = str(mc_norm.GetValue())
  print(mc_norm_str)

  df_sim = df_sim.Filter('nll>=1&&nmu==2&&ll_m[0]>80&&ll_m[0]<100')
  df_dat = df_dat.Filter('nll>=1&&nmu==2&&ll_m[0]>80&&ll_m[0]<100')
  df_sim = df_sim.Define('lead_mu_pt','mu_pt[0]')
  df_sim = df_sim.Define('subl_mu_pt','mu_pt[1]')
  df_dat = df_dat.Define('lead_mu_pt','mu_pt[0]')
  df_dat = df_dat.Define('subl_mu_pt','mu_pt[1]')
  df_sim = df_sim.Define('w_lumi_year','w_lumi*6077.22*60000.0/'+mc_norm_str)
  df_sim = df_sim.Define('w_lumi_year_singlemu','w_lumi_year*w_trig_singlemu')
  df_sim = df_sim.Define('w_lumi_year_doublemu','w_lumi_year*w_trig_doublemu')
  df_sim = df_sim.Define('w_lumi_year_trig','w_lumi_year*w_trig')

  leadmu_dat_ptr = {}
  sublmu_dat_ptr = {}
  leadmu_sim_ptrs = {}
  sublmu_sim_ptrs = {}

  for trig, trig_string in [('trig_single_mu','simu'),
                            ('trig_double_mu','dimu'),
                            ('trig_single_mu||trig_double_mu','ormu')]:
    leadmu_sim_ptrs[trig_string] = []
    sublmu_sim_ptrs[trig_string] = []
    df_dat_trig = df_dat.Filter(trig)
    df_sim_trig = df_sim.Filter(trig)
    leadmu_dat_ptr[trig_string] = df_dat_trig.Histo1D(
        ('leadmu_hist_'+trig_string,'Data;Lead Muon p_{T} [GeV]',len(LEAD_BINS)-1,
        array('d',LEAD_BINS)),'lead_mu_pt')
    sublmu_dat_ptr[trig_string] = df_dat_trig.Histo1D(
        ('sublmu_hist_'+trig_string,'Data;Sublead Muon p_{T} [GeV]',len(SUBL_BINS)-1,
        array('d',SUBL_BINS)),'subl_mu_pt')
    for weight, w_string, sim_name in [('w_lumi_year','lumi','MC (nominal)'),
        ('w_lumi_year_singlemu','simu','MC (w_mu_trig)'),
        ('w_lumi_year_doublemu','dimu','MC (w_dimu_trig)'),
        ('w_lumi_year_trig','ormu','MC (w_combmu_trig)')]:
      leadmu_sim_ptrs[trig_string].append(df_sim_trig.Histo1D(
          ('leadmu_hist_{}_{}'.format(trig_string,w_string),
          sim_name+';Lead Muon p_{T}[GeV]',len(LEAD_BINS)-1,array('d',LEAD_BINS)),
          'lead_mu_pt',weight))
      sublmu_sim_ptrs[trig_string].append(df_sim_trig.Histo1D(
          ('sublmu_hist_{}_{}'.format(trig_string,w_string),
          sim_name+';Sublead Muon p_{T}[GeV]',len(SUBL_BINS)-1,array('d',SUBL_BINS)),
          'subl_mu_pt',weight))

  #write histograms to file for ease of use
  out_file = ROOT.TFile('temp_trighists.root','RECREATE')
  for trig_string in ['simu','dimu','ormu']:
    leadmu_dat_ptr[trig_string].Write()
    sublmu_dat_ptr[trig_string].Write()
    for iweight in range(4):
      leadmu_sim_ptrs[trig_string][iweight].Write()
      sublmu_sim_ptrs[trig_string][iweight].Write()
  out_file.Close()

  #generate sumw2
  for dat_ptr, sim_ptrs in [(leadmu_dat_ptr, leadmu_sim_ptrs),
                            (sublmu_dat_ptr, sublmu_sim_ptrs)]:
    for trig_string in ['simu','dimu','ormu']:
      dat_ptr[trig_string].Sumw2()
      for ibin in range(4):
        sim_ptrs[trig_string][ibin].Sumw2()

  #make plots
  gc.collect()
  gc.disable()
  for dat_ptr, sim_ptrs, leg in [
          (leadmu_dat_ptr, leadmu_sim_ptrs,'leadmu'),
          (sublmu_dat_ptr, sublmu_sim_ptrs,'sublmu')]:
    for trig_string in ['simu','dimu','ormu']:
      plot = RplPlot()
      plot.lumi_data = [(60,13)]
      plot.y_title = 'Events/bin'
      plot.plot_points(dat_ptr[trig_string].GetPtr(),ROOT.kBlack)
      plot.plot_outline(sim_ptrs[trig_string][0].GetPtr())
      plot.plot_outline(sim_ptrs[trig_string][1].GetPtr())
      plot.plot_outline(sim_ptrs[trig_string][2].GetPtr())
      plot.plot_outline(sim_ptrs[trig_string][3].GetPtr())
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_lumi'.format(leg,trig_string),False)
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_simu'.format(leg,trig_string),False)
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_dimu'.format(leg,trig_string),False)
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_ormu'.format(leg,trig_string),False)
      plot.draw('plots/trigappcomp_{}_{}.pdf'.format(leg,trig_string))
    #for trig_string, trig_desc in [('simu','Single muon trigger'),
    #                               ('dimu','Dimuon trigger')]:
    for trig_string, trig_desc in [('dimu','Dimuon trigger')]:
      plot = RplPlot()
      dat_hist = dat_ptr['ormu'].Clone('dathist')
      dat_hist.Divide(dat_ptr[trig_string].GetPtr())
      sim_hist_lumi = sim_ptrs['ormu'][0].Clone('simhist_lumi')
      sim_hist_lumi.Divide(sim_ptrs[trig_string][2].GetPtr())
      sim_hist_ormu = sim_ptrs['ormu'][3].Clone('simhist_ormu')
      sim_hist_ormu.Divide(sim_ptrs[trig_string][2].GetPtr())
      plot.lumi_data = [(60,13)]
      plot.y_title = 'OR of triggers/'+trig_desc
      plot.plot_points(dat_hist,ROOT.kBlack)
      plot.plot_outline(sim_hist_lumi)
      plot.plot_outline(sim_hist_ormu)
      plot.add_ratio('dathist','simhist_lumi')
      plot.add_ratio('dathist','simhist_ormu')
      plot.draw('plots/trigappcomp_ratio_{}_{}.pdf'.format(leg,trig_string))
  gc.enable()
