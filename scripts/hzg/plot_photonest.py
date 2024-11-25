#!/usr/bin/env python3
"""@package docstring
Make plots related to photon weights
"""

from root_plot_lib import RplPlot
import ROOT

cpp_snippet = """
TH2D* lo_ratio_hist = static_cast<TH2D*>(gDirectory->Get("lo_ratio_hist"));
TH2D* hi_ratio_hist = static_cast<TH2D*>(gDirectory->Get("hi_ratio_hist"));
"""

if __name__ == '__main__':

  #what to plot
  variables = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','event_rho','ph_ESsigma','ph_esEnergyOverRawE','ph_mva94XV2','ph_et','ph_sc_abseta']
  xtitle = ['Photon R_{9}','Photon s_{4}','Photon supercluster #eta width','Photon supercluster #phi width','Photon #sigma_{i#eta i#eta}','Photon #sigma_{i#eta i#phi}','Photon I_{abs}(photon) [GeV]','Photon I_{abs}(charged, PV) [GeV]','Photon I_{abs}(charged, worst vertex) [GeV]','Photon supercluster E_{raw} [GeV]','Photon supercluster #eta','Event #rho [GeV]','Photon preshower #sigma_{eff}','Photon E_{ES}/E_{raw}','Photon Fall17v2 ID score','Photon p_{T} [GeV]','Photon supercluster |#eta|']
  ranges = [(0.4,1.0),(0.4,0.975),(0.004,0.0225),(0.005,0.1),(0.006,0.03),(-0.002,0.002),(0.0,6.0),(0.0,6.0),(0.0,25.0),(0.0,250.0),(-2.5,2.5),(0.0,45.0),(0.0,10.0),(0.0,0.16),(-0.6,1.0),(0.0,100.0),(0.0,2.5)]
  is_log = [False, False, False, False, False, True, False, True, True, False, False, False, True, True, False, False, False]

  #ROOT.EnableImplicitMT()
  df = ROOT.RDataFrame('tree','/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_data_2018_new.root').Filter('ph_et<20')
  lo_hist_ptrs = []
  md_hist_ptrs = []
  hi_hist_ptrs = []

  hist_ptr = df.Filter('ph_et>15&&ph_et<17.5&&ph_sc_abseta<0.25').Histo1D(('test_mll_hist','',80,50.0,130.0),'pair_mass')
  hist = hist_ptr.GetValue()
  out_file = ROOT.TFile('test_hist_file.root','RECREATE')
  hist.Write()
  out_file.Close()
  exit()

  #calculate weights
  print('Calculating weights.')
  df_lo = df.Filter('pair_mass<81')
  df_md = df.Filter('pair_mass>81&&pair_mass<101')
  df_hi = df.Filter('pair_mass>101')
  lo_pt_hist_ptr = df_lo.Histo2D(('lo_pt_hist','',5,15.0,20.0,20,0,2.5),'ph_et','ph_sc_abseta')
  md_pt_hist_ptr = df_md.Histo2D(('md_pt_hist','',5,15.0,20.0,20,0,2.5),'ph_et','ph_sc_abseta')
  hi_pt_hist_ptr = df_hi.Histo2D(('hi_pt_hist','',5,15.0,20.0,20,0,2.5),'ph_et','ph_sc_abseta')
  lo_pt_hist = lo_pt_hist_ptr.GetValue()
  md_pt_hist = md_pt_hist_ptr.GetValue()
  hi_pt_hist = hi_pt_hist_ptr.GetValue()
  lo_pt_hist.Scale(1.0/lo_pt_hist.Integral())
  md_pt_hist.Scale(1.0/md_pt_hist.Integral())
  hi_pt_hist.Scale(1.0/hi_pt_hist.Integral())
  lo_ratio_hist = md_pt_hist.Clone('lo_ratio_hist')
  hi_ratio_hist = md_pt_hist.Clone('hi_ratio_hist')
  lo_ratio_hist.Divide(lo_pt_hist)
  hi_ratio_hist.Divide(hi_pt_hist)
  ROOT.gInterpreter.Declare(cpp_snippet)
  df_lo = df_lo.Define('w_pteta','lo_ratio_hist->GetBinContent(lo_ratio_hist->FindBin(ph_et,ph_sc_abseta))')
  df_hi = df_hi.Define('w_pteta','hi_ratio_hist->GetBinContent(hi_ratio_hist->FindBin(ph_et,ph_sc_abseta))')
  
  #book histograms
  print('Generating histograms.')
  for ivar in range(len(variables)):
    lo_hist_ptrs.append(df_lo.Histo1D((
        'low_'+variables[ivar]+'_hist','m_{ll}<81 GeV; '+xtitle[ivar]+';',30,
        ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_pteta'))
    md_hist_ptrs.append(df_md.Histo1D((
        'mid_'+variables[ivar]+'_hist','81<m_{ll}<101 GeV; '+xtitle[ivar]+';',30,
        ranges[ivar][0],ranges[ivar][1]),variables[ivar]))
    hi_hist_ptrs.append(df_hi.Histo1D((
        'hi_'+variables[ivar]+'_hist','m_{ll}>101 GeV; '+xtitle[ivar]+';',30,
        ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_pteta'))

  #draw output
  for ivar in range(len(variables)):
    lo_hist = lo_hist_ptrs[ivar].GetValue()
    md_hist = md_hist_ptrs[ivar].GetValue()
    hi_hist = hi_hist_ptrs[ivar].GetValue()
    lo_hist.Scale(1.0/lo_hist.Integral())
    md_hist.Scale(1.0/md_hist.Integral())
    hi_hist.Scale(1.0/hi_hist.Integral())

    plot = RplPlot()
    plot.lumi_data = [(60,13)]
    plot.y_title = '% Events/bin'
    plot.log_y = is_log[ivar]
    plot.plot_outline(lo_hist)
    plot.plot_outline(md_hist)
    plot.plot_outline(hi_hist)
    plot.draw('plots/'+variables[ivar]+'_mllcomp.pdf')

  #ph_data_file = ROOT.TFile('/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_data_2018_new.root','READ')
  #ph_data_tree = ph_data_file.Get('tree')
  #ph_data_tree.Draw('pair_mass>>data_mll_hist(80,50,130)','ph_et>15&&ph_et<20')
  #hist_data = ROOT.gDirectory.Get('data_mll_hist')
  #hist_data.SetTitle('Data;m_{ll} [GeV];% Events/bin')
  ##hist_data.Scale(1.00/hist_data.Integral())

  #ph_simu_file = ROOT.TFile('/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_mc_2018.root','READ')
  #ph_simu_tree = ph_simu_file.Get('tree')
  #ph_simu_tree.Draw('pair_mass>>simu_mll_hist(80,50,130)','ph_et>15&&ph_et<20')
  #hist_simu = ROOT.gDirectory.Get('simu_mll_hist')
  #hist_simu.SetTitle('MC;m_{ll} [GeV];% Events/bin')
  #hist_simu.Scale(0.602993388*hist_data.Integral()/hist_simu.Integral())

  #print('Integral difference in (81,101):'+str(hist_data.Integral(31,51)-hist_simu.Integral(31,51)))

  #plot = RplPlot()
  #plot.lumi_data = [(60,13)]
  #plot.plot_outline(hist_simu)
  #plot.plot_outline(hist_data)
  #plot.draw('plots/pair_mass_test.pdf')
