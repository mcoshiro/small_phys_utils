'''
Script that validates photon DNN-based SFs
'''

from argparse import ArgumentParser
from array import array
import ROOT

ROOT.gInterpreter.AddIncludePath('inc/')
ROOT.gInterpreter.ProcessLine('#include "nn_model_phid2018.hpp"')
ROOT.gInterpreter.Declare("""

const nn_model_phid2018 dnn;

float get_sf(float var1, float var2, float var3, float var4, float var5, 
             float var6, float var7, float var8, float var9, float var10,
             float var11, float var12, float var13, float var14, float var15,
             float var16) {
  return dnn.get_sf({var1, var2, var3, var4, var5, var6, var7, var8, var9, 
                     var10, var11, var12, var13, var14, var15, var16});
}

""")

if __name__=='__main__':

  #parse arguments
  #argument_parser = ArgumentParser(prog='validate_sfs',
  #    description='Generates some validation histograms')
  #argument_parser.add_argument('-n','--num_filename')
  #argument_parser.add_argument('-d','--den_filename')
  #argument_parser.add_argument('-o','--output_filename',default='pteta_hist.hpp')
  #args = argument_parser.parse_args()

  df_simu = ROOT.RDataFrame('tree','/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_mc_2018.root').Filter('pair_mass>81&&pair_mass<101')
  #df_simu = ROOT.RDataFrame('tree','/homes/oshiro/temp_cms26notworkingsoihadtocopyntupleshere/photonidskim_mc_2018.root').Filter('pair_mass>81&&pair_mass<101')
  df_simu = df_simu.Define('sf','get_sf(ph_r9,ph_s4,ph_sc_etaWidth,ph_sc_phiWidth,ph_sieie,ph_sieip,ph_phoIso,ph_chIso,ph_chWorIso,ph_sc_rawEnergy,ph_sc_eta,event_rho,ph_ESsigma,ph_esEnergyOverRawE,ph_et,ph_sc_abseta)')
  df_simu = df_simu.Define('rwsf','w*sf')
  df_data = ROOT.RDataFrame('tree','/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_data_2018.root').Filter('pair_mass>81&&pair_mass<101')
  #df_data = ROOT.RDataFrame('tree','/homes/oshiro/temp_cms26notworkingsoihadtocopyntupleshere/photonidskim_data_2018.root').Filter('pair_mass>81&&pair_mass<101')

  variables = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','event_rho','ph_ESsigma','ph_esEnergyOverRawE','ph_mva94XV2','ph_et','ph_sc_abseta']
  xtitle = ['Photon R_{9}','Photon s_{4}','Photon supercluster #eta width','Photon supercluster #phi width','Photon #sigma_{i#eta i#eta}','Photon #sigma_{i#eta i#phi}','Photon I_{abs}(photon) [GeV]','Photon I_{abs}(charged, PV) [GeV]','Photon I_{abs}(charged, worst vertex) [GeV]','Photon supercluster E_{raw} [GeV]','Photon supercluster #eta','Event #rho [GeV]','Photon preshower #sigma_{eff}','Photon E_{ES}/E_{raw}','Photon Fall17v2 ID score','Photon p_{T} [GeV]','Photon supercluster |#eta|']
  #ranges = [(0.2,1.0),(0.3,1.0),(0.0025,0.0225),(0.0,0.1),(0.005,0.03),(-0.002,0.002),(0.0,10.0),(0.0,7.5),(0.0,25.0),(0.0,250.0),(-2.5,2.5),(0.0,45.0),(0.0,10.0),(0.0,0.2),(-1.0,1.0),(0.0,100.0),(0.0,2.5)]
  ranges = [(0.4,1.0),(0.4,0.975),(0.004,0.0225),(0.005,0.1),(0.006,0.03),(-0.002,0.002),(0.0,6.0),(0.0,6.0),(0.0,25.0),(0.0,250.0),(-2.5,2.5),(0.0,45.0),(0.0,10.0),(0.0,0.16),(-0.6,1.0),(0.0,100.0),(0.0,2.5)]
  is_log = [False, False, False, False, False, True, False, True, True, False, False, False, True, True, False, False, False]
  mcog_hist_ptrs = []
  mcsf_hist_ptrs = []
  data_hist_ptrs = []
  for ivar in range(len(variables)):
    mcog_hist_ptrs.append(df_simu.Histo1D((variables[ivar]+'_hist','',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w'))
    mcsf_hist_ptrs.append(df_simu.Histo1D((variables[ivar]+'_hist','',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'rwsf'))
    data_hist_ptrs.append(df_data.Histo1D((variables[ivar]+'_hist','',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar]))

  weight_histo = df_simu.Histo1D(('sf_hist','',30,0,2),'sf')
  can_weight = ROOT.TCanvas('c_sf','c',600,600)
  weight_histo.GetValue().Draw()
  can_weight.SaveAs('plots/validate_hist_sf.pdf')

  #save plot
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetErrorX(0)
  for ivar in range(len(variables)):
    can = ROOT.TCanvas('c_'+variables[ivar],'c',600,600)
    can.SetMargin(0,0,0,0)
    can.SetFillStyle(4000)
    mcog_hist = mcog_hist_ptrs[ivar].GetValue()
    mcsf_hist = mcsf_hist_ptrs[ivar].GetValue()
    data_hist = data_hist_ptrs[ivar].GetValue()
    ##dummy test
    #mcog_hist = ROOT.TH1D('ph_r9_hist_0','',30,0.0,1.0)
    #mcsf_hist = ROOT.TH1D('ph_r9_hist_1','',30,0.0,1.0)
    #data_hist = ROOT.TH1D('ph_r9_hist_2','',30,0.0,1.0)
    #rng = ROOT.TRandom3()
    #for i in range(500):
    #  mcog_hist.Fill(rng.Gaus(0.5,0.1))
    #  mcsf_hist.Fill(rng.Gaus(0.5,0.1))
    #  data_hist.Fill(rng.Gaus(0.5,0.1))
    ##/dummy test
    mcog_hist.Scale(1.0/mcog_hist.Integral())
    mcsf_hist.Scale(1.0/mcsf_hist.Integral())
    data_hist.Scale(1.0/data_hist.Integral())

    #top histogram settings
    mcog_hist.SetLineColor(ROOT.kBlue)
    mcog_hist.SetLineWidth(2)
    mcsf_hist.SetLineColor(ROOT.kRed)
    mcsf_hist.SetLineWidth(2)
    data_hist.SetMarkerStyle(ROOT.kFullCircle)
    data_hist.SetMarkerColor(ROOT.kBlack)
    data_hist.SetLineColor(ROOT.kBlack)

    #top plot settings
    mcog_hist.SetLabelSize(0,'x')
    mcog_hist.SetLabelSize(0.035,'y')
    mcog_hist.GetYaxis().SetTitle('% "Photons"/bin')
    mcog_hist.GetYaxis().SetNdivisions(606)
    if is_log[ivar]:
      mcog_hist.SetMaximum(10.0)
    else:
      mcog_hist.SetMinimum(0.0)
      mcog_hist.SetMaximum(max(mcog_hist.GetMaximum(),mcsf_hist.GetMaximum(),data_hist.GetMaximum())*1.5)

    #generate bottom histograms
    ratio_mcsf = data_hist.Clone()
    ratio_mcog = data_hist.Clone()
    ratio_mcsf.Divide(mcsf_hist)
    ratio_mcog.Divide(mcog_hist)

    #bottom histogram settings
    ratio_mcsf.SetLineColor(ROOT.kRed)
    ratio_mcsf.SetLineWidth(2)
    ratio_mcsf.SetMarkerColor(ROOT.kRed)
    ratio_mcog.SetLineColor(ROOT.kBlue)
    ratio_mcog.SetMarkerColor(ROOT.kBlue)
    ratio_mcog.SetLineWidth(2)

    #bottom plot settings
    ratio_mcsf.GetYaxis().SetTitle('Data/MC')
    ratio_mcsf.GetYaxis().SetNdivisions(606)
    ratio_mcsf.GetXaxis().SetTitle(xtitle[ivar])
    ratio_mcsf.SetLabelSize(0.035,'y')
    ratio_mcsf.SetMinimum(0.75)
    ratio_mcsf.SetMaximum(1.25)

    #draw everything
    top_pad = ROOT.TPad('top_pad_'+variables[ivar],'',0.0,0.0,1.0,1.0)
    bot_pad = ROOT.TPad('bot_pad_'+variables[ivar],'',0.0,0.0,1.0,1.0)
    #top_pad.SetTicks(1,1)
    #bot_pad.SetTicks(1,1)
    top_pad.SetMargin(0.19,0.06,0.3,0.05)
    bot_pad.SetMargin(0.19,0.06,0.1,0.7)
    top_pad.SetFillStyle(4000)
    bot_pad.SetFillStyle(4000)
    if is_log[ivar]:
      top_pad.SetLogy()
    top_pad.Draw()
    top_pad.cd()
    mcog_hist.Draw('HIST')
    mcsf_hist.Draw('HIST SAME')
    data_hist.Draw('P SAME')
    leg = ROOT.TLegend(0.24,0.78,0.4,0.93)
    leg.AddEntry(mcog_hist,'MC (No SFs)','f')
    leg.AddEntry(mcsf_hist,'MC (With SFs)','f')
    leg.AddEntry(data_hist,'Data','ep')
    leg.SetBorderSize(0)
    leg.Draw('SAME')
    label = ROOT.TLatex()
    label.SetTextSize(0.035)
    label.SetNDC(ROOT.kTRUE)
    label.SetTextAlign(11)
    label.DrawLatex(0.20,0.96,"#font[62]{CMS} #scale[0.8]{#font[52]{Work in Progress}}")
    label.SetTextAlign(31)
    label.SetTextSize(0.03)
    label.DrawLatex(0.93,0.96,"#font[42]{60 fb^{-1} (13 TeV)}")
    top_pad.Modified()
    can.cd()
    bot_pad.Draw('SAME')
    bot_pad.cd()
    ratio_mcsf.Draw()
    ratio_mcog.Draw('SAME')
    line = ROOT.TLine(ranges[ivar][0],1.0,ranges[ivar][1],1.0)
    line.SetNDC(ROOT.kFALSE)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw('SAME')
    bot_pad.Modified()
    can.Draw()
    can.SaveAs('plots/validate_hist_'+variables[ivar]+'.pdf')

