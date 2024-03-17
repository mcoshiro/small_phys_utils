'''
Script that validates photon DNN-based SFs
'''

from argparse import ArgumentParser
from array import array
from root_plot_lib import RplPlot
import ROOT

if __name__=='__main__':

  #parse arguments
  argument_parser = ArgumentParser(prog='validate_sfs',
      description='Generates validation histograms')
  argument_parser.add_argument('-m','--mc_filename')
  argument_parser.add_argument('-d','--data_filename')
  argument_parser.add_argument('-r','--region',choices=['barrel','endcap'])
  argument_parser.add_argument('-n','--nn_filename')
  args = argument_parser.parse_args()

  dir_loc = args.nn_filename.rfind('/')
  class_name = args.nn_filename[:-4]
  if (dir_loc != -1):
    class_name = class_name[dir_loc+1:]
  ROOT.gInterpreter.ProcessLine('#include "'+args.nn_filename+'"')
  ROOT.gInterpreter.Declare("""
  const """+class_name+""" dnn;

  float get_sf(std::vector<float> vars) {
    float dnn_output = dnn.evaluate(vars);
    return dnn_output/(1.0-dnn_output);
  }
  """)

  weight_definition = 'get_sf({ph_r9,ph_s4,ph_sc_etaWidth,ph_sc_phiWidth,ph_sieie,ph_sieip,ph_phoIso,ph_chIso,ph_chWorIso,ph_sc_rawEnergy,ph_sc_eta,ph_et,ph_sc_abseta,event_rho})'
  #dnn_definition = 'dnn.evaluate({ph_r9,ph_s4,ph_sc_etaWidth,ph_sc_phiWidth,ph_sieie,ph_sieip,ph_phoIso,ph_chIso,ph_chWorIso,ph_sc_rawEnergy,ph_sc_eta,ph_et,ph_sc_abseta,event_rho})'
  selection_string = '((event%10)>6)&&ph_sc_abseta<1.5'
  if args.region == 'endcap':
    weight_definition = 'get_sf({ph_r9,ph_s4,ph_sc_etaWidth,ph_sc_phiWidth,ph_sieie,ph_sieip,ph_phoIso,ph_chIso,ph_chWorIso,ph_sc_rawEnergy,ph_sc_eta,ph_ESsigma,ph_esEnergyOverRawE,ph_et,ph_sc_abseta,event_rho})'
    #dnn_definition = 'dnn.evaluate({ph_r9,ph_s4,ph_sc_etaWidth,ph_sc_phiWidth,ph_sieie,ph_sieip,ph_phoIso,ph_chIso,ph_chWorIso,ph_sc_rawEnergy,ph_sc_eta,ph_ESsigma,ph_esEnergyOverRawE,ph_et,ph_sc_abseta,event_rho})'
    selection_string = '((event%10)>6)&&ph_sc_abseta>1.5'

  df_simu = ROOT.RDataFrame('tree',args.mc_filename).Filter(selection_string)
  df_simu = df_simu.Define('w_photon',weight_definition)
  df_simu = df_simu.Define('w_all','w_pre*w_photon')
  df_data = ROOT.RDataFrame('tree',args.data_filename).Filter(selection_string)

  variables = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','event_rho','ph_mva94XV2','ph_et','ph_sc_abseta']
  xtitle = ['Photon R_{9}','Photon s_{4}','Photon supercluster #eta width','Photon supercluster #phi width','Photon #sigma_{i#eta i#eta}','Photon #sigma_{i#eta i#phi}','Photon I_{abs}(photon) [GeV]','Photon I_{abs}(charged, PV) [GeV]','Photon I_{abs}(charged, worst vertex) [GeV]','Photon supercluster E_{raw} [GeV]','Photon supercluster #eta','Event #rho [GeV]','Photon Fall17v2 ID score','Photon p_{T} [GeV]','Photon supercluster |#eta|']
  ranges = [(0.4,1.0),(0.4,0.975),(0.004,0.0225),(0.005,0.1),(0.006,0.03),(-0.002,0.002),(0.0,6.0),(0.0,6.0),(0.0,25.0),(0.0,250.0),(-2.5,2.5),(0.0,45.0),(-0.6,1.0),(0.0,100.0),(0.0,2.5)]
  is_log = [False, False, False, False, False, True, False, True, True, False, False, False, False, False, False]
  if args.region == 'endcap':
    variables += ['ph_ESsigma','ph_esEnergyOverRawE']
    xtitle += ['Photon preshower #sigma_{eff}','Photon E_{ES}/E_{raw}']
    ranges += [(0.0,10.0),(0.0,0.16)]
    is_log += [True, True]


  mcog_hist_ptrs = []
  mcsf_hist_ptrs = []
  data_hist_ptrs = []
  for ivar in range(len(variables)):
    mcog_hist_ptrs.append(df_simu.Histo1D((variables[ivar]+'_mcog_hist','Simulation (no photon weights)',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_pre'))
    mcsf_hist_ptrs.append(df_simu.Histo1D((variables[ivar]+'_mcsf_hist','Simulation (DNN photon weights)',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_all'))
    data_hist_ptrs.append(df_data.Histo1D((variables[ivar]+'_data_hist','Data',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_bkg'))
  weight_histo = df_simu.Histo1D(('sf_hist',';Photon weight;Events/bin',30,0,2),'w_photon')

  #save plots
  weight_plot = RplPlot()
  weight_plot.lumi_data = [(60,13)]
  weight_plot.plot_outline(weight_histo.GetValue())
  weight_plot.draw('plots/validate_sf_'+args.region+'.pdf')

  for ivar in range(len(variables)):
    mcog_hist = mcog_hist_ptrs[ivar].GetValue()
    mcsf_hist = mcsf_hist_ptrs[ivar].GetValue()
    data_hist = data_hist_ptrs[ivar].GetValue()
    mcog_hist.Scale(1.0/mcog_hist.Integral())
    mcsf_hist.Scale(1.0/mcsf_hist.Integral())
    data_hist.Scale(1.0/data_hist.Integral())

    var_plot = RplPlot()
    var_plot.lumi_data = [(60,13)]
    var_plot.x_title = xtitle[ivar]
    var_plot.y_title = '% Events/bin'
    var_plot.log_y = is_log[ivar]
    var_plot.y_max_lower = 1.125
    var_plot.y_min_lower = 0.875
    var_plot.plot_outline(mcog_hist)
    var_plot.plot_outline(mcsf_hist)
    var_plot.plot_outline(data_hist)
    var_plot.add_ratio(variables[ivar]+'_data_hist',variables[ivar]+'_mcog_hist')
    var_plot.add_ratio(variables[ivar]+'_data_hist',variables[ivar]+'_mcsf_hist')
    var_plot.draw('plots/validate_'+variables[ivar]+'_reweight_'+args.region+'.pdf')

    #can = ROOT.TCanvas('c_'+variables[ivar],'c',600,600)
    #can.SetMargin(0,0,0,0)
    #can.SetFillStyle(4000)

    ##top histogram settings
    #mcog_hist.SetLineColor(ROOT.kBlue)
    #mcog_hist.SetLineWidth(2)
    #mcsf_hist.SetLineColor(ROOT.kRed)
    #mcsf_hist.SetLineWidth(2)
    #data_hist.SetMarkerStyle(ROOT.kFullCircle)
    #data_hist.SetMarkerColor(ROOT.kBlack)
    #data_hist.SetLineColor(ROOT.kBlack)

    ##top plot settings
    #mcog_hist.SetLabelSize(0,'x')
    #mcog_hist.SetLabelSize(0.035,'y')
    #mcog_hist.GetYaxis().SetTitle('% "Photons"/bin')
    #mcog_hist.GetYaxis().SetNdivisions(606)
    #if is_log[ivar]:
    #  mcog_hist.SetMaximum(10.0)
    #else:
    #  mcog_hist.SetMinimum(0.0)
    #  mcog_hist.SetMaximum(max(mcog_hist.GetMaximum(),mcsf_hist.GetMaximum(),data_hist.GetMaximum())*1.5)

    ##generate bottom histograms
    #ratio_mcsf = data_hist.Clone()
    #ratio_mcog = data_hist.Clone()
    #ratio_mcsf.Divide(mcsf_hist)
    #ratio_mcog.Divide(mcog_hist)

    ##bottom histogram settings
    #ratio_mcsf.SetLineColor(ROOT.kRed)
    #ratio_mcsf.SetLineWidth(2)
    #ratio_mcsf.SetMarkerColor(ROOT.kRed)
    #ratio_mcog.SetLineColor(ROOT.kBlue)
    #ratio_mcog.SetMarkerColor(ROOT.kBlue)
    #ratio_mcog.SetLineWidth(2)

    ##bottom plot settings
    #ratio_mcsf.GetYaxis().SetTitle('Data/MC')
    #ratio_mcsf.GetYaxis().SetNdivisions(606)
    #ratio_mcsf.GetXaxis().SetTitle(xtitle[ivar])
    #ratio_mcsf.SetLabelSize(0.035,'y')
    #ratio_mcsf.SetMinimum(0.75)
    #ratio_mcsf.SetMaximum(1.25)

    ##draw everything
    #top_pad = ROOT.TPad('top_pad_'+variables[ivar],'',0.0,0.0,1.0,1.0)
    #bot_pad = ROOT.TPad('bot_pad_'+variables[ivar],'',0.0,0.0,1.0,1.0)
    ##top_pad.SetTicks(1,1)
    ##bot_pad.SetTicks(1,1)
    #top_pad.SetMargin(0.19,0.06,0.3,0.05)
    #bot_pad.SetMargin(0.19,0.06,0.1,0.7)
    #top_pad.SetFillStyle(4000)
    #bot_pad.SetFillStyle(4000)
    #if is_log[ivar]:
    #  top_pad.SetLogy()
    #top_pad.Draw()
    #top_pad.cd()
    #mcog_hist.Draw('HIST')
    #mcsf_hist.Draw('HIST SAME')
    #data_hist.Draw('P SAME')
    #leg = ROOT.TLegend(0.24,0.78,0.4,0.93)
    #leg.AddEntry(mcog_hist,'MC (No SFs)','f')
    #leg.AddEntry(mcsf_hist,'MC (With SFs)','f')
    #leg.AddEntry(data_hist,'Data','ep')
    #leg.SetBorderSize(0)
    #leg.Draw('SAME')
    #label = ROOT.TLatex()
    #label.SetTextSize(0.035)
    #label.SetNDC(ROOT.kTRUE)
    #label.SetTextAlign(11)
    #label.DrawLatex(0.20,0.96,"#font[62]{CMS} #scale[0.8]{#font[52]{Work in Progress}}")
    #label.SetTextAlign(31)
    #label.SetTextSize(0.03)
    #label.DrawLatex(0.93,0.96,"#font[42]{60 fb^{-1} (13 TeV)}")
    #top_pad.Modified()
    #can.cd()
    #bot_pad.Draw('SAME')
    #bot_pad.cd()
    #ratio_mcsf.Draw()
    #ratio_mcog.Draw('SAME')
    #line = ROOT.TLine(ranges[ivar][0],1.0,ranges[ivar][1],1.0)
    #line.SetNDC(ROOT.kFALSE)
    #line.SetLineStyle(2)
    #line.SetLineColor(ROOT.kBlack)
    #line.SetLineWidth(2)
    #line.Draw('SAME')
    #bot_pad.Modified()
    #can.Draw()
    #can.SaveAs('plots/validate_hist_'+variables[ivar]+'.pdf')

