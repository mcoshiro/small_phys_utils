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
  args = argument_parser.parse_args()

  ROOT.gInterpreter.ProcessLine('#include <stdexcept>')

  selection_string = 'ph_sc_abseta<1.5'
  if args.region == 'endcap':
    selection_string = 'ph_sc_abseta>1.5'

  df_simu = ROOT.RDataFrame('tree',args.mc_filename).Filter(selection_string)
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
  mcog_hists = []
  data_hists = []
  mcsf_hists = []

  #calculate SFs
  for ivar in range(len(variables)):
    mcog_hist_ptrs.append(df_simu.Histo1D((variables[ivar]+'_mcog_hist','Simulation (no photon weights)',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_pre'))
    data_hist_ptrs.append(df_data.Histo1D((variables[ivar]+'_data_hist','Data',30,ranges[ivar][0],ranges[ivar][1]),variables[ivar],'w_bkg'))

  for ivar in range(len(variables)):
    mcog_hists.append(mcog_hist_ptrs[ivar].GetValue())
    data_hists.append(data_hist_ptrs[ivar].GetValue())
    data_hists[ivar].Divide(mcog_hists[ivar])
    ROOT.gDirectory.Add(data_hists[ivar])

    ROOT.gInterpreter.Declare("""
    TH1D* """+variables[ivar]+"""_data_hist = static_cast<TH1D*>(gDirectory->Get("""+'"'+variables[ivar]+"""_data_hist"));
    float get_sf_"""+variables[ivar]+"""(float value) {
      if ("""+variables[ivar]+"""_data_hist == nullptr) {
        throw std::runtime_error("Couldn't find hist for """+variables[ivar]+""" ");
      }
      return """+variables[ivar]+"""_data_hist->GetBinContent("""+variables[ivar]+"""_data_hist->FindBin(value));
    }
    """)

    df_simu = df_simu.Define('w_'+variables[ivar],
        'get_sf_'+variables[ivar]+'('+variables[ivar]+')*w_pre')

  #make output plots
  for ivar in range(len(variables)):
    mcsf_hist_ptrs.append(df_simu.Histo1D(('idhist_'+variables[ivar],'Simulation (weights for '+variables[ivar]+');IDMVA',30,-0.6,1.0),'ph_mva94XV2','w_'+variables[ivar]))
  mcog_hist_ptr = df_simu.Histo1D(('idhist_mcog','Simulation (unweighted);IDMVA',30,-0.6,1.0),'ph_mva94XV2','w_pre')
  data_hist_ptr = df_data.Histo1D(('idhist_data','Data;IDMVA',30,-0.6,1.0),'ph_mva94XV2','w_bkg')

  for ivar in range(len(variables)):
    mcsf_hists.append(mcsf_hist_ptrs[ivar].GetValue())
    mcsf_hists[ivar].Scale(1.0/mcsf_hists[ivar].Integral())
  mcog_hist = mcog_hist_ptr.GetValue()
  mcog_hist.Scale(1.0/mcog_hist.Integral())
  data_hist = data_hist_ptr.GetValue()
  data_hist.Scale(1.0/data_hist.Integral())

  #save plots
  for ivar in range(len(variables)):
    var_plot = RplPlot()
    var_plot.lumi_data = [(60,13)]
    var_plot.y_title = '% Events/bin'
    var_plot.y_max_lower = 1.125
    var_plot.y_min_lower = 0.875
    var_plot.plot_outline(mcog_hist)
    var_plot.plot_outline(mcsf_hists[ivar])
    var_plot.plot_outline(data_hist)
    var_plot.add_ratio('idhist_data','idhist_mcog')
    var_plot.add_ratio('idhist_data','idhist_'+variables[ivar])
    var_plot.draw('plots/check_'+variables[ivar]+'_reweight_'+args.region+'.pdf')

