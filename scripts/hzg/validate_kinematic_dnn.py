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
  argument_parser.add_argument('-n','--nn_filename')
  argument_parser.add_argument('-t','--tag',default='')
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
    float sf = dnn_output/(1.0-dnn_output);
    if (sf > 5.0) return 5.0;
    return sf;
  }
  """)

  weight_definition = ('get_sf({photon_idmva0, photon_mht_dphi, '
                       'njet0, ll_pt0, '
                       'photon_pt0, llphoton_pt0, mht0, ht0})')
  selection_string = '(rng>0.6)'
  era = 'run3'

  df_simu = ROOT.RDataFrame('tree',args.mc_filename).Filter(selection_string)
  df_simu = df_simu.Define('w_rw',weight_definition)
  df_simu = df_simu.Define('w_combined','w_rw*w_param')
  df_data = ROOT.RDataFrame('tree',args.data_filename).Filter(selection_string)
  df_simu = df_simu.Define('mht_over_ht','mht0/ht0')
  df_data = df_data.Define('mht_over_ht','mht0/ht0')

  bins = [('1','inclusive'),
          ('njet0>0.9&&njet0<1.01','njet1'),
          ('njet0>1.9','dijet')]

  variables = [
      ('photon_idmva0','Photon IDMVA',(0.0,1.0),False),
      ('photon_mht_dphi','#Delta #phi(#gamma,H_{T}^{miss})',(0.0,3.142),False),
      ('njet0','N_{jet}',(-0.5,5.5),False),
      ('ll_pt0','Dilepton p_{T} [GeV]',(0.0,150.0),False),
      ('photon_pt0','Photon p_{T} [GeV]',(0.0,80.0),False),
      ('llphoton_pt0','ll#gamma p_{T} [GeV]',(0.0,150.0),False),
      ('mht0','H_{T}^{miss} [GeV]',(0.0,100.0),False),
      ('ht0','H_{T} [GeV]',(0.0,500.0),False),
      ('mht_over_ht','H_{T}^{miss}/H_{T}',(0.0,1.0),False)]
  nbins = 25

  mcog_hist_ptrs = []
  mcsf_hist_ptrs = []
  data_hist_ptrs = []
  weight_histo = df_simu.Histo1D(
      ('sf_hist','; DNN weight;Events/bin',nbins,0,2.5),'w_rw')

  for target_bin in bins:
    mcog_hist_ptrs.append([])
    mcsf_hist_ptrs.append([])
    data_hist_ptrs.append([])
    df_simu_bin = df_simu.Filter(target_bin[0])
    df_data_bin = df_data.Filter(target_bin[0])
    for var_name, xitlte, var_range, is_log in variables:
      mcog_hist_ptrs[-1].append(df_simu_bin.Histo1D((var_name+'_mcog_hist',
          'MC (nonreweighted)',nbins,var_range[0],
          var_range[1]),var_name,'w_param'))
      mcsf_hist_ptrs[-1].append(df_simu_bin.Histo1D((var_name+'_mcsf_hist',
          'MC (DNN weights)',nbins,var_range[0],
          var_range[1]),var_name,'w_combined'))
      data_hist_ptrs[-1].append(df_data_bin.Histo1D((var_name+'_data_hist',
          'Data',nbins,var_range[0],var_range[1]),var_name,'w_param'))

  lumi_data = [(138,13)]
  if era=='run3':
    lumi_data = [(62,13.6)]
  elif era=='run2+run3':
    lumi_data = [(138,13),(62,13.6)]

  #save plots
  weight_plot = RplPlot()
  weight_plot.lumi_data = lumi_data
  weight_plot.plot_outline(weight_histo.GetValue())
  weight_plot.draw(f'plots/validate_w_dnn_reweight{args.tag}.pdf')

  ibin = 0
  for target_bin in bins:
    for ivar in range(len(variables)):
      var_name = variables[ivar][0]
      xtitle = variables[ivar][1]
      var_range = variables[ivar][2]
      is_log = variables[ivar][3]
      mcog_hist = mcog_hist_ptrs[ibin][ivar].GetValue()
      mcsf_hist = mcsf_hist_ptrs[ibin][ivar].GetValue()
      data_hist = data_hist_ptrs[ibin][ivar].GetValue()
      mcog_hist.Scale(1.0/mcog_hist.Integral())
      mcsf_hist.Scale(1.0/mcsf_hist.Integral())
      data_hist.Scale(1.0/data_hist.Integral())

      var_plot = RplPlot()
      var_plot.lumi_data = lumi_data
      var_plot.x_title = xtitle
      var_plot.y_title = '% Events/bin'
      var_plot.log_y = is_log
      #var_plot.y_max_lower = 1.125
      #var_plot.y_min_lower = 0.875
      var_plot.y_max_lower = 1.24
      var_plot.y_min_lower = 0.76
      var_plot.plot_outline(mcog_hist)
      var_plot.plot_outline(mcsf_hist)
      var_plot.plot_outline(data_hist)
      var_plot.add_ratio(var_name+'_data_hist',var_name+'_mcog_hist')
      var_plot.add_ratio(var_name+'_data_hist',var_name+'_mcsf_hist')
      var_plot.draw('plots/validate_{}_reweight_{}{}.pdf'.format(var_name,
          target_bin[1],args.tag))
    ibin += 1

