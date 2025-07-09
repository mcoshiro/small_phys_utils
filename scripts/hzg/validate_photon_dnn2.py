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

  weight_definition = ('get_sf({ph_pt,ph_abseta,'
                      +'ph_idmva,static_cast<float>(ph_res)})')
  #weight_definition = ('get_sf({ph_idmva,static_cast<float>(ph_res)})')
  selection_string = '((rng%10)<2)'
  bins = [('1','inclusive'),
          ('ph_pt<25&&ph_abseta<1.5','pt15to25_eta0to1p5'),
          ('ph_pt<25&&ph_abseta>1.5','pt15to25_eta1p5to2p5'),
          ('ph_pt>25&&ph_abseta<1.5','pt25toInf_eta0to1p5'),
          ('ph_pt>25&&ph_abseta>1.5','pt25toInf_eta1p5to2p5')]

  df_simu = ROOT.RDataFrame('tree',args.mc_filename).Filter(selection_string)
  df_simu = df_simu.Define('w_photon',weight_definition)
  df_simu = df_simu.Define('w_combined','w_all*w_photon')
  df_simu = df_simu.Define('w_allbinned','w_all*w_binned')
  df_data = ROOT.RDataFrame('tree',args.data_filename).Filter(selection_string)
  df_simu = df_simu.Define('ph_res_inv','1.0/ph_res')
  df_data = df_data.Define('ph_res_inv','1.0/ph_res')

  #debug = df_simu.Display(['ph_pt','ph_abseta','ph_idmva','ph_res','w_binned','w_photon'],100)
  #debug.Print()
  #exit()

  variables = [
      ('ph_pt','Photon p_{T} [GeV]',(15.0,100.0),False),
      ('ph_abseta','Photon |#eta|',(0.0,2.5),False),
      #('npv','N_{PV}',(0.0,80.0),False),
      ('ph_idmva','Photon IDMVA',(0.0,1.0),False),
      ('ph_res_inv','Photon E/#sigma_{E}',(0.0,120.0),False),
      ('ph_res','Photon #sigma_{E}/E',(0.0,0.08),False)]
  nbins = 40

  mcog_hist_ptrs = []
  mcbn_hist_ptrs = []
  mcsf_hist_ptrs = []
  data_hist_ptrs = []
  weight_histo = df_simu.Histo1D(
      ('sf_hist',';Photon weight;Events/bin',50,0,2.5),'w_photon')

  for target_bin in bins:
    mcog_hist_ptrs.append([])
    mcbn_hist_ptrs.append([])
    mcsf_hist_ptrs.append([])
    data_hist_ptrs.append([])
    df_simu_bin = df_simu.Filter(target_bin[0])
    df_data_bin = df_data.Filter(target_bin[0])
    for var_name, xitlte, var_range, is_log in variables:
      mcog_hist_ptrs[-1].append(df_simu_bin.Histo1D((var_name+'_mcog_hist',
          'MC (nonreweighted)',nbins,var_range[0],
          var_range[1]),var_name,'w_all'))
      mcbn_hist_ptrs[-1].append(df_simu_bin.Histo1D((var_name+'_mcbn_hist',
          'MC (Binned weights)',nbins,var_range[0],
          var_range[1]),var_name,'w_allbinned'))
      mcsf_hist_ptrs[-1].append(df_simu_bin.Histo1D((var_name+'_mcsf_hist',
          'MC (DNN weights)',nbins,var_range[0],
          var_range[1]),var_name,'w_combined'))
      data_hist_ptrs[-1].append(df_data_bin.Histo1D((var_name+'_data_hist',
          'Data',nbins,var_range[0],var_range[1]),var_name,'w_all'))

  #save plots
  weight_plot = RplPlot()
  weight_plot.lumi_data = [(62,13.6)]
  weight_plot.plot_outline(weight_histo.GetValue())
  weight_plot.draw('plots/validate_w_photon_reweight.pdf')

  ibin = 0
  for target_bin in bins:
    for ivar in range(len(variables)):
      var_name = variables[ivar][0]
      xtitle = variables[ivar][1]
      var_range = variables[ivar][2]
      is_log = variables[ivar][3]
      mcog_hist = mcog_hist_ptrs[ibin][ivar].GetValue()
      mcbn_hist = mcbn_hist_ptrs[ibin][ivar].GetValue()
      mcsf_hist = mcsf_hist_ptrs[ibin][ivar].GetValue()
      data_hist = data_hist_ptrs[ibin][ivar].GetValue()
      mcog_hist.Scale(1.0/mcog_hist.Integral())
      mcbn_hist.Scale(1.0/mcbn_hist.Integral())
      mcsf_hist.Scale(1.0/mcsf_hist.Integral())
      data_hist.Scale(1.0/data_hist.Integral())

      var_plot = RplPlot()
      var_plot.lumi_data = [(62,13.6)]
      var_plot.x_title = xtitle
      var_plot.y_title = '% Events/bin'
      var_plot.log_y = is_log
      #var_plot.y_max_lower = 1.125
      #var_plot.y_min_lower = 0.875
      var_plot.y_max_lower = 1.24
      var_plot.y_min_lower = 0.76
      var_plot.plot_outline(mcog_hist)
      var_plot.plot_outline(mcsf_hist)
      var_plot.plot_outline(mcbn_hist)
      var_plot.plot_outline(data_hist)
      var_plot.add_ratio(var_name+'_data_hist',var_name+'_mcog_hist')
      var_plot.add_ratio(var_name+'_data_hist',var_name+'_mcsf_hist')
      var_plot.add_ratio(var_name+'_data_hist',var_name+'_mcbn_hist')
      var_plot.draw('plots/validate_{}_reweight_{}{}.pdf'.format(var_name,
          target_bin[1],args.tag))
    ibin += 1

