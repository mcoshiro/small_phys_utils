#!/usr/bin/env python3
"""@package docstring
Generates photon preselection corrections using tag-and-probe
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from tnp_utils import do_tnp_fit, integrate_tgraph
from root_plot_lib import RplPlot
from math import sqrt, hypot
import ROOT
import json
import subprocess
import ctypes

#constants
pt_bins = [15.0,17.5,20.0,25.0,35.0,50.0,200.0]
#pt_bins = [35.0,50.0] #havent run this yet 2017
eta_bins = [0.0,0.8,1.5,2.0,2.5]

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

def make_plots(data_input_file, mc_input_file):
  '''Generate mll spectra to fit in the next step
  '''
  #setup dataframe and add tag/probe preselections
  df = ROOT.RDataFrame('tnpPhoIDs/fitter_tree',data_input_file)
  df = df.Filter('tag_Ele_pt>35&&tag_sc_abseta<2.17&&(tag_sc_abseta>1.566||tag_sc_abseta<1.4442)')
  df = df.Filter('ph_et>15&&ph_sc_abseta<2.5&&(ph_sc_abseta<1.4442||ph_sc_abseta>1.566)')
  df = df.Define('pass_preselection','((ph_sc_abseta<1.4442&&ph_mva94XV2>-0.4)||(ph_sc_abseta>1.566&&ph_mva94XV2>-0.58))')

  df_mc = ROOT.RDataFrame('tnpPhoIDs/fitter_tree',mc_input_file)
  df_mc = df_mc.Filter('tag_Ele_pt>35&&tag_sc_abseta<2.17&&(tag_sc_abseta>1.566||tag_sc_abseta<1.4442)&&mcTrue')
  df_mc = df_mc.Filter('ph_et>15&&ph_sc_abseta<2.5&&(ph_sc_abseta<1.4442||ph_sc_abseta>1.566)')
  df_mc = df_mc.Define('pass_preselection','((ph_sc_abseta<1.4442&&ph_mva94XV2>-0.4)||(ph_sc_abseta>1.566&&ph_mva94XV2>-0.58))')

  probe_fail_ptrs = []
  probe_pass_ptrs = []
  mc_probe_fail_ptrs = []
  mc_probe_pass_ptrs = []

  #book histograms
  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):

      #Filter for pt-eta bin and for pass/fail
      df_bin = df.Filter('ph_et>'+str(pt_bins[ipt])+'&&ph_et<='
                         +str(pt_bins[ipt+1])+'&&ph_sc_abseta>'
                         +str(eta_bins[ieta])+'&&ph_sc_abseta<'
                         +str(eta_bins[ieta+1]))
      df_pass = df_bin.Filter('pass_preselection')
      df_fail = df_bin.Filter('!pass_preselection')
      df_mc_bin = df_mc.Filter('ph_et>'+str(pt_bins[ipt])+'&&ph_et<='
                         +str(pt_bins[ipt+1])+'&&ph_sc_abseta>'
                         +str(eta_bins[ieta])+'&&ph_sc_abseta<'
                         +str(eta_bins[ieta+1]))
      df_mc_pass = df_mc_bin.Filter('pass_preselection')
      df_mc_fail = df_mc_bin.Filter('!pass_preselection')

      #book plots to fit
      bin_string = '_bin'+str(ipt)+'_'+str(ieta)
      probe_pass_ptrs.append(df_pass.Histo1D(
          ('llfitdata_pass'+bin_string,'',80,50.0,130.0),'pair_mass'))
      probe_fail_ptrs.append(df_fail.Histo1D(
          ('llfitdata_fail'+bin_string,'',80,50.0,130.0),'pair_mass'))
      mc_probe_pass_ptrs.append(df_mc_pass.Histo1D(
          ('llfitsimu_pass'+bin_string,'',80,50.0,130.0),'pair_mass'))
      mc_probe_fail_ptrs.append(df_mc_fail.Histo1D(
          ('llfitsimu_fail'+bin_string,'',80,50.0,130.0),'pair_mass'))

  #save mll histograms
  output_file = ROOT.TFile('temp_zpeak_hists.root','RECREATE')

  ipteta = 0
  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):
      pass_hist = probe_pass_ptrs[ipteta].GetValue()
      fail_hist = probe_fail_ptrs[ipteta].GetValue()
      mc_pass_hist = mc_probe_pass_ptrs[ipteta].GetValue()
      mc_fail_hist = mc_probe_fail_ptrs[ipteta].GetValue()
      pass_hist.Write()
      fail_hist.Write()
      mc_pass_hist.Write()
      mc_fail_hist.Write()
      ipteta += 1

  output_file.Close()

def temp_debug(mc_input_file):
  df = ROOT.RDataFrame('tnpPhoIDs/fitter_tree',mc_input_file)
  df = df.Filter('tag_Ele_pt>35&&tag_sc_abseta<2.17&&(tag_sc_abseta>1.566||tag_sc_abseta<1.4442)')
  df = df.Filter('ph_et>15&&ph_sc_abseta<2.5&&(ph_sc_abseta<1.4442||ph_sc_abseta>1.566)')
  df = df.Define('pass_preselection','((ph_sc_abseta<1.4442&&ph_mva94XV2>-0.4)||(ph_sc_abseta>1.566&&ph_mva94XV2>-0.58))')
  plot_one = df.Filter('mcTrue==1').Histo1D(('plot_one','',50,-1.0,1.0),'ph_mva94XV2')
  plot_two = df.Filter('mcTrue==1&&ph_et>15&&ph_et<20').Histo1D(('plot_two','',50,-1.0,1.0),'ph_mva94XV2')
  plot_three = df.Filter('mcTrue==1&&ph_et>15&&ph_et<20').Histo1D(('plot_three','',50,50.0,130.0),'pair_mass')
  pplot_one = RplPlot()
  pplot_one.lumi_data = [(60,13)]
  pplot_one.plot_outline(plot_one.GetValue())
  pplot_one.draw('plots/debug_one.pdf')
  pplot_two = RplPlot()
  pplot_two.lumi_data = [(60,13)]
  pplot_two.plot_outline(plot_two.GetValue())
  pplot_two.draw('plots/debug_two.pdf')
  pplot_three = RplPlot()
  pplot_three.lumi_data = [(60,13)]
  pplot_three.plot_outline(plot_three.GetValue())
  pplot_three.draw('plots/debug_three.pdf')

def generate_weights(eff_json_filename, corr_json_filename, prefit_filename=''):
  '''Use histograms from the previous step, and perform fits to the Z-peak to
  get efficiencies for data. Directly extracts MC efficiencies for comparison
  '''
  ##Calculate MC efficiencies
  #df = ROOT.RDataFrame('tnpPhoIDs/fitter_tree',mc_input_file)
  #df = df.Filter('tag_Ele_pt>35&&tag_sc_abseta<2.17&&(tag_sc_abseta>1.566||tag_sc_abseta<1.4442)&&mcTrue')
  #df = df.Filter('ph_et>15&&ph_sc_abseta<2.5&&(ph_sc_abseta<1.4442||ph_sc_abseta>1.566)')
  #df = df.Define('pass_preselection','((ph_sc_abseta<1.4442&&ph_mva94XV2>-0.4)||(ph_sc_abseta>1.566&&ph_mva94XV2>-0.58))')

  #pass_count_ptrs = []
  #fail_count_ptrs = []

  mc_efficiency = []
  mc_uncertainty = []

  #for ipt in range(len(pt_bins)-1):
  #  for ieta in range(len(eta_bins)-1):
  #    df_bin = df.Filter('ph_et>'+str(pt_bins[ipt])+'&&ph_et<='
  #                       +str(pt_bins[ipt+1])+'&&ph_sc_abseta>'
  #                       +str(eta_bins[ieta])+'&&ph_sc_abseta<'
  #                       +str(eta_bins[ieta+1]))
  #    pass_count_ptrs.append(df_bin.Filter('pass_preselection').Count())
  #    fail_count_ptrs.append(df_bin.Filter('!pass_preselection').Count())
  #
  #ipteta = 0
  #for ipt in range(len(pt_bins)-1):
  #  for ieta in range(len(eta_bins)-1):
  #    pass_count = float(pass_count_ptrs[ipteta].GetValue())
  #    fail_count = float(fail_count_ptrs[ipteta].GetValue())
  #    #asymptotic limit
  #    print('pass and fail counts in MC:')
  #    print(pass_count)
  #    print(fail_count)
  #    eff = pass_count/(pass_count+fail_count)
  #    unc = eff*sqrt(1.0/pass_count+1.0/fail_count)
  #    mc_efficiency.append(eff)
  #    mc_uncertainty.append(unc)
  #    ipteta += 1

  prefit_strings = []
  prefit_idx = 0
  if prefit_filename != '':
    with open(prefit_filename,'r') as prefit_file:
      for line in prefit_file:
        prefit_strings.append(line[:-1])
  else:
    for i in range((len(pt_bins)-1)*(len(eta_bins)-1)*2*4):
      prefit_strings.append('')

  #Do Z-peak fits to extract data uncertainties
  input_file = ROOT.TFile('temp_zpeak_hists.root','READ')

  mass_hist_idx = 0
  output_string = ''

  data_efficiency = []
  data_uncertainty = []

  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):
      bin_string = '_bin'+str(ipt)+'_'+str(ieta)

      fit_types = ['nom','altbkg','altsig']
      fit_funcs = ['CB+CMS','CB+logisticexp','Gauss+CMS']
      pass_yield = []
      fail_yield = []

      for ifit in range(3):
        pass_hist = input_file.Get('llfitdata_pass'+bin_string)
        fit_result_pass = do_tnp_fit(pass_hist,fit_funcs[ifit],
            'plots/photon_corr/presel_pass_'+fit_types[ifit]+'_tnp'+bin_string+'.pdf',
            'curve',True,prefit_strings[prefit_idx])
        prefit_idx += 1
        if (fit_result_pass[0] != 0):
          raise RuntimeError('Fit failed.')

        integral_pass_sb = integrate_tgraph(fit_result_pass[1],51,129)
        integral_pass_bk = integrate_tgraph(fit_result_pass[2],51,129)
        pass_yield.append(integral_pass_sb-integral_pass_bk)

        fail_hist = input_file.Get('llfitdata_fail'+bin_string)
        fit_result_fail = do_tnp_fit(fail_hist,fit_funcs[ifit],
            'plots/photon_corr/presel_fail_'+fit_types[ifit]+'_tnp'+bin_string+'.pdf',
            'curve',True,prefit_strings[prefit_idx])
        prefit_idx += 1
        if (fit_result_fail[0] != 0):
          raise RuntimeError('Fit failed.')

        integral_fail_sb = integrate_tgraph(fit_result_fail[1],51,129)
        integral_fail_bk = integrate_tgraph(fit_result_fail[2],51,129)
        fail_yield.append(integral_fail_sb-integral_fail_bk)

      #do MC
      mc_pass_hist = input_file.Get('llfitsimu_pass'+bin_string)
      mc_fit_result_pass = do_tnp_fit(mc_pass_hist,'CB+CMS',
          'plots/photon_corr/presel_mc_pass_'+fit_types[0]+'_tnp'+bin_string+'.pdf',
          'curve',True,prefit_strings[prefit_idx])
      prefit_idx += 1
      if (mc_fit_result_pass[0] != 0):
        raise RuntimeError('Fit failed.')

      integral_pass_sb = integrate_tgraph(mc_fit_result_pass[1],51,129)
      integral_pass_bk = integrate_tgraph(mc_fit_result_pass[2],51,129)
      mc_pass_yield = integral_pass_sb-integral_pass_bk

      mc_fail_hist = input_file.Get('llfitsimu_fail'+bin_string)
      mc_fit_result_fail = do_tnp_fit(mc_fail_hist,'CB+CMS',
          'plots/photon_corr/presel_mc_fail_'+fit_types[0]+'_tnp'+bin_string+'.pdf',
          'curve',True,prefit_strings[prefit_idx])
      prefit_idx += 1
      if (mc_fit_result_fail[0] != 0):
        raise RuntimeError('Fit failed.')

      integral_fail_sb = integrate_tgraph(mc_fit_result_fail[1],51,129)
      integral_fail_bk = integrate_tgraph(mc_fit_result_fail[2],51,129)
      mc_fail_yield = integral_fail_sb-integral_fail_bk

      #calculate yields and whatnot
      nominal_efficiency = pass_yield[0]/(pass_yield[0]+fail_yield[0])
      stat_uncertainty = hypot(sqrt(pass_yield[0])/(pass_yield[0]+fail_yield[0]),
                               hypot(sqrt(pass_yield[0]),sqrt(fail_yield[0]))/(pass_yield[0]+fail_yield[0])**2*pass_yield[0])
      altbkg_efficiency = pass_yield[1]/(pass_yield[1]+fail_yield[1])
      altsig_efficiency = pass_yield[2]/(pass_yield[2]+fail_yield[2])
      syst_uncertainty = 0.0
      if (nominal_efficiency != 0):
        syst_uncertainty = max(abs(nominal_efficiency-altbkg_efficiency)/nominal_efficiency,
                              abs(nominal_efficiency-altsig_efficiency)/nominal_efficiency)
        syst_uncertainty *= nominal_efficiency
      mc_eff = mc_pass_yield/(mc_pass_yield+mc_fail_yield)
      mc_unc = hypot(sqrt(mc_pass_yield)/(mc_pass_yield+mc_fail_yield),
                     hypot(sqrt(mc_pass_yield),sqrt(mc_fail_yield))/(mc_pass_yield+mc_fail_yield)**2*mc_pass_yield)

      data_efficiency.append(min(nominal_efficiency,1.0))
      data_uncertainty.append(hypot(stat_uncertainty,syst_uncertainty))
      mc_efficiency.append(mc_eff)
      mc_uncertainty.append(mc_unc)

  input_file.Close()

  #calculate scale factors
  pass_sfs_nom = []
  fail_sfs_nom = []
  pass_sfs_up = []
  fail_sfs_up = []
  pass_sfs_down = []
  fail_sfs_down = []

  ipteta = 0
  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):
      data_eff = data_efficiency[ipteta]
      data_unc = data_uncertainty[ipteta]
      simu_eff = mc_efficiency[ipteta]
      simu_unc = mc_uncertainty[ipteta]
      
      pass_sf = 1.0
      pass_unc = 0.0
      fail_sf = 1.0
      fail_unc = 0.0
      if (simu_eff != 0):
        pass_sf = data_eff/simu_eff
      if (simu_eff != 1):
        fail_sf = (1.0-data_eff)/(1.0-simu_eff)
      pass_unc = sqrt((data_unc*simu_eff)**2+(simu_unc*data_eff)**2)
      fail_unc = sqrt((data_unc*(1.0-simu_eff))**2+(simu_unc*(1.0-data_eff))**2)
      pass_sfs_nom.append(pass_sf)
      fail_sfs_nom.append(fail_sf)
      pass_sfs_up.append(pass_sf+pass_unc)
      fail_sfs_up.append(fail_sf-fail_unc)
      pass_sfs_down.append(pass_sf-pass_unc)
      fail_sfs_down.append(fail_sf+fail_unc)
      ipteta += 1

  #package corrections in a JSON file
  simu_eff_set = schemav2.Correction(
      name='photon_preselection_mc_eff',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='eff', type='real', description='MC efficiency'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=mc_efficiency,
          flow='clamp',
          ),
      )
  data_eff_set = schemav2.Correction(
      name='photon_preselection_data_eff',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='eff', type='real', description='Data efficiency'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=data_efficiency,
          flow='clamp',
          ),
      )
  simu_unc_set = schemav2.Correction(
      name='photon_preselection_mc_unc',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='eff', type='real', description='MC uncertainty'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=mc_uncertainty,
          flow='clamp',
          ),
      )
  data_unc_set = schemav2.Correction(
      name='photon_preselection_data_unc',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='eff', type='real', description='Data uncertainty'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=data_uncertainty,
          flow='clamp',
          ),
      )

  #
  sf_nom_pass_set = schemav2.Correction(
      name='photon_preselection_sf_pass_nom',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='sf', type='real', description='Nominal weight'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=pass_sfs_nom,
          flow='clamp',
          ),
      )
  sf_nom_fail_set = schemav2.Correction(
      name='photon_preselection_sf_fail_nom',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='sf', type='real', description='Nominal weight'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=fail_sfs_nom,
          flow='clamp',
          ),
      )
  sf_up_pass_set = schemav2.Correction(
      name='photon_preselection_sf_pass_up',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='sf', type='real', description='Weight upward variation'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=pass_sfs_up,
          flow='clamp',
          ),
      )
  sf_up_fail_set = schemav2.Correction(
      name='photon_preselection_sf_fail_up',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='sf', type='real', description='Weight downward variation'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=fail_sfs_up,
          flow='clamp',
          ),
      )
  sf_down_pass_set = schemav2.Correction(
      name='photon_preselection_sf_pass_down',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='sf', type='real', description='Weight downward variation'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=pass_sfs_down,
          flow='clamp',
          ),
      )
  sf_down_fail_set = schemav2.Correction(
      name='photon_preselection_sf_fail_down',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta')],
      output=schemav2.Variable(name='sf', type='real', description='Weight downward variation'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta'],
          edges=[pt_bins,eta_bins],
          content=fail_sfs_down,
          flow='clamp',
          ),
      )

  with open(eff_json_filename,'w') as output_file:
    output_file.write(fix_correctionlib_json(
      [data_eff_set.json(exclude_unset=True),
       simu_eff_set.json(exclude_unset=True),
       data_unc_set.json(exclude_unset=True),
       simu_unc_set.json(exclude_unset=True)]))

  with open(corr_json_filename,'w') as output_file:
    output_file.write(fix_correctionlib_json(
      [sf_nom_pass_set.json(exclude_unset=True),
       sf_nom_fail_set.json(exclude_unset=True),
       sf_up_pass_set.json(exclude_unset=True),
       sf_up_fail_set.json(exclude_unset=True),
       sf_down_pass_set.json(exclude_unset=True),
       sf_down_fail_set.json(exclude_unset=True)]))

if __name__ == '__main__':

  argument_parser = ArgumentParser(prog='fit_zpeak',
  description='Script to generate plots and weights for photon correction background subtraction')
  argument_parser.add_argument('-s','--steps',default='make_plots,generate_weights,clean')
  argument_parser.add_argument('-d','--data_input_file')
  argument_parser.add_argument('-m','--mc_input_file')
  argument_parser.add_argument('-e','--eff_json_file')
  argument_parser.add_argument('-c','--corr_json_file')

  args = argument_parser.parse_args()
  steps = args.steps.split(',')

  if 'make_plots' in steps:
    make_plots(args.data_input_file, args.mc_input_file)

  if 'generate_weights' in steps:
    #temp_debug(args.mc_input_file)
    generate_weights(args.eff_json_file, args.corr_json_file)

  if 'do_precalculated_fits' in steps:
    generate_weights(args.eff_json_file, args.corr_json_file, 'temp_tnp.txt')

  if 'clean' in steps:
    subprocess.run('rm temp_zpeak_hists.root'.split())


