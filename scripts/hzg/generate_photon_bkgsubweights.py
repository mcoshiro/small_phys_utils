#!/usr/bin/env python3
"""@package docstring
Fit Z peaks
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from tnp_utils import fit_cmsshape_cb, integrate_tgraph
from root_plot_lib import RplPlot
import ROOT
import math
import subprocess
import ctypes

#constants
variables = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','event_rho','ph_ESsigma','ph_esEnergyOverRawE','ph_mva94XV2','ph_et','ph_sc_abseta']
xtitle = ['Photon R_{9}','Photon s_{4}','Photon supercluster #eta width','Photon supercluster #phi width','Photon #sigma_{i#eta i#eta}','Photon #sigma_{i#eta i#phi}','Photon I_{abs}(photon) [GeV]','Photon I_{abs}(charged, PV) [GeV]','Photon I_{abs}(charged, worst vertex) [GeV]','Photon supercluster E_{raw} [GeV]','Photon supercluster #eta','Event #rho [GeV]','Photon preshower #sigma_{eff}','Photon E_{ES}/E_{raw}','Photon Fall17v2 ID score','Photon p_{T} [GeV]','Photon supercluster |#eta|']
ranges = [(0.4,1.0),(0.4,0.975),(0.004,0.0225),(0.005,0.1),(0.006,0.03),(-0.002,0.002),(0.0,6.0),(0.0,6.0),(0.0,25.0),(0.0,250.0),(-2.5,2.5),(0.0,45.0),(0.0,10.0),(0.0,0.16),(-0.6,1.0),(0.0,100.0),(0.0,2.5)]
is_log = [False, False, False, False, False, True, False, True, True, False, False, False, True, True, False, False, False]
pt_bins = [15.0,17.5,20.0,22.5,25.0]
eta_bins = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5]
mll_bins = [50.0,81.0,101.0,130.0]

def fix_correctionlib_json(json_text):
  '''Fixes the format of correctionlib json created using corr.json, since 
  it is not properly formatted by default
  '''
  corr = json.loads(json_text)
  json_dict = {
    'schema_version' : 2,
    'description' : '',
    'corrections' : [corr]
  }
  return json.dumps(json_dict,indent=2)

def make_plots(input_file):
  '''Generate eta-pt binned plots of photon ID variables to show relative 
  indepedence of mll, and mll spectra to fit in the next step
  '''
  df = ROOT.RDataFrame('tree',input_file)

  mass_hist_ptrs = []
  lo_hist_ptrs = []
  md_hist_ptrs = []
  hi_hist_ptrs = []

  #book histograms
  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):
      df_bin = df.Filter('ph_et>'+str(pt_bins[ipt])+'&&ph_et<='
                         +str(pt_bins[ipt+1])+'&&ph_sc_abseta>'
                         +str(eta_bins[ieta])+'&&ph_sc_abseta<'
                         +str(eta_bins[ieta+1]))
      bin_string = '_bin'+str(ipt)+'_'+str(ieta)
      mass_hist_ptrs.append(df_bin.Histo1D(('pair_mass_hist'+bin_string,'',80,50.0,130.0),'pair_mass'))
      df_lo = df_bin.Filter('pair_mass<81')
      df_md = df_bin.Filter('pair_mass>81&&pair_mass<101')
      df_hi = df_bin.Filter('pair_mass>101')
      for ivar in range(len(variables)):
        lo_hist_ptrs.append(df_lo.Histo1D((
            'low_'+variables[ivar]+'_hist','m_{ll}<81 GeV; '+xtitle[ivar]+';',30,
            ranges[ivar][0],ranges[ivar][1]),variables[ivar]))
        md_hist_ptrs.append(df_md.Histo1D((
            'mid_'+variables[ivar]+'_hist','81<m_{ll}<101 GeV; '+xtitle[ivar]+';',30,
            ranges[ivar][0],ranges[ivar][1]),variables[ivar]))
        hi_hist_ptrs.append(df_hi.Histo1D((
            'hi_'+variables[ivar]+'_hist','m_{ll}>101 GeV; '+xtitle[ivar]+';',30,
            ranges[ivar][0],ranges[ivar][1]),variables[ivar]))

  #draw output and save mll histograms
  mass_hist_idx = 0
  var_hist_idx = 0
  code_string = ''
  output_file = ROOT.TFile('temp_zpeak_hists.root','RECREATE')

  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):
      bin_string = '_bin'+str(ipt)+'_'+str(ieta)

      mass_hist = mass_hist_ptrs[mass_hist_idx].GetValue()
      mass_hist.Write()

      for ivar in range(len(variables)):
        lo_hist = lo_hist_ptrs[var_hist_idx].GetValue()
        md_hist = md_hist_ptrs[var_hist_idx].GetValue()
        hi_hist = hi_hist_ptrs[var_hist_idx].GetValue()
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
        plot.draw('plots/photon_corr/'+variables[ivar]+bin_string+'_mllcomp.pdf')
        var_hist_idx += 1
      mass_hist_idx += 1

  output_file.Close()

def generate_weights(json_filename):
  '''Use histograms from the previous step, and perform fits to the Z-peak to
  get background subtraction weights
  '''
  input_file = ROOT.TFile('temp_zpeak_hists.root','READ')

  mass_hist_idx = 0
  output_string = ''
  weight_content = []

  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(eta_bins)-1):
      bin_string = '_bin'+str(ipt)+'_'+str(ieta)

      mass_hist = input_file.Get('pair_mass_hist'+bin_string)
      fit_result = fit_cmsshape_cb(mass_hist,'curve',True,'plots/photon_corr/')
      if (fit_result[0] != 0):
        raise RuntimeError('Fit failed.')

      integral_losb = integrate_tgraph(fit_result[1],51,81)
      integral_lobk = integrate_tgraph(fit_result[2],51,81)
      integral_mdsb = integrate_tgraph(fit_result[1],81,101)
      integral_mdbk = integrate_tgraph(fit_result[2],81,101)
      integral_hisb = integrate_tgraph(fit_result[1],101,129)
      integral_hibk = integrate_tgraph(fit_result[2],101,129)

      output_string += 'For pt-eta bin '+str(ipt)+'-'+str(ieta)+': '
      output_string += 'm<81: S+B='+str(integral_losb)
      output_string += ', B='+str(integral_lobk)+'\n'
      output_string += '81<m<101: S+B='+str(integral_mdsb)
      output_string += ', B='+str(integral_mdbk)+'\n'
      output_string += 'm>101: S+B='+str(integral_hisb)
      output_string += ', B='+str(integral_hibk)+'\n'

      weight_content.append(-1.0*integral_mdbk/integral_losb/2.0)
      weight_content.append(1.0)
      weight_content.append(-1.0*integral_mdbk/integral_hisb/2.0)

      mass_hist_idx += 1

  weight_set = schemav2.Correction(
      name='bkg_subtraction_weight',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='"Photon" pt'),
              schemav2.Variable(name='abseta', type='real', description='"Photon" absolute eta'),
              schemav2.Variable(name='mll', type='real', description='Di"photon" mass')],
      output=schemav2.Variable(name='weight', type='real', description='Event weight'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta','mll'],
          edges=[pt_bins,eta_bins,mll_bins],
          content=weight_content,
          flow='clamp',
          ),
      )

  with open(json_filename,'w') as output_file:
    output_file.write(weight_set.json(exclude_unset=True))
    output_file.write(fix_correctionlib_json(weight_set.json(exclude_unset=True)))
  print(output_string)
  input_file.Close()

if __name__ == '__main__':

  argument_parser = ArgumentParser(prog='fit_zpeak',
  description='Script to generate plots and weights for photon correction background subtraction')
  argument_parser.add_argument('-s','--steps',default='make_plots,generate_weights,clean')
  argument_parser.add_argument('-i','--input_file',default='/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_data_2018_new.root')
  argument_parser.add_argument('-j','--json_file',default='json/bkg_weights.json')

  args = argument_parser.parse_args()
  steps = args.steps.split(',')

  if 'make_plots' in steps:
    make_plots(args.input_file)

  if 'generate_weights' in steps:
    generate_weights(args.json_file)

  if 'clean' in steps:
    subprocess.run('rm temp_zpeak_hists.root'.split())


