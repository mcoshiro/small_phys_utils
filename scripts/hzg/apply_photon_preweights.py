'''
Script that makes ratio of histograms
'''

from argparse import ArgumentParser
from array import array
from root_plot_lib import RplPlot
import ROOT
import correctionlib
import correctionlib.convert
import subprocess
import skim_and_slim
import json
import numpy as np

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

def th2_get_bin_center(hist, x_idx, y_idx):
  '''Returns (float,float) indicating center
  of bins
    hist   TH2  ROOT 2D histogram
    x_idx  int  index of x bin
    y_idx  int  index of y bin
  '''
  x_lo = hist.GetXaxis().GetBinLowEdge(x_idx)
  y_lo = hist.GetYaxis().GetBinLowEdge(y_idx)
  x_hi = hist.GetXaxis().GetBinUpEdge(x_idx)
  y_hi = hist.GetYaxis().GetBinUpEdge(y_idx)
  return ((x_lo+x_hi)/2.0, (y_lo+y_hi)/2.0)

def generate_weights_file_correctionlib(numerator_filename, denominator_filename, weights_tag):
  '''Generates ROOT file needed for reweighting
  '''
  print('Generating pt-eta weights')
  pt_bins_list = [15,16.25,17.5,18.75,20,22.5,25,27.5,30,32,34,36,38,40,42,44,46,48,50,52.5,55,60,65,70,75,100,125,150,400]
  abseta_bins_list = [d*0.25 for d in range(11)]
  abseta_bins_list = [0.0,0.167,0.333,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.167,2.333,2.5]
  pt_bins = array('f', pt_bins_list)
  abseta_bins = array('f', abseta_bins_list)

  df_num = ROOT.RDataFrame('tree',numerator_filename)
  df_den = ROOT.RDataFrame('tree',denominator_filename)
  #assumes background subtraction weights already applied
  hist_num_ptr = df_num.Histo2D(
      ('num_hist','',len(pt_bins)-1,pt_bins,len(abseta_bins)-1,abseta_bins),'ph_et','ph_sc_abseta','w_bkg')
  hist_den_ptr = df_den.Histo2D(
      ('den_hist','',len(pt_bins)-1,pt_bins,len(abseta_bins)-1,abseta_bins),'ph_et','ph_sc_abseta')
  rho_num_ptr = df_num.Histo1D(
      ('rho_num_hist','',45,0.0,45.0),'event_rho','w_bkg')
  rho_den_ptr = df_den.Histo1D(
      ('rho_den_hist','',45,0.0,45.0),'event_rho')
  ess_num_ptr = df_num.Histo1D(
      ('ess_num_hist','',45,0.0,45.0),'ph_ESsigma','w_bkg')
  ess_den_ptr = df_den.Histo1D(
      ('ess_den_hist','',45,0.0,45.0),'ph_ESsigma')

  hist_num = hist_num_ptr.GetValue()
  hist_den = hist_den_ptr.GetValue()
  hist_num.Scale(1.0/hist_num.Integral())
  hist_den.Scale(1.0/hist_den.Integral())
  hist_ratio = hist_num.Clone('pteta_ratio')
  hist_ratio.SetTitle('pteta_weights;pt;eta')
  hist_ratio.Divide(hist_den)
  ROOT.gStyle.SetOptStat(0)
  can = ROOT.TCanvas('c','c',600,400)
  can.SetLogx()
  hist_ratio.Draw('colz')
  can.SaveAs('plots/pteta_ratio_hist.pdf')

  rho_num = rho_num_ptr.GetValue()
  rho_den = rho_den_ptr.GetValue()
  rho_num.Scale(1.0/rho_num.Integral())
  rho_den.Scale(1.0/rho_den.Integral())
  rho_ratio = rho_num.Clone('rho_ratio')
  rho_ratio.SetTitle('rho_weights;rho;')
  rho_ratio.Divide(rho_den)
  can.SetLogx(False)
  rho_ratio.Draw()
  can.SaveAs('plots/rho_ratio_hist.pdf')

  ess_num = ess_num_ptr.GetValue()
  ess_den = ess_den_ptr.GetValue()
  ess_num.Scale(1.0/ess_num.Integral())
  ess_den.Scale(1.0/ess_den.Integral())
  ess_ratio = ess_num.Clone('ess_ratio')
  ess_ratio.SetTitle('ess_weights;ess;')
  ess_ratio.Divide(ess_den)
  ess_ratio.Draw()
  can.SaveAs('plots/ess_ratio_hist.pdf')

  weights_filename = 'json/pteta_weights_'+weights_tag+'.json'
  rho_weights_filename = 'json/rho_weights_'+weights_tag+'.json'
  ess_weights_filename = 'json/ess_weights_'+weights_tag+'.json'

  ##METHOD 2: polynomial fit
  #nbins_x = hist_ratio.GetNbinsX()
  #nbins_y = hist_ratio.GetNbinsY()
  #x_values = []
  #y_values = []
  #z_values = []
  #w_values = []
  #for ix in range(1,nbins_x+1):
  #  for iy in range(1,nbins_y+1):
  #    x_center, y_center = th2_get_bin_center(hist_ratio,ix,iy)
  #    x_values.append(x_center)
  #    y_values.append(y_center)
  #    z_values.append(hist_ratio.GetBinContent(ix,iy))
  #    w_values.append(1.0/hist_ratio.GetBinErrorLow(ix,iy))

  #corr_poly, corr_fit = correctionlib.convert.ndpolyfit(
  #    points=[np.array(x_values),np.array(y_values)],
  #    values=np.array(z_values),
  #    weights=np.array(w_values),
  #    varnames=['pt','abseta'],
  #    degree=(8,3))
  #corr_poly.name = 'pteta_ratio'

  #with open(weights_filename,'w') as output_file:
  #  output_file.write(fix_correctionlib_json(corr_poly.json(exclude_unset=True)))

  #METHOD 1: directly use bin content
  hist_file = ROOT.TFile('temp.root','RECREATE')
  hist_ratio.Write()
  hist_file.Close()
  corr = correctionlib.convert.from_uproot_THx('temp.root:pteta_ratio')
  corr.description = 'pt-eta reweighting factors for photon DNN training'
  corr.data.flow = 'clamp'
  with open(weights_filename,'w') as output_file:
    output_file.write(fix_correctionlib_json(corr.json(exclude_unset=True)))

  hist_file = ROOT.TFile('temp.root','RECREATE')
  rho_ratio.Write()
  hist_file.Close()
  corr = correctionlib.convert.from_uproot_THx('temp.root:rho_ratio')
  corr.description = 'rho reweighting factors for photon DNN training'
  corr.data.flow = 'clamp'
  with open(rho_weights_filename,'w') as output_file:
    output_file.write(fix_correctionlib_json(corr.json(exclude_unset=True)))

  hist_file = ROOT.TFile('temp.root','RECREATE')
  ess_ratio.Write()
  hist_file.Close()
  corr = correctionlib.convert.from_uproot_THx('temp.root:ess_ratio')
  corr.description = 'ess reweighting factors for photon DNN training'
  corr.data.flow = 'clamp'
  with open(ess_weights_filename,'w') as output_file:
    output_file.write(fix_correctionlib_json(corr.json(exclude_unset=True)))

  subprocess.run('rm temp.root'.split())

def apply_weights_correctionlib(input_filename, output_filename, json_tag):
  '''Applies weights in weight_filename to tree
  '''
  print('JITting C++')
  ROOT.gInterpreter.AddIncludePath('inc/')
  ROOT.gInterpreter.ProcessLine('#include "correction_wrapper.hpp"')
  ROOT.gSystem.Load('libSmallPhysUtils.so')
  ROOT.gInterpreter.Declare("""
  
  using std::vector;
  
  const CorrectionWrapper corr_pteta("""+'"json/pteta_weights_'+json_tag+'.json"'+""","pteta_ratio");
  
  float get_w_pteta(float ph_et, float ph_sc_abseta) {
    vector<double> eval_args;
    eval_args.push_back(ph_et);
    eval_args.push_back(ph_sc_abseta);
    return corr_pteta.evaluate(eval_args);
  }

  const CorrectionWrapper corr_rho("""+'"json/rho_weights_'+json_tag+'.json"'+""","rho_ratio");
  
  float get_w_rho(float rho) {
    vector<double> eval_args;
    eval_args.push_back(rho);
    return corr_rho.evaluate(eval_args);
  }
  
  const CorrectionWrapper corr_ess("""+'"json/ess_weights_'+json_tag+'.json"'+""","ess_ratio");
  
  float get_w_ess(float ess) {
    vector<double> eval_args;
    eval_args.push_back(ess);
    return corr_ess.evaluate(eval_args);
  }
  
  """)
  print('Applying pt-eta weights')
  defines = [('w_pre','get_w_pteta(ph_et,ph_sc_abseta)*get_w_rho(event_rho)*get_w_ess(ph_ESsigma)')]
  skim_and_slim.write_ntuples(
      [input_filename],
      [],
      output_filename,
      defines,
      'tree')

def calculate_weights(numerator_filename, denominator_filename, temp_filename):
  '''Reweights distributions in MC to match data
  Assumes background subtraction weights have already been calculated for data
  '''
  pt_bins_list = [15,16.25,17.5,18.75,20,22.5,25,27.5,30,32,34,36,38,40,42,44,46,48,50,52.5,55,60,65,70,75,100,125,150,400]
  abseta_bins_list = [d*0.25 for d in range(11)]
  abseta_bins_list = [0.0,0.167,0.333,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.167,2.333,2.5]
  pt_bins = array('f', pt_bins_list)
  abseta_bins = array('f', abseta_bins_list)

  df_num = ROOT.RDataFrame('tree',numerator_filename)
  df_den = ROOT.RDataFrame('tree',denominator_filename)
  ess_num_ptr = df_num.Histo1D(
      ('ess_num_hist','',30,0.0,10.0),'ph_ESsigma','w_bkg')
  ess_den_ptr = df_den.Histo1D(
      ('ess_den_hist','',30,0.0,10.0),'ph_ESsigma')
  rho_num_ptr = df_num.Histo1D(
      ('rho_num_hist','',45,0.0,45.0),'event_rho','w_bkg')
  rho_den_ptr = df_den.Histo1D(
      ('rho_den_hist','',45,0.0,45.0),'event_rho')

  rho_num = rho_num_ptr.GetValue()
  rho_den = rho_den_ptr.GetValue()
  rho_num.Scale(1.0/rho_num.Integral())
  rho_den.Scale(1.0/rho_den.Integral())
  rho_ratio = rho_num.Clone('rho_ratio')
  rho_ratio.SetTitle('rho_weights;rho;')
  rho_ratio.Divide(rho_den)
  ROOT.gStyle.SetOptStat(0)
  can = ROOT.TCanvas('c','c',600,400)
  rho_ratio.Draw()
  can.SaveAs('plots/rho_ratio_hist.pdf')

  ess_num = ess_num_ptr.GetValue()
  ess_den = ess_den_ptr.GetValue()
  ess_num.Scale(1.0/ess_num.Integral())
  ess_den.Scale(1.0/ess_den.Integral())
  ess_ratio = ess_num.Clone('ess_ratio')
  ess_ratio.SetTitle('ess_weights;ess;')
  ess_ratio.Divide(ess_den)
  ess_ratio.Draw()
  can.SaveAs('plots/ess_ratio_hist.pdf')

  ROOT.gDirectory.Add(rho_ratio)
  ROOT.gDirectory.Add(ess_ratio)

  ROOT.gInterpreter.Declare("""
  TH1D* rho_ratio = static_cast<TH1D*>(gDirectory->Get("rho_ratio"));
  TH1D* ess_ratio = static_cast<TH1D*>(gDirectory->Get("ess_ratio"));
  """)
  df_den = df_den.Define('w_rho','rho_ratio->GetBinContent(rho_ratio->FindBin(event_rho))')
  df_den = df_den.Define('w_ess','ess_ratio->GetBinContent(ess_ratio->FindBin(ph_ESsigma))')
  df_den = df_den.Define('w_rho_ess','w_rho*w_ess')

  hist_num_ptr = df_num.Histo2D(
      ('num_hist','',len(pt_bins)-1,pt_bins,len(abseta_bins)-1,abseta_bins),'ph_et','ph_sc_abseta','w_bkg')
  hist_den_ptr = df_den.Histo2D(
      ('den_hist','',len(pt_bins)-1,pt_bins,len(abseta_bins)-1,abseta_bins),'ph_et','ph_sc_abseta','w_rho_ess')

  hist_num = hist_num_ptr.GetValue()
  hist_den = hist_den_ptr.GetValue()
  hist_num.Scale(1.0/hist_num.Integral())
  hist_den.Scale(1.0/hist_den.Integral())
  hist_ratio = hist_num.Clone('pteta_ratio')
  hist_ratio.SetTitle('pteta_weights;pt;eta')
  hist_ratio.Divide(hist_den)
  can.SetLogx()
  hist_ratio.Draw('colz')
  can.SaveAs('plots/pteta_ratio_hist.pdf')

  ROOT.gDirectory.Add(hist_ratio)

  ROOT.gInterpreter.Declare("""
  TH2D* pteta_ratio = static_cast<TH2D*>(gDirectory->Get("pteta_ratio"));

  float get_prew(float ph_pt, float ph_abseta, float rho, float ph_ESsigma) {
    float w_rho = rho_ratio->GetBinContent(rho_ratio->FindBin(rho));
    float w_ess = ess_ratio->GetBinContent(ess_ratio->FindBin(ph_ESsigma));
    float w_pteta = pteta_ratio->GetBinContent(pteta_ratio->FindBin(ph_pt, ph_abseta));
    return w_rho*w_ess*w_pteta;
  }
  """)

  defines = [('w_pre','get_prew(ph_et,ph_sc_abseta,event_rho,ph_ESsigma)')]
  skim_and_slim.write_ntuples(
      [denominator_filename],
      [],
      temp_filename,
      defines,
      'tree')

def validate_weights(data_filename, weighted_mc_filename):
  '''Make plots to show effect of pt-eta reweighting
  '''
  print('Validating pt-eta weights')
  df_simu = ROOT.RDataFrame('tree',weighted_mc_filename)
  df_data = ROOT.RDataFrame('tree',data_filename)
  hist_pt_simu_og_ptr = df_simu.Histo1D(('hist_pt_simu_og','Simulation unweighted;p_{T} [GeV];Events/bin',30,15,75),'ph_et')
  hist_pt_simu_wg_ptr = df_simu.Histo1D(('hist_pt_simu_wg','Simulation reweighted;p_{T} [GeV];Events/bin',30,15,75),'ph_et','w_pre')
  hist_pt_data_og_ptr = df_data.Histo1D(('hist_pt_data_og','Data;p_{T} [GeV];Events/bin',30,15,75),'ph_et','w_bkg')
  hist_eta_simu_og_ptr = df_simu.Histo1D(('hist_eta_simu_og','Simulation unweighted;#eta;Events/bin',30,-2.5,2.5),'ph_sc_eta')
  hist_eta_simu_wg_ptr = df_simu.Histo1D(('hist_eta_simu_wg','Simulation reweighted;#eta;Events/bin',30,-2.5,2.5),'ph_sc_eta','w_pre')
  hist_eta_data_og_ptr = df_data.Histo1D(('hist_eta_data_og','Data;#eta;Events/bin',30,-2.5,2.5),'ph_sc_eta','w_bkg')
  hist_rho_simu_og_ptr = df_simu.Histo1D(('hist_rho_simu_og','Simulation unweighted;#rho [GeV];Events/bin',30,0.0,45.0),'event_rho')
  hist_rho_simu_wg_ptr = df_simu.Histo1D(('hist_rho_simu_wg','Simulation reweighted;#rho [GeV];Events/bin',30,0.0,45.0),'event_rho','w_pre')
  hist_rho_data_og_ptr = df_data.Histo1D(('hist_rho_data_og','Data;#rho [GeV];Events/bin',30,0.0,45.0),'event_rho','w_bkg')
  hist_ess_simu_og_ptr = df_simu.Histo1D(('hist_ess_simu_og','Simulation unweighted;ES #sigma_{eff};Events/bin',30,0.0,10.0),'ph_ESsigma')
  hist_ess_simu_wg_ptr = df_simu.Histo1D(('hist_ess_simu_wg','Simulation reweighted;ES #sigma_{eff};Events/bin',30,0.0,10.0),'ph_ESsigma','w_pre')
  hist_ess_data_og_ptr = df_data.Histo1D(('hist_ess_data_og','Data;ES #sigma_{eff};Events/bin',30,0.0,10.0),'ph_ESsigma','w_bkg')
  hist_pt_simu_og = hist_pt_simu_og_ptr.GetValue()
  hist_pt_simu_wg = hist_pt_simu_wg_ptr.GetValue()
  hist_pt_data_og = hist_pt_data_og_ptr.GetValue()
  hist_eta_simu_og = hist_eta_simu_og_ptr.GetValue()
  hist_eta_simu_wg = hist_eta_simu_wg_ptr.GetValue()
  hist_eta_data_og = hist_eta_data_og_ptr.GetValue()
  hist_rho_simu_og = hist_rho_simu_og_ptr.GetValue()
  hist_rho_simu_wg = hist_rho_simu_wg_ptr.GetValue()
  hist_rho_data_og = hist_rho_data_og_ptr.GetValue()
  hist_ess_simu_og = hist_ess_simu_og_ptr.GetValue()
  hist_ess_simu_wg = hist_ess_simu_wg_ptr.GetValue()
  hist_ess_data_og = hist_ess_data_og_ptr.GetValue()
  hist_pt_simu_og.Scale(hist_pt_data_og.Integral()/hist_pt_simu_og.Integral())
  hist_pt_simu_wg.Scale(hist_pt_data_og.Integral()/hist_pt_simu_wg.Integral())
  hist_eta_simu_og.Scale(hist_eta_data_og.Integral()/hist_eta_simu_og.Integral())
  hist_eta_simu_wg.Scale(hist_eta_data_og.Integral()/hist_eta_simu_wg.Integral())
  hist_rho_simu_og.Scale(hist_rho_data_og.Integral()/hist_rho_simu_og.Integral())
  hist_rho_simu_wg.Scale(hist_rho_data_og.Integral()/hist_rho_simu_wg.Integral())
  hist_ess_simu_og.Scale(hist_ess_data_og.Integral()/hist_ess_simu_og.Integral())
  hist_ess_simu_wg.Scale(hist_ess_data_og.Integral()/hist_ess_simu_wg.Integral())
  pt_plot = RplPlot()
  pt_plot.lumi_data = [(60,13)]
  pt_plot.y_title = 'Events/bin'
  pt_plot.plot_outline(hist_pt_simu_og)
  pt_plot.plot_outline(hist_pt_simu_wg)
  pt_plot.plot_outline(hist_pt_data_og)
  pt_plot.add_ratio('hist_pt_data_og','hist_pt_simu_og')
  pt_plot.add_ratio('hist_pt_data_og','hist_pt_simu_wg')
  pt_plot.draw('plots/validate_pt_reweight.pdf')
  eta_plot = RplPlot()
  eta_plot.lumi_data = [(60,13)]
  eta_plot.y_title = 'Events/bin'
  eta_plot.plot_outline(hist_eta_simu_og)
  eta_plot.plot_outline(hist_eta_simu_wg)
  eta_plot.plot_outline(hist_eta_data_og)
  eta_plot.add_ratio('hist_eta_data_og','hist_eta_simu_og')
  eta_plot.add_ratio('hist_eta_data_og','hist_eta_simu_wg')
  eta_plot.draw('plots/validate_eta_reweight.pdf')
  rho_plot = RplPlot()
  rho_plot.lumi_data = [(60,13)]
  rho_plot.y_title = 'Events/bin'
  rho_plot.plot_outline(hist_rho_simu_og)
  rho_plot.plot_outline(hist_rho_simu_wg)
  rho_plot.plot_outline(hist_rho_data_og)
  rho_plot.add_ratio('hist_rho_data_og','hist_rho_simu_og')
  rho_plot.add_ratio('hist_rho_data_og','hist_rho_simu_wg')
  rho_plot.draw('plots/validate_rho_reweight.pdf')
  ess_plot = RplPlot()
  ess_plot.lumi_data = [(60,13)]
  ess_plot.y_title = 'Events/bin'
  ess_plot.plot_outline(hist_ess_simu_og)
  ess_plot.plot_outline(hist_ess_simu_wg)
  ess_plot.plot_outline(hist_ess_data_og)
  ess_plot.add_ratio('hist_ess_data_og','hist_ess_simu_og')
  ess_plot.add_ratio('hist_ess_data_og','hist_ess_simu_wg')
  ess_plot.draw('plots/validate_ess_reweight.pdf')

if __name__=='__main__':

  #parse arguments
  argument_parser = ArgumentParser(prog='add_pteta_weights',
      description='Applies pt-eta reweighting to an MC file')
  argument_parser.add_argument('-m','--mc_filename',default='/net/cms26/cms26r0/oshiro/tnp_tuples/shuffled_photonidskim_mc_2018.root')
  argument_parser.add_argument('-d','--data_filename',default='/net/cms26/cms26r0/oshiro/tnp_tuples/shuffled_photonidskim_data_2018.root')
  argument_parser.add_argument('-j','--json_tag',default='2018')
  argument_parser.add_argument('-s','--step',default='generate,validate')
  args = argument_parser.parse_args()
  steps = args.step.split(',')

  ROOT.EnableImplicitMT()

  if (args.mc_filename[-5:] != '.root'):
    raise ValueError('MC Input file extension should be root')
  temp_filename = args.mc_filename[:-5]+'_temp.root'

  if 'generate_correctionlib' in steps: #old version
    generate_weights_file_correctionlib(args.data_filename, args.mc_filename, args.json_tag)
    apply_weights_correctionlib(args.mc_filename, temp_filename, args.json_tag)
    subprocess.run(('rm '+args.mc_filename).split())
    subprocess.run(('mv '+temp_filename+' '+args.mc_filename).split())

  if 'generate' in steps: #new version
    calculate_weights(args.data_filename, args.mc_filename, temp_filename)
    subprocess.run(('rm '+args.mc_filename).split())
    subprocess.run(('mv '+temp_filename+' '+args.mc_filename).split())

  if 'validate' in steps:
    validate_weights(args.data_filename, args.mc_filename)

