#!/usr/bin/env python3
"""@package docstring
Generates photon preselection corrections using tag-and-probe
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from math import sqrt
from root_plot_lib import RplPlot
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from tnp_utils import do_tnp_fit, integrate_tgraph
import array
import ctypes
import json
import numpy
import ROOT
import subprocess

#constants
PT_BINS = [15.0,17.5,20.0,25.0,30.0,35.0,45.0,500.0]
ETA_BINS = [0.0,0.8,1.5,2.5]
#MVA_BINS = [-1.0,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0]
MVA_BINS = [-1.0,-0.2,0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.925,0.95,0.975,1.0]
#MVA_BINS_BARREL = [-1.0,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0]
#MVA_BINS_ENDCAP = [-1.0,-0.58,-0.5,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0]

PT_BINS_COARSE = [15.0,20.0,25.0,35.0,100.0]
ETA_BINS_COARSE = [0.0,1.5,2.5]
MVA_BINS_COARSE = [-1.0,-0.5,0.0,0.5,1.0]

#ROOT JIT'ed C++ definitions
ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using ROOT::Math::PtEtaPhiMVector;
using std::vector;

int get_probe_photon(RVec<float> photon_pt, RVec<float> photon_drmin) {
  int ph_idx = -1;
  float max_photon_pt = 0;
  for (unsigned iph = 0; iph < photon_pt.size(); iph++ ) {
    if (photon_drmin[iph]>0.03 && photon_pt[iph]>15) {
      if (photon_pt[iph] > max_photon_pt) {
        ph_idx = static_cast<int>(iph);
        max_photon_pt = photon_pt[iph];
      }
    }
  }
  return ph_idx;
}

float get_mllph(RVec<float> photon_pt, RVec<float> photon_eta, 
                RVec<float> photon_phi, int photon_idx,
                RVec<float> mu_pt, RVec<float> mu_eta,
                RVec<float> mu_phi, RVec<bool> mu_sig) {
  vector<PtEtaPhiMVector> p;
  p.push_back(PtEtaPhiMVector(photon_pt[photon_idx],photon_eta[photon_idx],
                              photon_phi[photon_idx],0.0));
  for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
    if (mu_sig[imu]) {
      p.push_back(PtEtaPhiMVector(mu_pt[imu],mu_eta[imu],
                                  mu_phi[imu],0.106));
    }
  }
  return (p[0]+p[1]+p[2]).M();
}

const vector<float> MVA_BINS_BARREL = {-1.0,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0};
const vector<float> MVA_BINS_ENDCAP = {-1.0,-0.58,-0.5,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0};

const float* MVA_BINS_BARREL_PTR = &MVA_BINS_BARREL[0];
const float* MVA_BINS_ENDCAP_PTR = &MVA_BINS_ENDCAP[0];

TH1D* cast_tobject_to_th1d(TObject* object)
{ return static_cast<TH1D*>(object); }

""")

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


def book_hist1d_zpeak_bkgsub(dataframe, hist_desc, branch_name, weight_name = ''):
  '''Return pointers to histograms used to extract yield after nonresonant background subtraction'''
  peak_ptr = None
  subt_ptr = None
  if weight_name == '':
    peak_ptr = dataframe.Filter('mllph>80&&mllph<98').Histo1D(hist_desc,branch_name)
    subt_ptr = dataframe.Filter('(mllph>78&&mllph<80)||(mllph>98&&mllph<100)').Histo1D(hist_desc,branch_name)
  else:
    peak_ptr = dataframe.Filter('mllph>80&&mllph<98').Histo1D(hist_desc,branch_name,weight_name)
    subt_ptr = dataframe.Filter('(mllph>78&&mllph<80)||(mllph>98&&mllph<100)').Histo1D(hist_desc,branch_name,weight_name)
  return (peak_ptr, subt_ptr)


def get_bkgsub_hist(bkgsub_hist_ptrs):
  '''Returns histogram with background subtracted'''
  tot_hist = bkgsub_hist_ptrs[0].GetValue()
  sub_hist = bkgsub_hist_ptrs[1].GetValue()
  tot_hist.Add(sub_hist, -4.0/18.0)
  return tot_hist

def get_hist_content(hist):
  '''Return bin content of TH1D histogram as python list, excluding under/overflow'''
  content = []
  for ibin in range(1,hist.GetNbinsX()+1):
    content.append(hist.GetBinContent(ibin))
  return content

def get_points_from_hist(hist, x_condition):
  '''Return lists of x and y points from ROOT TH1D. 
  Condition allows one to provide a lambda, which can be used to blind part of 
  the histogram.'''
  x_list = []
  y_list = []
  for ibin in range(1,hist.GetNbinsX()+1):
    x_lo = hist.GetBinLowEdge(ibin)
    x_wd = hist.GetBinWidth(ibin)
    x = x_lo+x_wd/2.0
    y = hist.GetBinContent(ibin)
    if (x_condition(x)):
      x_list.append(x)
      y_list.append(y)
  return x_list, y_list

def get_uncs_from_hist(hist, x_condition):
  '''Return lists of y variances from ROOT TH1D
  Condition allows one to provide a lambda, which can be used to blind part of 
  the histogram.'''
  u_list = []
  for ibin in range(1,hist.GetNbinsX()+1):
    x_lo = hist.GetBinLowEdge(ibin)
    x_wd = hist.GetBinWidth(ibin)
    x = x_lo+x_wd/2.0
    if (x_condition(x)):
      u_list.append(hist.GetBinError(ibin)**2)
  return u_list

def list_to_1darray(lst):
  '''Convert list into 1D array (list of lists)'''
  return [[item] for item in lst]

def fill_th1d(h, content):
  """Fill TH1D h with content from python list"""
  for ibin in range(len(content)):
    if (content[ibin] >= 0):
      h.SetBinContent(ibin+1, content[ibin])

def double_bins(h, nbins_original, lo, hi):
  """Returns TH1D with twice as many bins and content filled by 
  interpolation"""
  h2 = ROOT.TH1D('','',nbins_original*2, lo, hi)
  for ibin in range(1,nbins_original+1):
    content1 = h.GetBinContent(ibin)
    content2 = h.GetBinContent(ibin+1)
    if (ibin == nbins_original):
      content2 = content1
    h2.SetBinContent(ibin*2-1, content1)
    h2.SetBinContent(ibin*2, (content1+content2)/2.0)
  return h2

def vec_sum(lista, listb):
  """Vector sum of two lists"""
  return [lista[i]+listb[i] for i in range(len(lista))]

def scalar_mul(scalar, lst):
  """Scalar multiplication of scalar and list"""
  return [scalar*c for c in lst]

def gpr_fit_zpeak_yield(hist, output_basename):
  '''Fit histogram with Gaussian Process Regression, save hists for sanity checks, and
  return yields as a tuple (nominal, up, down)'''
  #Build templates from GPR
  mllg = ROOT.RooRealVar('mllg','m_{ll#gamma} [GeV]',50.0,110.0)
  data = ROOT.RooDataHist('','',ROOT.RooArgList(mllg),hist)

  #do some pre-guessing to get better fits
  content_centre = hist.GetBinContent(hist.FindBin(90.5))
  content_upref = hist.GetBinContent(hist.FindBin(94.5))
  content_reldiff = 1.0
  if content_centre != 0.0:
    content_reldiff = (content_centre-content_upref)/content_centre
  blind_radius = 8
  sig_fraction_guess = 0.8
  if (content_reldiff < 0.5):
    blind_radius = 5
    sig_fraction_guess = 0.2
  elif (content_reldiff < 0.1):
    blind_radius = 3
    sig_fraction_guess = 0.05

  x_train, y_train = get_points_from_hist(hist, lambda mll : (mll < 90-blind_radius or mll > 90+blind_radius))
  #u_train = numpy.array([y for y in y_train]) #variance, not stddev
  u_train = numpy.array(get_uncs_from_hist(hist, lambda mll : (mll < 90-blind_radius or mll > 90+blind_radius))) #variance, not stddev
  x_train_gpr = list_to_1darray(x_train)
  kernel = 1 * RBF(length_scale=0.1, length_scale_bounds=(1e-6, 1e2))
  gaussian_process = GaussianProcessRegressor(kernel=kernel, alpha=u_train, n_restarts_optimizer=30)
  gaussian_process.fit(x_train_gpr, y_train)
  #fit_test, unc_test = gaussian_process.predict(x_test_gpr,return_std=True)

  x_template = [[50.125+0.25*i] for i in range(240)]
  template_th1d = ROOT.TH1D('template_th1d','',240,50.0,110.0)
  template_th1d_up = ROOT.TH1D('template_th1d','',240,50.0,110.0)
  template_th1d_dn = ROOT.TH1D('template_th1d','',240,50.0,110.0)
  y_template, y_unc = gaussian_process.predict(x_template, return_std=True)
  y_syst = []
  for i in range(240):
    x = 50.125+0.25*i
    if x < 85 or x > 97:
      y_syst.append(0.0)
    else:
      #print(y_template[i])
      #print(y_unc[i])
      #y_syst.append(sqrt(y_unc[i]**2-y_template[i]))
      y_syst.append(y_unc[i])
  fill_th1d(template_th1d, y_template)
  fill_th1d(template_th1d_up, vec_sum(y_template,y_syst))
  fill_th1d(template_th1d_dn, vec_sum(y_template,scalar_mul(-1.0,y_syst)))
  template_th1d_2 = double_bins(template_th1d,240,50.0,110.0)
  template_th1d_up_2 = double_bins(template_th1d_up,240,50.0,110.0)
  template_th1d_dn_2 = double_bins(template_th1d_dn,240,50.0,110.0)

  #debug
  can_debug = ROOT.TCanvas()
  template_th1d_2.Draw()
  can_debug.SaveAs(output_basename+'_debugtemplate.pdf')
  #/debug

  template_roodatahist = ROOT.RooDataHist('','',ROOT.RooArgList(mllg),template_th1d_2)
  template_up_roodatahist = ROOT.RooDataHist('','',ROOT.RooArgList(mllg),template_th1d_up_2)
  template_dn_roodatahist = ROOT.RooDataHist('','',ROOT.RooArgList(mllg),template_th1d_dn_2)
  shape_b = []
  shape_b.append(ROOT.RooHistPdf('shape_b','',ROOT.RooArgSet(mllg),template_roodatahist))
  shape_b.append(ROOT.RooHistPdf('shape_b_up','',ROOT.RooArgSet(mllg),template_up_roodatahist))
  shape_b.append(ROOT.RooHistPdf('shape_b_dn','',ROOT.RooArgSet(mllg),template_dn_roodatahist))

  sig_mu = ROOT.RooRealVar('sig_mu','CB #mu',87.0,93.0)
  sig_sigma = ROOT.RooRealVar('sig_sigma','CB #sigma',1.0,5.0)
  sig_alphal = ROOT.RooRealVar('sig_alphal','CB #alpha_{L}',0.5,9.9)
  sig_alphar = ROOT.RooRealVar('sig_alphar','CB #alpha_{R}',0.5,9.9)
  sig_nl = ROOT.RooRealVar('sig_nl','CB n_{L}',0.5,9.9)
  sig_nr = ROOT.RooRealVar('sig_nr','CB n_{R}',0.5,9.9)
  sig_mu.setVal(90.5)
  sig_sigma.setVal(2.8)
  sig_alphal.setVal(3.0)
  sig_alphar.setVal(3.0)
  sig_nl.setVal(2.5)
  sig_nr.setVal(2.5)
  shape_s = ROOT.RooCrystalBall('signal','',mllg,sig_mu,sig_sigma,sig_alphal,sig_nl,sig_alphar,sig_nr)
  norm_s = ROOT.RooRealVar('norm_s','Signal norm',0.0,1.0) 
  norm_s.setVal(sig_fraction_guess)

  #Do fits
  file_suffix = ['_nom','_up','_dn']
  outputs = []
  for i in range(1):
    shape_sb = ROOT.RooAddPdf('pdf_sb','',ROOT.RooArgList(shape_s,shape_b[i]),ROOT.RooArgList(norm_s))
    shape_sb.fitTo(data)

    plot_postfit = mllg.frame()
    data.plotOn(plot_postfit)
    shape_sb.plotOn(plot_postfit,ROOT.RooFit.Name('postfit_sb'))
    sig_norm_temp = norm_s.getValV()
    bak_norm_temp = 1.0-sig_norm_temp
    norm_s.setVal(0.0)
    shape_sb.plotOn(plot_postfit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
    norm_s.setVal(sig_norm_temp)
    postfit_canvas = ROOT.TCanvas()
    plot_postfit.Draw()
    postfit_canvas.SaveAs(output_basename+file_suffix[i]+'.pdf')
    outputs.append(sig_norm_temp*data.sumEntries())
  outputs.append(0.0)
  outputs.append(0.0)

  return outputs

def generate_hists(year):
  '''Generate weights to correct photons using a simple procedure'''

  
  simu_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
  data_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  if (year == '2018'):
    simu_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    data_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif (year == '2017'):
    simu_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2017/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    data_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2017/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif (year == '2016'):
    simu_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    data_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif (year == '2016APV'):
    simu_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    data_filename = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2016APV/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  else:
    raise ValueError('Unsupported year')

  #load files, filter, and define needed columns
  df_simu = ROOT.RDataFrame('tree',simu_filename)
  df_data = ROOT.RDataFrame('tree',data_filename)

  df_simu = df_simu.Filter('use_event&&nmu==2&&trig_double_mu&&nel==0')
  df_data = df_data.Filter('nmu==2&&trig_double_mu&&nel==0')
  df_simu = df_simu.Define('probe_photon_idx','get_probe_photon(photon_pt,photon_drmin)')
  df_data = df_data.Define('probe_photon_idx','get_probe_photon(photon_pt,photon_drmin)')
  df_simu = df_simu.Filter('probe_photon_idx != -1')
  df_data = df_data.Filter('probe_photon_idx != -1')
  df_simu = df_simu.Define('probe_photon_pt','photon_pt[probe_photon_idx]')
  df_data = df_data.Define('probe_photon_pt','photon_pt[probe_photon_idx]')
  df_simu = df_simu.Define('probe_photon_abseta','fabs(photon_origin_eta[probe_photon_idx])')
  df_data = df_data.Define('probe_photon_abseta','fabs(photon_origin_eta[probe_photon_idx])')
  df_simu = df_simu.Define('probe_photon_idmva','photon_idmva[probe_photon_idx]')
  df_data = df_data.Define('probe_photon_idmva','photon_idmva[probe_photon_idx]')
  df_simu = df_simu.Define('mllph','get_mllph(photon_pt,photon_eta,photon_phi,probe_photon_idx,mu_pt,mu_eta,mu_phi,mu_sig)')
  df_data = df_data.Define('mllph','get_mllph(photon_pt,photon_eta,photon_phi,probe_photon_idx,mu_pt,mu_eta,mu_phi,mu_sig)')
  df_simu = df_simu.Define('w_lumi_pu','w_lumi*w_pu')
  #df_simu = df_simu.Filter('(ll_m[0]+mllph)<175')
  #df_data = df_data.Filter('(ll_m[0]+mllph)<175')

  ##book mllph histograms for reference
  #mass_simu_hist_ptrs = []
  #mass_data_hist_ptrs = []

  #for ipt in range(len(PT_BINS_COARSE)-1):
  #  for ieta in range(len(ETA_BINS_COARSE)-1):
  #    for imva in range(len(MVA_BINS_COARSE)-1):
  #      df_simu_bin = df_simu.Filter('probe_photon_pt>'+str(PT_BINS_COARSE[ipt])+
  #                                   '&&probe_photon_pt<'+str(PT_BINS_COARSE[ipt+1])+
  #                                   '&&probe_photon_abseta>'+str(ETA_BINS_COARSE[ieta])+
  #                                   '&&probe_photon_abseta<'+str(ETA_BINS_COARSE[ieta+1])+
  #                                   '&&probe_photon_idmva>'+str(MVA_BINS_COARSE[imva])+
  #                                   '&&probe_photon_idmva<'+str(MVA_BINS_COARSE[imva+1]))
  #      df_data_bin = df_data.Filter('probe_photon_pt>'+str(PT_BINS_COARSE[ipt])+
  #                                   '&&probe_photon_pt<'+str(PT_BINS_COARSE[ipt+1])+
  #                                   '&&probe_photon_abseta>'+str(ETA_BINS_COARSE[ieta])+
  #                                   '&&probe_photon_abseta<'+str(ETA_BINS_COARSE[ieta+1])+
  #                                   '&&probe_photon_idmva>'+str(MVA_BINS_COARSE[imva])+
  #                                   '&&probe_photon_idmva<'+str(MVA_BINS_COARSE[imva+1]))
  #      mass_simu_hist_ptrs.append(df_simu_bin.Histo1D(('simu_hist','Z#gamma;m_{ll#gamma} [GeV];',60,50.0,110.0),'mllph','w_lumi_pu'))
  #      mass_data_hist_ptrs.append(df_data_bin.Histo1D(('data_hist','Data;m_{ll#gamma} [GeV];',60,50.0,110.0),'mllph'))

  #book histograms to generate weights
  simu_hist_ptrs = []
  data_hist_ptrs = []

  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      for imva in range(len(MVA_BINS)-1):
        df_simu_bin = df_simu.Filter('probe_photon_pt>'+str(PT_BINS[ipt])+
                                     '&&probe_photon_pt<'+str(PT_BINS[ipt+1])+
                                     '&&probe_photon_abseta>'+str(ETA_BINS[ieta])+
                                     '&&probe_photon_abseta<'+str(ETA_BINS[ieta+1])+
                                     '&&probe_photon_idmva>'+str(MVA_BINS[imva])+
                                     '&&probe_photon_idmva<'+str(MVA_BINS[imva+1]))
        df_data_bin = df_data.Filter('probe_photon_pt>'+str(PT_BINS[ipt])+
                                     '&&probe_photon_pt<'+str(PT_BINS[ipt+1])+
                                     '&&probe_photon_abseta>'+str(ETA_BINS[ieta])+
                                     '&&probe_photon_abseta<'+str(ETA_BINS[ieta+1])+
                                     '&&probe_photon_idmva>'+str(MVA_BINS[imva])+
                                     '&&probe_photon_idmva<'+str(MVA_BINS[imva+1]))
        #binning = ROOT.MVA_BINS_BARREL_PTR
        #nbins = len(MVA_BINS_BARREL)-1
        #if (ieta > 1):
        #  binning = ROOT.MVA_BINS_ENDCAP_PTR
        #  nbins = len(MVA_BINS_ENDCAP)-1
        #simu_hist_ptrs.append(book_hist1d_zpeak_bkgsub(df_simu_bin,('simu_hist','Z#gamma;m_{ll#gamma} [GeV];',nbins,binning),'probe_photon_idmva','w_lumi_pu'))
        #data_hist_ptrs.append(book_hist1d_zpeak_bkgsub(df_data_bin,('data_hist','Data;m_{ll#gamma} [GeV];',nbins,binning),'probe_photon_idmva'))
        simu_hist_ptrs.append(df_simu_bin.Histo1D(('simu_hist','Z#gamma;m_{ll#gamma} [GeV];',60,50.0,110.0),'mllph','w_lumi_pu'))
        data_hist_ptrs.append(df_data_bin.Histo1D(('data_hist','Data;m_{ll#gamma} [GeV];',60,50.0,110.0),'mllph'))

  ##make mllph output plots for reference
  #mass_idx = 0
  #for ipt in range(len(PT_BINS_COARSE)-1):
  #  for ieta in range(len(ETA_BINS_COARSE)-1):
  #    for imva in range(len(MVA_BINS_COARSE)-1):
  #      simu_hist = mass_simu_hist_ptrs[mass_idx].GetValue()
  #      data_hist = mass_data_hist_ptrs[mass_idx].GetValue()
  #      simu_hist.Scale(data_hist.Integral()/simu_hist.Integral())
  #      mass_plot = RplPlot()
  #      mass_plot.title_type = 'cms private work'
  #      mass_plot.plot_filled(simu_hist, ROOT.TColor.GetColor('#5790fc'))
  #      mass_plot.plot_points(data_hist, ROOT.kBlack)
  #      mass_plot.add_ratio('data_hist','simu_hist',True)
  #      mass_plot.draw('plots/check_bkgphwgt_mllph_'+str(ipt)+'_'+str(ieta)+'_'+str(imva)+'.pdf')
  #      mass_idx += 1

  ##generate weights and plots
  #weights = []
  #pteta_idx = 0
  #for ipt in range(len(PT_BINS)-1):
  #  for ieta in range(len(ETA_BINS)-1): 
  #    simu_hist = get_bkgsub_hist(simu_hist_ptrs[pteta_idx])
  #    data_hist = get_bkgsub_hist(data_hist_ptrs[pteta_idx])
  #    data_hist.Scale(1.0/data_hist.Integral())
  #    simu_hist.Scale(1.0/simu_hist.Integral())
  #    mva_plot = RplPlot()
  #    mva_plot.title_type = 'cms private work'
  #    mva_plot.y_title = '% Events/bin'
  #    mva_plot.plot_filled(simu_hist, ROOT.TColor.GetColor('#5790fc'))
  #    mva_plot.plot_points(data_hist, ROOT.kBlack)
  #    mva_plot.add_ratio('data_hist','simu_hist',True)
  #    mva_plot.y_title_lower = 'Scale factor'
  #    mva_plot.draw('plots/bkgphwgt_idmva_'+str(ipt)+'_'+str(ieta)+'.pdf')
  #    pteta_idx += 1
  #    sf_hist = data_hist.Clone('sf_hist')
  #    sf_hist.Divide(simu_hist)
  #    weights.append(get_hist_content(sf_hist))
  #print(weights)

  #write histograms to file
  output_file = ROOT.TFile('temp_llphpeak_hists.root','RECREATE')
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      for imva in range(len(MVA_BINS)-1):
        simu_hist_ptrs[ibin].GetValue().Write('simu_hist_bin'+str(ibin))
        data_hist_ptrs[ibin].GetValue().Write('data_hist_bin'+str(ibin))
        ibin += 1
  output_file.Close()

def fit_hists(year):
  """Fit histograms to generate weights"""

  if not year in ['2016APV','2016','2017','2018']:
    raise ValueError('Unsupported year')

  #fit hists
  input_file = ROOT.TFile('temp_llphpeak_hists.root','READ')
  simu_yields = []
  data_yields = []
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    simu_yields.append([])
    data_yields.append([])
    for ieta in range(len(ETA_BINS)-1):
      simu_yields[ipt].append([])
      data_yields[ipt].append([])
      for imva in range(len(MVA_BINS)-1):
        #if not (ieta == 0 and ipt == 0):
        #  ibin += 1
        #  continue
        simu_hist = ROOT.cast_tobject_to_th1d(input_file.Get('simu_hist_bin'+str(ibin)))
        data_hist = ROOT.cast_tobject_to_th1d(input_file.Get('data_hist_bin'+str(ibin)))
        bin_name = 'bin'+str(ipt)+'_'+str(ieta)+'_'+str(imva)
        simu_yields[ipt][ieta].append(gpr_fit_zpeak_yield(simu_hist, 'plots/simu_tnpfit_'+bin_name))
        data_yields[ipt][ieta].append(gpr_fit_zpeak_yield(data_hist, 'plots/data_tnpfit_'+bin_name))
        ibin += 1
  input_file.Close()
  print('simu yields:')
  print(simu_yields)
  print('data yields:')
  print(data_yields)

  #calculate SFs
  simu_total_yield = [0.0,0.0,0.0]
  data_total_yield = [0.0,0.0,0.0]
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      for imva in range(len(MVA_BINS)-1):
        simu_total_yield[0] += simu_yields[ipt][ieta][imva][0]
        data_total_yield[0] += data_yields[ipt][ieta][imva][0]
        simu_total_yield[1] += simu_yields[ipt][ieta][imva][1]
        data_total_yield[1] += data_yields[ipt][ieta][imva][1]
        simu_total_yield[2] += simu_yields[ipt][ieta][imva][2]
        data_total_yield[2] += data_yields[ipt][ieta][imva][2]
  sfs = []
  sfs_unc = []
  sfs_up = []
  sfs_dn = []
  for ipt in range(len(PT_BINS)-1):
    #sfs.append([])
    for ieta in range(len(ETA_BINS)-1):
      #sfs[ipt].append([])
      x_list = []
      y_list_data = []
      y_list_simu = []
      ex_list = []
      ey_list_data = []
      ey_list_simu = []
      for imva in range(len(MVA_BINS)-1):
        if (data_total_yield[0] != 0.0 and simu_total_yield[0] != 0.0 and simu_yields[ipt][ieta][imva][0] != 0.0):
          sfs.append(data_yields[ipt][ieta][imva][0]/data_total_yield[0]/(simu_yields[ipt][ieta][imva][0]/simu_total_yield[0]))
          if (data_yields[ipt][ieta][imva][0] != 0.0):
            data_rel_unc = 1.0/sqrt(data_yields[ipt][ieta][imva][0])
            data_sum_rel_unc = 1.0/sqrt(data_total_yield[0])
            simu_rel_unc = sqrt(0.2835/simu_yields[ipt][ieta][imva][0])
            simu_sum_rel_unc = sqrt(0.2835/simu_total_yield[0])
            sfs_unc.append(sqrt(data_rel_unc**2+data_sum_rel_unc**2+simu_rel_unc**2+simu_sum_rel_unc**2)*sfs[-1])
          else:
            sfs_unc.append(0.0)
        else:
          sfs.append(1.0)
          sfs_unc.append(0.0)
        #if (data_total_yield[1] != 0.0 and simu_total_yield[1] != 0.0):
        #  sfs_up.append(data_yields[ipt][ieta][imva][1]/data_total_yield[1]/(simu_yields[ipt][ieta][imva][1]/simu_total_yield[1]))
        #else:
        #  sfs_up.append(1.0)
        #if (data_total_yield[2] != 0.0 and simu_total_yield[2] != 0.0):
        #  sfs_dn.append(data_yields[ipt][ieta][imva][2]/data_total_yield[2]/(simu_yields[ipt][ieta][imva][2]/simu_total_yield[2]))
        #else:
        #  sfs_dn.append(1.0)
        x_width = MVA_BINS[imva+1]-MVA_BINS[imva]
        if data_total_yield[0] != 0.0:
          y_list_data.append(data_yields[ipt][ieta][imva][0]/data_total_yield[0]/x_width)
          ey_list_data.append(sqrt(data_yields[ipt][ieta][imva][0])/data_total_yield[0]/x_width)
        else:
          y_list_data.append(0.0)
          ey_list_data.append(0.0)
        if simu_total_yield[0] != 0.0:
          y_list_simu.append(simu_yields[ipt][ieta][imva][0]/simu_total_yield[0]/x_width)
          ey_list_simu.append(sqrt(0.2835*simu_yields[ipt][ieta][imva][0])/simu_total_yield[0]/x_width)
        else:
          y_list_simu.append(0.0)
          ey_list_simu.append(0.0)
        x_list.append((MVA_BINS[imva+1]+MVA_BINS[imva])/2.0)
        ex_list.append((MVA_BINS[imva+1]-MVA_BINS[imva])/2.0)
        ibin += 1
      x_vals = array.array('d',x_list)
      y_vals_data = array.array('d',y_list_data)
      y_vals_simu = array.array('d',y_list_simu)
      ex_vals = array.array('d',ex_list)
      ey_vals_data = array.array('d',ey_list_data)
      ey_vals_simu = array.array('d',ey_list_simu)
      data_graph = ROOT.TGraphErrors(len(x_list),x_vals,y_vals_data,ex_vals,ey_vals_data)
      data_graph.SetTitle('data')
      simu_graph = ROOT.TGraphErrors(len(x_list),x_vals,y_vals_simu,ex_vals,ey_vals_simu)
      data_graph.SetTitle('simulation')
      sf_plot = RplPlot()
      sf_plot.lumi_data = [(60,13)]
      sf_plot.plot_graph(data_graph)
      sf_plot.plot_graph(simu_graph)
      sf_plot.x_title = 'Photon MVA ID'
      sf_plot.y_title = '% Events/ID score'
      sf_plot.draw('plots/simplemvashapes_'+year+'_bin'+str(ipt)+'_'+str(ieta)+'.pdf')

  #make graph and save SFs
  ibin = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      x_list = []
      y_list = []
      ex_list = []
      ey_list = []
      for imva in range(len(MVA_BINS)-1):
        x_list.append((MVA_BINS[imva+1]+MVA_BINS[imva])/2.0)
        ex_list.append((MVA_BINS[imva+1]-MVA_BINS[imva])/2.0)
        y_list.append(sfs[ibin])
        ey_list.append(sfs_unc[ibin])
        ibin += 1
      x_vals = array.array('d',x_list)
      y_vals = array.array('d',y_list)
      ex_vals = array.array('d',ex_list)
      ey_vals = array.array('d',ey_list)
      sf_graph = ROOT.TGraphErrors(len(x_list),x_vals,y_vals,ex_vals,ey_vals)
      sf_plot = RplPlot()
      sf_plot.lumi_data = [(60,13)]
      sf_plot.plot_graph(sf_graph)
      sf_plot.x_title = 'Photon MVA ID'
      sf_plot.y_title = 'Scale Factor'
      sf_plot.draw('plots/simplesf_'+year+'_bin'+str(ipt)+'_'+str(ieta)+'.pdf')

  sf_set = schemav2.Correction(
      name='photon_idmva_sf',
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
              schemav2.Variable(name='abseta', type='real', description='Photon absolute eta'),
              schemav2.Variable(name='mva', type='real', description='Photon Fall17v2 IDMVA')],
      output=schemav2.Variable(name='sf', type='real', description='Data-MC SF'),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','abseta','mva'],
          edges=[PT_BINS,ETA_BINS,MVA_BINS],
          content=sfs,
          flow='clamp',
          ),
      )
  #sf_up_set = schemav2.Correction(
  #    name='photon_idmva_sf_up',
  #    version=1,
  #    inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
  #            schemav2.Variable(name='abseta', type='real', description='Photon absolute eta'),
  #            schemav2.Variable(name='mva', type='real', description='Photon Fall17v2 IDMVA')],
  #    output=schemav2.Variable(name='sf_up', type='real', description='Data-MC SF up variation'),
  #    data=schemav2.MultiBinning(
  #        nodetype='multibinning',
  #        inputs=['pt','abseta','mva'],
  #        edges=[PT_BINS,ETA_BINS,MVA_BINS],
  #        content=sfs,
  #        flow='clamp',
  #        ),
  #    )
  #sf_dn_set = schemav2.Correction(
  #    name='photon_idmva_sf_dn',
  #    version=1,
  #    inputs=[schemav2.Variable(name='pt', type='real', description='Photon pt'),
  #            schemav2.Variable(name='abseta', type='real', description='Photon absolute eta'),
  #            schemav2.Variable(name='mva', type='real', description='Photon Fall17v2 IDMVA')],
  #    output=schemav2.Variable(name='sf_dn', type='real', description='Data-MC SF dn variation'),
  #    data=schemav2.MultiBinning(
  #        nodetype='multibinning',
  #        inputs=['pt','abseta','mva'],
  #        edges=[PT_BINS,ETA_BINS,MVA_BINS],
  #        content=sfs,
  #        flow='clamp',
  #        ),
  #    )

  with open('json/simple_photon_weights'+year+'.json','w') as output_file:
    output_file.write(fix_correctionlib_json(
      [sf_set.json(exclude_unset=True)]))
       #sf_up_set.json(exclude_unset=True),
       #sf_dn_set.json(exclude_unset=True)]))

def write_weights():
  """Save photon weights as correctionlib file"""

if __name__ == '__main__':

  year = '2017'
  print('Generating histograms')
  generate_hists(year)
  print('Fitting histograms')
  fit_hists(year)
