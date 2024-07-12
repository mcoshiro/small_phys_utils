#!/usr/bin/env python3
"""@package docstring
Generates correction scale factors for the IDMVA distribution of fake photons
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from math import sqrt
from root_plot_lib import RplPlot
import array
import ctypes
import json
import ROOT
import subprocess

#constants
PT_BINS = [15.0,17.5,20.0,25.0,30.0,500.0]
ETA_BINS = [0.0,0.8,1.5,2.5]
MVA_BINS = [-1.0,-0.4,-0.2,0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
year = '2017'

#ROOT JIT'ed C++ definitions
ROOT.gInterpreter.AddIncludePath('inc/')
ROOT.gInterpreter.ProcessLine('#include "correction_wrapper.hpp"')
ROOT.gSystem.Load('libSmallPhysUtils.so')
ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using ROOT::Math::PtEtaPhiMVector;
using std::vector;

const CorrectionWrapper corr_real("json/simple_photon_weights"""+year+""".json","photon_idmva_sf");
const CorrectionWrapper corr_fake("json/simple_fakephoton_weights"""+year+""".json","fakephoton_idmva_sf");

const double mvabins_cpp[12] = {-1.0,-0.4,-0.2,0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

float get_realweight(float ph_et, float ph_sc_abseta, float idmva) {
  vector<double> eval_args;
  eval_args.push_back(ph_et);
  eval_args.push_back(ph_sc_abseta);
  eval_args.push_back(idmva);
  return corr_real.evaluate(eval_args);
}

float get_fakeweight(float ph_et, float ph_sc_abseta, float idmva) {
  vector<double> eval_args;
  eval_args.push_back(ph_et);
  eval_args.push_back(ph_sc_abseta);
  eval_args.push_back(idmva);
  return corr_fake.evaluate(eval_args);
}

//returns H->Zy signal photon: highest pT photon passing WP80, or highest pT photon if none passing WP80
int get_probe_photon(RVec<float> photon_pt, RVec<float> photon_drmin, 
                     RVec<bool> photon_elveto, RVec<bool> photon_id80, 
                     RVec<int> photon_elidx, RVec<bool> el_sig) {
  int ph_idx = -1;
  bool found_wp80 = false;
  float max_photon_pt = 0;
  for (unsigned iph = 0; iph < photon_pt.size(); iph++ ) {
    if (photon_drmin[iph]>0.03 && photon_pt[iph]>15 && photon_elveto[iph]) {
      bool pass_overlap = true;
      if (photon_elidx[iph] != -1) {
        pass_overlap = el_sig[photon_elidx[iph]];
      }
      if (pass_overlap) {
        if (photon_id80[iph] && !found_wp80) {
          ph_idx = static_cast<int>(iph);
          max_photon_pt = photon_pt[iph];
          found_wp80 = true;
        }
        else if (photon_id80[iph] && found_wp80) {
          if (photon_pt[iph] > max_photon_pt) {
            ph_idx = static_cast<int>(iph);
            max_photon_pt = photon_pt[iph];
          }
        }
        else if (!found_wp80) {
          if (photon_pt[iph] > max_photon_pt) {
            ph_idx = static_cast<int>(iph);
            max_photon_pt = photon_pt[iph];
          }
        }
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

//currently 2018
bool pass_trig_cuts(bool trig_single_el, bool trig_double_el, 
                    bool trig_single_mu, bool trig_double_mu, int nel, 
                    int nmu, RVec<float> el_pt, RVec<float> mu_pt) {
  if (trig_single_el && nel >= 1)
    return (el_pt[0]>35);
  else if (trig_double_el && nel >= 2)
    return (el_pt[0]>25 && el_pt[1]>15);
  else if (trig_single_mu && nmu >= 1)
    return (mu_pt[0]>25);
  else if (trig_double_mu && nmu >= 2)
    return (mu_pt[0]>20 && mu_pt[1]>10);
  return false;
}

const vector<float> MVA_BINS_BARREL = {-1.0,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0};
const vector<float> MVA_BINS_ENDCAP = {-1.0,-0.58,-0.5,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0};

const float* MVA_BINS_BARREL_PTR = &MVA_BINS_BARREL[0];
const float* MVA_BINS_ENDCAP_PTR = &MVA_BINS_ENDCAP[0];

TH1D* cast_tobject_to_th1d(TObject* object)
{ return static_cast<TH1D*>(object); }

""")

def sum_quadrature(values):
  '''Return sum of values in list in quadrature'''
  quad_sum = 0
  for value in values:
    quad_sum += value**2
  return sqrt(quad_sum)

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

def generate_weights(year):
  '''Generate scale factors for correcting fake photon IDMVA distribution'''

  file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
  file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
  file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  if year == '2018':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif year == '2017':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2017/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2017/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2017/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif year == '2016':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif year == '2016APV':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016APV/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016APV/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016APV/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  else:
    raise ValueError('Unsupported year')

  #load files, filter, and define needed columns
  #selections: overlap removal, 2 muons, dimuon trigger, 1 reco photon,
  #            80<mll<100, mllph 100-120 or 130-180
  #nb: compared with baseline, no mll+mllph, no ptph/mllph, no trigger plateau, no e channel
  df_real = ROOT.RDataFrame('tree',file_real)
  df_fake = ROOT.RDataFrame('tree',file_fake)
  df_data = ROOT.RDataFrame('tree',file_data)

  df_real = df_real.Filter('use_event&&nmu==2&&trig_double_mu&&nel==0')
  df_fake = df_fake.Filter('use_event&&nmu==2&&trig_double_mu&&nel==0')
  df_data = df_data.Filter('nmu==2&&trig_double_mu&&nel==0')

  df_real = df_real.Filter('ll_m[0]>80&&ll_m[0]<100')
  df_fake = df_fake.Filter('ll_m[0]>80&&ll_m[0]<100')
  df_data = df_data.Filter('ll_m[0]>80&&ll_m[0]<100')

  df_real = df_real.Define('probe_photon_idx','get_probe_photon(photon_pt,photon_drmin,photon_elveto,photon_id80,photon_elidx,el_sig)')
  df_fake = df_fake.Define('probe_photon_idx','get_probe_photon(photon_pt,photon_drmin,photon_elveto,photon_id80,photon_elidx,el_sig)')
  df_data = df_data.Define('probe_photon_idx','get_probe_photon(photon_pt,photon_drmin,photon_elveto,photon_id80,photon_elidx,el_sig)')

  df_real = df_real.Filter('probe_photon_idx != -1')
  df_fake = df_fake.Filter('probe_photon_idx != -1')
  df_data = df_data.Filter('probe_photon_idx != -1')

  df_real = df_real.Define('probe_photon_pt','photon_pt[probe_photon_idx]')
  df_fake = df_fake.Define('probe_photon_pt','photon_pt[probe_photon_idx]')
  df_data = df_data.Define('probe_photon_pt','photon_pt[probe_photon_idx]')

  df_real = df_real.Define('probe_photon_abseta','fabs(photon_origin_eta[probe_photon_idx])')
  df_fake = df_fake.Define('probe_photon_abseta','fabs(photon_origin_eta[probe_photon_idx])')
  df_data = df_data.Define('probe_photon_abseta','fabs(photon_origin_eta[probe_photon_idx])')

  df_real = df_real.Define('probe_photon_idmva','photon_idmva[probe_photon_idx]')
  df_fake = df_fake.Define('probe_photon_idmva','photon_idmva[probe_photon_idx]')
  df_data = df_data.Define('probe_photon_idmva','photon_idmva[probe_photon_idx]')

  df_real = df_real.Define('mllph','get_mllph(photon_pt,photon_eta,photon_phi,probe_photon_idx,mu_pt,mu_eta,mu_phi,mu_sig)')
  df_fake = df_fake.Define('mllph','get_mllph(photon_pt,photon_eta,photon_phi,probe_photon_idx,mu_pt,mu_eta,mu_phi,mu_sig)')
  df_data = df_data.Define('mllph','get_mllph(photon_pt,photon_eta,photon_phi,probe_photon_idx,mu_pt,mu_eta,mu_phi,mu_sig)')

  df_real = df_real.Filter('(mllph>100&mllph<120)||(mllph>130&&mllph<180)')
  df_fake = df_fake.Filter('(mllph>100&mllph<120)||(mllph>130&&mllph<180)')
  df_data = df_data.Filter('(mllph>100&mllph<120)||(mllph>130&&mllph<180)')

  df_real = df_real.Define('w_ph_lumi_year_pu','get_realweight(probe_photon_pt,probe_photon_abseta,probe_photon_idmva)*w_lumi*w_pu*60.0')
  df_fake = df_fake.Define('w_lumi_year_pu','w_lumi*w_pu*60.0')

  #book histograms to generate weights
  real_hist_ptrs = []
  fake_hist_ptrs = []
  data_hist_ptrs = []

  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      df_real_bin = df_real.Filter('probe_photon_pt>'+str(PT_BINS[ipt])+
                                   '&&probe_photon_pt<'+str(PT_BINS[ipt+1])+
                                   '&&probe_photon_abseta>'+str(ETA_BINS[ieta])+
                                   '&&probe_photon_abseta<'+str(ETA_BINS[ieta+1]))
      df_fake_bin = df_fake.Filter('probe_photon_pt>'+str(PT_BINS[ipt])+
                                   '&&probe_photon_pt<'+str(PT_BINS[ipt+1])+
                                   '&&probe_photon_abseta>'+str(ETA_BINS[ieta])+
                                   '&&probe_photon_abseta<'+str(ETA_BINS[ieta+1]))
      df_data_bin = df_data.Filter('probe_photon_pt>'+str(PT_BINS[ipt])+
                                   '&&probe_photon_pt<'+str(PT_BINS[ipt+1])+
                                   '&&probe_photon_abseta>'+str(ETA_BINS[ieta])+
                                   '&&probe_photon_abseta<'+str(ETA_BINS[ieta+1]))
      real_hist_ptrs.append(df_real_bin.Histo1D(('real_hist','Z#gamma;Photon IDMVA;',11,ROOT.mvabins_cpp),'probe_photon_idmva','w_ph_lumi_year_pu'))
      fake_hist_ptrs.append(df_fake_bin.Histo1D(('fake_hist','Z+fake photon;Photon IDMVA;',11,ROOT.mvabins_cpp),'probe_photon_idmva','w_lumi_year_pu'))
      data_hist_ptrs.append(df_data_bin.Histo1D(('data_hist','Data;Photon IDMVA;',11,ROOT.mvabins_cpp),'probe_photon_idmva'))

  print('Processing data')
  #draw IDMVA plots used to derive weights
  plot_idx = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      real_hist_copy = real_hist_ptrs[plot_idx].GetValue().Clone()
      fake_hist_copy = fake_hist_ptrs[plot_idx].GetValue().Clone()
      data_hist_copy = data_hist_ptrs[plot_idx].GetValue().Clone()
      overall_scale = data_hist_copy.Integral()/(real_hist_copy.Integral()+fake_hist_copy.Integral())
      real_hist_copy.Scale(overall_scale)
      fake_hist_copy.Scale(overall_scale)
      for imva in range(len(MVA_BINS)-1):
        bin_width = MVA_BINS[imva+1]-MVA_BINS[imva]
        real_hist_copy.SetBinContent(imva+1, real_hist_copy.GetBinContent(imva+1)/bin_width)
        fake_hist_copy.SetBinContent(imva+1, fake_hist_copy.GetBinContent(imva+1)/bin_width)
        data_hist_copy.SetBinContent(imva+1, data_hist_copy.GetBinContent(imva+1)/bin_width)
      fake_hist_copy.Add(real_hist_copy)
      fakeratio_plot = RplPlot()
      fakeratio_plot.title_type = 'cms private work'
      fakeratio_plot.y_title = 'Events/IDMVA'
      fakeratio_plot.lumi_data = [(60,13)]
      fakeratio_plot.plot_filled(fake_hist_copy, ROOT.TColor.GetColor('#f89c20'))
      fakeratio_plot.plot_filled(real_hist_copy, ROOT.TColor.GetColor('#5790fc'))
      fakeratio_plot.plot_points(data_hist_copy, ROOT.kBlack)
      fakeratio_plot.add_ratio('data_hist','fake_hist',True)
      fakeratio_plot.draw('plots/fakephoton_weightderivation_'+year+'_'+str(ipt)+'_'+str(ieta)+'.pdf')
      plot_idx += 1

  print('Calculating weights')
  #calculate weights
  sfs = []
  sfs_unc = []
  plot_idx = 0
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      real_hist = real_hist_ptrs[plot_idx].GetValue()
      fake_hist = fake_hist_ptrs[plot_idx].GetValue()
      data_hist = data_hist_ptrs[plot_idx].GetValue()
      overall_scale = data_hist.Integral()/(real_hist.Integral()+fake_hist.Integral())
      print('Overall scale factor for bin '+str(ipt)+', '+str(ieta)+': '+str(overall_scale))
      real_hist.Scale(overall_scale)
      fake_hist.Scale(overall_scale)
      data_hist.Add(real_hist, -1.0)
      data_total_yield = data_hist.Integral()
      fake_total_yield = fake_hist.Integral()
      for imva in range(len(MVA_BINS)-1):
        data_bin_yield = data_hist.GetBinContent(imva+1)
        fake_bin_yield = fake_hist.GetBinContent(imva+1)
        if data_bin_yield > 0.0 and data_total_yield > 0.0 and fake_bin_yield > 0.0:
          sfs.append(data_bin_yield/data_total_yield*fake_total_yield/fake_bin_yield)
          if fake_total_yield > 0.0:
            sfs_unc.append(sum_quadrature([1.0/sqrt(data_total_yield),
                                           1.0/sqrt(fake_total_yield),
                                           1.0/sqrt(data_bin_yield),
                                           1.0/sqrt(fake_bin_yield)])*sfs[-1])
          else:
            sfs_unc.append(sfs[-1]) #could do something fancier, but this is conservative
        else:
          sfs.append(1.0)
          sfs_unc.append(1.0)
      plot_idx += 1

  print('Writing output')
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
      sf_plot.y_title = 'Fake Photon Scale Factor'
      sf_plot.draw('plots/simplesffake_'+year+'_bin'+str(ipt)+'_'+str(ieta)+'.pdf')

  sf_set = schemav2.Correction(
      name='fakephoton_idmva_sf',
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

  with open('json/simple_fakephoton_weights'+year+'.json','w') as output_file:
    output_file.write(fix_correctionlib_json(
      [sf_set.json(exclude_unset=True)]))

def validate_weights(year):
  """Generate some histograms to demonstrate weight (non)closure"""

  file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
  file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
  file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  if year == '2018':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif year == '2017':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2017/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2017/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2017/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif year == '2016':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  elif year == '2016APV':
    file_real = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016APV/mc/merged_zgmc_ll/merged_pico_ll_ZGToLLG_01J_5f_lowMLL*.root'
    file_fake = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016APV/mc/merged_zgmc_ll/merged_pico_ll_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root'
    file_data = '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2016APV/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root'
  else:
    raise ValueError('Unsupported year')

  #load files, filter, and define needed columns
  #selections: overlap removal, 2 leptons + photons, pT cuts and plateau cuts
  #            80<mll<100, mllph 100-120 or 130-180
  #            mll+mllph < , ptph/mllph < 
  df_real = ROOT.RDataFrame('tree',file_real)
  df_fake = ROOT.RDataFrame('tree',file_fake)
  df_data = ROOT.RDataFrame('tree',file_data)

  df_real = df_real.Filter('use_event&&nll>=1&&nphoton>=1&&pass_trig_cuts(trig_single_el,trig_double_el,trig_single_mu,trig_double_mu,nel,nmu,el_pt,mu_pt)')
  df_fake = df_fake.Filter('use_event&&nll>=1&&nphoton>=1&&pass_trig_cuts(trig_single_el,trig_double_el,trig_single_mu,trig_double_mu,nel,nmu,el_pt,mu_pt)')
  df_data = df_data.Filter('nlep>=2&&nphoton>=1&&pass_trig_cuts(trig_single_el,trig_double_el,trig_single_mu,trig_double_mu,nel,nmu,el_pt,mu_pt)')

  df_real = df_real.Filter('photon_id80[0]')
  df_fake = df_fake.Filter('photon_id80[0]')
  df_data = df_data.Filter('photon_id80[0]')

  df_real = df_real.Filter('(ll_m[0]>80&&ll_m[0]<100)&&((llphoton_m[0]>100&&llphoton_m[0]<120)||(llphoton_m[0]>130&&llphoton_m[0]<180))')
  df_fake = df_fake.Filter('(ll_m[0]>80&&ll_m[0]<100)&&((llphoton_m[0]>100&&llphoton_m[0]<120)||(llphoton_m[0]>130&&llphoton_m[0]<180))')
  df_data = df_data.Filter('(ll_m[0]>80&&ll_m[0]<100)&&((llphoton_m[0]>100&&llphoton_m[0]<120)||(llphoton_m[0]>130&&llphoton_m[0]<180))')

  df_real = df_real.Filter('((ll_m[0]+llphoton_m[0])>185)&&((photon_pt[0]/llphoton_m[0])>(15.0/110.0))')
  df_fake = df_fake.Filter('((ll_m[0]+llphoton_m[0])>185)&&((photon_pt[0]/llphoton_m[0])>(15.0/110.0))')
  df_data = df_data.Filter('((ll_m[0]+llphoton_m[0])>185)&&((photon_pt[0]/llphoton_m[0])>(15.0/110.0))')

  df_real = df_real.Define('lead_photon_pt','photon_pt[0]')
  df_fake = df_fake.Define('lead_photon_pt','photon_pt[0]')
  df_data = df_data.Define('lead_photon_pt','photon_pt[0]')

  df_real = df_real.Define('lead_photon_eta','photon_eta[0]')
  df_fake = df_fake.Define('lead_photon_eta','photon_eta[0]')
  df_data = df_data.Define('lead_photon_eta','photon_eta[0]')

  df_real = df_real.Define('lead_photon_absorigineta','fabs(photon_origin_eta[0])')
  df_fake = df_fake.Define('lead_photon_absorigineta','fabs(photon_origin_eta[0])')
  df_data = df_data.Define('lead_photon_absorigineta','fabs(photon_origin_eta[0])')

  #df_real = df_real.Filter('lead_photon_absorigineta<2.5&&lead_photon_absorigineta>1.5')
  #df_fake = df_fake.Filter('lead_photon_absorigineta<2.5&&lead_photon_absorigineta>1.5')
  #df_data = df_data.Filter('lead_photon_absorigineta<2.5&&lead_photon_absorigineta>1.5')

  df_real = df_real.Define('lead_photon_idmva','photon_idmva[0]')
  df_fake = df_fake.Define('lead_photon_idmva','photon_idmva[0]')
  df_data = df_data.Define('lead_photon_idmva','photon_idmva[0]')

  df_real = df_real.Define('mllph','llphoton_m[0]')
  df_fake = df_fake.Define('mllph','llphoton_m[0]')
  df_data = df_data.Define('mllph','llphoton_m[0]')

  df_real = df_real.Define('w_lumi_year_pu','w_lumi*w_pu*60.0')
  df_fake = df_fake.Define('w_lumi_year_pu','w_lumi*w_pu*60.0')

  df_real = df_real.Define('w_ph_lumi_year_pu','get_realweight(lead_photon_pt,lead_photon_absorigineta,lead_photon_idmva)*w_lumi_year_pu')
  df_fake = df_fake.Define('w_ph_lumi_year_pu','get_fakeweight(lead_photon_pt,lead_photon_absorigineta,lead_photon_idmva)*w_lumi_year_pu')

  #make histograms
  print('validating weights')
  print('booking histograms')
  var_names = ('lead_photon_pt','lead_photon_eta','lead_photon_idmva','mllph')
  var_titles = ('Photon p_{T} [GeV]','Photon #eta','Photon IDMVA','m_{ll#gamma} [GeV]')
  var_bins = ((80,15.0,100.0),(40,-2.5,2.5),(40,-1.0,1.0),(40,100,180))

  #convention: lists of lists with outer list corresponding to variable, inner list corresponding to [real, fake, data]
  hist_ptrs_lumi = []
  hist_ptrs_realsf = []
  hist_ptrs_allsf = []

  for ivar in range(len(var_names)):
    hist_ptrs_lumi.append([df_real.Histo1D(('real_hist','Z#gamma;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar],'w_lumi_year_pu'),
                           df_fake.Histo1D(('fake_hist','Z+fake photon;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar],'w_lumi_year_pu'),
                           df_data.Histo1D(('data_hist','Data;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar])])
    hist_ptrs_realsf.append([df_real.Histo1D(('real_hist','Z#gamma;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar],'w_ph_lumi_year_pu'),
                             df_fake.Histo1D(('fake_hist','Z+fake photon;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar],'w_lumi_year_pu'),
                             df_data.Histo1D(('data_hist','Data;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar])])
    hist_ptrs_allsf.append([df_real.Histo1D(('real_hist','Z#gamma;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar],'w_ph_lumi_year_pu'),
                            df_fake.Histo1D(('fake_hist','Z+fake photon;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar],'w_ph_lumi_year_pu'),
                            df_data.Histo1D(('data_hist','Data;'+var_titles[ivar]+';',var_bins[ivar][0],var_bins[ivar][1],var_bins[ivar][2]),var_names[ivar])])

  #make nice plots
  print('drawing output')
  hist_ptrs = (hist_ptrs_lumi, hist_ptrs_realsf, hist_ptrs_allsf)
  hist_type_names = ('wlumi','wrealphoton','wallphoton')

  for ivar in range(len(var_names)):
    for hist_ptr_type in range(3):
      real_hist = hist_ptrs[hist_ptr_type][ivar][0].GetValue()
      fake_hist = hist_ptrs[hist_ptr_type][ivar][1].GetValue()
      data_hist = hist_ptrs[hist_ptr_type][ivar][2].GetValue()
      overall_scale = data_hist.Integral()/(real_hist.Integral()+fake_hist.Integral())
      real_hist.Scale(overall_scale)
      fake_hist.Scale(overall_scale)
      fake_hist.Add(real_hist)
      fakeratio_plot = RplPlot()
      fakeratio_plot.title_type = 'cms private work'
      fakeratio_plot.y_title = 'Events/bin'
      fakeratio_plot.lumi_data = [(60,13)]
      fakeratio_plot.plot_filled(fake_hist, ROOT.TColor.GetColor('#f89c20'))
      fakeratio_plot.plot_filled(real_hist, ROOT.TColor.GetColor('#5790fc'))
      fakeratio_plot.plot_points(data_hist, ROOT.kBlack)
      fakeratio_plot.add_ratio('data_hist','fake_hist',True)
      fakeratio_plot.draw('plots/validatephoton_'+year+'_'+var_names[ivar]+'_'+hist_type_names[hist_ptr_type]+'.pdf')

if __name__ == '__main__':
  #generate_weights(year)
  validate_weights(year)

