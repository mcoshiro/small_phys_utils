#!/usr/bin/env python3
'''
Script that generates slimmed n-tuples for photon Fall17v2 IDMVA corrections
'''
import subprocess
import ROOT
import skim_and_slim

ROOT.gInterpreter.ProcessLine('#include "nn_model_phid2018_old.hpp"')
ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;

const nn_model_phid2018 dnn;

float pass_ele_fall17v2_wp90(float mva_score, float el_pt, float el_sc_abseta) {
  float cut = -1.0;
  if (el_pt < 10) {
    if (el_sc_abseta < 0.8)
      cut = 2.84704783417 - TMath::Exp(-1.0*el_pt / 3.32529515837) * 9.38050947827;
    else if (el_sc_abseta < 1.4446) 
      cut = 2.03833922005 - TMath::Exp(-1.0*el_pt / 1.93288758682) * 15.364588247;
    cut = 1.82704158461 - TMath::Exp(-1.0*el_pt / 1.89796754399) * 19.1236071158;
  }
  else {
    if (el_sc_abseta < 0.8)
      cut = 6.12931925263 - TMath::Exp(-1.0*el_pt / 13.281753835) * 8.71138432196;
    else if (el_sc_abseta < 1.4446) 
      cut = 5.26289004857 - TMath::Exp(-1.0*el_pt / 13.2154971491) * 8.0997882835;
    cut = 4.37338792902 - TMath::Exp(-1.0*el_pt / 14.0776094696) * 8.48513324496;
  } 
  return cut;
}

bool pass_ele_fall17v2_wpl(float mva_score, float el_pt, float el_sc_abseta) {
  float cut = -1.0;
  if (el_pt < 10) {
    if (el_sc_abseta < 0.8)
      cut = 0.700642584415;
    else if (el_sc_abseta < 1.4446) 
      cut = 0.739335420875;
    cut = 1.45390456109;
  }
  else {
    if (el_sc_abseta < 0.8)
      cut = -0.146270871164;
    else if (el_sc_abseta < 1.4446) 
      cut = -0.0315850882679;
    cut = -0.0321841194737;
  } 
  return (mva_score > cut);
}

float get_rw(float ph_et, float ph_sc_abseta) {
  return 5.2200224/dnn.get_hist_interp(ph_et, ph_sc_abseta);
}

""")

if __name__=='__main__':
  ROOT.EnableImplicitMT()
  cuts = ['tag_Ele_pt>35&&tag_sc_abseta<2.17&&(tag_sc_abseta>1.566||tag_sc_abseta<1.4442)',
          'ph_et>=15&&ph_sc_abseta<2.5&&((ph_sc_abseta<1.4442&&ph_mva94XV2>-0.4)||(ph_sc_abseta>1.566&&ph_mva94XV2>-0.58))']
  defines_mc = [('tag_ele_fall17v2_wp90','pass_ele_fall17v2_wp90(tag_Ele_IsoMVA94XV2,tag_Ele_pt,tag_sc_abseta)'),
                ('tag_ele_fall17v2_wpl','pass_ele_fall17v2_wpl(tag_Ele_IsoMVA94XV2,tag_Ele_pt,tag_sc_abseta)'),
                ('ph_esEnergyOverRawE','static_cast<float>((ph_preshower_energy_plane1+ph_preshower_energy_plane2)/ph_sc_rawEnergy)'),
                ('w','get_rw(ph_et, ph_sc_abseta)')]
  defines_data = [('tag_ele_fall17v2_wpl','pass_ele_fall17v2_wpl(tag_Ele_IsoMVA94XV2,tag_Ele_pt,tag_sc_abseta)'),
                  ('ph_esEnergyOverRawE','static_cast<float>((ph_preshower_energy_plane1+ph_preshower_energy_plane2)/ph_sc_rawEnergy)'),
                  ('w','0')]
  branches = ('pair_mass','ph_et','ph_eta','ph_r9','ph_s4','ph_sc_etaWidth',
              'ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso',
              'ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','ph_sc_abseta','event_rho',
              'ph_ESsigma','ph_esEnergyOverRawE','ph_mva94XV2','w','tag_Ele_IsoMVA94XV2')
  simu_files_2018 = ['/net/cms26/cms26r0/oshiro/tnp_tuples/DY_LO2018.root']
  data_files_2018 = ['/net/cms26/cms26r0/oshiro/tnp_tuples/Run2018A.root',
                     '/net/cms26/cms26r0/oshiro/tnp_tuples/Run2018B.root',
                     '/net/cms26/cms26r0/oshiro/tnp_tuples/Run2018C.root',
                     '/net/cms26/cms26r0/oshiro/tnp_tuples/Run2018D.root']
  skim_and_slim.write_ntuples(
      data_files_2018,
      cuts,
      'photonidskim_data_2018_new.root',
      defines_data,
      'tnpPhoIDs/fitter_tree',
      branches,
      '/net/cms26/cms26r0/oshiro/tnp_tuples/')

