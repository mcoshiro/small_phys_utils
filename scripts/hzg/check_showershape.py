#!/usr/bin/env python3
"""@package docstring
Generates photon preselection corrections using tag-and-probe
"""

from argparse import ArgumentParser
from correctionlib import schemav2 
from tnp_utils import do_tnp_fit, integrate_tgraph
from root_plot_lib import RplPlot
from math import sqrt
import ROOT
import json
import subprocess
import ctypes

ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using std::vector;
using ROOT::Math::PtEtaPhiMVector;

const float PI = 3.141593;

float delta_phi(float phi1, float phi2){
  float dphi = fmod(fabs(phi2-phi1), 2.0*PI);
  return dphi>PI ? 2.0*PI-dphi : dphi;
}

float delta_r(float eta1, float eta2, float phi1, float phi2) {
  float deta = eta1-eta2;
  float dphi = delta_phi(phi1, phi2);
  return sqrt(deta*deta+dphi*dphi);
}

int get_el_photon_idx(RVec<float> photon_pt, RVec<bool> photon_isScEtaEB,
                           RVec<bool> photon_isScEtaEE, RVec<float> photon_idmva,
                           RVec<float> photon_eta, RVec<float> photon_phi,
                           RVec<float> el_sig, RVec<float> el_eta, 
                           RVec<float> el_phi, RVec<int> el_charge) {
  //get highest pT photon passing criteria except eVeto, overlap removal, and minLepDR cut
  //w.r.t. postive electrons (positrons) only
  //in order to get "photon" that is electron for tag-and-probe
  int best_photon_idx = -1;
  float max_photon_pt = 0;
  for (unsigned iph = 0; iph < photon_pt.size(); iph++) {
    if (((photon_isScEtaEB[iph] && photon_idmva[iph] > -0.4) || 
         (photon_isScEtaEE[iph] && photon_idmva[iph] > -0.58)) &&
         photon_pt[iph] > 15) {
      if (photon_pt[iph] > max_photon_pt) {
        //check dR from positive electrons
        float min_dr = 999;
        for (unsigned iel = 0; iel < el_sig.size(); iel++) {
          if (el_sig[iel] && (el_charge[iel]==1)) {
            float dr = delta_r(photon_eta[iph], el_eta[iel], photon_phi[iph], el_phi[iel]);
            if (dr < min_dr) min_dr = dr;
          }
        }
        if (min_dr > 0.4) {
          max_photon_pt = photon_pt[iph];
          best_photon_idx = static_cast<int>(iph);
        }
      }
    }
  }
  return best_photon_idx;
}

float get_el_z_mass(RVec<float> el_pt, RVec<float> el_eta, RVec<float> el_phi,
                    RVec<bool> el_sig, RVec<int> el_charge,
                    int el_photon_idx, RVec<float> photon_pt, 
                    RVec<float> photon_eta, RVec<float> photon_phi) {
  if (el_photon_idx == -1) return -999.0;
  int el_idx = -1;
  for (unsigned iel = 0; iel < el_pt.size(); iel++) {
    if (el_sig[iel] && el_charge[iel]==1) {
      el_idx = static_cast<int>(iel);
      break;
    }
  }
  if (el_idx == -1) return -999.0;
  PtEtaPhiMVector el(el_pt[el_idx], el_eta[el_idx], el_phi[el_idx], 0.000511);
  PtEtaPhiMVector ph(photon_pt[el_photon_idx], photon_eta[el_photon_idx], 
                     photon_phi[el_photon_idx], 0.0);
  return (el+ph).M();
}

float get_el_phton_idmva(int el_photon_idx, RVec<float> photon_idmva) {
  if (el_photon_idx == -1) return -2.0;
  return photon_idmva[el_photon_idx];
}


//Return index of highest pT photon passing old signal criteria except electron veto and dR w.r.t. negative lepton
int nano_get_electron_photon_idx(RVec<float> Photon_pt, RVec<bool> Photon_isScEtaEB,
                                 RVec<bool> Photon_isScEtaEE, RVec<float> Photon_mvaID,
                                 RVec<float> Photon_eta, RVec<float> Photon_phi,
                                 RVec<bool> Electron_signal, RVec<int> Electron_charge,
                                 RVec<float> Electron_eta, RVec<float> Electron_phi) {

  int ph_idx = -1;
  float max_photon_pt = 0;
  
  for (unsigned iph = 0; iph < Photon_pt.size(); iph++) {
    if ((Photon_pt[iph]>15) &&
        ((Photon_isScEtaEB[iph] && Photon_mvaID[iph] > -0.4) ||
         (Photon_isScEtaEE[iph] && Photon_mvaID[iph] > -0.58))) {
      bool fail_dr = false;
      for (unsigned iel = 0; iel < Electron_signal.size(); iel++) {
        if ((Electron_charge[iel] == 1) && (Electron_signal[iel])) {
          if (delta_r(Electron_eta[iel], Photon_eta[iph], 
                      Electron_phi[iel], Photon_phi[iph]) < 0.3) {
            fail_dr = true;
          }
        }
      } //electron loop
      if (!fail_dr && (Photon_pt[iph] > max_photon_pt)) {
        max_photon_pt = Photon_pt[iph];
        ph_idx = static_cast<int>(iph);
      }
    }
  } //photon loop

  return ph_idx;
}

//conditions for passing TnP preselection for NanoAOD
bool pass_tnp(RVec<int> Electron_charge, RVec<bool> Electron_signal,
              RVec<float> Electron_eta, RVec<float> Electron_phi,
              RVec<int> TrigObj_id, RVec<int> TrigObj_filterBits,
              RVec<float> TrigObj_eta, RVec<float> TrigObj_phi,
              int LeadPhoton_index) {

  //flag indicating signal positron that passes trigger
  bool pass_tag = false;

  for (unsigned iel = 0; iel < Electron_charge.size(); iel++) {
    if ((Electron_charge[iel] == 1) && (Electron_signal[iel])) {
      for (unsigned itrig = 0; itrig < TrigObj_id.size(); itrig++) {
        if (TrigObj_id[itrig]==11) {
          if (delta_r(Electron_eta[iel], TrigObj_eta[itrig], 
                      Electron_phi[iel], TrigObj_phi[itrig])<0.2) {
            if ((TrigObj_filterBits[itrig] & 0x2) != 0) {
              pass_tag = true;
            }
          }
        }
      } //trigobj loop
    }
  } //electron loop

  return pass_tag && (LeadPhoton_index != -1);

}

//mvaID for lead photon
float LeadPhoton_mvaID(RVec<float> Photon_mvaID, int LeadPhoton_index) {
  if (LeadPhoton_index == -1) return -2.0;
  return Photon_mvaID[LeadPhoton_index];
}

RVec<bool> get_Electron_signal(RVec<float> Electron_pt, RVec<float> Electron_eta, 
                           RVec<float> Electron_deltaEtaSC, RVec<float> Electron_dxy, 
                           RVec<float> Electron_dz, RVec<bool> Electron_mvaFall17V2Iso_WPL) {
  RVec<bool> el_sig;
  for (unsigned iel = 0; iel < Electron_pt.size(); iel++) {
    el_sig.push_back((Electron_pt[iel] > 35) && 
                     ((Electron_eta[iel]+Electron_deltaEtaSC[iel])<2.5) && 
                     (Electron_dxy[iel] < 0.5) && 
                     (Electron_dz[iel] < 1.0) && 
                     (Electron_mvaFall17V2Iso_WPL[iel]));
  }
  return el_sig;
}

""")

nano_photon_idx_args = ['Photon_pt','Photon_isScEtaEB','Photon_isScEtaEE','Photon_mvaID','Photon_eta','Photon_phi','Electron_signal','Electron_charge','Electron_eta','Electron_phi']
nano_pass_tnp_args = ['Electron_charge','Electron_signal','Electron_eta','Electron_phi','TrigObj_id','TrigObj_filterBits','TrigObj_eta','TrigObj_phi','LeadPhoton_index']
nano_el_sig_args = ['Electron_pt','Electron_eta','Electron_deltaEtaSC','Electron_dxy','Electron_dz','Electron_mvaFall17V2Iso_WPL']

nano_photon_idx_string = 'nano_get_electron_photon_idx(Photon_pt,Photon_isScEtaEB,Photon_isScEtaEE,Photon_mvaID,Photon_eta,Photon_phi,Electron_signal,Electron_charge,Electron_eta,Electron_phi)'
nano_pass_tnp_string = 'pass_tnp(Electron_charge,Electron_signal,Electron_eta,Electron_phi,TrigObj_id,TrigObj_filterBits,TrigObj_eta,TrigObj_phi,LeadPhoton_index)'
#nano_el_signal_string = '(Electron_pt > 35) && ((Electron_eta+Electron_deltaEtaSC)<2.5) && (Electron_dxy < 0.5) && (Electron_dz < 1.0) && (Electron_mvaFall17V2Iso_WPL)'
nano_el_signal_string = 'get_Electron_signal(Electron_pt,Electron_eta,Electron_deltaEtaSC, Electron_dxy,Electron_dz,Electron_mvaFall17V2Iso_WPL)'

def debug_plots(): 
  '''temp debugging
  '''

  ROOT.EnableImplicitMT()

  #preselection is two signal leptons
  zg_filename = '/net/cms11/cms11r0/pico/NanoAODv9UCSB1/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_77.root'

  #preselection is tag_pt>35,tag_eta<2.17 not in gap, probe_pt>15, probe_eta<2.5, HIG-19-014 ID
  tnp_filename = '/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root'
  dy_filename = '/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8*.root'

  df_el = ROOT.RDataFrame('Events',dy_filename)
  print('DEBUG 0')
  df_el = df_el.Define('nel','nElectron')
  print('DEBUG 1')
  df_el = df_el.Define('Electron_signal',nano_el_signal_string)
  print('DEBUG 2')
  df_el = df_el.Define('LeadPhoton_index',nano_photon_idx_string)
  print('DEBUG 3')
  df_el = df_el.Define('Pass_tnp',nano_pass_tnp_string)
  print('DEBUG 4')
  df_el = df_el.Ftiler('Pass_tnp')
  print('DEBUG 5')
  df_el = df_el.Define('LeadPhoton_mvaID','LeadPhoton_mvaID(Photon_mvaID,LeadPhoton_index)')
  #df_el = df_el.Filter('ph_et>30')

  df_ph = ROOT.RDataFrame('tree',zg_filename)
  df_ph = df_ph.Filter('nmu==2&&trig_double_mu')
  df_ph = df_ph.Filter('mu_pt[0]>20&&mu_pt[1]>10')
  df_ph = df_ph.Define('lead_photon_idmva','photon_idmva')
  #df_ph = df_ph.Filter('photon_pt[0]>30')

  #make plots
  #el_hist_ptr = df_el.Histo1D(('el_hist','Electron;Photon IDMVA',30,-0.58,1.0),'ph_mva94XV2')
  el_hist_ptr = df_el.Histo1D(('el_hist','Electron;Photon IDMVA',30,-0.58,1.0),'LeadPhoton_mvaID','genWeight')
  ph_hist_ptr = df_ph.Histo1D(('ph_hist','Photon;Photon IDMVA',30,-0.58,1.0),'lead_photon_idmva','w_lumi')
  el_hist = el_hist_ptr.GetValue()
  ph_hist = ph_hist_ptr.GetValue()
  el_hist.Scale(1.0/el_hist.Integral())
  ph_hist.Scale(1.0/ph_hist.Integral())
  var_plot = RplPlot()
  var_plot.lumi_data = [(60,13)]
  var_plot.y_title = 'Events/bin'
  var_plot.plot_outline(el_hist)
  var_plot.plot_outline(ph_hist)
  var_plot.add_ratio('el_hist','ph_hist', True)
  var_plot.draw('plots/comp_phel.pdf')

if __name__ == '__main__':

  debug_plots()

  exit()

  argument_parser = ArgumentParser(prog='generate_photonid_dnn_syst',
  description='Script to validate DNN performance')
  argument_parser.add_argument('-s','--steps',default='make_plots,fit_plots,clean')
  argument_parser.add_argument('-d','--data_input_file')
  argument_parser.add_argument('-m','--mc_input_file')
  argument_parser.add_argument('-n','--nn_filename')
  argument_parser.add_argument('-r','--region',choices=['barrel','endcap'])

  args = argument_parser.parse_args()


