#!/usr/bin/env python3
'''
Script that generates slimmed n-tuples for photon Fall17v2 IDMVA corrections
'''
from argparse import ArgumentParser
from glob import glob
from skim_and_slim import write_ntuples
from subprocess import run as subprocessrun
import ROOT

#ROOT.gInterpreter.ProcessLine('#include "nn_model_phid2018_old.hpp"')
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

//returns preselection decision assuming leading muon is tag
//preselection is: nmu(pt>5 eta<2.4 looseId)=2
//                 tag pt>29 |dz|<1 |dxy|<0.5 sip3d<4 reliso03<0.35 trigger
//                 pair charge = 0
bool tag_muon_preselection(RVec<float> Muon_pt, RVec<float> Muon_eta, 
    RVec<float> Muon_phi, RVec<float> Muon_dxy, RVec<float> Muon_dz,
    RVec<float> Muon_sip3d, RVec<float> Muon_pfRelIso03_all, 
    RVec<bool> Muon_looseId, RVec<int> Muon_charge, RVec<float> TrigObj_pt, 
    RVec<float> TrigObj_eta, RVec<float> TrigObj_phi, RVec<int> TrigObj_id, 
    RVec<int> TrigObj_filterBits, unsigned tag_idx) {

  unsigned imu_sig = 0;
  int charge = 0;
  float tag_eta = 0.0, tag_phi = 0.0;

  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == tag_idx) {
        if (Muon_pt[imu] < 29) return false;
        if (Muon_dz[imu] > 1.0) return false;
        if (Muon_dxy[imu] > 0.5) return false;
        if (Muon_sip3d[imu] > 4.0) return false;
        if (Muon_pfRelIso03_all[imu] > 0.35) return false;
        tag_eta = Muon_eta[imu];
        tag_phi = Muon_phi[imu];
      }
      charge += Muon_charge[imu];
      imu_sig++;
    }
  }
  if (imu_sig != 2) return false;
  if (charge != 0) return false;
  //make sure tag triggered
  bool tag_triggered = false;
  for (unsigned itrig = 0; itrig<TrigObj_pt.size(); itrig++) {
    if (TrigObj_id[itrig]==13) {
      if ((TrigObj_filterBits[itrig] & 0xA) == 0xA) { //IsoMuXX
        if (delta_r(TrigObj_eta[itrig], TrigObj_phi[itrig], tag_eta, tag_phi)<0.1) {
          tag_triggered = true;
        }
      }
    }
  }
  return tag_triggered;
}

//get probe muon pt
float probe_muon_pt(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                    RVec<bool> Muon_looseId, unsigned probe_idx) {
  unsigned imu_sig = 0;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == probe_idx) {
        return Muon_pt[imu];
      }
      imu_sig++;
    }
  }
  return -999;
}

//get probe muon eta
float  probe_muon_eta(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                      RVec<bool> Muon_looseId, unsigned probe_idx) {
  unsigned imu_sig = 0;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == probe_idx) {
        return Muon_eta[imu];
      }
      imu_sig++;
    }
  }
  return -999;
}

//get probe muon reliso
float probe_muon_reliso(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                        RVec<bool> Muon_looseId, RVec<float> Muon_pfRelIso03_all,
                        unsigned probe_idx) {
  unsigned imu_sig = 0;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == probe_idx) {
        return Muon_pfRelIso03_all[imu];
      }
      imu_sig++;
    }
  }
  return -999;
}

//get probe muon dxy
float probe_muon_dxy(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                    RVec<bool> Muon_looseId, RVec<float> Muon_dxy,
                    unsigned probe_idx) {
  unsigned imu_sig = 0;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == probe_idx) {
        return Muon_dxy[imu];
      }
      imu_sig++;
    }
  }
  return -999;
}

//get probe muon dz
float probe_muon_dz(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                    RVec<bool> Muon_looseId, RVec<float> Muon_dz,
                    unsigned probe_idx) {
  unsigned imu_sig = 0;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == probe_idx) {
        return Muon_dz[imu];
      }
      imu_sig++;
    }
  }
  return -999;
}

//get probe muon sip3d
float probe_muon_sip3d(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                       RVec<bool> Muon_looseId, RVec<float> Muon_sip3d,
                       unsigned probe_idx) {
  unsigned imu_sig = 0;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
      //is signal muon
      if (imu_sig == probe_idx) {
        return Muon_sip3d[imu];
      }
      imu_sig++;
    }
  }
  return -999;
}

//get pair mass
float dimu_mass(RVec<float> Muon_pt, RVec<float> Muon_eta, 
                RVec<float> Muon_phi) {
                //, RVec<bool> Muon_looseId) {
  unsigned imu_sig = 0;
  std::vector<PtEtaPhiMVector> four_vecs;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    //if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4 && Muon_looseId[imu]) {
    if (Muon_pt[imu]>5.0 && fabs(Muon_eta[imu])<2.4) {
      //is signal muon
      four_vecs.push_back(PtEtaPhiMVector(Muon_pt[imu],Muon_eta[imu],Muon_phi[imu],0.106));
      imu_sig++;
    }
  }
  if (four_vecs.size() != 2) return -999;
  return (four_vecs[0]+four_vecs[1]).M();
}

//preselection requirements
bool tag_photon_preselection(RVec<float> Muon_pt, RVec<float> Muon_eta,
    RVec<float> Muon_phi, RVec<float> Muon_dxy, RVec<float> Muon_dz,
    RVec<float> Muon_sip3d, RVec<float> Muon_pfRelIso03_all,
    RVec<bool> Muon_looseId, RVec<int> Muon_charge, RVec<float> Photon_pt,
    RVec<float> Photon_eta, RVec<float>Photon_phi,) {
  return true;
}

""")

if __name__=='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='generate_photonid_skim',
      description='Generates a skim for deriving muon isolation correction from input file (-i) and outputs it to output file (-o) requiring the tag muon pass a trigger (-t)')
  argument_parser.add_argument('-i','--input_filename')
  argument_parser.add_argument('-o','--output_filename')
  args = argument_parser.parse_args()

  tag = 'lead'
  tag_idx_string = '0'
  probe_idx_string = '1'
  if (tag=='sublead'):
    tag_idx_string = '1'
    probe_idx_string = '0'

  ROOT.EnableImplicitMT()

  cuts = ['tag_photon_preselection(Muon_pt,Muon_eta,Muon_phi,Muon_dxy,Muon_dz,Muon_sip3d,Muon_pfRelIso03_all,Muon_looseId,Muon_charge,Photon_pt,Photon_eta,Photon_phi,)']
      
      tag_muon_preselection(Muon_pt, Muon_eta, Muon_phi, Muon_dxy, Muon_dz, Muon_sip3d, Muon_pfRelIso03_all, Muon_looseId, Muon_charge, TrigObj_pt, TrigObj_eta, TrigObj_phi, TrigObj_id, TrigObj_filterBits, '+tag_idx_string+')']
  defines = [('pair_mass','dimu_mass(Muon_pt, Muon_eta, Muon_phi, Muon_looseId)'),
             ('mu_pt','probe_muon_pt(Muon_pt, Muon_eta, Muon_looseId, '+probe_idx_string+')'),
             ('mu_eta','probe_muon_eta(Muon_pt, Muon_eta, Muon_looseId, '+probe_idx_string+')'),
             ('mu_reliso','probe_muon_reliso(Muon_pt, Muon_eta, Muon_looseId, Muon_pfRelIso03_all, '+probe_idx_string+')'),
             ('mu_dxy','probe_muon_dxy(Muon_pt, Muon_eta, Muon_looseId, Muon_dxy, '+probe_idx_string+')'),
             ('mu_dz','probe_muon_dz(Muon_pt, Muon_eta, Muon_looseId, Muon_dz, '+probe_idx_string+')'),
             ('mu_sip3d','probe_muon_sip3d(Muon_pt, Muon_eta, Muon_looseId, Muon_sip3d, '+probe_idx_string+')')]
  branches = ('pair_mass','mu_pt','mu_eta','mu_reliso','mu_dxy','mu_dz','mu_sip3d')

  input_files = glob(args.input_filename)

  write_ntuples(
      input_files,
      cuts,
      args.output_filename,
      defines,
      'tree',
      branches)

