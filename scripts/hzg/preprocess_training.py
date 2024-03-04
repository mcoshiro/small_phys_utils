#!/usr/bin/env python3
'''
Script that generates slimmed n-tuples for H->Zgamma MVA training
'''
import subprocess
import ROOT
import skim_and_slim

ROOT.gInterpreter.AddIncludePath('inc/')
ROOT.gInterpreter.ProcessLine('#include "hzg_mem.hpp"')
#ROOT.gSystem.AddLinkedLibs('-Llib -lSmallPhysUtils')
ROOT.gSystem.Load('libSmallPhysUtils.so')
ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;

//gSystem->AddIncludePath("-Iinc");
//gSystem->AddLinkedLibs("-Llib -lSmallPhysUtils");

float get_correct_cosTheta(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  TLorentzVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetPtEtaPhiM(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetPtEtaPhiM(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetPtEtaPhiM(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetPtEtaPhiM(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return get_cosTheta_alt2(p_lm, p_lp, p_gam);
}

float get_correct_costheta(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  TLorentzVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetPtEtaPhiM(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetPtEtaPhiM(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetPtEtaPhiM(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetPtEtaPhiM(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return get_costheta_alt2(p_lm, p_lp, p_gam);
}

float get_correct_phi(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  TLorentzVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetPtEtaPhiM(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetPtEtaPhiM(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetPtEtaPhiM(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetPtEtaPhiM(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetPtEtaPhiM(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return get_phi_alt2(p_lm, p_lp, p_gam);
}

float get_mela_disc(float cosTheta, float costheta, float phi) {
  return mela_discriminant(cosTheta, costheta, phi);
}

float get_dr(float eta1, float phi1, float eta2, float phi2) {
  const double PI = 3.1415;
  double dphi = fmod(fabs(phi2-phi1), 2.*PI);
  dphi = dphi>PI ? 2.*PI-dphi : dphi;
  double deta = fabs(eta1-eta2);
  return sqrt(deta*deta+dphi*dphi);
}

float get_max_dr(RVec<float> photon_eta, RVec<float> photon_phi, 
    RVec<float> el_eta, RVec<float> el_phi, RVec<float> mu_eta,
    RVec<float> mu_phi, RVec<int> ll_lepid, RVec<int> ll_i1,
    RVec<int> ll_i2) {
  float dr1, dr2;
  if (ll_lepid[0]==11) {
    dr1 = get_dr(photon_eta[0],photon_phi[0],el_eta[ll_i1[0]],el_phi[ll_i1[0]]);
    dr2 = get_dr(photon_eta[0],photon_phi[0],el_eta[ll_i2[0]],el_phi[ll_i2[0]]);
    return dr1 > dr2 ? dr1 : dr2;
  }
  dr1 = get_dr(photon_eta[0],photon_phi[0],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]]);
  dr2 = get_dr(photon_eta[0],photon_phi[0],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]]);
  return dr1 > dr2 ? dr1 : dr2;
}

float get_l1_rapidity(RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> mu_pt, RVec<float> mu_eta, RVec<int> ll_lepid, 
    RVec<int> ll_i1, RVec<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i1[0]] : el_eta[ll_i2[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i1[0]] : mu_eta[ll_i2[0]];
}

float get_l2_rapidity(RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> mu_pt, RVec<float> mu_eta, RVec<int> ll_lepid, 
    RVec<int> ll_i1, RVec<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i2[0]] : el_eta[ll_i1[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i2[0]] : mu_eta[ll_i1[0]];
}

float assoc_lepton_pdgid(RVec<int> ll_i1, RVec<int> ll_i2, RVec<int> ll_lepid,
    RVec<float> el_pt, RVec<bool> el_sig, RVec<float> mu_pt, RVec<bool> mu_sig) {
  float min_pt = 999;
  bool z_is_el = false;
  float pdgid = 0;
  unsigned int z_lep_i1 = static_cast<unsigned int>(ll_i1[0]);
  unsigned int z_lep_i2 = static_cast<unsigned int>(ll_i2[0]);
  if (ll_lepid[0]==11)
    z_is_el = true;
  for (unsigned iel = 0; iel < el_pt.size(); iel++) {
    if (!el_sig[iel]) continue;
    if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
    if (el_pt[iel]<min_pt) {
      min_pt = el_pt[iel];
      pdgid = 11;
    }
  }
  for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
    if (mu_sig[imu]) continue;
    if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
    if (mu_pt[imu]<min_pt) {
      min_pt = mu_pt[imu];
      pdgid = 13;
    }
  }
  return pdgid;
}

float assoc_lepton_quality(RVec<int> ll_i1, RVec<int> ll_i2, RVec<int> ll_lepid,
    RVec<float> el_pt, RVec<bool> el_sig, RVec<float> mu_pt, RVec<bool> mu_sig,
    RVec<bool> el_id80, RVec<bool> el_id90, RVec<bool> el_idLoose,
    RVec<bool> mu_tightid, RVec<bool> mu_mediumid, RVec<bool> mu_highptid, RVec<bool> mu_id) {
    float min_pt = 999;
    bool z_is_el = false;
    float quality = 0;
    unsigned int z_lep_i1 = static_cast<unsigned int>(ll_i1[0]);
    unsigned int z_lep_i2 = static_cast<unsigned int>(ll_i2[0]);
    if (ll_lepid[0]==11)
      z_is_el = true;
    for (unsigned iel = 0; iel < el_pt.size(); iel++) {
      if (!el_sig[iel]) continue;
      if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
      if (el_pt[iel]<min_pt) {
        min_pt = el_pt[iel];
        if (el_id80[iel]) quality = 3;
        else if (el_id90[iel]) quality = 2;
        else if (el_idLoose[iel]) quality = 1;
        else quality = 0;
      }
    }
    for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
      if (!mu_sig[imu]) continue;
      if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
      if (mu_pt[imu]<min_pt) {
        min_pt = mu_pt[imu];
        if (mu_tightid[imu]) quality = 3;
        else if (mu_mediumid[imu] 
            || (mu_highptid[imu] && mu_pt[imu]>200)) quality = 2;
        else if (mu_id[imu]) quality = 1;
        else quality = 0;
      }
    }
    return quality;
}

float assoc_lep_miniso(RVec<int> ll_i1, RVec<int> ll_i2, RVec<int> ll_lepid,
    RVec<float> el_pt, RVec<bool> el_sig, RVec<float> mu_pt, RVec<bool> mu_sig,
    RVec<float> el_miniso, RVec<float> mu_miniso) {
  float min_pt = 999;
  bool z_is_el = false;
  float miniso = 0;
  unsigned int z_lep_i1 = static_cast<unsigned int>(ll_i1[0]);
  unsigned int z_lep_i2 = static_cast<unsigned int>(ll_i2[0]);
  if (ll_lepid[0]==11)
    z_is_el = true;
  for (unsigned iel = 0; iel < el_pt.size(); iel++) {
    if (!el_sig[iel]) continue;
    if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
    if (el_pt[iel]<min_pt) {
      min_pt = el_pt[iel];
      miniso = el_miniso[iel];
    }
  }
  for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
    if (!mu_sig[imu]) continue;
    if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
    if (mu_pt[imu]<min_pt) {
      min_pt = mu_pt[imu];
      miniso = mu_miniso[imu];
    }
  }
  return miniso;
}

float assoc_lep_pt(RVec<int> ll_i1, RVec<int> ll_i2, RVec<int> ll_lepid,
  RVec<float> el_pt, RVec<bool> el_sig, RVec<float> mu_pt, RVec<bool> mu_sig) {
  float min_pt = 999;
  bool z_is_el = false;
  unsigned int z_lep_i1 = static_cast<unsigned int>(ll_i1[0]);
  unsigned int z_lep_i2 = static_cast<unsigned int>(ll_i2[0]);
  if (ll_lepid[0]==11)
    z_is_el = true;
  for (unsigned iel = 0; iel < el_pt.size(); iel++) {
    if (!el_sig[iel]) continue;
    if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
    if (el_pt[iel]<min_pt) {
      min_pt = el_pt[iel];
    }
  }
  for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
    if (!mu_sig[imu]) continue;
    if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
    if (mu_pt[imu]<min_pt) {
      min_pt = mu_pt[imu];
    }
  }
  return min_pt;
}

float min_lepton_quality(RVec<bool> el_sig, RVec<bool> mu_sig, RVec<float> mu_pt,
    RVec<bool> el_id80, RVec<bool> el_id90, RVec<bool> el_idLoose,
    RVec<bool> mu_tightid, RVec<bool> mu_mediumid, RVec<bool> mu_highptid, RVec<bool> mu_id) {
    float quality = 3;
    for (unsigned iel = 0; iel < el_sig.size(); iel++) {
      if (!el_sig[iel]) continue;
      float this_quality = 0;
      if (el_id80[iel]) this_quality = 3;
      else if (el_id90[iel]) this_quality = 2;
      else if (el_idLoose[iel]) this_quality = 1;
      if (this_quality < quality) quality = this_quality;
    }
    for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
      if (!mu_sig[imu]) continue;
      float this_quality = 0;
      if (mu_tightid[imu]) this_quality = 3;
      else if (mu_mediumid[imu] 
          || (mu_highptid[imu] && mu_pt[imu]>200)) this_quality = 2;
      else if (mu_id[imu]) this_quality = 1;
      if (this_quality < quality) quality = this_quality;
    }
    return quality;
}

float max_lep_miniso(RVec<bool> el_sig, RVec<bool> mu_sig,
    RVec<float> el_miniso, RVec<float> mu_miniso) {
  float miniso = 0;
  for (unsigned iel = 0; iel < el_sig.size(); iel++) {
    if (!el_sig[iel]) continue;
    if (el_miniso[iel] > miniso)
      miniso = el_miniso[iel];
  }
  for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
    if (!mu_sig[imu]) continue;
    if (mu_miniso[imu] > miniso)
      miniso = mu_miniso[imu];
  }
  return miniso;
}

float min_lep_pt(RVec<float> el_pt, RVec<bool> el_sig, RVec<float> mu_pt, RVec<bool> mu_sig) {
  float min_pt = 999;
  for (unsigned iel = 0; iel < el_pt.size(); iel++) {
    if (!el_sig[iel]) continue;
    if (el_pt[iel]<min_pt) {
      min_pt = el_pt[iel];
    }
  }
  for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
    if (!mu_sig[imu]) continue;
    if (mu_pt[imu]<min_pt) {
      min_pt = mu_pt[imu];
    }
  }
  return min_pt;
}

float closest_top_mass(RVec<bool> jet_isgood, RVec<float> jet_deepflav, RVec<float> jet_pt,
  RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_m) {
  //first find 2 highest b-tag score jets
  float lead_btag_score = -1;
  float subl_btag_score = -1;
  unsigned bidx1 = -1;
  unsigned bidx2 = -1;
  TLorentzVector j1, j2, w;
  for (unsigned ijet = 0; ijet < jet_pt.size(); ijet++) {
    if (!jet_isgood[ijet]) continue;
    float this_btag_score = jet_deepflav[ijet];
    if (this_btag_score > lead_btag_score) {
      subl_btag_score = lead_btag_score;
      lead_btag_score = this_btag_score;
      bidx2 = bidx1;
      bidx1 = ijet;
    }
    else if (this_btag_score > subl_btag_score) {
      subl_btag_score = this_btag_score;
      bidx2 = ijet;
    }
  }
  //from other jets, find pair closest to W
  float closest_w = -999;
  for (unsigned ijet1 = 0; ijet1 < jet_pt.size(); ijet1++) {
    if (!jet_isgood[ijet1]) continue;
    if (ijet1==bidx1 || ijet1==bidx2) continue;
    j1.SetPtEtaPhiM(jet_pt[ijet1], jet_eta[ijet1], 
                    jet_phi[ijet1], jet_m[ijet1]);
    for (unsigned ijet2 = ijet1+1; ijet2 < jet_pt.size(); ijet2++) {
      if (!jet_isgood[ijet2]) continue;
      if (ijet2==bidx1 || ijet2==bidx2) continue;
      j2.SetPtEtaPhiM(jet_pt[ijet2], jet_eta[ijet2], 
                      jet_phi[ijet2], jet_m[ijet2]);
      float this_mass = (j1+j2).M();
      if (fabs(this_mass-80.4) < fabs(closest_w-80.4)) {
        closest_w = this_mass;
        w = j1+j2;
      }
    }
  }
  //refit W??
  //combine W with one of the b-tagged jets and take value closest to top
  j1.SetPtEtaPhiM(jet_pt[bidx1], jet_eta[bidx1], 
                  jet_phi[bidx1], jet_m[bidx1]);
  j2.SetPtEtaPhiM(jet_pt[bidx2], jet_eta[bidx2], 
                  jet_phi[bidx2], jet_m[bidx2]);
  float top1_m = (j1+w).M();
  float top2_m = (j2+w).M();
  if (fabs(top1_m-172.7)<fabs(top2_m-172.7)) return top1_m;
  return top2_m;
}

float assoc_lep_mt(int nlep, RVec<int> ll_i1, RVec<int> ll_i2, RVec<int> ll_lepid, RVec<bool> el_sig,
    RVec<float> el_pt, RVec<float> el_phi, RVec<bool> mu_sig, RVec<float> mu_pt, RVec<float> mu_phi, float met, float met_phi) {
  if (nlep<3) return 0;
  float max_pt = 0;
  float max_phi = 0;
  bool z_is_el = false;
  unsigned int z_lep_i1 = static_cast<unsigned int>(ll_i1[0]);
  unsigned int z_lep_i2 = static_cast<unsigned int>(ll_i2[0]);
  if (ll_lepid[0]==11)
    z_is_el = true;
  for (unsigned iel = 0; iel < el_pt.size(); iel++) {
    if (!el_sig[iel]) continue;
    if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
    if (el_sig[iel] && el_pt[iel]>max_pt) {
      max_pt = el_pt[iel];
      max_phi = el_phi[iel];
    }
  }
  for (unsigned imu = 0; imu < mu_pt.size(); imu++) {
    if (!mu_sig[imu]) continue;
    if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
    if (mu_sig[imu] && mu_pt[imu]>max_pt) {
      max_pt = mu_pt[imu];
      max_phi = mu_phi[imu];
    }
  }
  return sqrt(2.0*met*max_pt*(1.0-cos(met_phi-max_phi)));
}

float dphi_h_met(RVec<float> llphoton_phi, float met_phi) {
  if (llphoton_phi.size()==0) return 0;
  float dphi = abs(llphoton_phi[0]-met_phi);
  if (fabs(llphoton_phi[0]+2.0*3.1415-met_phi)<dphi) 
    dphi = fabs(llphoton_phi[0]+2.0*3.1415-met_phi);
  if (fabs(llphoton_phi[0]-2.0*3.1415-met_phi)<dphi) 
    dphi = fabs(llphoton_phi[0]-2.0*3.1415-met_phi);
  return dphi;
}

float h_pz(RVec<float> llphoton_pt, RVec<float> llphoton_eta) {
  return llphoton_pt[0]*sinh(llphoton_eta[0]);
}


template<class coords1, class coords2>
double cosangle(ROOT::Math::LorentzVector<coords1> lv1, ROOT::Math::LorentzVector<coords2> lv2) {
  ROOT::Math::XYZVector v1 = lv1.Vect();
  ROOT::Math::XYZVector v2 = lv2.Vect();
  return v1.Dot(v2)/v1.R()/v2.R();
}

template<class coords1, class coords2, class coords3>
double get_phi(ROOT::Math::LorentzVector<coords1> lv, ROOT::Math::LorentzVector<coords2> lv_zref, ROOT::Math::LorentzVector<coords3> lv_xref) {
  ROOT::Math::XYZVector v_target = lv.Vect();
  ROOT::Math::XYZVector v_xref = lv_xref.Vect();
  ROOT::Math::XYZVector z_unit = lv_zref.Vect()/lv_zref.Vect().R();
  ROOT::Math::XYZVector y_unit = lv_zref.Vect().Cross(lv_xref.Vect());
  y_unit = y_unit/y_unit.R();
  ROOT::Math::XYZVector x_unit = y_unit.Cross(z_unit);
  v_target = v_target - v_target.Dot(z_unit)*z_unit;
  float abs_x_angle = TMath::ACos(v_target.Dot(x_unit)/v_target.R());
  float abs_y_angle = TMath::ACos(v_target.Dot(y_unit)/v_target.R());
  if (abs_y_angle < TMath::Pi()/2.0) return abs_x_angle;
  return -1.0*abs_x_angle;
}

//returns 4-vector of +z gluon in lab frame assuming the usual convention
template<class coords1, class coords2, class coords3>
ROOT::Math::PxPyPzEVector poszgluon(ROOT::Math::LorentzVector<coords1> const & p_neglep, ROOT::Math::LorentzVector<coords2> const & p_poslep, ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //set up boost to "zgtransverse" frame
  ROOT::Math::PtEtaPhiMVector p_zph = (p_poslep+p_neglep+p_photon);
  ROOT::Math::XYZVector boost_vec_lab_to_zgtrans = p_zph.BoostToCM();
  boost_vec_lab_to_zgtrans.SetZ(0.0);
  ROOT::Math::Boost boost_lab_to_zgtrans(boost_vec_lab_to_zgtrans);
  ROOT::Math::Boost boost_zgtrans_to_lab(boost_lab_to_zgtrans.Inverse());

  //calculate +z gluon in "zgtransverse" frame and boost back
  ROOT::Math::PtEtaPhiMVector p_zph_zgtrans(boost_lab_to_zgtrans(p_zph));
  double partz_zgtrans = p_zph_zgtrans.Pz()+p_zph_zgtrans.E();
  ROOT::Math::PxPyPzEVector p_g1(0.0,0.0,partz_zgtrans,partz_zgtrans); 
  return boost_zgtrans_to_lab(p_g1);
}

//returns cosTheta of Z in lly CM frame with +z gluon as z reference and lab +x as X reference
template<class coords1, class coords2, class coords3>
float coscaptheta(ROOT::Math::LorentzVector<coords1> const & p_neglep, ROOT::Math::LorentzVector<coords2> const & p_poslep, ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::PtEtaPhiMVector p_zph = (p_poslep+p_neglep+p_photon);
  ROOT::Math::Boost boost_lab_to_zgcm(p_zph.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PxPyPzEVector p_g1_cm = boost_lab_to_zgcm(poszgluon(p_neglep,p_poslep,p_photon));
  ROOT::Math::PtEtaPhiMVector p_z_cm = boost_lab_to_zgcm(p_z);
  return cosangle(p_z_cm, p_g1_cm);
}

//returns Phi of Z in lly CM frame with +z gluon as z reference and lab +x as X reference
template<class coords1, class coords2, class coords3>
float capphi(ROOT::Math::LorentzVector<coords1> const & p_neglep, ROOT::Math::LorentzVector<coords2> const & p_poslep, ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PxPyPzEVector p_xlab(1.0,0.0,0.0,0.0); 
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::PtEtaPhiMVector p_zph = (p_poslep+p_neglep+p_photon);
  ROOT::Math::Boost boost_lab_to_zgcm(p_zph.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PxPyPzEVector p_g1_cm = boost_lab_to_zgcm(poszgluon(p_neglep,p_poslep,p_photon));
  ROOT::Math::PxPyPzEVector p_xlab_cm = boost_lab_to_zgcm(p_xlab);
  ROOT::Math::PtEtaPhiMVector p_z_cm = boost_lab_to_zgcm(p_z);
  return get_phi(p_z_cm, p_g1_cm, p_xlab_cm);
}

//returns costheta of positive lepton in ll CM frame with photon as z reference and +z gluon as x reference
template<class coords1, class coords2, class coords3>
float coslowertheta(ROOT::Math::LorentzVector<coords1> const & p_neglep, ROOT::Math::LorentzVector<coords2> const & p_poslep, ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::Boost boost_lab_to_zcm(p_z.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PtEtaPhiMVector p_poslep_z = boost_lab_to_zcm(p_poslep);
  ROOT::Math::PtEtaPhiMVector p_photon_z = boost_lab_to_zcm(p_photon);
  return cosangle(p_poslep_z, p_photon_z);
}

//returns phi of positive lepton in ll CM frame with photon as z reference and +z gluon as x reference
template<class coords1, class coords2, class coords3>
float lowerphi(ROOT::Math::LorentzVector<coords1> const & p_neglep, ROOT::Math::LorentzVector<coords2> const & p_poslep, ROOT::Math::LorentzVector<coords3> const & p_photon) {
  //initialize necessary 4-vectors and boosts
  ROOT::Math::PtEtaPhiMVector p_z = (p_poslep+p_neglep);
  ROOT::Math::Boost boost_lab_to_zcm(p_z.BoostToCM());

  //boost to CM frame and calculate
  ROOT::Math::PxPyPzEVector p_g1_z = boost_lab_to_zcm(poszgluon(p_neglep,p_poslep,p_photon));
  ROOT::Math::PtEtaPhiMVector p_poslep_z = boost_lab_to_zcm(p_poslep);
  ROOT::Math::PtEtaPhiMVector p_photon_z = boost_lab_to_zcm(p_photon);
  return get_phi(p_poslep_z, p_photon_z, p_g1_z);
}

float get_new_cosTheta(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  ROOT::Math::PtEtaPhiMVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetCoordinates(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return coscaptheta(p_lm, p_lp, p_gam);
}

float get_new_Phi(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  ROOT::Math::PtEtaPhiMVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetCoordinates(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return capphi(p_lm, p_lp, p_gam);
}

float get_new_costheta(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  ROOT::Math::PtEtaPhiMVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetCoordinates(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return coslowertheta(p_lm, p_lp, p_gam);
}

float get_new_phi(RVec<float> photon_pt, RVec<float> photon_eta, 
    RVec<float> photon_phi, RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> el_phi, RVec<int> el_charge, RVec<float> mu_pt, 
    RVec<float> mu_eta, RVec<float> mu_phi, RVec<int> mu_charge, 
    RVec<int> ll_lepid, RVec<int> ll_i1, RVec<int> ll_i2) {
  ROOT::Math::PtEtaPhiMVector p_lm, p_lp, p_gam;
  if (ll_lepid[0]==11) {
    if (el_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0);
    }
  }
  else {
    if (mu_charge[ll_i1[0]]==-1) {
      p_lm.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lp.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
    else {
      p_lp.SetCoordinates(mu_pt[ll_i1[0]],mu_eta[ll_i1[0]],mu_phi[ll_i1[0]],0);
      p_lm.SetCoordinates(mu_pt[ll_i2[0]],mu_eta[ll_i2[0]],mu_phi[ll_i2[0]],0);
    }
  }
  p_gam.SetCoordinates(photon_pt[0],photon_eta[0],photon_phi[0],0);
  return lowerphi(p_lm, p_lp, p_gam);
}

bool trig_lep_cuts(float year, int nel, int nmu, bool trig_single_el, bool trig_single_mu,
    bool trig_double_el, bool trig_double_mu, RVec<float> el_pt, 
    RVec<float> mu_pt) {
  float el_cut = 35;
  float mu_cut = 25;
  if (year+0.1 < 2017)
    el_cut = 30;
  else if (year+0.1 > 2017 && year+0.1 < 2018)
    mu_cut = 28;
  if (nel>=1)
    if (trig_single_el && el_pt[0]>el_cut) return true;
  if (nel>=2)
    if (trig_double_el && el_pt[0]>25 && el_pt[1]>15) return true;
  if (nmu>=1)
    if (trig_single_mu && mu_pt[0]>mu_cut) return true;
  if (nmu>=2)
    if (trig_double_mu && mu_pt[0]>20 && mu_pt[1]>10) return true;
  return false;
}

float get_lead_jet_pt(RVec<float> jet_pt, RVec<bool> jet_isgood) {
  float max_pt = 0;
  for (unsigned ijet = 0; ijet < jet_pt.size(); ijet++) {
    if (jet_isgood[ijet])
      if (jet_pt[ijet] > max_pt)
        max_pt = jet_pt[ijet];
  }
  return max_pt;
}

float get_sublead_jet_pt(RVec<float> jet_pt, RVec<bool> jet_isgood) {
  float lead_pt = 0;
  float subl_pt = -1;
  for (unsigned ijet = 0; ijet < jet_pt.size(); ijet++) {
    if (jet_isgood[ijet]) {
      if (jet_pt[ijet] > lead_pt) {
        subl_pt = lead_pt;
        lead_pt = jet_pt[ijet];
      }
      else if (jet_pt[ijet] > subl_pt) {
        subl_pt = jet_pt[ijet];
      }
    }
  }
  return subl_pt;
}

""")

def produce_vbf_samples():
  ROOT.EnableImplicitMT()
  #note trigger and lepton pt cuts are year specific, see below
  cuts = ['use_event',
          'nlep==2&&nphoton>=1&&njet>=2&&nbm==0',
          'll_m[0]>81&&ll_m[0]<101&&ll_charge[0]==0',
          '(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&((llphoton_m[0]+ll_m[0])>185)',
          'photon_id80[0]',
          'llphoton_m[0]>120&&llphoton_m[0]<130']
  #note lumi is year specific, see below
  defines = [('photon_mva','photon_idmva[0]'),
             ('min_dR','photon_drmin[0]'),
             ('max_dR','get_max_dr(photon_eta,photon_phi,el_eta,el_phi,mu_eta,mu_phi,ll_lepid,ll_i1,ll_i2)'),
             ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
             ('cosTheta','llphoton_cosTheta[0]'),
             ('costheta','llphoton_costheta[0]'),
             ('phi','llphoton_phi[0]'),
             ('photon_res','static_cast<float>(photon_energyErr[0]/(photon_pt[0]*TMath::CosH(photon_eta[0])))'),
             ('photon_rapidity','photon_eta[0]'),
             ('l1_rapidity','get_l1_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('l2_rapidity','get_l2_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             #above are the original kinematic mva variables
             ('detajj','dijet_deta'),
             ('dphizgjj','llphoton_dijet_dphi[0]'),
             ('zgjj_balance','llphoton_dijet_balance[0]'),
             ('ptt','llphoton_pTt2[0]'), #N.B. not same as run 2 pTt
             ('dphijj','dijet_dphi'), 
             ('zeppenfeld','photon_zeppenfeld[0]'), 
             ('ptj2','get_sublead_jet_pt(jet_pt,jet_isgood)'), 
             ('ptj1','get_lead_jet_pt(jet_pt,jet_isgood)'), 
             ('drgj','photon_jet_mindr[0]'), 
             #the above (except pTt) are run 2 dijet variables
             ('mllg','llphoton_m[0]')]
  branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta',
              'phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity',
              'detajj','dphizgjj','zgjj_balance','ptt','dphijj','zeppenfeld',
              'ptj2','ptj1','drgj','mllg','w_lumi_year')
  #make n-tuples for 3 years, then merge
  years = ['2016APV','2016','2017','2018']
  lumis = ['19.51','16.80','41.48','59.83']
  year_float = ['2016.0','2016.5','2017.0','2018.0']
  base_name = 'dijet'
  for iyear in range(len(years)):
    cuts = cuts + ['trig_lep_cuts('+year_float[iyear]+', nel, nmu, trig_single_el, trig_single_mu,trig_double_el, trig_double_mu, el_pt, mu_pt)']
    defines_year = defines + [('w_lumi_year','float w_lumi_year=w_lumi*'+lumis[iyear]+'; return w_lumi_year;')]
    skim_and_slim.write_ntuples(
        ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*GluGluHToZG_ZToLL_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*VBFHToZG_ZToLL_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*WminusH_HToZG_WToAll_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*WplusH_HToZG_WToAll_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ZH_HToZG_ZToAll_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ttHToZG_M125*'],
        cuts,
        base_name+'_sig_'+years[iyear]+'.root',
        defines_year,
        'tree',
        branches)
    skim_and_slim.write_ntuples(
        ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ZGToLLG_01J_5f_lowMLL*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*TTTo2L2Nu*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*TGJets*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ZGamma2JToGamma2L2J_EWK*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*_WW_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*_WZ_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*_ZZ_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*WZG_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ZZG_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ttWJets_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*ttZJets_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/'+years[iyear]+'/mc/skim_llg/*WGToLNuG*.root'],
        cuts,
        base_name+'_bak_'+years[iyear]+'.root',
        defines_year,
        'tree',
        branches)
  #merge n-tuples and shuffle
  print('Merging and shuffling output files')
  for sigbak in ['sig','bak']:
    merged_tree = ROOT.TChain("tree", "tree")
    for year in years:
      merged_tree.Add('ntuples/'+base_name+'_'+sigbak+'_'+year+'.root')
    out_file = ROOT.TFile('ntuples/'+base_name+'_'+sigbak+'.root', "recreate", "", 0)
    merged_tree.Merge(out_file,0,'fast keep')
    out_file.Close()
    subprocess.run(('/data1/jbkim/Linux/el7_v1/bin/python3.9 scripts/shuffle_tree.py ntuples/'+base_name+'_'+sigbak+'.root tree').split())
  #clean up
  for sigbak in ['sig','bak']:
    for year in years:
      subprocess.run(('rm ntuples/'+base_name+'_'+sigbak+'_'+year+'.root').split())
    subprocess.run(('rm ntuples/'+base_name+'_'+sigbak+'.root').split())

def produce_ggh_samples():
  ROOT.EnableImplicitMT()
  #cuts = ['HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
  #    'nlep>=2&&nphoton>=1',
  #    'll_m.size()>0&&llphoton_m.size()>0',
  #    'stitch_dy||(type/1000!=6)',
  #    '(ll_m[0]>50)&&(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&((llphoton_m[0]+ll_m[0])>185)&&(photon_drmin[0]>0.4)',
  #    '(ll_lepid[0]==11&&el_pt[ll_i1[0]]>25&&el_pt[ll_i2[0]]>15)||(ll_lepid[0]==13&&mu_pt[ll_i1[0]]>20&&mu_pt[ll_i2[0]]>10)',
  #    #'photon_idmva[0]>0.5',
  #    #'photon_id80[0]',
  #    #'(ll_m[0]>80&&ll_m[0]<100)',
  #    'llphoton_m[0]>120&&llphoton_m[0]<130']
  #    #'llphoton_m[0]>100&&llphoton_m[0]<160']
  #    #'nlep==2&&nbdfm==0&&met<95&&njet>=2&&dijet_m>750']
  #    #'photon_pt[0]<35']
  cuts = ['(trig_single_el||trig_double_el||trig_single_mu||trig_double_mu)',
      'nlep>=2&&nphoton>=1',
      '(type!=1600&&use_event)||(type==1600&&photon_pflavor[0]!=1)||(type==16000)',
      '(ll_m[0]>50)&&(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&((llphoton_m[0]+ll_m[0])>185)',
      '(ll_lepid[0]==11&&el_pt[0]>25&&el_pt[1]>15)||(ll_lepid[0]==13&&mu_pt[0]>20&&mu_pt[1]>10)',
      'photon_id80[0]',
      'llphoton_m[0]>120&&llphoton_m[0]<130']
  defines = [('photon_mva','photon_idmva[0]'),
             ('min_dR','photon_drmin[0]'),
             ('max_dR','get_max_dr(photon_eta,photon_phi,el_eta,el_phi,mu_eta,mu_phi,ll_lepid,ll_i1,ll_i2)'),
             ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
             ('cosTheta','llphoton_cosTheta[0]'),
             ('costheta','llphoton_costheta[0]'),
             ('phi','llphoton_phi[0]'),
             ('photon_res','static_cast<float>(photon_pterr[0]/(photon_pt[0]*TMath::CosH(photon_eta[0])))'),
             ('photon_rapidity','photon_eta[0]'),
             ('l1_rapidity','get_l1_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('l2_rapidity','get_l2_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             #above are the original 11
             #('pt_llg','llphoton_pt[0]'),
             #('eta_llg','llphoton_eta[0]'),
             #('phi_llg','llphoton_phi[0]'),
             #('cosTheta','get_new_cosTheta(photon_pt,photon_eta,photon_phi,el_pt,el_eta,el_phi,el_charge,mu_pt,mu_eta,mu_phi,mu_charge,ll_lepid,ll_i1,ll_i2)'),
             #('Phi','get_new_Phi(photon_pt,photon_eta,photon_phi,el_pt,el_eta,el_phi,el_charge,mu_pt,mu_eta,mu_phi,mu_charge,ll_lepid,ll_i1,ll_i2)'),
             #('costheta','get_new_costheta(photon_pt,photon_eta,photon_phi,el_pt,el_eta,el_phi,el_charge,mu_pt,mu_eta,mu_phi,mu_charge,ll_lepid,ll_i1,ll_i2)'),
             #('phi','get_new_phi(photon_pt,photon_eta,photon_phi,el_pt,el_eta,el_phi,el_charge,mu_pt,mu_eta,mu_phi,mu_charge,ll_lepid,ll_i1,ll_i2)'),
             #('photon_ptransverse','photon_pt[0]'),
             #('mll','ll_m[0]'),
             ('mllg','llphoton_m[0]')]
             #('higgsdr','llphoton_dr[0]'),
             #('z_pt_mass','ll_pt[0]/llphoton_m[0]'),
             #('mela','get_mela_disc(cosTheta,costheta,phi)'),
             #('photon_pt_mass','photon_pt[0]/llphoton_m[0]'),
             #('ph_irel','photon_reliso[0]'),
             #('ph_r9','photon_r9[0]'),
             #('ph_sieie','photon_sieie[0]'),
             #('ph_hoe','photon_hoe[0]')]
  branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta',
      'phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity',
      'mllg','w_lumi_year')
  #defines = [('photon_mva','photon_idmva[0]'),
  #           ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
  #           ('min_dR','photon_drmin[0]'),
  #           ('min_lepton_quality','min_lepton_quality(el_sig,mu_sig,mu_pt,el_id80,el_id90,el_idLoose,mu_tightid,mu_mediumid,mu_highptid,mu_id)'),
  #           ('assoc_lepton_quality','assoc_lepton_quality(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig,el_id80,el_id90,el_idLoose,mu_tightid,mu_mediumid,mu_highptid,mu_id)'),

  #           ('max_lepton_miniso','max_lep_miniso(el_sig,mu_sig,el_miniso,mu_miniso)'),
  #           ('assoc_lepton_miniso','assoc_lep_miniso(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig,el_miniso,mu_miniso)'),

  #           ('min_lepton_pt','min_lep_pt(el_pt,el_sig,mu_pt,mu_sig)'),
  #           ('assoc_lepton_pt','assoc_lep_pt(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig)'),
  #           ('assoc_lepton_pdgid','assoc_lepton_pdgid(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig)'),
  #           ('njet_f','static_cast<float>(njet)'),
  #           ('nbdfl_f','static_cast<float>(nbdfl)'),
  #           ('nbdfm_f','static_cast<float>(nbdfm)'),
  #           ('nbdft_f','static_cast<float>(nbdft)'),
  #           ('nlep_f','static_cast<float>(nlep)'),
  #           ('mllg','llphoton_m[0]')]
  ##also: ht, met
  #branches = ('w_lumi_year','mllg','photon_mva','pt_mass','min_dR',
  #    'min_lepton_quality','assoc_lepton_quality','max_lepton_miniso',
  #    'assoc_lepton_miniso','min_lepton_pt','assoc_lepton_pt',
  #    'assoc_lepton_pdgid','ht','met','njet_f','nbdfl_f','nbdfm_f','nbdft_f','nlep_f')
  #defines = [('photon_mva','photon_idmva[0]'),
  #           ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
  #           ('min_dR','photon_drmin[0]'),
  #           ('njet_f','static_cast<float>(njet)'),
  #           ('nbdfl_f','static_cast<float>(nbdfl)'),
  #           ('nbdfm_f','static_cast<float>(nbdfm)'),
  #           ('nbdft_f','static_cast<float>(nbdft)'),
  #           ('top_mass','closest_top_mass(jet_isgood,jet_deepflav,jet_pt,jet_eta,jet_phi,jet_m)'),
  #           ('mllg','llphoton_m[0]')]
  #branches = ('w_lumi_year','mllg','photon_mva','pt_mass','min_dR',
  #    'ht','met','njet_f','nbdfl_f','nbdfm_f','nbdft_f','top_mass')
  #defines = [('photon_mva','photon_idmva[0]'),
  #           ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
  #           ('min_dR','photon_drmin[0]'),
  #           ('min_lepton_quality','min_lepton_quality(el_sig,mu_sig,mu_pt,el_id80,el_id90,el_idLoose,mu_tightid,mu_mediumid,mu_highptid,mu_id)'),
  #           ('max_lepton_miniso','max_lep_miniso(el_sig,mu_sig,el_miniso,mu_miniso)'),
  #           ('min_lepton_pt','min_lep_pt(el_pt,el_sig,mu_pt,mu_sig)'),
  #           ('assoc_lepton_quality','assoc_lepton_quality(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig,el_id80,el_id90,el_idLoose,mu_tightid,mu_mediumid,mu_highptid,mu_id)'),
  #           ('assoc_lepton_miniso','assoc_lep_miniso(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig,el_miniso,mu_miniso)'),
  #           ('assoc_lepton_pt','assoc_lep_pt(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig)'),
  #           ('assoc_lepton_pdgid','assoc_lepton_pdgid(ll_i1,ll_i2,ll_lepid,el_pt,el_sig,mu_pt,mu_sig)'),
  #           ('assoc_lepton_mt','assoc_lep_mt(nlep,ll_i1,ll_i2,ll_lepid,el_sig,el_pt,el_phi,mu_sig,mu_pt,mu_phi,met,met_phi)'),
  #           ('njet_f','static_cast<float>(njet)'),
  #           ('nlep_f','static_cast<float>(nlep)'),
  #           ('mllg','llphoton_m[0]')]
  #branches = ('w_lumi_year','mllg','photon_mva','pt_mass','min_dR',
  #    'min_lepton_quality','max_lepton_miniso','min_lepton_pt',
  #    'assoc_lepton_quality','assoc_lepton_miniso','assoc_lepton_pt',
  #    'assoc_lepton_pdgid','assoc_lepton_mt','ht','met','njet_f','nlep_f')
  #defines = [('photon_mva','photon_idmva[0]'),
  #           ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
  #           ('min_dR','photon_drmin[0]'),
  #           ('dphi_h_met','dphi_h_met(llphoton_phi, met_phi)'),
  #           ('njet_f','static_cast<float>(njet)'),
  #           ('nlep_f','static_cast<float>(nlep)'),
  #           ('mllg','llphoton_m[0]')]
  #branches = ('w_lumi_year','mllg','photon_mva','pt_mass','min_dR',
  #            'dphi_h_met','met','njet_f','nlep_f')
  #also: ht, met
  #make n-tuples for 3 years, then merge
  years = ['2016APV','2016','2017','2018']
  lumis = ['19.51','16.80','41.48','59.83']
  #lumis = ['89.4','102.2','147.1'] #scaled to run 3
  base_name = 'phidcomp_post'
  for iyear in range(len(years)):
    defines_year = defines + [('w_lumi_year','float w_lumi_year=w_lumi*'+lumis[iyear]+'; return w_lumi_year;')]
    skim_and_slim.write_ntuples(
        #['/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*HToZG*M-125*',
        #'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*ttHToZG_M125*'],
        #['/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*GluGluHToZG*M-125*'],
        ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*GluGluHToZG_ZToLL_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*VBFHToZG_ZToLL_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*WminusH_HToZG_WToAll_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*WplusH_HToZG_WToAll_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*ZH_HToZG_ZToAll_M-125*',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*ttHToZG_M125*'],
        cuts,
        base_name+'_sig_'+years[iyear]+'.root',
        defines_year,
        'tree',
        branches)
    skim_and_slim.write_ntuples(
        ##['/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM*.root',
        #['/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root',
        #'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*ZGToLLG_01J_5f_Tune*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*EWKZ2Jets*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*TGJets*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*ttZJets*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*ttWJets*.root',
        #'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*TTTo2L2Nu*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*_WW_Tune*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*_WZ_Tune*.root',
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*_ZZ_Tune*.root'],
        ##'/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*ZGToLLG_01J_5f_Tune*.root'],
        ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*DYJetsToLL_M-50_TuneCP5_13TeV-amcatnlo*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*ZGToLLG_01J_5f_lowMLL*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*TTTo2L2Nu*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*TGJets*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*_WW_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*_WZ_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*_ZZ_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*_WWG_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*_WZG_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*EWKZ2Jets*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*ttWJets_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*ttZJets_Tune*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v0/'+years[iyear]+'/mc/skim_llg/*WGToLNuG*.root'],
        cuts,
        base_name+'_bak_'+years[iyear]+'.root',
        defines_year,
        'tree',
        branches)
    #only fake photon background
    #skim_and_slim.write_ntuples(
    #    ['/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/'+years[iyear]+'/mc/skim_llg/*DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM*.root'],
    #    cuts+['photon_pflavor[0]!=1'],
    #    base_name+'_spr_'+years[iyear]+'.root',
    #    defines_year,
    #    'tree',
    #    branches)
  #merge n-tuples and shuffle
  print('Merging and shuffling output files')
  for sigbak in ['sig','bak']:
  #for sigbak in ['sig','bak','spr']:
    merged_tree = ROOT.TChain("tree", "tree")
    for year in years:
      merged_tree.Add('ntuples/'+base_name+'_'+sigbak+'_'+year+'.root')
    out_file = ROOT.TFile('ntuples/'+base_name+'_'+sigbak+'.root', "recreate", "", 0)
    merged_tree.Merge(out_file,0,'fast keep')
    out_file.Close()
    subprocess.run(('/data1/jbkim/Linux/el7_v1/bin/python3.9 scripts/shuffle_tree.py ntuples/'+base_name+'_'+sigbak+'.root tree').split())
  #clean up
  for sigbak in ['sig','bak']:
  #for sigbak in ['sig','bak','spr']:
    for year in years:
      subprocess.run(('rm ntuples/'+base_name+'_'+sigbak+'_'+year+'.root').split())
    subprocess.run(('rm ntuples/'+base_name+'_'+sigbak+'.root').split())

if __name__=='__main__':
  produce_vbf_samples()
