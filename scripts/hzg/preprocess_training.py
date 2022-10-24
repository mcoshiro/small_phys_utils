#!/usr/bin/env python3
'''
Script that generates slimmed n-tuples for H->Zgamma MVA training
'''
import ROOT
import skim_and_slim

ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;


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

""")

if __name__=='__main__':
  ROOT.EnableImplicitMT()
  cuts = ['HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
      'll_m.size()>0&&llphoton_m.size()>0',
      'stitch_dy||(type/1000!=6)',
      '(ll_m[0]>50)&&(photon_pt[0]/llphoton_m[0]>=15.0/110.0)&&((llphoton_m[0]+ll_m[0])>185)&&(photon_drmin[0]>0.4)',
      '(ll_lepid[0]==11&&el_pt[ll_i1[0]]>25&&el_pt[ll_i2[0]]>15)||(ll_lepid[0]==13&&mu_pt[ll_i1[0]]>20&&mu_pt[ll_i2[0]]>10)',
      'photon_idmva[0]>0.5',
      'llphoton_m[0]>120&&llphoton_m[0]<130']
  #defines = [('higgsdr','llphoton_dr[0]'),('higgspt','llphoton_pt[0]'),('zpt','ll_pt[0]'),('phpt','photon_pt[0]')]
  #branches = ('higgsdr','higgspt','zpt','phpt')
  defines = [('photon_mva','photon_idmva[0]'),
             ('min_dR','photon_drmin[0]'),
             ('max_dR','get_max_dr(photon_eta,photon_phi,el_eta,el_phi,mu_eta,mu_phi,ll_lepid,ll_i1,ll_i2)'),
             ('pt_mass','llphoton_pt[0]/llphoton_m[0]'),
             ('cosTheta','llphoton_cosTheta[0]'),
             ('costheta','llphoton_costheta[0]'),
             ('phi','llphoton_phi[0]'),
             ('photon_res','photon_pterr[0]/photon_pt[0]'),
             ('photon_rapidity','photon_eta[0]'),
             ('l1_rapidity','get_l1_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('l2_rapidity','get_l2_rapidity(el_pt,el_eta,mu_pt,mu_eta,ll_lepid,ll_i1,ll_i2)'),
             ('photon_ptransverse','photon_pt[0]'),
             #('decorr_photon_pt','photon_pt[0]-0.207*llphoton_m[0]'),
             ('photon_pt_mass','photon_pt[0]/llphoton_m[0]')]
  branches = ('photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta',
      'phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity','photon_ptransverse','photon_pt_mass','weight')
  #define drmax, pt_mass, first index
  #make n-tuples
  write_ntuples(['/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal/2017/signal/skim_llg/*.root'],
      cuts,
      'train_kinbdt_idmva_nomasscut_sig.root',
      defines,
      'tree',
      branches)
  write_ntuples(['/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/2016/mc/skim_llg/pico_llg_DYJetsToLL*.root','/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/2016/mc/skim_llg/pico_llg_ZGToLLG*.root'],
      cuts,
      'train_kinbdt_idmva_nomasscut_bak.root',
      defines,
      'tree',
      branches)

  
