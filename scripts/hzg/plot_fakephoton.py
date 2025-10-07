#!/usr/bin/env python3
"""@package docstring
Make plots for truth matching
"""

from array import array
from root_plot_lib import RplPlot
import ROOT

ROOT.gInterpreter.Declare("""

template <class C>
using RVec = ROOT::VecOps::RVec<C>;

float get_sign(float a) {
  if (a > 0) return 1.0;
  else if (a==0) return 0.0;
  return -1.0;
}

float bound(float val, float lower, float upper) {
  if (val < lower)
    return lower;
  else if (val > upper)
    return upper;
  return val;
}

float dr(float eta1, float phi1, float eta2, float phi2) {
  float dphi = fmod(fabs(phi2-phi1), 2*M_PI);
  dphi = dphi > M_PI ? 2*M_PI-dphi : dphi;
  return hypot(eta2-eta1, dphi);
}

bool get_photon_ishard(RVec<float> photon_eta, RVec<float> photon_phi, 
    RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi, 
    RVec<int> mc_statusflag, RVec<int> mc_id, float dr_cone) {
  if (photon_eta.size()==0) return false;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    if (mc_id[imc] == 6 || mc_id[imc] == -6 || mc_id[imc] == 23 
        || mc_id[imc] == 24 || mc_id[imc] == -24 || mc_id[imc] == 25)
      continue;
    if ((mc_statusflag[imc] & 0x100) && (mc_pt[imc] > 15)) {
      if (dr(photon_eta[0], photon_phi[0], mc_eta[imc], mc_phi[imc])<dr_cone) {
        return true;
      }
    }
  }
  return false;
}

""")

def main():
  """Make plots for checking photon matching
  """
  ROOT.EnableImplicitMT()

  df = ROOT.RDataFrame('tree','/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/merged_zgmc_llg/merged_pico_llg_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_zgmc_llg_nfiles_693.root')

  df = df.Define('photon_ishard0p6','get_photon_ishard(photon_eta, photon_phi,'
    'mc_pt, mc_eta, mc_phi, mc_statusflag, mc_id, 0.6);')
  df = df.Define('photon_ishard0p8','get_photon_ishard(photon_eta, photon_phi,'
    'mc_pt, mc_eta, mc_phi, mc_statusflag, mc_id, 0.8);')
  df = df.Define('lead_photon_eta','photon_eta[0]')
  df = df.Define('lead_photon_hardprocess','photon_hardprocess[0]')
      
  pheta_0p4_hist_ptr = df.Filter('!lead_photon_hardprocess').Histo1D(
        ('pheta_0p4_hist',
        f'#Delta R=0.4;Photon #eta;Fraction events',40,-2.5,2.5),
        'lead_photon_eta','weight')
  pheta_0p6_hist_ptr = df.Filter('!photon_ishard0p6').Histo1D(
        ('pheta_0p6_hist',
        f'#Delta R=0.6;Photon #eta;Fraction events',40,-2.5,2.5),
        'lead_photon_eta','weight')
  pheta_0p8_hist_ptr = df.Filter('!photon_ishard0p8').Histo1D(
        ('pheta_0p8_hist',
        f'#Delta R=0.8;Photon #eta;Fraction events',40,-2.5,2.5),
        'lead_photon_eta','weight')

  pheta_0p4_hist = pheta_0p4_hist_ptr.GetValue()
  pheta_0p6_hist = pheta_0p6_hist_ptr.GetValue()
  pheta_0p8_hist = pheta_0p8_hist_ptr.GetValue()
  pheta_0p4_hist.Scale(1.0/pheta_0p4_hist.Integral())
  pheta_0p6_hist.Scale(1.0/pheta_0p6_hist.Integral())
  pheta_0p8_hist.Scale(1.0/pheta_0p8_hist.Integral())

  plot = RplPlot()
  plot.lumi_data = [(27,13.6)]
  plot.y_title = '% Events/bin'
  plot.plot_outline(pheta_0p4_hist)
  plot.plot_outline(pheta_0p6_hist)
  plot.plot_outline(pheta_0p8_hist)
  plot.add_ratio('pheta_0p6_hist','pheta_0p4_hist',True)
  plot.add_ratio('pheta_0p8_hist','pheta_0p4_hist',True)
  plot.y_title_lower = 'Variation/nom.'
  plot.draw('plots/pheta_matching_comp.pdf')

if __name__ == '__main__':
  main()
