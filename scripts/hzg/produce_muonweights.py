#!/usr/bin/env python3
"""@package docstring
Generates our format muon weights from H->ZZ
"""

from argparse import ArgumentParser
from array import array
from correctionlib import schemav2 
from math import sqrt, hypot
import json
import ROOT

#constants
ETA_BINS = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4]
PT_BINS = [5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 200.0]
PT_BINS_2022 = [5.0, 12.0, 17.0, 22.0, 27.0, 32.0, 40.0, 50.0, 60.0, 80.0, 200.0]

#ROOT JIT'ed C++ definitions
ROOT.gInterpreter.Declare('''
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using std::vector;

float dr(float eta1, float phi1, float eta2, float phi2) {
  float dphi = fmod(fabs(phi2-phi1), 2.0*M_PI);
  if (dphi > M_PI) {
    dphi = 2*M_PI-dphi;
  }
  return sqrt((eta2-eta1)*(eta2-eta1)+dphi*dphi);
}

RVec<bool> mc_muonreco(RVec<bool> mu_sig, RVec<float> mu_eta, 
                         RVec<float> mu_phi, RVec<float> mc_eta, 
                         RVec<float> mc_phi) {
  vector<bool> found_muon;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    bool match = false;
    for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
      if (mu_sig[imu]) {
        if (dr(mu_eta[imu], mu_phi[imu], mc_eta[imc], mc_phi[imc])<0.2) {
          match = true;
        }
      }
    }
    found_muon.push_back(match);
  }
  return found_muon;
}

''')

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

def generate_mchists(pico_filename,pt_bins):
  '''Generate MC efficiency ratio plot

  pico_filename     filename of Z->mumu pico(s) for appropriate era
  pt_bins           list of bin boundaries (floats)
  '''
  #get MC efficiencies from picos
  df = ROOT.RDataFrame('tree',pico_filename)
  #is prompt and is first copy
  status_str = '((mc_statusflag&0x1)!=0)&&((mc_statusflag&0x1000)!=0)'
  df = df.Define('mc_mureco','mc_muonreco(mu_sig,mu_eta,mu_phi,mc_eta,mc_phi)')
  df = df.Define('mcmu_pt','mc_pt[abs(mc_id)==13&&' + status_str+']')
  df = df.Define('mcmu_eta','mc_eta[abs(mc_id)==13&&' + status_str+']')
  df = df.Define('mcmu_reco','mc_mureco[abs(mc_id)==13&&' +status_str+']')
  df = df.Define('mcmatchmu_pt','mcmu_pt[mcmu_reco]')
  df = df.Define('mcmatchmu_eta','mcmu_eta[mcmu_reco]')
  den_hist_ptr = df.Histo2D(('den_hist','den_hist',len(ETA_BINS)-1,
                            array('d',ETA_BINS),len(pt_bins)-1,
                            array('d',pt_bins)),'mcmu_eta','mcmu_pt','w_lumi')
  num_hist_ptr = df.Histo2D(('num_hist','num_hist',len(ETA_BINS)-1,
                            array('d',ETA_BINS),len(pt_bins)-1,
                            array('d',pt_bins)),'mcmatchmu_eta','mcmatchmu_pt',
                            'w_lumi')
  print('Processing DY samples, this may take some time.')
  ratio_hist = num_hist_ptr.GetPtr()
  den_hist = den_hist_ptr.GetPtr()
  ratio_hist.Divide(den_hist)
  return ratio_hist

def make_correction(name, values, pt_bins):
  '''Generates correctionlib corrections with given name and values
  '''
  return schemav2.Correction(
      name=name,
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='Muon pt'),
              schemav2.Variable(name='eta', type='real', description='Muon eta')],
      output=schemav2.Variable(name='sf', type='real', description=name),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','eta'],
          edges=[pt_bins,ETA_BINS],
          content=values,
          flow='clamp',
          ),
      )

def generate_correctionlib(hzzroot_filename, pico_filename, out_filename, year):
  '''Generate correctionlib file with muon scale factors

  hzzroot_filename  filename of H->ZZ derived SFs
  pico_filename     filename of Z->mumu pico(s) for appropriate era
  '''

  #get inputs
  hzzfile = ROOT.TFile(hzzroot_filename,'READ')
  sf_hist = hzzfile.Get('FINAL')
  unc_hist = hzzfile.Get('ERROR')
  sf_hist.SetDirectory(ROOT.nullptr)
  unc_hist.SetDirectory(ROOT.nullptr)
  hzzfile.Close()
  pt_bins = PT_BINS
  if year=='2022':
    pt_bins = PT_BINS_2022
  mceff_hist = generate_mchists(pico_filename,pt_bins)

  #calculate SFs and uncertainties
  pass_sf = []
  pass_unc = []
  fail_sf = []
  fail_unc = []
  for ipt in range(len(pt_bins)-1):
    for ieta in range(len(ETA_BINS)-1):
      pass_sf.append(sf_hist.GetBinContent(ieta+1,ipt+1))
      pass_unc.append(unc_hist.GetBinContent(ieta+1,ipt+1))
      mc_eff = mceff_hist.GetBinContent(ieta+1,ipt+1)
      mc_unc = mceff_hist.GetBinError(ieta+1,ipt+1)
      mc_relunc = 1.0
      if (mc_eff > 0.0):
        mc_relunc = mc_unc/mc_eff
      data_eff = mc_eff*pass_sf[-1]
      data_unc = data_eff*hypot(mc_relunc,pass_unc[-1]/pass_sf[-1])
      data_relunc = 1.0
      if (data_eff > 0.0):
        data_relunc = data_unc/data_eff
      fail_sf.append(1.0)
      if (mc_eff < 1.0):
        fail_sf[-1] = (1.0-data_eff)/(1.0-mc_eff)
      fail_unc.append(fail_sf[-1]*hypot(mc_relunc,data_relunc))

  #generate output
  json_pass_sf = make_correction('sf_pass',pass_sf,pt_bins)
  json_pass_unc = make_correction('unc_pass',pass_unc,pt_bins)
  json_fail_sf = make_correction('sf_fail',fail_sf,pt_bins)
  json_fail_unc = make_correction('unc_fail',fail_unc,pt_bins)

  with open(out_filename,'w') as out_file:
    out_file.write(fix_correctionlib_json([
        json_pass_sf.json(exclude_unset=True),
        json_pass_unc.json(exclude_unset=True),
        json_fail_sf.json(exclude_unset=True),
        json_fail_unc.json(exclude_unset=True)
        ]))

if __name__ == '__main__':
  parser = ArgumentParser(prog='process_muonweights',
                          description='Generates muon correctionlib files')
  parser.add_argument('-z','--hzzroot_filename')
  parser.add_argument('-y','--year',default='2016')
  parser.add_argument('-d','--dypico_filename')
  parser.add_argument('-o','--output_filename')
  args  = parser.parse_args()

  ROOT.EnableImplicitMT()

  generate_correctionlib(args.hzzroot_filename, args.dypico_filename, 
                         args.output_filename, args.year)
