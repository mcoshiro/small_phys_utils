#!/usr/bin/env python3
"""@package docstring
Generates MC photon efficiencies
"""

from argparse import ArgumentParser
from array import array
from correctionlib import schemav2 
from math import sqrt, hypot
import json
import ROOT

#constants
ETA_BINS = [-2.4,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.4]
PT_BINS = [15.0,20.0,35.0,50.0,100.0,200.0,500.0]

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

RVec<bool> found_ph(RVec<bool> ph_sig, RVec<float> ph_eta, 
                    RVec<float> ph_phi, RVec<float> mc_eta, 
                    RVec<float> mc_phi) {
  vector<bool> found_photon;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    bool match = false;
    for (unsigned iph = 0; iph < ph_sig.size(); iph++) {
      if (ph_sig[iph]) {
        if (dr(ph_eta[iph], ph_phi[iph], mc_eta[imc], mc_phi[imc])<0.2) {
          match = true;
        }
      }
    }
    found_photon.push_back(match);
  }
  return found_photon;
}

//checks there is a separation of at least 0.3 from the nearest (prompt) lepton
RVec<bool> mc_leptonsep(RVec<float> mc_eta, RVec<float> mc_phi, RVec<int> mc_id,
                        RVec<int> mc_statusflag) {
  vector<bool> no_nearby_lepton;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    bool nolepton = true;
    for (unsigned imc2 = 0; imc2 < mc_eta.size(); imc2++) {
      if ((mc_statusflag[imc2]&0x1)!=0&&(mc_statusflag[imc2]&0x1000)!=0) {
        if (abs(mc_id[imc2])==11 || abs(mc_id[imc2])==13) {
          if (dr(mc_eta[imc], mc_phi[imc], mc_eta[imc2], mc_phi[imc2]) < 0.3) {
            nolepton = false;
          }
        }
      }
    }
    no_nearby_lepton.push_back(nolepton);
  }
  return no_nearby_lepton;
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

def generate_mchists(pico_filename):
  '''Generate MC efficiency ratio plot

  pico_filename     filename of ZG pico(s) for appropriate era
  '''
  #get MC efficiencies from picos
  df = ROOT.RDataFrame('tree',pico_filename)
  #is prompt, stable, and first copy photon not near lepton
  selection_str = ('(mc_id==22)&&((mc_statusflag&0x1)!=0)'
                  +'&&((mc_statusflag&0x1000)!=0)'
                  +'&&(mc_status==1)&&(mc_lepveto)')
  df = df.Define('mc_phreco',
                 'found_ph(photon_sig,photon_eta,photon_phi,mc_eta,mc_phi)')
  df = df.Define('mc_lepveto',
                 'mc_leptonsep(mc_eta,mc_phi,mc_id,mc_statusflag)')
  df = df.Define('mcph_pt','mc_pt['+selection_str+']')
  df = df.Define('mcph_eta','mc_eta['+selection_str+']')
  df = df.Define('mcph_reco','mc_phreco['+selection_str+']')
  df = df.Define('mcmatchph_pt','mcph_pt[mcph_reco]')
  df = df.Define('mcmatchph_eta','mcph_eta[mcph_reco]')
  den_hist_ptr = df.Histo2D(('den_hist','den_hist',len(ETA_BINS)-1,
                            array('d',ETA_BINS),len(PT_BINS)-1,
                            array('d',PT_BINS)),'mcph_eta','mcph_pt','w_lumi')
  num_hist_ptr = df.Histo2D(('num_hist','num_hist',len(ETA_BINS)-1,
                            array('d',ETA_BINS),len(PT_BINS)-1,
                            array('d',PT_BINS)),'mcmatchph_eta','mcmatchph_pt',
                            'w_lumi')
  print('Processing DY samples, this may take some time.')
  ratio_hist = num_hist_ptr.GetPtr()
  den_hist = den_hist_ptr.GetPtr()
  ratio_hist.Divide(den_hist)
  return ratio_hist

def make_correction(name, values):
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
          edges=[PT_BINS,ETA_BINS],
          content=values,
          flow='clamp',
          ),
      )

def generate_correctionlib(pico_filename, out_filename):
  '''Generate correctionlib file with muon scale factors

  pico_filename  filename of Z->mumu pico(s) for appropriate era
  out_filename   output filename
  '''

  #get inputs
  mceff_hist = generate_mchists(pico_filename)

  #calculate effs and generate output
  eff_mc = []
  unc_mc = []
  for ipt in range(len(PT_BINS)-1):
    for ieta in range(len(ETA_BINS)-1):
      eff_mc.append(mceff_hist.GetBinContent(ieta+1,ipt+1))
      unc_mc.append(mceff_hist.GetBinError(ieta+1,ipt+1))

  json_eff_mc = make_correction('effmc',eff_mc)
  json_unc_mc = make_correction('systmc',unc_mc)

  with open(out_filename,'w') as out_file:
    out_file.write(fix_correctionlib_json([
        json_eff_mc.json(exclude_unset=True),
        json_unc_mc.json(exclude_unset=True),
        ]))

if __name__ == '__main__':
  parser = ArgumentParser(prog='process_muonweights',
                          description='Generates muon correctionlib files')
  parser.add_argument('-i','--input_filename')
  parser.add_argument('-o','--output_filename')
  args  = parser.parse_args()

  ROOT.EnableImplicitMT()

  generate_correctionlib(args.input_filename, args.output_filename)
