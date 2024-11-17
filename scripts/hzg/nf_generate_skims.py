#!/usr/bin/env python3
'''
Script that generates slimmed n-tuples for photon Fall17v2 IDMVA corrections
'''
from argparse import ArgumentParser
from glob import glob
from skim_and_slim import write_ntuples
from subprocess import run as subprocessrun
import ROOT

#Define helper functions in C++
ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;

//get leading pT photon index in pico
int get_photon_idx_lead(RVec<float> photon_pt) {
  int lead_idx = -1;
  float max_pt = -1;
  for (unsigned iph = 0; iph < photon_pt.size(); iph++) {
    if (photon_pt[iph] > max_pt) {
      lead_idx = static_cast<int>(iph);
      max_pt = photon_pt[iph];
    }
  }
  return lead_idx;
}

//basic cuts to apply to photon with leading pT
bool lead_photon_cuts(RVec<float> photon_pt, RVec<bool> photon_isScEtaEB,
                      RVec<bool> photon_isScEtaEE, RVec<float> photon_idmva) {
  int lead_idx = get_photon_idx_lead(photon_pt);
  return (photon_pt[lead_idx]>15&&((photon_isScEtaEB[lead_idx]&&photon_idmva[lead_idx]>-0.4)||(photon_isScEtaEE[lead_idx]&&photon_idmva[lead_idx]>-0.58)));
}

""")

if __name__=='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='generate_photonid_skim',
      description='Generates a skim for deriving photon ID correction from input file (-i) and outputs it to output file (-o)')
  argument_parser.add_argument('-i','--input_filename')
  argument_parser.add_argument('-f','--input_format')
  argument_parser.add_argument('-o','--output_filename')
  args = argument_parser.parse_args()

  if not args.input_format in ['pico','EGMtnp']:
    raise ValueError('Unsupported format')

  #setup skim for EGM TnP file format
  cuts = ['tag_Ele_pt>35&&tag_sc_abseta<2.17&&(tag_sc_abseta>1.566||tag_sc_abseta<1.4442)',
          'ph_et>=15&&ph_sc_abseta<2.5&&((ph_sc_abseta<1.4442&&ph_mva94XV2>-0.4)||(ph_sc_abseta>1.566&&ph_mva94XV2>-0.58))']
  defines = [('ph_esEnergyOverRawE','static_cast<float>((ph_preshower_energy_plane1+ph_preshower_energy_plane2)/ph_sc_rawEnergy)')]
  branches = ('pair_mass','ph_et','ph_eta','ph_r9','ph_s4','ph_sc_etaWidth',
              'ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso',
              'ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','ph_sc_abseta',
              'event_rho','ph_ESsigma','ph_esEnergyOverRawE','ph_mva94XV2',
              'event')
  tree_name = 'tnpPhoIDs/fitter_tree',

  #setup skim for Pico format. More fomats can be added in a similar fashion
  if args.input_format == 'pico':
    #select mumugamma events with photon passing HIG-19-014 selection
    cuts = ['nmu==2&&photon_pt.size()>0',
            'lead_photon_cuts(photon_pt, photon_isScEtaEB, photon_isScEtaEE, photon_idmva)']
    defines = [('idx_lead','get_photon_idx_lead(photon_pt)'),
               ('ph_et','photon_pt[idx_lead]'),
               ('ph_eta','photon_eta[idx_lead]'),
               ('ph_r9','photon_r9[idx_lead]'),
               ('ph_s4','photon_s4[idx_lead]'),
               ('ph_sc_etaWidth','photon_etawidth[idx_lead]'),
               ('ph_sc_phiWidth','photon_phiwidth[idx_lead]'),
               ('ph_sieie','photon_sieie[idx_lead]'),
               ('ph_sieip','photon_sieip[idx_lead]'),
               ('ph_phoIso','photon_phiso[idx_lead]'),
               ('ph_chIso','photon_chiso[idx_lead]'),
               ('ph_chWorIso','photon_chiso_worst[idx_lead]'),
               ('ph_sc_rawEnergy','photon_energy_raw[idx_lead]'),
               ('ph_sc_eta','photon_origin_eta[idx_lead]'),
               ('ph_sc_abseta','fabs(photon_origin_eta[idx_lead])'),
               ('event_rho','rho'),
               ('ph_ESsigma','photon_essigmarr[idx_lead]'),
               ('ph_esEnergyOverRawE','photon_esoversc[idx_lead]'),
               ('ph_mva94XV2','photon_idmva[idx_lead]')]
    branches = ('ph_et','ph_eta','ph_r9','ph_s4','ph_sc_etaWidth',
                'ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso',
                'ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta','ph_sc_abseta',
                'event_rho','ph_ESsigma','ph_esEnergyOverRawE','ph_mva94XV2')
    tree_name = 'tree',

  #Get input files and run skim
  input_files = glob(args.input_filename)

  ROOT.EnableImplicitMT()
  write_ntuples(
      input_files,
      cuts,
      args.output_filename,
      defines,
      tree_name,
      branches)

