'''Script for evaluating effects of requiring trigger object matching

06.30.2025 MO
'''

from argparse import ArgumentParser
from array import array
import ROOT

#JIT C++
ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
const float PI = 3.14159;

bool get_match_singleel(RVec<int> TrigObj_id, RVec<int> TrigObj_filterBits, 
                        RVec<float> TrigObj_pt, float pt_threshold) {
  for (unsigned itrig = 0; itrig < TrigObj_id.size(); itrig++) {
    if (TrigObj_id[itrig] == 11 && (TrigObj_filterBits[itrig] && 0x2 != 0) 
        && TrigObj_pt[itrig] > pt_threshold) {
      return true;
    }
  }
  return false;
}

bool get_match_singlemu(RVec<int> TrigObj_id, RVec<int> TrigObj_filterBits, 
                        RVec<float> TrigObj_pt, float pt_threshold) {
  for (unsigned itrig = 0; itrig < TrigObj_id.size(); itrig++) {
    if (TrigObj_id[itrig] == 13 && (TrigObj_filterBits[itrig] && 0x2 != 0) 
        && TrigObj_pt[itrig] > pt_threshold) {
      return true;
    }
  }
  return false;
}

bool get_match_doubleel(RVec<int> TrigObj_id, RVec<int> TrigObj_filterBits, 
                        RVec<float> TrigObj_pt) {
  int upper = 0;
  int lower = 0;
  for (unsigned itrig = 0; itrig < TrigObj_id.size(); itrig++) {
    if (TrigObj_id[itrig] == 11 && (TrigObj_filterBits[itrig] && 0x1 != 0)) {
      if (TrigObj_pt[itrig] > 23) {
        upper++;
        lower++;
      }
      else if (TrigObj_pt[itrig] > 12) {
        lower++;
      }
    }
  }
  return (upper>=1 && lower>=2);
}

bool get_match_doublemu(RVec<int> TrigObj_id, RVec<int> TrigObj_filterBits, 
                        RVec<float> TrigObj_pt) {
  int upper = 0;
  int lower = 0;
  for (unsigned itrig = 0; itrig < TrigObj_id.size(); itrig++) {
    if (TrigObj_id[itrig] == 13 && (TrigObj_filterBits[itrig] && 0x1 != 0)) {
      if (TrigObj_pt[itrig] > 17) {
        upper++;
        lower++;
      }
      else if (TrigObj_pt[itrig] > 8) {
        lower++;
      }
    }
  }
  return (upper>=1 && lower>=2);
}

float delta_phi(float phi1, float phi2){
  float dphi = fmod(fabs(phi2-phi1), 2.0*PI);
  return dphi>PI ? 2.0*PI-dphi : dphi;
}

float delta_r(float eta1, float eta2, float phi1, float phi2) {
  float deta = eta1-eta2;
  float dphi = delta_phi(phi1, phi2);
  return sqrt(deta*deta+dphi*dphi);
}

float ConvertMVA(float mva_mini) {
  // 2.0 / (1.0 + exp(-2.0 * response)) - 1)
  float mva_nano = 2.0 / (1.0 + exp(-2.0 * mva_mini)) - 1;
  return mva_nano;
}

bool HzzId_WP2022(float pt, float etasc, float hzzmvaid) {
  //2022 WPs for 2018 ID training taken from https://indico.cern.ch/event/1429005/contributions/6039535/attachments/2891374/5077286/240712_H4lrun3_Approval.pdf
  if (pt < 10.0f) {
    if (fabs(etasc) < 0.8f) {
     return (hzzmvaid > ConvertMVA(1.6339));
    }
    else if (fabs(etasc) < 1.479f) {
      return (hzzmvaid > ConvertMVA(1.5499));
    }
    else {
      return (hzzmvaid > ConvertMVA(2.0629));
    }
  }
  else {
    if (fabs(etasc) < 0.8f) {
      return (hzzmvaid > ConvertMVA(0.3685));
    }
    else if (fabs(etasc) < 1.479f) {
      return (hzzmvaid > ConvertMVA(0.2662));
    }
    else {
      return (hzzmvaid > ConvertMVA(-0.5444));
    }
  }
}

RVec<bool> get_Muon_signal(RVec<float> Muon_pt, RVec<float> Muon_eta, 
    RVec<float> Muon_dxy, RVec<float> Muon_dz, RVec<bool> Muon_looseId, 
    RVec<float> Muon_pfRelIso03_all, RVec<float> Muon_sip3d) {
  RVec<bool> mu_sig;
  for (unsigned imu = 0; imu < Muon_pt.size(); imu++) {
    mu_sig.push_back((Muon_pt[imu] > 5) && 
                     (fabs(Muon_eta[imu])<2.4) && 
                     (fabs(Muon_dxy[imu]) < 0.5) && 
                     (fabs(Muon_dz[imu]) < 1.0) && 
                     (Muon_looseId[imu]) &&
                     (Muon_pfRelIso03_all[imu] < 0.35) &&
                     (Muon_sip3d[imu] < 4.0));
  }
  return mu_sig;
}

RVec<bool> get_Electron_signal_run2(RVec<float> Electron_pt, 
    RVec<float> Electron_eta, RVec<float> Electron_deltaEtaSC, 
    RVec<float> Electron_dxy, RVec<float> Electron_dz, 
    RVec<bool> Electron_mvaFall17V2Iso_WPL) {
  RVec<bool> el_sig;
  for (unsigned iel = 0; iel < Electron_pt.size(); iel++) {
    el_sig.push_back((Electron_pt[iel] > 7) && 
                     (fabs(Electron_eta[iel]+Electron_deltaEtaSC[iel])<2.5) && 
                     (fabs(Electron_dxy[iel]) < 0.5) && 
                     (fabs(Electron_dz[iel]) < 1.0) && 
                     (Electron_mvaFall17V2Iso_WPL[iel]));
  }
  return el_sig;
}

RVec<bool> get_Electron_signal_run3(RVec<float> Electron_pt, 
    RVec<float> Electron_eta, RVec<float> Electron_deltaEtaSC, 
    RVec<float> Electron_dxy, RVec<float> Electron_dz, 
    RVec<bool> Electron_mvaHZZIso) {
  RVec<bool> el_sig;
  for (unsigned iel = 0; iel < Electron_pt.size(); iel++) {
    float etasc = Electron_eta[iel]+Electron_deltaEtaSC[iel];
    el_sig.push_back((Electron_pt[iel] > 7) && 
                     (fabs(etasc)<2.5) && 
                     (fabs(Electron_dxy[iel]) < 0.5) && 
                     (fabs(Electron_dz[iel]) < 1.0) && 
                     HzzId_WP2022(Electron_pt[iel], etasc, 
                                  Electron_mvaHZZIso[iel]));
  }
  return el_sig;
}

int reduce_bool_sum(RVec<bool> vec_entry) {
  int sum = 0;
  for (unsigned i = 0; i < vec_entry.size(); i++) {
    if (vec_entry[i]) sum++;
  }
  return sum;
}

//golden json loader
class GoldenJsonLoader {
public:
  GoldenJsonLoader();
  GoldenJsonLoader(GoldenJsonLoader &&) = default;
  GoldenJsonLoader& operator=(GoldenJsonLoader &&) = default;
  ~GoldenJsonLoader() = default;
  bool pass_json(unsigned int run, unsigned int luminosityBlock);
private:
  std::vector<std::vector<int>> VVRunLumi;
};

GoldenJsonLoader::GoldenJsonLoader() {
  std::vector<std::string> golden_filenames = {
    "json/golden_Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "json/golden_Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
    "json/golden_Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
    "json/Cert_Collisions2022_355100_362760_Golden.json",
    "json/Cert_Collisions2023_366442_370790_Golden.json"};
  for (std::string golden_filename : golden_filenames) {
    std::ifstream orgJSON;
    orgJSON.open(golden_filename.c_str());
    std::vector<int> VRunLumi;
    if(orgJSON.is_open()){
      char inChar;
      int inInt;
      std::string str;
      while(!orgJSON.eof()){
        char next = orgJSON.peek();
        if( next == '1' || next == '2' || next == '3' ||
            next == '4' || next == '5' || next == '6' ||
            next == '7' || next == '8' || next == '9' || 
            next == '0'){     
          orgJSON >> inInt;
          VRunLumi.push_back(inInt);        
        }
        else if(next == ' '){
          getline(orgJSON,str,' ');
        }
        else{
          orgJSON>>inChar;
        }
      }
    }//check if the file opened.
    else{
      std::cout<<"Invalid JSON File:"<<golden_filename<<"!\\n";
    }
    orgJSON.close();
    if(VRunLumi.size() == 0){
      std::cout<<"No Lumiblock found in JSON file\\n";
    }
    for(unsigned int i = 0; i+2 < VRunLumi.size();){
      if(VRunLumi[i] > 130000){
        std::vector<int> RunLumi;
        RunLumi.push_back(VRunLumi[i]);
        while(VRunLumi[i+1] < 130000 && i+1 < VRunLumi.size()){
          RunLumi.push_back(VRunLumi[i+1]);
          ++i;
        }
        VVRunLumi.push_back(RunLumi);
        ++i;
      }
    }
  }
}

bool GoldenJsonLoader::pass_json(unsigned int run, 
                                 unsigned int luminosityBlock) {
  bool answer = false;
  if(run < 120000){
    answer = true;
  }
  else{
    for(unsigned int i = 0; i < VVRunLumi.size();++i){
      if(run == VVRunLumi[i][0]){
        for(unsigned int j = 1; j+1 < VVRunLumi[i].size();j=j+2){
          if(luminosityBlock >= VVRunLumi[i][j] && luminosityBlock <= VVRunLumi[i][j+1]){
            answer = true;
          }
        }
      }
    }
  }
  return answer;
}


""")

def get_triggers(year):
  '''Gets trigger names for given year

  Args:
    years: data-taking period

  Returns:
    names of single electron, dielectron, single muon, and dimuon triggers
  '''
  if year=='2016' or year=='2016APV':
    return ['HLT_Ele27_WPTight_Gsf',
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
            'HLT_IsoMu24||HLT_IsoTkMu24',
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL||'
            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL||' 
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ||'
            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ']
  elif year=='2017':
    return ['HLT_Ele32_WPTight_Gsf_L1DoubleEG',
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
            'HLT_IsoMu27',
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||'
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8']
  elif year=='2018':
    return ['HLT_Ele32_WPTight_Gsf',
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
            'HLT_IsoMu24',
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8']
  else:
    return ['HLT_Ele30_WPTight_Gsf',
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
            'HLT_IsoMu24',
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8']

def get_run(year):
  '''Returns if year is part of run 2 or run 3 campaign

  Args:
    year: data-taking period

  Returns:
    data-taking campaign
  '''
  if (year=='2016APV' or year=='2016' or year=='2017' or year=='2018'):
    return 'run2'
  return 'run3'

def get_singlelep_threshold(year):
  '''Returns pt threshold for lowest unprescaled single lepton triggers

  Args:
    year: data-taking period

  Returns:
    electron and muon trigger thresholds
  '''
  if (year=='2016' or year=='2016APV'):
    return ['27','24']
  elif (year=='2017'):
    return ['32','27']
  elif (year=='2018'):
    return ['32','24']
  else:
    return ['30','24']

if __name__=='__main__':

  #parse arguments
  argument_parser = ArgumentParser(prog='trigger_study',
      description='Script to check trigger efficiencies.')
  argument_parser.add_argument('-y','--year')
  argument_parser.add_argument('-m','--mc',action='store_true')
  argument_parser.add_argument('-i','--input_file', help='NanoAOD input')
  args = argument_parser.parse_args()

  trig_names = get_triggers(args.year)
  year_threshold = get_singlelep_threshold(args.year)
  run = get_run(args.year)

  ROOT.EnableImplicitMT()
  ROOT.gROOT.ProcessLine('GoldenJsonLoader golden_json_loader;')

  df = ROOT.RDataFrame('Events',args.input_file)
  if (run == 'run2'):
    df = df.Define('Electron_signal','get_Electron_signal_run2(Electron_pt, '
                   'Electron_eta, Electron_deltaEtaSC, Electron_dxy, '
                   'Electron_dz, Electron_mvaFall17V2Iso_WPL)')
  else:
    df = df.Define('Electron_signal','get_Electron_signal_run3(Electron_pt, '
                   'Electron_eta, Electron_deltaEtaSC, Electron_dxy, '
                   'Electron_dz, Electron_mvaHZZIso)')
  df = df.Define('Muon_signal','get_Muon_signal(Muon_pt, Muon_eta, Muon_dxy, '
                 'Muon_dz, Muon_looseId, Muon_pfRelIso03_all, Muon_sip3d)')
  df = df.Define('nSignalElectron','reduce_bool_sum(Electron_signal)')
  df = df.Define('nSignalMuon','reduce_bool_sum(Muon_signal)')
  df = df.Define('match_singleel','get_match_singleel(TrigObj_id, '
                 'TrigObj_filterBits, TrigObj_pt, '+year_threshold[0]+')')
  df = df.Define('match_singlemu','get_match_singlemu(TrigObj_id, '
                 'TrigObj_filterBits, TrigObj_pt, '+year_threshold[1]+')')
  df = df.Define('match_doubleel','get_match_doubleel(TrigObj_id, '
                 'TrigObj_filterBits, TrigObj_pt)')
  df = df.Define('match_doublemu','get_match_doublemu(TrigObj_id, '
                 'TrigObj_filterBits, TrigObj_pt)')
  df = df.Define('trig_singleel',trig_names[0])
  df = df.Define('trig_doubleel',trig_names[1])
  df = df.Define('trig_singlemu',trig_names[2])
  df = df.Define('trig_doublemu',trig_names[3])
  if args.mc:
    df = df.Define('w_lumi','Generator_weight>0 ? 1.0 : -1.0')
  else:
    df = df.Filter('golden_json_loader.pass_json(run, luminosityBlock)')
    df = df.Define('w_lumi','1.0')

  trig_cats = ['singleel','doubleel','singlemu','doublemu']
  meas_denom = ['nSignalElectron>=2','nSignalElectron>=2','nSignalMuon>=2',
                'nSignalMuon>=2', '(nSignalElectron>=2||nSignalMuon>=2)']
  meas_numr1 = ['trig_' + trig_cat for trig_cat in trig_cats]
  meas_numr2 = ['trig_' + trig_cat + '&&match_' + trig_cat for trig_cat 
                in trig_cats]
  trig_cats.append('or')
  meas_numr1.append('(trig_singleel||trig_doubleel||trig_singlemu||'
                    'trig_doublemu)')
  meas_numr2.append('((trig_singleel&&match_singleel)||'
                    '(trig_doubleel&&match_doubleel)||'
                    '(trig_singlemu&&match_singlemu)||'
                    '(trig_doublemu&&match_doublemu))')
  denom_ptrs = []
  numr1_ptrs = []
  numr2_ptrs = []

  for icat in range(len(trig_cats)):
    denom_ptrs.append(df.Filter(meas_denom[icat]).Sum('w_lumi'))
    numr1_ptrs.append(df.Filter(meas_denom[icat]+'&&'+meas_numr1[icat]).Sum(
        'w_lumi'))
    numr2_ptrs.append(df.Filter(meas_denom[icat]+'&&'+meas_numr2[icat]).Sum(
        'w_lumi'))

  for icat in range(len(trig_cats)):
    denom = denom_ptrs[icat].GetValue()
    if (denom==0):
      print('In category '+trig_cats[icat]+' no denominator events')
      continue
    trig_eff = numr1_ptrs[icat].GetValue()/denom
    match_eff = numr2_ptrs[icat].GetValue()/denom
    print('In category '+trig_cats[icat]+' trig eff: '+str(trig_eff))
    print('In category '+trig_cats[icat]+' match eff: '+str(match_eff))
  
