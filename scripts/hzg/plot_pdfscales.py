#!/usr/bin/env python3
"""@package docstring
Make plots related to differential PDF and perturbative truncation 
uncertainties
"""

from array import array
from root_plot_lib import RplPlot
import ROOT

ROOT.gInterpreter.Declare("""

template <class C>
using RVec = ROOT::VecOps::RVec<C>;

float get_GenHiggs_pt(RVec<int> GenPart_pdgId, RVec<float> GenPart_pt, RVec<int> GenPart_statusFlags) {
  float GenHiggs_pt = -999;
  for (unsigned imc = 0; imc < GenPart_pdgId.size(); imc++) {
    //status flags hardprocess, lastcopy
    if (GenPart_pdgId[imc]==25 && (GenPart_statusFlags[imc] & 0x1080)) {
      GenHiggs_pt = GenPart_pt[imc];
    }
  }
  return GenHiggs_pt;
}

float get_GenHiggs_eta(RVec<int> GenPart_pdgId, RVec<float> GenPart_eta, RVec<int> GenPart_statusFlags) {
  float GenHiggs_eta = -999;
  for (unsigned imc = 0; imc < GenPart_pdgId.size(); imc++) {
    //status flags hardprocess, lastcopy
    if (GenPart_pdgId[imc]==25 && (GenPart_statusFlags[imc] & 0x1080)) {
      GenHiggs_eta = GenPart_eta[imc];
    }
  }
  return GenHiggs_eta;
}

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

""")

def get_sample(sample: str, era: str) -> ROOT.RDataFrame:
  """Get dataframe with sample

  Args:
    sample: ggf, vbf, vh, or tth
    era: 2018 or 2022EE
  """
  if not sample in ['ggf','vbf','vh','tth']:
    raise RuntimeError('Invalid sample')
  if not era in ['2018','2022EE']:
    raise RuntimeError('Invalid era')
  df = None
  if sample=='ggf':
    if era=='2018':
      df = ROOT.RDataFrame('Events','/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/GluGluHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8*.root')
    elif era=='2022EE':
      df = ROOT.RDataFrame('Events','/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8*.root')
  elif sample=='vbf':
    if era=='2018':
      df = ROOT.RDataFrame('Events','/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/VBFHToZG_ZToLL_M-125_TuneCP5_withDipoleRecoil_13TeV-powheg-pythia8*.root')
    elif era=='2022EE':
      df = ROOT.RDataFrame('Events','/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/VBFHtoZG_Zto2L_M-125_TuneCP5_withDipoleRecoil_13p6TeV_powheg-pythia8*.root')
  elif sample=='vh':
    #note that the JHUGEN-less samples do not have variations and will crash
    if era=='2018':
      df = ROOT.RDataFrame('Events',[
        '/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/WplusH_HToZG_WToAll_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/WminusH_HToZG_WToAll_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018*.root',
        '/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/ZH_HToZG_ZToAll_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018*.root'])
    elif era=='2022EE':
      df = ROOT.RDataFrame('Events',[
        '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/WminusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8*.root',
        '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/WplusH_HtoZG_WtoAll_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8*.root',
        '/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/ZH_HtoZG_ZtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8*.root'])
  elif sample=='tth':
    if era=='2018':
      df = ROOT.RDataFrame('Events','/net/cms11/cms11r0/pico/NanoAODv9/nano/2018/mc/ttHToZG_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8*.root')
    elif era=='2022EE':
      df = ROOT.RDataFrame('Events','/net/cms11/cms11r0/pico/NanoAODv12/nano/2022EE/mc/ttHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8*.root')
  return df

def make_murf_pdf_plots(sample: str, era: str):
  """Make plots for scale and pdf variations of given sample and era

  Args:
    sample: ggf, vbf, vh, or tth
    era: 2018 or 2022EE
  """
  ROOT.EnableImplicitMT()
  df = get_sample(sample, era)

  df = df.Define('GenHiggs_pt','get_GenHiggs_pt(GenPart_pdgId, GenPart_pt, '
                 +'GenPart_statusFlags)')
  df = df.Define('GenHiggs_eta','get_GenHiggs_eta(GenPart_pdgId, GenPart_eta, '
                 +'GenPart_statusFlags)')
  df = df.Define('genWeightSign','get_sign(genWeight)')

  print('Calculating normalizations.')
  murf_idxs = [0, 1, 3, 5, 7, 8]
  murf_names = [
      '#Lambda_{r} down #Lambda_{f} down', 
      '#Lambda_{r} down #Lambda_{f} nom.', 
      '#Lambda_{r} nom. #Lambda_{f} down', 
      '#Lambda_{r} nom. #Lambda_{f} up', 
      '#Lambda_{r} up #Lambda_{f} nom.', 
      '#Lambda_{r} up #Lambda_{f} up']
  murf_norm_ptrs = []
  murf_norm = []
  pdf_norm_ptrs = []
  pdf_norm = []
  #murf_norm_ptrs.append(df.Sum(f'genWeight'))
  for murf_idx in murf_idxs:
    df = df.Define(f'LHEScaleWeight{murf_idx}',
                   f'bound(LHEScaleWeight[{murf_idx}],-20,20)')
    murf_norm_ptrs.append(df.Sum(f'LHEScaleWeight{murf_idx}'))
  for pdf_idx in range(102):
    df = df.Define(f'LHEPdfWeight{pdf_idx}',
                   f'bound(LHEPdfWeight[{pdf_idx}],-20,20)')
    pdf_norm_ptrs.append(df.Sum(f'LHEPdfWeight{pdf_idx}'))
  imurf = 0
  for murf_norm_ptr in murf_norm_ptrs:
    murf_norm.append(murf_norm_ptr.GetValue())
    df = df.Define(f'norm_LHEScaleWeight{imurf}',
        f'genWeightSign*LHEScaleWeight{murf_idxs[imurf]}/{murf_norm[-1]}')
    imurf += 1
  ipdf = 0
  for pdf_norm_ptr in pdf_norm_ptrs:
    pdf_norm.append(pdf_norm_ptr.GetValue())
    df = df.Define(f'norm_LHEPdfWeight{ipdf}',
        f'genWeightSign*LHEPdfWeight{ipdf}/{pdf_norm[-1]}')
    ipdf += 1

  print('Generating histograms')
  plot_types = [('GenHiggs_pt','Gen Higgs p_{T} [GeV]',40,0,150), 
                ('GenHiggs_eta','Gen Higgs #eta',40,-6,6)]
  if sample in ['vbf', 'vh']:
    plot_types = [('GenHiggs_pt','Gen Higgs p_{T} [GeV]',40,0,350), 
                  ('GenHiggs_eta','Gen Higgs #eta',40,-6,6)]
  if sample=='tth':
    plot_types = [('GenHiggs_pt','Gen Higgs p_{T} [GeV]',40,0,400), 
                  ('GenHiggs_eta','Gen Higgs #eta',40,-6,6)]
  #lists [indexed by type] of lists [indexed by variation] of plots
  murf_hist_ptrs = []
  pdf_hist_ptrs = []
  for plot_type in plot_types:
    plot_var = plot_type[0]
    plot_desc = plot_type[1]
    murf_hist_ptrs.append([])
    pdf_hist_ptrs.append([])
    murf_hist_ptrs[-1].append(df.Histo1D((f'murf_nom_{plot_var}',
        f'nominal;{plot_desc};Fraction events',plot_type[2],plot_type[3],
        plot_type[4]),plot_var,'genWeightSign'))

  for imurf in range(len(murf_idxs)):
    iplot = 0
    for plot_type in plot_types:
      plot_var = plot_type[0]
      plot_desc = plot_type[1]
      murf_hist_ptrs[iplot].append(df.Histo1D((f'murf_var_{plot_var}{imurf}',
          f'{murf_names[imurf]};{plot_desc};Fraction events',
          plot_type[2],plot_type[3],plot_type[4]),plot_var,
          f'norm_LHEScaleWeight{imurf}'))
      iplot += 1
  for ipdf in range(102):
    iplot = 0
    for plot_type in plot_types:
      plot_var = plot_type[0]
      plot_desc = plot_type[1]
      pdf_hist_ptrs[iplot].append(df.Histo1D((f'pdf_var_{plot_var}{ipdf}',
          f'variations (101 plots);{plot_desc};Fraction events',
          plot_type[2],plot_type[3],plot_type[4]),plot_var,
          f'norm_LHEPdfWeight{ipdf}'))
      iplot += 1

  lumi_data = [(60,13)]
  if era=='2022EE':
    lumi_data = [(27,13.6)]

  murf_hists = []
  for iplot in range(len(plot_types)):
    plot_var = plot_types[iplot][0]
    murf_hists.append([])
    for imurf in range(len(murf_hist_ptrs[iplot])):
      murf_hists[-1].append(murf_hist_ptrs[iplot][imurf].GetValue())
      murf_hists[-1][-1].Scale(1.0/murf_hists[-1][-1].Integral())
    plot = RplPlot()
    plot.lumi_data = lumi_data
    plot.y_title = '% Events/bin'
    for imurf in range(len(murf_hist_ptrs[iplot])-1,-1,-1):
      plot.plot_outline(murf_hists[iplot][imurf])
    for imurf in range(len(murf_hist_ptrs[iplot])-1):
      plot.add_ratio(f'murf_var_{plot_var}{imurf}',f'murf_nom_{plot_var}',True)
    plot.y_title_lower = 'Variation/nom.'
    plot.draw(f'plots/{plot_var}_murf_{sample}_{era}.pdf')

  pdf_hists = []
  color_blue = ROOT.TColor.GetColor('#3f90da')
  color_orange = ROOT.TColor.GetColor('#ffa90e')
  for iplot in range(len(plot_types)):
    plot_var = plot_types[iplot][0]
    pdf_hists.append([])
    for ipdf in range(len(pdf_hist_ptrs[iplot])):
      pdf_hists[-1].append(pdf_hist_ptrs[iplot][ipdf].GetValue())
      pdf_hists[-1][-1].Scale(1.0/pdf_hists[-1][-1].Integral())
    plot = RplPlot()
    plot.lumi_data = lumi_data
    plot.y_title = '% Events/bin'
    for ipdf in range(len(pdf_hist_ptrs[iplot])-1,-1,-1):
      plot.plot_outline(pdf_hists[iplot][ipdf], color_orange)
      if ipdf != 0:
        plot.hist_legend[-1] = False
    plot.plot_outline(murf_hists[iplot][0], color_blue)
    for ipdf in range(len(pdf_hist_ptrs[iplot])-1):
      plot.add_ratio(f'pdf_var_{plot_var}{ipdf}',f'murf_nom_{plot_var}',True)
    plot.y_title_lower = 'Variation/nom.'
    plot.legend_ncolumns = 1
    plot.y_min_lower = 0.9
    plot.y_max_lower = 1.1
    plot.draw(f'plots/{plot_var}_pdf_{sample}_{era}.pdf')

def print_cpp_vector_initializer(name: str, init_list: list[float]):
  '''Prints const C++ vector initializer

  Args:
    name: variable name
    init_list: initializer values
  '''
  print(f'const static vector<float> {name} = {{',end='')
  for i in range(len(init_list)):
    if i != 0:
      print(',',end='')
    print(f'{init_list[i]:.6f}',end='')
  print('};')

def get_vh_variations():
  """Gets scale weights as a function of VH pT and prints C++ evaluator
  """
  ROOT.EnableImplicitMT()
  df = get_sample('vh', '2018')
  df = df.Define('GenHiggs_pt','get_GenHiggs_pt(GenPart_pdgId, GenPart_pt, '
                 +'GenPart_statusFlags)')
  df = df.Define('GenHiggs_eta','get_GenHiggs_eta(GenPart_pdgId, GenPart_eta, '
                 +'GenPart_statusFlags)')
  df = df.Define('genWeightSign','get_sign(genWeight)')
  pt_bins = [i*20.0 for i in range(18)]+[360.0]

  print('Calculating normalizations.')
  df = df.Define(f'LHEScaleWeight_up',f'LHEScaleWeight[3]')
  df = df.Define(f'LHEScaleWeight_dn',f'LHEScaleWeight[8]')
  weight_nm_norm_ptr = df.Sum(f'genWeightSign')
  weight_up_norm_ptr = df.Sum(f'LHEScaleWeight_up')
  weight_dn_norm_ptr = df.Sum(f'LHEScaleWeight_dn')
  weight_up_norm = weight_nm_norm_ptr.GetValue()/weight_up_norm_ptr.GetValue()
  weight_dn_norm = weight_nm_norm_ptr.GetValue()/weight_dn_norm_ptr.GetValue()
  df = df.Define(f'norm_LHEScaleWeight_up',
                 f'LHEScaleWeight[3]/{weight_up_norm}')
  df = df.Define(f'norm_LHEScaleWeight_dn',
                 f'LHEScaleWeight[8]/{weight_dn_norm}')

  print('Calculating weights.')
  hist_nm_ptr = df.Histo1D(('hist_nm',
      'Nominal;Gen Higgs p_{T} [GeV];Fraction events',len(pt_bins)-1,
      array('d',pt_bins)),'GenHiggs_pt','genWeightSign')
  hist_up_ptr = df.Histo1D(('hist_up',
      'Scale up;Gen Higgs p_{T} [GeV];Fraction events',len(pt_bins)-1,
      array('d',pt_bins)),'GenHiggs_pt','norm_LHEScaleWeight_up')
  hist_dn_ptr = df.Histo1D(('hist_dn',
      'Scale down;Gen Higgs p_{T} [GeV];Fraction events',len(pt_bins)-1,
      array('d',pt_bins)),'GenHiggs_pt','norm_LHEScaleWeight_dn')
  hist_nm = hist_nm_ptr.GetValue()
  hist_up = hist_up_ptr.GetValue()
  hist_dn = hist_dn_ptr.GetValue()
  hist_nm.Scale(1.0/hist_nm.Integral())
  hist_up.Scale(1.0/hist_up.Integral())
  hist_dn.Scale(1.0/hist_dn.Integral())
  hist_up.Divide(hist_nm)
  hist_dn.Divide(hist_nm)
  sfs_up = []
  sfs_dn = []
  for ibin in range(len(pt_bins)-1):
    sfs_up.append(hist_up.GetBinContent(ibin+1))
    sfs_dn.append(hist_dn.GetBinContent(ibin+1))

  #Print evaluator
  print_cpp_vector_initializer('pt_bin_edges',pt_bins)
  print_cpp_vector_initializer('sfs_up',sfs_up)
  print_cpp_vector_initializer('sfs_dn',sfs_dn)
  print('float sf = 1.0')
  print('for (unsigned ibin = 0; ibin < pt_bin_edges.size()-1; ibin++) {')
  print('  if (gen_higgs_pt >= pt_bin_edges[ibin] ')
  print('      && gen_higgs_pt < pt_bin_edges[ibin+1]) {')
  print('    if (variation_up)')
  print('      sf = sfs_up[ibin];')
  print('    else')
  print('      sf = sfs_dn[ibin];')
  print('  }')
  print('}')

if __name__ == '__main__':
  #for sample in ['ggf','vbf','vh','tth']:
  #  make_murf_pdf_plots(sample, '2022EE')
  get_vh_variations()
