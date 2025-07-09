"""@package docstring
Library of Python utilities for working with TMVA BDTs

"""
import make_simple_datacard
import os
import ROOT
import root_utils
import subprocess
from functools import reduce

class BdtOptions:
  '''Simple class to hold BDT options
  '''
  def __init__(self):
    self.nTrain_Signal = 0
    self.nTrain_Background = 0
    self.SplitMode = 'Block'
    self.NTrees = 850
    self.MinNodeSize = 2.5
    self.MaxDepth = 3
    self.AdaBoostBeta = 0.5
    self.Shrinkage = 0.1
    self.BoostType = 'AdaBoost'

  def ptat_string(self):
    '''Return option to be given to DataLoader::PrepareTrainingAndTestTree
    '''
    if (self.SplitMode=='Block'):
      print('                         : ',end='')
      print('Using split mode Block, make sure your samples are pre-shuffled!')
    ret = 'nTrain_Signal='+str(self.nTrain_Signal)
    ret += ':nTrain_Background='+str(self.nTrain_Background)
    ret += ':SplitMode='+self.SplitMode
    ret += ':NormMode=NumEvents:!V'
    return ret

  def bm_string(self):
    '''Return option to be given to Factory::BookMethod
    '''
    ret = '!H:!V:NTrees='+str(self.NTrees)
    ret += ':MinNodeSize='+str(self.MinNodeSize)
    ret += '%:MaxDepth='+str(self.MaxDepth)
    ret += ':BoostType='+self.BoostType
    if (self.BoostType=='AdaBoost'):
      ret += ':AdaBoostBeta='+str(self.AdaBoostBeta)
    elif (self.BoostType=='Grad'):
      ret += ':Shrinkage='+str(self.Shrinkage)
    else:
      raise TypeError('Unknown BDT BoostType')
    ret += ':UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20'
    return ret

class BdtSummary:
  '''Class to hold summary of BDT performance
  '''
  def __init__(self):
    self.auc_train = 0
    self.auc_tests = 0
    self.csi_train = 0
    self.csi_tests = 0
    self.bn4_signif = 0
    self.bn1_signif = 0
    self.sig_ks = 0
    self.bak_ks = 0

  def print_summary(self):
    '''Print digest of BDT results
    '''
    print('AUC (train/test): '+str(self.auc_train)+'/'+str(self.auc_tests))
    print('CSI (train/test): '+str(self.csi_train)+'/'+str(self.csi_tests))
    print('4 bin sig (train+test): '+str(self.bn4_signif)+' ('+str(self.bn4_signif/self.bn1_signif-1.0)+'%impovement)')
    print('K-S p-value (sig/bak): '+str(self.sig_ks)+'/'+str(self.bak_ks))

def train_bdt(base_name, variables, weight, options, output_suffix='', 
              cut='1'):
  '''Function that trains a BDT

  @params
  base_name - filename of training samples, signal should be 
              base_name+_sig.root and background base_name+_bak.root
  variables - list of variable names to consider for BDT training
  weight - name of weight branch in file
  options - BDT options class
  output_suffix - optional tag for output file
  cut - selection applied to signal and background
  '''
  ROOT.TMVA.Tools.Instance()
  sig_input_file = ROOT.TFile('ntuples/'+base_name+'_sig.root','READ')
  bak_input_file = ROOT.TFile('ntuples/'+base_name+'_bak.root','READ')
  output_file = ROOT.TFile('ntuples/output_'+base_name+output_suffix+'.root','RECREATE') 
  bdt_factory = ROOT.TMVA.Factory(base_name+output_suffix,output_file,
      '!V:ROC:!Correlations:!Silent:Color:!DrawProgressBar:AnalysisType=Classification')
  bdt_loader = ROOT.TMVA.DataLoader('dataset')
  variable_types = [root_utils.get_root_type(cpp_type) for cpp_type in 
      root_utils.get_column_types(variables,'ntuples/'+base_name+'_sig.root','tree')]
  for ivar in range(len(variables)):
    bdt_loader.AddVariable(variables[ivar],variable_types[ivar])
  bdt_loader.SetBackgroundWeightExpression(weight)
  bdt_loader.SetSignalWeightExpression(weight)
  bdt_loader.AddSignalTree(sig_input_file.tree)
  bdt_loader.AddBackgroundTree(bak_input_file.tree)
  cut_s = ROOT.TCut(cut)
  cut_b = ROOT.TCut(cut)
  bdt_loader.PrepareTrainingAndTestTree(cut_s,cut_b,options.ptat_string())
  bdt_factory.BookMethod(bdt_loader,ROOT.TMVA.Types.kBDT,'BDT',options.bm_string())
  bdt_factory.TrainAllMethods()
  bdt_factory.TestAllMethods()
  bdt_factory.EvaluateAllMethods()
  output_file.Close()

#TODO: maybe change this so it can be combined with evaluate_bdt to avoid duplication
def optimize_binning(base_name, variables, weight, options, output_suffix='', mva_base_name='',nbins_max=4,min_evts=1.5,lumi_multiplier=1.0):
  '''Prints optimized binning for BDT given input ROOT file

  @params
  base_name - filename of training samples, see train_bdt
  variables - list of variables used in BDT
  weight - name of weight branch
  options - BdtOptions with options used for training
  output_suffix - tag for output file
  mva_base_name - filename used for evaluation
  nbins_max - maximum number of bins
  min_evts - minimum number of events per category
  lumi_multiplier - multiplier for luminosity
  '''
  if (mva_base_name == ''):
    mva_base_name = base_name
  variables_aug = variables+[weight]
  variable_types = root_utils.get_column_types(variables_aug,'ntuples/'+base_name+'_sig.root','tree')
  columns = [(variables_aug[i], variable_types[i]) for i in range(len(variables_aug))]
  out_file = open('bdt_evaluate_macro.cxx','w')
  write_event_loop_boilerplate(out_file, base_name, columns, weight, output_suffix, mva_base_name)
  out_file.write('  //Set up histograms\n')
  out_file.write('  TH1D hist_sig("hist_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_bak("hist_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  //Set up histograms\n')
  out_file.write('\n')
  out_file.write('  //Perform event loops\n')
  out_file.write('  for (long ievt = 0; ievt < sig_tree->GetEntries(); ievt++) {\n')
  out_file.write('    sig_tree->GetEntry(ievt);\n')
  out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    hist_sig.Fill(BDT,'+weight+');\n')
  out_file.write('  }\n')
  out_file.write('  for (long ievt = 0; ievt < bak_tree->GetEntries(); ievt++) {\n')
  out_file.write('    bak_tree->GetEntry(ievt);\n')
  out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    hist_bak.Fill(BDT,'+weight+');\n')
  out_file.write('  }\n')
  out_file.write('\n')
  out_file.write('  //get ouput\n')
  out_file.write('  get_roc_auc(&hist_sig, &hist_bak);\n')
  out_file.write('  binning_optimizer(&hist_sig, &hist_bak,'+str(nbins_max+1)+','+str(min_evts)+','+str(lumi_multiplier)+');\n')
  out_file.write('}\n')
  out_file.close()
  #run ROOT macro and parse output to get statistics and yields
  process_result = subprocess.run('root -l -q bdt_evaluate_macro.cxx'.split())
  subprocess.run('rm bdt_evaluate_macro.cxx'.split())

def make_bdt_friend_ttree(base_name, variables, weight, options, output_suffix='', mva_base_name=''):
  '''Makes a friend TTree with BDT score
  
  Prints optimized binning for BDT given input ROOT file

  @params
  base_name - filename of training samples, see train_bdt
  variables - list of variables used in BDT
  weight - name of weight branch
  options - BdtOptions with options used for training
  output_suffix - tag for output file
  mva_base_name - filename used for evaluation
  '''
  if (mva_base_name == ''):
    mva_base_name = base_name
  variables_aug = variables+[weight]
  variable_types = root_utils.get_column_types(variables_aug,'ntuples/'+base_name+'_sig.root','tree')
  columns = [(variables_aug[i], variable_types[i]) for i in range(len(variables_aug))]
  out_file = open('bdt_evaluate_macro.cxx','w')
  write_event_loop_boilerplate(out_file, base_name, columns, weight, output_suffix, mva_base_name)
  out_file.write('  Float_t BDT = 0;\n')
  out_file.write('\n')
  out_file.write('  //Write signal event loop\n')
  out_file.write('  TFile* sig_out_file = new TFile("ntuples/'+base_name+'_sig_bdtscore.root","RECREATE");\n')
  out_file.write('  TTree * sig_out_tree = new TTree("tree","tree");\n')
  out_file.write('  sig_out_tree->Branch("bdt_score",&BDT,"bdt_score/F");\n')
  out_file.write('  for (long ievt = 0; ievt < sig_tree->GetEntries(); ievt++) {\n')
  out_file.write('    sig_tree->GetEntry(ievt);\n')
  out_file.write('    BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    sig_out_tree->Fill();\n')
  out_file.write('  }\n')
  out_file.write('  sig_out_tree->Write();\n')
  out_file.write('  sig_out_file->Close();\n')
  out_file.write('\n')
  out_file.write('  TFile* bak_out_file = new TFile("ntuples/'+base_name+'_bak_bdtscore.root","RECREATE");\n')
  out_file.write('  TTree * bak_out_tree = new TTree("tree","tree");\n')
  out_file.write('  bak_out_tree->Branch("bdt_score",&BDT,"bdt_score/F");\n')
  out_file.write('  //Write background event loop\n')
  out_file.write('  for (long ievt = 0; ievt < bak_tree->GetEntries(); ievt++) {\n')
  out_file.write('    bak_tree->GetEntry(ievt);\n')
  out_file.write('    BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    bak_out_tree->Fill();\n')
  out_file.write('  }\n')
  out_file.write('  bak_out_tree->Write();\n')
  out_file.write('  bak_out_file->Close();\n')
  out_file.write('}\n')
  out_file.close()
  #run ROOT macro and parse output to get statistics and yields
  process_result = subprocess.run('root -l -q bdt_evaluate_macro.cxx'.split())
  subprocess.run('rm bdt_evaluate_macro.cxx'.split())

def evaluate_bdt(base_name, variables, weight, options, output_suffix=''):
  '''Function that evaluates various metrics related to a BDT
  Note: ASSUMES BLOCK SPLITTING IS USED

  @params
  base_name - filename of training samples, see train_bdt
  variables - list of variables used in BDT
  weight - name of weight branch
  options - BdtOptions with options used for training
  output_suffix - tag for output file

  @returns
  a BdtSummary of performance statistics
  '''
  variables_aug = variables+[weight]
  variable_types = root_utils.get_column_types(variables_aug,'ntuples/'+base_name+'_sig.root','tree')
  columns = [(variables_aug[i], variable_types[i]) for i in range(len(variables_aug))]
  out_file = open('bdt_evaluate_macro.cxx','w')
  write_event_loop_boilerplate(out_file, base_name, columns, weight, output_suffix)
  out_file.write('  //Set up histograms\n')
  out_file.write('  TH1D hist_sig("hist_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_bak("hist_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_train_sig("hist_train_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_train_bak("hist_train_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_tests_sig("hist_tests_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_tests_bak("hist_tests_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  //Set up histograms\n')
  out_file.write('\n')
  out_file.write('  //Perform event loops\n')
  if (options.nTrain_Signal == 0):
    out_file.write('  long ntrain = (sig_tree->GetEntries())/2;\n')
  else:
    out_file.write('  long ntrain = '+str(options.nTrain_Signal)+';\n')
  out_file.write('  for (long ievt = 0; ievt < sig_tree->GetEntries(); ievt++) {\n')
  out_file.write('    sig_tree->GetEntry(ievt);\n')
  out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    if (ievt < ntrain) {\n')
  out_file.write('      hist_train_sig.Fill(BDT,'+weight+');\n')
  out_file.write('    }\n')
  out_file.write('    else {\n')
  out_file.write('      hist_tests_sig.Fill(BDT,'+weight+');\n')
  out_file.write('    }\n')
  out_file.write('  }\n')
  if (options.nTrain_Background == 0):
    out_file.write('  ntrain = (bak_tree->GetEntries())/2;\n')
  else:
    out_file.write('  ntrain = '+str(options.nTrain_Background)+';\n')
  out_file.write('  for (long ievt = 0; ievt < bak_tree->GetEntries(); ievt++) {\n')
  out_file.write('    bak_tree->GetEntry(ievt);\n')
  out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    if (ievt < ntrain) {\n')
  out_file.write('      hist_train_bak.Fill(BDT,'+weight+');\n')
  out_file.write('    }\n')
  out_file.write('    else {\n')
  out_file.write('      hist_tests_bak.Fill(BDT,'+weight+');\n')
  out_file.write('    }\n')
  out_file.write('  }\n')
  out_file.write('\n')
  out_file.write('  //add histograms\n')
  out_file.write('  hist_sig.Add(&hist_train_sig);\n')
  out_file.write('  hist_sig.Add(&hist_tests_sig);\n')
  out_file.write('  hist_bak.Add(&hist_train_bak);\n')
  out_file.write('  hist_bak.Add(&hist_tests_bak);\n')
  out_file.write('\n')
  out_file.write('  //get ouput\n')
  out_file.write('  get_roc_auc(&hist_sig, &hist_bak);\n')
  out_file.write('  get_roc_auc(&hist_train_sig, &hist_train_bak);\n')
  out_file.write('  get_roc_auc(&hist_tests_sig, &hist_tests_bak);\n')
  out_file.write('  std::cout << "Sig KS p: " << hist_tests_sig.KolmogorovTest(&hist_train_sig) << "\\n";\n')
  out_file.write('  std::cout << "Bak KS p: " << hist_tests_bak.KolmogorovTest(&hist_train_bak) << "\\n";\n')
  out_file.write('  std::cout << "Binning optimization for test+train\\n";\n')
  out_file.write('  binning_optimizer(&hist_sig, &hist_bak,4,1.5,1.0);\n')
  out_file.write('  std::cout << "Binning optimization for test+train with 1 throw\\n";\n')
  out_file.write('  binning_optimizer(&hist_sig, &hist_bak,5,1.5,1.0,1);\n')
  out_file.write('  //add factor of 2 since we are using train/test set only\n')
  out_file.write('  std::cout << "Binning optimization for train\\n";\n')
  out_file.write('  binning_optimizer(&hist_train_sig, &hist_train_bak,4,1.5,2.0);\n')
  out_file.write('  std::cout << "Binning optimization for tests\\n";\n')
  out_file.write('  binning_optimizer(&hist_tests_sig, &hist_tests_bak,4,1.5,2.0);\n')
  out_file.write('  std::cout << "Binning optimization for train with 1 throw\\n";\n')
  out_file.write('  binning_optimizer(&hist_train_sig, &hist_train_bak,5,1.5,2.0,1);\n')
  out_file.write('  std::cout << "Binning optimization for tests with 1 throw\\n";\n')
  out_file.write('  binning_optimizer(&hist_tests_sig, &hist_tests_bak,5,1.5,2.0,1);\n')
  out_file.write('  cut_optimizer(&hist_tests_sig, &hist_tests_bak,0.5,false,2.0);\n')
  out_file.write('}\n')
  out_file.close()
  #run ROOT macro and parse output to get statistics and yields
  bdt_summary = BdtSummary()
  process_result = subprocess.run('root -l -q bdt_evaluate_macro.cxx'.split(),capture_output=True)
  print(process_result.stdout)
  print(process_result.stderr)
  evaluate_lines = (process_result.stdout.decode('utf-8')).split('\n')
  auc_counter = 0
  csi_counter = 0
  for line in evaluate_lines:
    if line[:12] == 'ROC AUC is: ':
      if (auc_counter==1):
        bdt_summary.auc_train = float(line[12:])
      elif (auc_counter==2):
        bdt_summary.auc_tests = float(line[12:])
      auc_counter += 1
    if line[:12] == 'ROC CSI is: ':
      if (csi_counter==1):
        bdt_summary.csi_train = float(line[12:])
      elif (csi_counter==2):
        bdt_summary.csi_tests = float(line[12:])
      csi_counter += 1
    if line[:10] == 'Sig KS p: ':
      bdt_summary.sig_ks = float(line[10:])
    if line[:10] == 'Bak KS p: ':
      bdt_summary.bak_ks = float(line[10:])
  #get yields, make datacard, and run through combine
  bin1_yields = get_yields_from_binning_optimizer(evaluate_lines, 1)
  bin4_yields = get_yields_from_binning_optimizer(evaluate_lines, 4, 3)
  bin4_cuts = get_cuts_from_binning_optimizer(evaluate_lines, 4, 3)
  print(bin4_cuts)
  make_simple_datacard.make_simple_datacard('temp_datacard.txt',bin1_yields[0],bin1_yields[1])
  evaluate_lines = ((subprocess.run('./scripts/combine_anyenv.py temp_datacard.txt -M Significance'.split(),capture_output=True)).stdout.decode('utf-8')).split('\n')
  for line in evaluate_lines:
    if line[:14] == 'Significance: ':
      bdt_summary.bn1_signif = float(line[14:])
  subprocess.run('rm temp_datacard.txt'.split())
  make_simple_datacard.make_simple_datacard('temp_datacard.txt',bin4_yields[0],bin4_yields[1])
  evaluate_lines = ((subprocess.run('./scripts/combine_anyenv.py temp_datacard.txt -M Significance'.split(),capture_output=True)).stdout.decode('utf-8')).split('\n')
  for line in evaluate_lines:
    if line[:14] == 'Significance: ':
      bdt_summary.bn4_signif = float(line[14:])
  subprocess.run('rm temp_datacard.txt'.split())
  #make variable plots
  print('Making BDT distributions')
  subprocess.run('rm bdt_evaluate_macro.cxx'.split())
  out_file = open('bdt_evaluate_macro.cxx','w')
  write_event_loop_boilerplate(out_file, base_name, columns, weight, output_suffix)
  out_file.write('  ')
  out_file.write('  //set ranges\n')
  for column in columns:
    if (column[0] != weight):
      out_file.write('  '+root_utils.get_cpp_type(column[1])+' max_'+column[0]+' = -999999999;\n');
      out_file.write('  '+root_utils.get_cpp_type(column[1])+' min_'+column[0]+' = 999999999;\n');
  out_file.write('  long sample_nevts = sig_tree->GetEntries();\n')
  out_file.write('  sample_nevts = sample_nevts<1000 ? sample_nevts : 1000;\n')
  out_file.write('  for (long ievt = 0; ievt < sample_nevts; ievt++) {\n')
  out_file.write('    sig_tree->GetEntry(ievt);\n')
  for column in columns:
    if (column[0] != weight):
      out_file.write('    if ('+column[0]+'>max_'+column[0]+')\n')
      out_file.write('      max_'+column[0]+' = '+column[0]+';\n')
      out_file.write('    if ('+column[0]+'<min_'+column[0]+')\n')
      out_file.write('      min_'+column[0]+' = '+column[0]+';\n')
  out_file.write('  }\n')
  out_file.write('  ')
  out_file.write('  //Set up histograms\n')
  for column in columns:
    if (column[0] != weight):
      for mva_bin in range(4):
        out_file.write('  TH1D hist_sig_'+column[0]+str(mva_bin)+'("hist_sig_'+column[0]+str(mva_bin)+'","'+column[0]+'",30,min_'+column[0]+',max_'+column[0]+');\n')
        out_file.write('  TH1D hist_bak_'+column[0]+str(mva_bin)+'("hist_bak_'+column[0]+str(mva_bin)+'","'+column[0]+'",30,min_'+column[0]+',max_'+column[0]+');\n')
  out_file.write('\n')
  out_file.write('  //Perform event loops\n')
  for category in ['sig','bak']:
    out_file.write('  for (long ievt = 0; ievt < '+category+'_tree->GetEntries(); ievt++) {\n')
    out_file.write('    '+category+'_tree->GetEntry(ievt);\n')
    out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
    for ibin in range(len(bin4_cuts)-1):
      out_file.write('    if (BDT > '+str(bin4_cuts[ibin])+' and BDT < '+str(bin4_cuts[ibin+1])+') {\n')
      for column in columns:
        if (column[0] != weight):
          out_file.write('      hist_'+category+'_'+column[0]+str(ibin)+'.Fill('+column[0]+','+weight+');\n')
      out_file.write('    }\n')
    out_file.write('  }\n')
  out_file.write('\n')
  out_file.write('  //draw histograms and write to file\n')
  mva_colors = ['kRed','kGreen','kBlue','kViolet']
  for category in ['sig','bak']:
    for column in columns:
      if (column[0] != weight):
        out_file.write('  TCanvas can_'+category+'_'+column[0]+';\n')
        for mva_bin in range(4):
          hist_name = 'hist_'+category+'_'+column[0]+str(mva_bin)
          out_file.write('  '+hist_name+'.SetLineColor('+mva_colors[mva_bin]+');\n')
          out_file.write('  '+hist_name+'.Scale(1.0/'+hist_name+'.Integral());\n')
          if (mva_bin==0):
            out_file.write('  '+hist_name+'.GetYaxis()->SetRangeUser(0.0,3.5*'+hist_name+'.GetMaximum());\n')
            out_file.write('  '+hist_name+'.Draw();\n')
          else:
            out_file.write('  '+hist_name+'.Draw("SAME");\n')
        out_file.write('  can_'+category+'_'+column[0]+'.Print("dataset/plots/'+base_name+output_suffix+'_'+category+'_'+column[0]+'.pdf");\n')
  out_file.write('}\n')
  out_file.close()
  subprocess.run('root -l -q bdt_evaluate_macro.cxx'.split())
  subprocess.run('rm bdt_evaluate_macro.cxx'.split())
  return bdt_summary

def evaluate_bdt_new(base_name, variables, weight, options, cut_sr, 
                     output_suffix=''):
  '''Function that evaluates various metrics related to a BDT
  Note: ASSUMES BLOCK SPLITTING IS USED

  @params
  base_name - filename of training samples, see train_bdt
  variables - list of variables used in BDT
  weight - name of weight branch
  options - BdtOptions with options used for training
  cut_sr - selection for signal region (inverted for sideband)
  output_suffix - tag for output file

  @returns
  a BdtSummary of performance statistics
  '''
  #JIT C++ evaluator so we can use it through RDataFrame
  cpp_snippet = 'TMVA::Reader bdt_reader;\n'
  cpp_snippet += f'vector<float> bdt_args({len(variables)},0.0);\n'
  for ivar in range(len(variables)):
    cpp_snippet += f'bdt_reader.AddVariable("{variables[ivar]}",'
    cpp_snippet += f'&bdt_args[{ivar}]);\n'
  cpp_snippet += 'bdt_reader.BookMVA("BDT","dataset/weights/'
  cpp_snippet += f'{base_name}{output_suffix}_BDT.weights.xml"'
  cpp_snippet += ');\n'
  ROOT.gInterpreter.ProcessLine(cpp_snippet)
  cpp_snippet = 'float get_bdt_score(vector<float> args) {\n'
  cpp_snippet += '  for (unsigned iarg = 0; iarg<bdt_args.size(); iarg++) {\n'
  cpp_snippet += '    bdt_args[iarg] = args[iarg];\n'
  cpp_snippet += '  }\n'
  cpp_snippet += '  return bdt_reader.EvaluateMVA("BDT");\n'
  cpp_snippet += '}\n'
  ROOT.gInterpreter.Declare(cpp_snippet)

  #do not multithread to get right rdfentry_
  #get significance
  df_sig = ROOT.RDataFrame('tree',f'ntuples/{base_name}_sig.root')
  df_bak = ROOT.RDataFrame('tree',f'ntuples/{base_name}_bak.root')
  df_dat = ROOT.RDataFrame('tree',f'ntuples/{base_name}_dat.root')
  df_sig = df_sig.Define('bdt_score',
      'get_bdt_score({'+','.join([var for var in variables])+'})')
  df_bak = df_bak.Define('bdt_score',
      'get_bdt_score({'+','.join([var for var in variables])+'})')
  df_dat = df_dat.Define('bdt_score',
      'get_bdt_score({'+','.join([var for var in variables])+'})')
  df_sig_sr = df_sig.Filter(cut_sr)
  df_bak_sr = df_bak.Filter(cut_sr)
  df_bak_cr = df_bak.Filter(f'!({cut_sr})')
  df_dat_cr = df_dat.Filter(f'!({cut_sr})')
  ntrain_sig = options.nTrain_Signal
  ntrain_bak = options.nTrain_Background
  if (ntrain_sig==0):
    ntrain_sig = df_sig_sr.Count().GetValue()//2.0
  if (ntrain_bak==0):
    ntrain_bak = df_bak_sr.Count().GetValue()//2.0
  df_bak_sr_train = df_bak_sr.Filter(f'rdfentry_<={ntrain_bak}')
  df_bak_sr_tests = df_bak_sr.Filter(f'rdfentry_>{ntrain_bak}')
  df_sig_sr_train = df_sig_sr.Filter(f'rdfentry_<={ntrain_sig}')
  df_sig_sr_tests = df_sig_sr.Filter(f'rdfentry_>{ntrain_sig}')
  hist_model = ('','',110,-1.1,1.1)
  hist_bak_sr_train_ptr = df_bak_sr_train.Histo1D(hist_model,'bdt_score',
                                                  weight)
  hist_sig_sr_train_ptr = df_sig_sr_train.Histo1D(hist_model,'bdt_score',
                                                  weight)
  ROOT.gInterpreter.ProcessLine('.L src/distribution_analyzer.cpp+')
  n_bdt_bins = 4
  bin_boundaries_vec = ROOT.binning_optimizer(hist_sig_sr_train_ptr.GetPtr(), 
      hist_bak_sr_train_ptr.GetPtr(),n_bdt_bins,1.5,2.0,0)
  bin_boundaries = [-1.0]
  for ibin in range(bin_boundaries_vec.size()):
    bin_boundaries.append(-1.1+2.2/110.0*float(bin_boundaries_vec.at(ibin)))
  bin_boundaries.append(1.0)
  print(bin_boundaries)
  sig_train_yield_ptrs = []
  sig_tests_yield_ptrs = []
  bak_train_yield_ptrs = []
  bak_tests_yield_ptrs = []
  bak_cr_yield_ptrs = []
  dat_yield_ptrs = []
  for ibin in range(n_bdt_bins):
    loedge = bin_boundaries[ibin]
    hiedge = bin_boundaries[ibin+1]
    cut_string = f'bdt_score>={loedge}&&bdt_score<{hiedge}'
    sig_train_yield_ptrs.append(df_sig_sr_train.Filter(cut_string).Sum(weight))
    sig_tests_yield_ptrs.append(df_sig_sr_tests.Filter(cut_string).Sum(weight))
    bak_train_yield_ptrs.append(df_bak_sr_train.Filter(cut_string).Sum(weight))
    bak_tests_yield_ptrs.append(df_bak_sr_tests.Filter(cut_string).Sum(weight))
    bak_cr_yield_ptrs.append(df_bak_cr.Filter(cut_string).Sum(weight))
    dat_yield_ptrs.append(df_dat_cr.Filter(cut_string).Sum(weight))
  train_signif = 0.0
  tests_signif = 0.0
  tests_adj_signif = 0.0
  for ibin in range(n_bdt_bins):
    train_signif += ((sig_train_yield_ptrs[ibin].GetValue()
        /((bak_train_yield_ptrs[ibin].GetValue())**0.5))**2)
    tests_signif += ((sig_tests_yield_ptrs[ibin].GetValue()
        /((bak_tests_yield_ptrs[ibin].GetValue())**0.5))**2)
    bak_cr_ratio = (dat_yield_ptrs[ibin].GetValue()
                    /(bak_cr_yield_ptrs[ibin].GetValue()))
    print(f'Sideband scale factor: {bak_cr_ratio}')
    tests_adj_signif += ((sig_tests_yield_ptrs[ibin].GetValue()
        /(bak_cr_ratio*bak_tests_yield_ptrs[ibin].GetValue())**0.5)**2)
  train_signif = (train_signif**0.5)*(2.0**0.5)
  tests_signif = (tests_signif**0.5)*(2.0**0.5)
  tests_adj_signif = (tests_adj_signif**0.5)*(2.0**0.5)
  print(f'Train signif: {train_signif}')
  print(f'Tests signif: {tests_signif}')
  print(f'Adjs. signif: {tests_adj_signif}')

def write_event_loop_boilerplate(out_file, base_name, columns, weight, output_suffix, mva_base_name=''):
  '''Helper function for writing event loop boilerplate

  @params
  out_file - file to write output to
  base_name - base filename
  columns - list of tuples of column names and types
  weight - name of weight column
  output_suffix - tag for BDT
  '''
  if (mva_base_name==''):
    mva_base_name = base_name
  out_file.write('#include "distribution_analyzer.hpp"\n')
  out_file.write('\n')
  out_file.write('void bdt_evaluate_macro() {\n')
  out_file.write('  gSystem->AddIncludePath("-Iinc");')
  out_file.write('  gSystem->AddLinkedLibs("-Llib -lSmallPhysUtils");')
  out_file.write('  //Get trees from file\n')
  out_file.write('  TFile in_file_sig("ntuples/'+base_name+'_sig.root","READ");\n')
  out_file.write('  TFile in_file_bak("ntuples/'+base_name+'_bak.root","READ");\n')
  #out_file.write('  TDirectoryFile* dataset = static_cast<TDirectoryFile*>(in_file.Get("dataset"));\n')
  out_file.write('  TTree* sig_tree = static_cast<TTree*>(in_file_sig.Get("tree"));\n')
  out_file.write('  TTree* bak_tree = static_cast<TTree*>(in_file_bak.Get("tree"));\n')
  out_file.write('\n')
  out_file.write('  //Set up branches and BDT\n')
  out_file.write('  TMVA::Reader bdt_reader;\n')
  for column in columns:
    out_file.write('  '+root_utils.get_cpp_type(column[1])+' '+column[0]+' = 0;\n');
    if column[1]=='Float_t' or column[1]=='Int_t':
      out_file.write('  sig_tree->SetBranchAddress("'+column[0]+'", &'+column[0]+');\n');
      out_file.write('  bak_tree->SetBranchAddress("'+column[0]+'", &'+column[0]+');\n');
      if (column[0] != weight):
        out_file.write('  bdt_reader.AddVariable("'+column[0]+'", &'+column[0]+');\n')
  out_file.write('  bdt_reader.BookMVA("BDT","dataset/weights/'+mva_base_name
      +output_suffix+'_BDT.weights.xml");\n')
  out_file.write('\n')

def get_yields_from_binning_optimizer(captured_text, nbins, nth=1):
  '''Helper function to extract yields from binning_optimizer print-out

  @params
  captured_text - stdout from bdt_evaluation_macro
  nbins - number of bins to consider
  nth - int, find nth instance of bins
  '''
  yields_line_index = captured_text.index('With '+str(nbins)+' bins: ')+1
  if (nth>1):
    for i in range(nth-1):
      yields_line_index = captured_text.index('With '+str(nbins)+' bins: ',yields_line_index+1)+1
  yields_line = captured_text[yields_line_index]
  paren_idx = yields_line.index('(')
  sig_yields = []
  bak_yields = []
  for i in range(nbins):
    paren_idx = yields_line.index('(',paren_idx+1)
    comma_idx = yields_line.index(',',paren_idx+1)
    close_idx = yields_line.index(')',paren_idx+1)
    sig_yields.append(float(yields_line[paren_idx+1:comma_idx]))
    bak_yields.append(float(yields_line[comma_idx+1:close_idx]))
  return (sig_yields, bak_yields)

def get_cuts_from_binning_optimizer(captured_text, nbins, nth=1):
  '''Helper function to extract cuts from binning_optimizer print-out

  @params
  captured_text - stdout from bdt_evaluation_macro
  nbins - number of bins to consider
  '''
  yields_line_index = captured_text.index('With '+str(nbins)+' bins: ')+1
  if (nth>1):
    for i in range(nth-1):
      yields_line_index = captured_text.index('With '+str(nbins)+' bins: ',yields_line_index+1)+1
  yields_line = captured_text[yields_line_index]
  comma_idx = yields_line.index(':')
  cuts = []
  for i in range(nbins-1):
    paren_idx = yields_line.index('(',comma_idx+1)
    cuts.append(float(yields_line[comma_idx+2:paren_idx]))
    comma_idx = yields_line.index(',',comma_idx+1)
    comma_idx = yields_line.index(',',comma_idx+1)
  paren_idx = yields_line.index('(',comma_idx+1)
  cuts.append(float(yields_line[comma_idx+2:paren_idx]))
  cuts.append(1.1)
  return cuts

def clean_bdt(base_name, output_suffix=''):
  '''Function to clean up after BDT training

  @params
  base_name - filename of training samples, see train_bdt
  output_suffix - tag for output file
  '''
  subprocess.run(('rm ntuples/output_'+base_name+output_suffix+'.root').split())
  subprocess.run(('rm dataset/weights/'+base_name+output_suffix+'_BDT.class.C').split())
  subprocess.run(('rm dataset/weights/'+base_name+output_suffix+'_BDT.weights.xml').split())
