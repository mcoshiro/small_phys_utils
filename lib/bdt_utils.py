"""@package docstring
Library of Python utilities for working with TMVA BDTs

"""
import make_simple_datacard
import os
import ROOT
import root_utils
import subprocess

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
    ret += ':BoostType=AdaBoost:AdaBoostBeta='+str(self.AdaBoostBeta)
    ret += ':UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20'
    return ret

def train_bdt(base_name, variables, weight, options, output_suffix=''):
  '''Function that trains a BDT

  @params
  base_name - filename of training samples, signal should be 
              base_name+_sig.root and background base_name+_bak.root
  variables - list of variable names to consider for BDT training
  weight - name of weight branch in file
  options - BDT options class
  output_suffix - optional tag for output file
  '''
  ROOT.TMVA.Tools.Instance()
  sig_input_file = ROOT.TFile('ntuples/'+base_name+'_sig.root','READ')
  bak_input_file = ROOT.TFile('ntuples/'+base_name+'_bak.root','READ')
  output_file = ROOT.TFile('ntuples/output_'+base_name+output_suffix+'.root','RECREATE') 
  #TODO incorporate options into output filename
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
  cut_s = ROOT.TCut()
  cut_b = ROOT.TCut()
  bdt_loader.PrepareTrainingAndTestTree(cut_s,cut_b,options.ptat_string())
  bdt_factory.BookMethod(bdt_loader,ROOT.TMVA.Types.kBDT,'BDT',options.bm_string())
  bdt_factory.TrainAllMethods()
  bdt_factory.TestAllMethods()
  bdt_factory.EvaluateAllMethods()
  output_file.Close()

def evaluate_bdt(base_name, variables, weight, output_suffix=''):
  '''Function that evaluates various metrics related to a BDT

  @params
  base_name - filename of training samples, see train_bdt
  variables - list of variables used in BDT
  weight - name of weight branch
  output_suffix - tag for output file
  '''
  variables.append(weight)
  variable_types = root_utils.get_column_types(variables,'ntuples/'+base_name+'_sig.root','tree')
  columns = [(variables[i],variable_types[i]) for i in range(len(variables))]
  out_file = open('bdt_evaluate_macro.cxx','w')
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
  out_file.write('  bdt_reader.BookMVA("BDT","dataset/weights/'+base_name
      +output_suffix+'_BDT.weights.xml");\n')
  out_file.write('\n')
  out_file.write('  //Set up histograms\n')
  out_file.write('  TH1D hist_sig("hist_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_bak("hist_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_train_sig("hist_train_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_train_bak("hist_train_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_tests_sig("hist_tests_sig","BDT Output",110,-1.1,1.1);\n')
  out_file.write('  TH1D hist_tests_bak("hist_tests_bak","BDT Output",110,-1.1,1.1);\n')
  out_file.write('\n')
  out_file.write('  //Perform event loops\n')
  out_file.write('  long half_evts = (sig_tree->GetEntries())/2;\n')
  out_file.write('  for (long ievt = 0; ievt < sig_tree->GetEntries(); ievt++) {\n')
  out_file.write('    sig_tree->GetEntry(ievt);\n')
  out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    if (ievt < half_evts) {\n')
  out_file.write('      hist_train_sig.Fill(BDT,'+weight+');\n')
  out_file.write('    }\n')
  out_file.write('    else {\n')
  out_file.write('      hist_tests_sig.Fill(BDT,'+weight+');\n')
  out_file.write('    }\n')
  out_file.write('  }\n')
  out_file.write('  half_evts = (bak_tree->GetEntries())/2;\n')
  out_file.write('  for (long ievt = 0; ievt < bak_tree->GetEntries(); ievt++) {\n')
  out_file.write('    bak_tree->GetEntry(ievt);\n')
  out_file.write('    Float_t BDT = bdt_reader.EvaluateMVA("BDT");\n')
  out_file.write('    if (ievt < half_evts) {\n')
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
  out_file.write('  std::cout << "Signal KS p-value:" << hist_tests_sig.KolmogorovTest(&hist_train_sig) << "\\n";\n')
  out_file.write('  std::cout << "Background KS p-value:" << hist_tests_bak.KolmogorovTest(&hist_train_bak) << "\\n";\n')
  out_file.write('  binning_optimizer(&hist_sig, &hist_bak,5,1.5,138.0/3.0);\n')
  out_file.write('}\n')
  out_file.close()
  #run ROOT macro and parse output to get yields
  evaluate_lines = ((subprocess.run('root -l -q bdt_evaluate_macro.cxx'.split(),capture_output=True)).stdout.decode('utf-8')).split('\n')
  yields_line = evaluate_lines[evaluate_lines.index('With 4 bins: ')+1]
  paren_idx = yields_line.index('(')
  sig_yields = []
  bak_yields = []
  for i in range(4):
    paren_idx = yields_line.index('(',paren_idx+1)
    comma_idx = yields_line.index(',',paren_idx+1)
    close_idx = yields_line.index(')',paren_idx+1)
    sig_yields.append(float(yields_line[paren_idx+1:comma_idx]))
    bak_yields.append(float(yields_line[comma_idx+1:close_idx]))
  make_simple_datacard.make_simple_datacard('temp_datacard.txt',sig_yields,bak_yields)
  subprocess.run('./scripts/combine_anyenv.py temp_datacard.txt -M Significance'.split())
  subprocess.run('rm bdt_evaluate_macro.cxx'.split())
  subprocess.run('rm temp_datacard.txt'.split())


def clean_bdt(base_name, output_suffix=''):
  '''Function to clean up after BDT training

  @params
  base_name - filename of training samples, see train_bdt
  output_suffix - tag for output file
  '''
  subprocess.run(('rm ntuples/output_'+base_name+output_suffix+'.root').split())
  subprocess.run(('rm dataset/weights/'+base_name+output_suffix+'_BDT.class.C').split())
  subprocess.run(('rm dataset/weights/'+base_name+output_suffix+'_BDT.weights.xml').split())
