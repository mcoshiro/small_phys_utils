"""@package docstring
Script to try other TMVA methods for comparison with BDT

"""

import ROOT
import root_utils

base_name = 'shuffled_dnncomp'
output_suffix = '_pders'
variables = ['photon_mva','pt_mass','cosTheta']
weight = 'w_lumi_year'

if __name__=='__main__':
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
  cut_s = ROOT.TCut()
  cut_b = ROOT.TCut()
  bdt_loader.PrepareTrainingAndTestTree(cut_s,cut_b,'nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:NormMode=NumEvents:!V')
  bdt_factory.BookMethod(bdt_loader,ROOT.TMVA.Types.kPDERS,'PDERS','!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600')
  bdt_factory.TrainAllMethods()
  bdt_factory.TestAllMethods()
  bdt_factory.EvaluateAllMethods()
  output_file.Close()
