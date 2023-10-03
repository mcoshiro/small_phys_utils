#!/usr/bin/env python3 

"""@package docstring 
Script to convert h5 format Keras DNN into a format that TMVA can use
Currently, still retains various hardcoded references for density estimation project

"""

from argparse import ArgumentParser
from array import array
from tensorflow.keras import backend
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.optimizers import SGD
from ROOT import TMVA, TFile, TTree, TCut
from subprocess import run

if __name__=='__main__':
  #parse args
  parser = ArgumentParser(prog='convert_keras_to_tmva',
      description='Converts DNN from Keras H5 to TMVA weights file')
  parser.add_argument('input_filename')
  parser.add_argument('output_filename')
  #parser.add_argument('root_filename')
  #parser.add_argument('target_branchname')
  args = parser.parse_args()

  # Setup TMVA
  TMVA.Tools.Instance()
  TMVA.PyMethodBase.PyInitialize()
  output = TFile.Open(args.output_filename+'.root', 'RECREATE')
  factory = TMVA.Factory('TMVARegression', output,
      '!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression')
  
  #auto-generate dummy TTree
  #hardcoded for now
  tree = TTree('tree','dummy tree for Keras->TMVA conversion')
  mllg = array('f',[0])
  ttru = array('f',[0])
  tree.Branch('mllg',mllg,'mllg/F')
  tree.Branch('ttru',ttru,'ttru/F')
  #TMVA errors if not given at least 10 values to train on
  dummy_values = [0 for i in range(5)]+[1 for i in range(5)]+[0 for i in range(5)]+[1 for i in range(5)]
  for value in dummy_values:
    mllg[0] = value
    ttru[0] = value
    tree.Fill()
  dataloader = TMVA.DataLoader('dataset')
  dataloader.AddVariable('mllg')
  dataloader.AddTarget('ttru')
  dataloader.AddRegressionTree(tree, 1.0)
  dataloader.PrepareTrainingAndTestTree(TCut(''),
          'nTrain_Regression=10:SplitMode=Block:NormMode=NumEvents:!V')
   
  # Load model
  #hack: set training rate to 0 so that the TMVA training won't affect the model
  model = load_model(args.input_filename)
  backend.set_value(model.optimizer.learning_rate, 0.0)
  model.save('keras_model_learn0.h5')
   
  # Book methods and run TMVA
  # You will be very sorry when trying to evaluate the MVA if you don't set Verbose=0
  factory.BookMethod(dataloader, TMVA.Types.kPyKeras, 'PyKeras',
      'H:!V:Verbose=0:FilenameModel=keras_model_learn0.h5:FilenameTrainedModel=trainedModelRegression.h5:NumEpochs=1:BatchSize=32')
  factory.TrainAllMethods()
  factory.TestAllMethods()
  factory.EvaluateAllMethods()
  #run('rm TMVA_Regression_Keras.root'.split())
