#!/usr/bin/env python3
'''
Script that generates neural net used for calculating kinematic scale factors
'''

from argparse import ArgumentParser
from array import array
#import mathplotlib.pyplot
import uproot, numpy, pandas, sklearn.preprocessing, collections, os
from tensorflow import math, constant, where, boolean_mask, size, float32, clip_by_value
from tensorflow.keras.backend import epsilon
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.models import Sequential, Model, load_model
from tensorflow.keras.layers import BatchNormalization, Dropout, concatenate, Dense, Activation
from tensorflow.keras.losses import Loss, MeanSquaredError, BinaryCrossentropy
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam

#fix for weird bug that happens sometimes?
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

if __name__ =='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='generate_photon_correctiondnn',
      description=('Generate photon ID corrections from data file (-d) and '
                   +'simulation file (-m)'))
  argument_parser.add_argument('-d','--data_file')
  argument_parser.add_argument('-m','--mc_file')
  argument_parser.add_argument('-t','--tag', default='kin')
  args = argument_parser.parse_args()

  train_columns = ['photon_idmva0','photon_mht_dphi','njet','ll_pt0',
                   'photon_pt0','llphoton_pt0','mht0','ht0']

  #read files
  simu_file = uproot.open(args.mc_file)
  data_file = uproot.open(args.data_file)

  #make datasets using pandas
  all_columns = train_columns+['w_param','rng']
  simu_data_frame = simu_file['tree'].arrays(all_columns, library='pd')
  data_data_frame = data_file['tree'].arrays(all_columns, library='pd')

  #renormalize weights
  w_mc_scale = (data_data_frame['w_param'].sum()
                /simu_data_frame['w_param'].sum())
  simu_data_frame['w_param'] = (simu_data_frame['w_param']*w_mc_scale)

  #split datasets
  simu_train_data = simu_data_frame[simu_data_frame.rng%10>=4]
  data_train_data = data_data_frame[data_data_frame.rng%10>=4]
  simu_valid_data = simu_data_frame[(simu_data_frame.rng%10==2)
                                    |(simu_data_frame.rng%10==3)]
  data_valid_data = data_data_frame[(data_data_frame.rng%10==2)
                                    |(data_data_frame.rng%10==3)]
  #simu_tests_data = simu_data_frame[simu_data_frame.event%10>7]
  #data_tests_data = data_data_frame[data_data_frame.event%10>7]
  #validation done elsewhere
  inp_vars_train = numpy.concatenate(
      (simu_train_data[train_columns].to_numpy(), 
       data_train_data[train_columns].to_numpy()))
  inp_vars_valid = numpy.concatenate(
      (simu_valid_data[train_columns].to_numpy(), 
       data_valid_data[train_columns].to_numpy()))
  out_vars_train = numpy.concatenate(
      (simu_train_data[['is_data']].to_numpy(), 
       data_train_data[['is_data']].to_numpy()))
  out_vars_valid = numpy.concatenate(
      (simu_valid_data[['is_data']].to_numpy(), 
       data_valid_data[['is_data']].to_numpy()))
  wgt_vars_train = numpy.concatenate(
      (simu_train_data['w_param'].to_numpy(), 
       data_train_data['w_param'].to_numpy()))
  wgt_vars_valid = numpy.concatenate(
      (simu_valid_data['w_param'].to_numpy(), 
       data_valid_data['w_param'].to_numpy()))

  #process inputs to be gaussian-like with zero mean and unit norm
  #to save parameters for external application, apply normalization and
  #translation/scaling separately
  print('Cleaning input.')
  scaler_list = []
  
  scaler = sklearn.preprocessing.StandardScaler().fit(inp_vars_train)
  inp_vars_valid = scaler.transform(inp_vars_valid)
  inp_vars_train = scaler.transform(inp_vars_train)
  print('scaler_mean = [',end='')
  for i in range(len(scaler.mean_)):
    if i != 0:
      print(',',end='')
    print(scaler.mean_[i],end='')
    scaler_list.append(scaler.mean_[i])
  print(']')
  print('scaler_vari = [',end='')
  for i in range(len(scaler.var_)):
    if i != 0:
      print(',',end='')
    print(scaler.var_[i],end='')
    scaler_list.append(scaler.var_[i])
  print(']')
  scaler_file = open('json/scaler_file_'+args.tag,'wb')
  scaler_array = array('d',scaler_list)
  print(scaler_array)
  scaler_array.tofile(scaler_file)
  scaler_file.close()

  #build MVA model
  print('Building model.')
  n_input = len(train_columns)
  width = [10,10]
  lr = 0.0001
  dropout = 0.1 #currently commented out
  depth = len(width) 
  model = Sequential()
  model.add(Dense(units=width[0], 
            input_dim=n_input, 
            activation='elu'))
  for i in range(1, depth):
    #model.add(Dropout(dropout))
    model.add(Dense(units=width[i], 
              activation='elu'))
  #model.add(Dropout(dropout))
  model.add(Dense(1, activation='sigmoid'))

  #perform training
  print('Training model.')
  model.compile(loss=BinaryCrossentropy(),
                optimizer=Adam(learning_rate=lr, clipnorm=50.0),
                weighted_metrics=[])
  callbacks = [
      # if no decrease of the loss for several epochs, terminate training.
      EarlyStopping(verbose=True, patience=8, monitor='val_loss'),
      # Make sure that we're saving the model weights with the best val loss.
      ModelCheckpoint('json/nn_model_'+args.tag+'.h5', monitor='val_loss', 
                       verbose=True, save_best_only=True)]
  #shuffles by default
  modelMetricsHistory = model.fit(
      inp_vars_train,
      out_vars_train,
      sample_weight=wgt_vars_train,
      shuffle=True,
      epochs=100,
      batch_size=1024,
      validation_data=(inp_vars_valid,out_vars_valid,wgt_vars_valid),
      callbacks=callbacks,
      verbose=1)
  
