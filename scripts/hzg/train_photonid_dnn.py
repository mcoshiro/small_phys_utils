#!/usr/bin/env python3
'''
Script that generates neural net used for calculating photon scale factors
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
from tensorflow.keras.losses import Loss
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam

#fix for weird bug that happens sometimes?
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

#We use the hack where both the weight and true value (1=data, 0=MC) are encoded in y_true
#This is done by defining y_true = weight + (is_MC ? 500 : 0)
class weighted_bce_loss(Loss):
  def __init__(self):
    super().__init__()

  def call(self, y_true, y_pred):
    actual_true = where(math.less(y_true,400.0),constant(1.0),constant(0.0))
    weight = where(math.less(y_true,400.0),y_true,y_true-500.0)
    eps = epsilon()
    y_pred = clip_by_value(y_pred, eps, 1.0-eps)
    bce = actual_true*math.log(y_pred+eps)
    bce += (1.0-actual_true)*math.log(1.0-y_pred+eps)
    return math.reduce_mean(-1.0*weight*bce)

if __name__ =='__main__':
  #parse arguments
  argument_parser = ArgumentParser(prog='generate_photon_correctiondnn',
      description='Generate photon ID corrections from data file (-d) and simulation file (-m)')
  argument_parser.add_argument('-d','--data_file',default='/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_data_2018.root')
  argument_parser.add_argument('-m','--mc_file',default='/net/cms26/cms26r0/oshiro/tnp_tuples/photonidskim_simu_2018.root')
  argument_parser.add_argument('-r','--region',default='barrel',choices=['barrel','endcap'])
  argument_parser.add_argument('-t','--tag',default='barrel2018')
  args = argument_parser.parse_args()

  train_columns = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_rawEnergy','ph_sc_eta']+['ph_et','ph_sc_abseta','event_rho']
  if args.region=='endcap':
    #train_columns = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_eta','ph_ESsigma','ph_esEnergyOverRawE']+['ph_et','ph_sc_abseta','event_rho']
    train_columns = ['ph_r9','ph_s4','ph_sc_etaWidth','ph_sc_phiWidth','ph_sieie','ph_sieip','ph_phoIso','ph_chIso','ph_chWorIso','ph_sc_eta','ph_esEnergyOverRawE']+['ph_et','ph_sc_abseta','event_rho','ph_ESsigma']

  #read files
  simu_file = uproot.open(args.mc_file)
  data_file = uproot.open(args.data_file)

  #make datasets using pandas
  simu_data_frame = simu_file['tree'].arrays(train_columns+['w_pre','event'],library='pd')
  data_data_frame = data_file['tree'].arrays(train_columns+['w_bkg','event'],library='pd')
  w_mc_scale = data_data_frame['w_bkg'].sum()/simu_data_frame['w_pre'].sum()
  print('MC scale factor: '+str(w_mc_scale))
  simu_data_frame['w'] = (simu_data_frame['w_pre']*w_mc_scale+500.0)
  data_data_frame['w'] = data_data_frame['w_bkg']
  if args.region=='barrel':
    simu_data_frame = simu_data_frame[simu_data_frame.ph_sc_abseta<1.5]
    data_data_frame = data_data_frame[data_data_frame.ph_sc_abseta<1.5]
  else:
    simu_data_frame = simu_data_frame[simu_data_frame.ph_sc_abseta>1.5]
    data_data_frame = data_data_frame[data_data_frame.ph_sc_abseta>1.5]
  #simu_ntrain = int(round(simu_data_frame.shape[0] * 0.6))
  #data_ntrain = int(round(data_data_frame.shape[0] * 0.6))
  #simu_nvalid = int(round(simu_data_frame.shape[0] * 0.7))
  #data_nvalid = int(round(data_data_frame.shape[0] * 0.7))
  simu_train_data = simu_data_frame[simu_data_frame.event%10<=5]
  data_train_data = data_data_frame[data_data_frame.event%10<=5]
  simu_valid_data = simu_data_frame[(simu_data_frame.event%10==6)|(simu_data_frame.event%10==7)]
  data_valid_data = data_data_frame[(data_data_frame.event%10==6)|(data_data_frame.event%10==7)]
  #simu_tests_data = simu_data_frame[simu_nvalid:]
  #data_tests_data = data_data_frame[data_nvalid:]
  #for i in range(5):
  #  print('--simu '+str(i))
  #  print(simu_train_data.iloc[i:i+1])
  #  print('--data '+str(i))
  #  print(data_train_data.iloc[i:i+1])
  #exit()

  inp_vars_train = numpy.concatenate((simu_train_data[train_columns].to_numpy(), data_train_data[train_columns].to_numpy()))
  inp_vars_valid = numpy.concatenate((simu_valid_data[train_columns].to_numpy(), data_valid_data[train_columns].to_numpy()))
  out_vars_train = numpy.concatenate((simu_train_data[['w']].to_numpy(), data_train_data[['w']].to_numpy()))
  out_vars_valid = numpy.concatenate((simu_valid_data[['w']].to_numpy(), data_valid_data[['w']].to_numpy()))

  #scale input to have zero mean and unit norm
  print('Cleaning input.')
  scaler = sklearn.preprocessing.StandardScaler().fit(inp_vars_train)
  inp_vars_train = scaler.transform(inp_vars_train)
  scaler_list = []
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
  width = [15,15,15]
  dropout = 0.017
  depth = len(width) 
  model = Sequential()
  model.add(Dense(units=width[0], 
            input_dim=n_input, 
            #kernel_regularizer=L1L2(0.0,0.001),
            activation='elu'))
  #model.add(Dropout(dropout))
  for i in range(1,depth):
    model.add(Dense(units=width[i], 
              #kernel_regularizer=L1L2(0.0,0.001),
              activation='elu'))
  model.add(Dropout(dropout))
  model.add(Dense(1, activation='sigmoid'))

  #perform training
  print('Training model.')
  model.compile(loss=weighted_bce_loss(),optimizer=Adam(learning_rate=0.0001,clipnorm=10.0),metrics=[])
  callbacks = [
      # if we don't have a decrease of the loss for 7 epochs, terminate training.
      EarlyStopping (verbose=True, patience=7, monitor='val_loss'),
      # Always make sure that we're saving the model weights with the best val loss.
      ModelCheckpoint ('json/nn_model_'+args.tag+'.h5', monitor='val_loss', verbose=True, save_best_only=True)]
  #shuffles by default
  modelMetricsHistory = model.fit(
      inp_vars_train,
      out_vars_train,
      #class_weight={
      #0 : w_categ0,
      #1 : w_categ1},
      epochs=100,
      batch_size=1024,
      validation_data=(inp_vars_valid,out_vars_valid),
      callbacks=callbacks,
      verbose=1)
  
