#!/usr/bin/env python

#uproot documentation and tutorials
#https://masonproffitt.github.io/uproot-tutorial/03-trees/index.html
#https://masonproffitt.github.io/uproot-tutorial/03-trees/index.html
#https://uproot.readthedocs.io/en/latest/basic.html
#https://twiki.cern.ch/twiki/bin/view/Sandbox/MachineLearningTutorialSiimplified
#pandas
#https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html

from array import array
import uproot, numpy, pandas, matplotlib.pyplot, sklearn.preprocessing, collections
from tensorflow.keras.models import Sequential, Model, load_model
from tensorflow.keras.layers import BatchNormalization, Dropout, concatenate, Dense, Activation
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam

def update_vars(values, axes):
  '''Returns updated values of random variables following Gaussian 
  distribution centered at current values. Helper method for mcmc.

  values  list  random variable current values
  axes    list  random variable meta info [[lower, upper, sigma],...]
  '''
  updated_values = []
  for ivar in range(len(axes)):
    axis_range = axes[ivar][1]-axes[ivar][0]
    updated_value = values[ivar]+axes[ivar][2]*numpy.random.normal()
    #wrap around
    while updated_value > axes[ivar][1]:
      updated_value -= axis_range
    while updated_value < axes[ivar][0]:
      updated_value += axis_range
    updated_values.append(updated_value)
  return updated_values

def mcmc(axes, pdf, nsample):
  '''Python implementation of Metropolis-Hastings. May have to port to C(++) 
  if too slow. Returns numpy array of samples

  axes     list    random variable meta info [[lower, upper, sigma],...]
  pdf      method  method that returns value proportional to PDF
  nsample  int     size of sample to generate
  '''
  burn_in = 100
  prev_random_vars = [(axes[ivar][1]-axes[ivar][0])/2.0 for ivar in range(len(axes))]
  prev_pdf = pdf(prev_random_vars)
  isample = 0
  sample = numpy.array([])
  while isample < (burn_in+nsample):
    next_random_vars = update_vars(prev_random_vars, axes)
    next_pdf = pdf(next_random_vars)
    rnd = numpy.random.uniform()
    if (rnd < next_pdf/prev_pdf):
      prev_random_vars = next_random_vars
      prev_pdf = next_pdf
      if (isample >= burn_in):
        sample = numpy.append(sample, numpy.array(next_random_vars), 0)
      isample += 1
      #debug
      if (isample % 250 == 0):
        print(isample)
  return sample

#read file 
#estimation_file = uproot.open('shuffled_esttree_bak.root')
estimation_file = uproot.open('shuffled_dnnstudies_zg.root')
#uniformref_file = uproot.open('uniformref_file.root')

#make datasets using pandas
#column_names = ['plx','ply','plz','nlx','nly','nlz','phx','phy','phz']
#column_mins = [-200, -200, -500, -200, -200, -500, -150, -150, -300]
#column_maxs = [200, 200, 500, 200, 200, 500, 150, 150, 300]
#column_names = ['pt_llg','eta_llg','cosTheta','Phi','costheta','phi']
#column_mins = [0.0,-6.0,-1.0,-3.1416,-1.0,-3.1416]
#column_maxs = [225.0, 6.0, 1.0, 3.1416, 1.0, 3.1416]
additional_columns = ['w_lumi_year']
column_names = ['mllg']
column_mins = [100]
column_maxs = [160]
estimation_data_frame = estimation_file['tree'].arrays(column_names+additional_columns,library='pd')
estimation_data_frame = estimation_data_frame[estimation_data_frame['w_lumi_year']>=0]
#num_sample_uniform = estimation_data_frame.shape[0]
num_sample_uniform = 5000000
uniformref_data_frame = pandas.DataFrame()
for icolumn in range(len(column_names)):
  uniformref_data_frame[column_names[icolumn]] = numpy.random.uniform(column_mins[icolumn],column_maxs[icolumn],(num_sample_uniform))
estimation_data_frame = estimation_data_frame.assign(type_true=1.0)
uniformref_data_frame = uniformref_data_frame.assign(type_true=0.0)
ntrain_stop = int(round(estimation_data_frame.shape[0] * 0.5))
estimation_train_data = estimation_data_frame[:ntrain_stop]
estimation_tests_data = estimation_data_frame[ntrain_stop:]
uniformref_train_data = uniformref_data_frame[:ntrain_stop]
uniformref_tests_data = uniformref_data_frame[ntrain_stop:]
#free to concatenate since Keras' fit method automatically shuffles
inp_vars_all = numpy.concatenate((estimation_data_frame[column_names].to_numpy(), uniformref_data_frame[column_names].to_numpy()))
inp_vars_train = numpy.concatenate((estimation_train_data[column_names].to_numpy(), uniformref_train_data[column_names].to_numpy()))
inp_vars_tests = numpy.concatenate((estimation_tests_data[column_names].to_numpy(), uniformref_tests_data[column_names].to_numpy()))
out_vars_train = numpy.concatenate((estimation_train_data[['type_true']].to_numpy(), uniformref_train_data[['type_true']].to_numpy()))
out_vars_tests = numpy.concatenate((estimation_tests_data[['type_true']].to_numpy(), uniformref_tests_data[['type_true']].to_numpy()))
inp_vars_train_cat0 = uniformref_train_data[column_names].to_numpy()
inp_vars_train_cat1 = estimation_train_data[column_names].to_numpy()
inp_vars_tests_cat0 = uniformref_tests_data[column_names].to_numpy()
inp_vars_tests_cat1 = estimation_tests_data[column_names].to_numpy()

#check data is correct
#matplotlib.pyplot.hist2d(train_data[train_data['type_tru']==0]['x'], train_data[train_data['type_tru']==0]['y'], [10,10], range=[[0,10],[0,10]])
#matplotlib.pyplot.ylabel('y')
#matplotlib.pyplot.xlabel('x')
#matplotlib.pyplot.savefig('./hist_cat0test.pdf')
#matplotlib.pyplot.clf()
#matplotlib.pyplot.hist2d(train_data[train_data['type_tru']==1]['x'], train_data[train_data['type_tru']==1]['y'], [10,10], range=[[0,10],[0,10]])
#matplotlib.pyplot.ylabel('y')
#matplotlib.pyplot.xlabel('x')
#matplotlib.pyplot.savefig('./hist_cat1test.pdf')
#matplotlib.pyplot.clf()

#scale input and get weights for output if categories have different numbers of events
print('Cleaning input.')

scaler = sklearn.preprocessing.StandardScaler().fit(inp_vars_train)
inp_vars_train = scaler.transform(inp_vars_train)
inp_vars_tests = scaler.transform(inp_vars_tests)
inp_vars_train_cat0 = scaler.transform(inp_vars_train_cat0)
inp_vars_tests_cat0 = scaler.transform(inp_vars_tests_cat0)
inp_vars_train_cat1_og = inp_vars_train_cat1
inp_vars_tests_cat1_og = inp_vars_tests_cat1
inp_vars_train_cat1 = scaler.transform(inp_vars_train_cat1)
inp_vars_tests_cat1 = scaler.transform(inp_vars_tests_cat1)
inp_vars_all = scaler.transform(inp_vars_all)
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
scaler_file = open('scaler_file','wb')
scaler_array = array('d',scaler_list)
print(scaler_array)
scaler_array.tofile(scaler_file)
scaler_file.close()

w_categ0 = estimation_train_data.shape[0]/uniformref_train_data.shape[0]
#w_categ0 = 1.0
w_categ1 = 1.0

#matplotlib.pyplot.hist2d(inp_vars_train[::2,0], inp_vars_train[::2,1], [10,10], range=[[-2,2],[-2,2]])
#matplotlib.pyplot.ylabel('y')
#matplotlib.pyplot.xlabel('x')
#matplotlib.pyplot.savefig('./hist_cat0test.pdf')
#matplotlib.pyplot.clf()
#matplotlib.pyplot.hist2d(inp_vars_train[1::2,0], inp_vars_train[1::2,1], [10,10], range=[[-2,2],[-2,2]])
#matplotlib.pyplot.ylabel('y')
#matplotlib.pyplot.xlabel('x')
#matplotlib.pyplot.savefig('./hist_cat1test.pdf')
#matplotlib.pyplot.clf()

#build MVA model
print('Building model.')

n_input = len(column_names)
width = [30,30,30]
#dropout = 0.3
#dropout = 0.15
dropout = 0.1
depth = 3 
model = Sequential()
model.add(Dense(units=width[0], input_dim=n_input, activation='relu'))
model.add(Dropout(dropout))
for i in range(1,depth):
  model.add(Dense(units=width[i],activation='relu'))
  model.add(Dropout(dropout))
model.add(Dense(1, activation='sigmoid'))

#perform training
print('Training model.')

model.compile(loss='binary_crossentropy',optimizer='Adam',metrics=['accuracy'])
model.save('model.h5')
callbacks = [
    # if we don't have a decrease of the loss for 4 epochs, terminate training.
    EarlyStopping (verbose=True, patience=20, monitor='val_loss'),
    # Always make sure that we're saving the model weights with the best val loss.
    ModelCheckpoint ('./nn_model.h5', monitor='val_loss', verbose=True, save_best_only=True)]
modelMetricsHistory = model.fit(
    inp_vars_train,
    out_vars_train,
    class_weight={
    0 : w_categ0,
    1 : w_categ1},
    epochs=400,
    batch_size=2048,
    validation_split=0.2,
    callbacks=callbacks,
    verbose=1)

#model = load_model('nn_model.h5')
#for layer in model.layers:
#  print(layer.get_config(), layer.get_weights())

##evaluate performance
#print('Evaluating perfomance.')
#constant_density = 1.0/(160.0-100.0)
##constant_density = 1.0
#density_graph_x = numpy.linspace(101,159,116)
#density_graph_x_trans = scaler.transform(density_graph_x.reshape((-1,1)))
#density_graph_ratio = (model.predict(density_graph_x_trans, batch_size=2048))
#density_graph_ratio = density_graph_ratio.reshape((density_graph_ratio.shape[0],))
#density_graph_y = numpy.fromfunction( lambda x : (density_graph_ratio[x]*constant_density/(1.0-density_graph_ratio[x])), density_graph_ratio.shape, dtype=int)
#density_graph_y_norm = 0.5*numpy.sum(density_graph_y)
#print('PDF norm: '+str(density_graph_y_norm))
#density_graph_y = numpy.fromfunction( lambda x : (density_graph_y[x]/density_graph_y_norm), density_graph_y.shape, dtype=int)
##density_estimator = lambda x : model.predict(scaler.transform(numpy.array([x])),batch_size=2048)
##mc_sample = mcmc([[100,160,3.0]],density_estimator,10000)
##matplotlib.pyplot.hist(mc_sample, bins=60, range=[100,160], histtype='step', lw=2, alpha=0.5, label=['MC sample'], density=True)
#matplotlib.pyplot.hist(inp_vars_train_cat1_og, bins=60, range=[100,160], histtype='step', lw=2, alpha=0.5, label=['Training data'], density=True)
#matplotlib.pyplot.hist(inp_vars_tests_cat1_og, bins=60, range=[100,160], histtype='step', lw=2, alpha=0.5, label=['Test data'], density=True)
#matplotlib.pyplot.plot(density_graph_x, density_graph_y, color='red', label='Predicted density')
#matplotlib.pyplot.ylabel('Probability/% Entries')
#matplotlib.pyplot.xlabel('Higgs candidate mass [GeV]')
#matplotlib.pyplot.legend(loc='upper right')
#matplotlib.pyplot.xlim(100, 160)
#matplotlib.pyplot.savefig('./hist2.pdf')
#matplotlib.pyplot.clf()

#out_vars_train_cat0_pred = model.predict(inp_vars_train_cat0, batch_size=2048)
#out_vars_train_cat1_pred = model.predict(inp_vars_train_cat1, batch_size=2048)
#out_vars_tests_cat0_pred = model.predict(inp_vars_tests_cat0, batch_size=2048)
#out_vars_tests_cat1_pred = model.predict(inp_vars_tests_cat1, batch_size=2048)
#print('cat0')
#print(model.predict([[0.0, -1.0]], batch_size=2048)) #very cat0 like
#print('cat1')
#print(model.predict([[1.0, 1.0]], batch_size=2048)) #very cat1 like
#matplotlib.pyplot.hist(out_vars_train_cat0_pred, bins=20, range=[0,1], histtype='step', lw=2, alpha=0.5, label=['category 0 training'], density=False)
#matplotlib.pyplot.hist(out_vars_train_cat1_pred, bins=20, range=[0,1], histtype='step', lw=2, alpha=0.5, label=['category 1 training'], density=False)
#matplotlib.pyplot.ylabel('Entries')
#matplotlib.pyplot.xlabel('MVA discriminant')
#matplotlib.pyplot.legend(loc='upper center')
#matplotlib.pyplot.savefig('./hist2.pdf')
#matplotlib.pyplot.clf()
#matplotlib.pyplot.hist(out_vars_tests_cat0_pred, bins=20, range=[0,1], histtype='step', lw=2, alpha=0.5, label=['category 0 validation'], density=True)
#matplotlib.pyplot.hist(out_vars_tests_cat1_pred, bins=20, range=[0,1], histtype='step', lw=2, alpha=0.5, label=['category 1 validation'], density=True)
#matplotlib.pyplot.ylabel('Entries')
#matplotlib.pyplot.xlabel('MVA discriminant')
#matplotlib.pyplot.legend(loc='upper center')
#matplotlib.pyplot.savefig('./hist3.pdf')
#matplotlib.pyplot.clf()

#write output
#out_vars_all = model.predict(inp_vars_all, batch_size=2048)
#data_frame['nn_disc'] = out_vars_all.tolist()
#out_file = uproot.recreate('out_file_with_mldisc.root')
#out_file['out_tree'] = data_frame

#examples
#inp_vars_train_cat0 = train_data[train_data['type_tru']==0][['x','y']].to_numpy()

