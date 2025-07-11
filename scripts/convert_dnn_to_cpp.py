#!/usr/bin/env python

""" Script to convert TensorFlow DNN to C++ for faster single point evaluation 
(ex. for Metropolis-Hastings)

"""

from argparse import ArgumentParser
from array import array
from math import sqrt
from tensorflow.keras.models import load_model

if __name__=='__main__':

  #parse arguments and get set up
  argument_parser = ArgumentParser(prog='convert_dnn_to_cpp',
      description='Converts an h5 Keras model into a header-only C++ library.')
  argument_parser.add_argument('-i','--input_filename',default='model.h5')
  argument_parser.add_argument('-s','--scaler_filename',default='')
  argument_parser.add_argument('-t','--scaler_type',
      choices=['scale','normscale'],default='scale')
  argument_parser.add_argument('-v','--verbose',action='store_true')
  argument_parser.add_argument('-o','--output_filename',default='model.hpp')
  args = argument_parser.parse_args()
  class_name = args.output_filename[:-4]
  model = load_model(args.input_filename, compile=False)

  scaler_mean = []
  scaler_vari = []
  if (args.scaler_filename != ''):
    scaler_file = open(args.scaler_filename, 'rb')
    scaler_array = array('d')
    #scaler_array.fromstring(scaler_file.read())
    size_in_doubles = len(scaler_file.read())//8
    scaler_file.seek(0)
    scaler_array.fromfile(scaler_file,size_in_doubles)
    scaler_file.close()
    scaler_mean = []
    scaler_vari = []
    scaler_lamb = []
    if args.scaler_type=='normscale':
      ncolumns = len(scaler_array)//3
      scaler_lamb = scaler_array[:ncolumns].tolist()
      scaler_mean = scaler_array[ncolumns:2*ncolumns].tolist()
      scaler_vari = scaler_array[2*ncolumns:].tolist()
    else:
      ncolumns = len(scaler_array)//2
      scaler_mean = scaler_array[:ncolumns].tolist()
      scaler_vari = scaler_array[ncolumns:].tolist()
    if (args.verbose):
      print(len(scaler_array))
      print(scaler_array)
      print(scaler_mean)
      print(scaler_vari)
      if args.scaler_type=='normscale':
        print(scaler_lamb)

  #parse Keras model and generate relevant text
  n_neuron_layers = 0
  n_input = 0
  first_neuron_layer = True
  for layer in model.layers:
    if not ('dropout' in layer.get_config()['name']):
      n_neuron_layers += 1
      if first_neuron_layer:
        n_input = layer.get_weights()[0].shape[0]
        first_neuron_layer = False
  initialize_code = ''
  initialize_code += '  n_layer('+str(n_neuron_layers)+'),\n'
  initialize_code += '  n_input('+str(n_input)+'),\n'
  initialize_code += '  scale_lamb{'
  for iinput in range(len(scaler_lamb)):
    if (iinput != 0):
      initialize_code += ','
    initialize_code += str(scaler_lamb[iinput])
  initialize_code += '},\n'
  initialize_code += '  scale_mean{'
  for iinput in range(len(scaler_mean)):
    if (iinput != 0):
      initialize_code += ','
    initialize_code += str(scaler_mean[iinput])
  initialize_code += '},\n'
  initialize_code += '  scale_stdv{'
  for iinput in range(len(scaler_mean)):
    if (iinput != 0):
      initialize_code += ','
    initialize_code += str(sqrt(scaler_vari[iinput]))
  initialize_code += '},\n'
  initialize_code += '  n_unit{'
  first_neuron_layer = True
  for layer in model.layers:
    if not ('dropout' in layer.get_config()['name']):
      if first_neuron_layer:
        first_neuron_layer = False
      else:
        initialize_code += ','
      initialize_code += str(layer.get_config()['units'])
  initialize_code += '},\n'
  initialize_code += '  activation_type{'
  first_neuron_layer = True
  for layer in model.layers:
    if not ('dropout' in layer.get_config()['name']):
      if first_neuron_layer:
        first_neuron_layer = False
      else:
        initialize_code += ','
      initialize_code += 'ActivationType::'+str(layer.get_config()['activation'])
  initialize_code += '},\n'
  initialize_code += '  weight{'
  first_neuron_layer = True
  for layer in model.layers:
    if not ('dropout' in layer.get_config()['name']):
      if first_neuron_layer:
        first_neuron_layer = False
      else:
        initialize_code += ','
      initialize_code += '{'
      weights = layer.get_weights()
      for iunit in range(weights[0].shape[1]):
        if (iunit != 0):
          initialize_code += ','
        initialize_code += '{'
        for iinput in range(weights[0].shape[0]):
          if (iinput != 0):
            initialize_code += ','
          initialize_code += str(weights[0][iinput][iunit])
        initialize_code += '}'
      initialize_code += '}'
  initialize_code += '},\n'
  initialize_code += '  bias{'
  first_neuron_layer = True
  for layer in model.layers:
    if not ('dropout' in layer.get_config()['name']):
      if first_neuron_layer:
        first_neuron_layer = False
      else:
        initialize_code += ','
      initialize_code += '{'
      weights = layer.get_weights()
      for iunit in range(weights[0].shape[1]):
        if (iunit != 0):
          initialize_code += ','
        initialize_code += str(weights[1][iunit])
      initialize_code += '}'
  initialize_code += '} {}\n'
  
  #write output to file
  output_code = ''
  output_code += '//this file was auto-generated by convert_dnn_to_cpp.py\n'
  output_code += '#ifndef H_'+class_name.upper()+'\n'
  output_code += '#define H_'+class_name.upper()+'\n'
  output_code += '\n'
  output_code += '#include <cmath>\n'
  output_code += '#include <vector>\n'
  output_code += '\n'
  output_code += '/*!\\class '+class_name+'\n'
  output_code += ' class containing DNN such as one produced by TensorFlow\n'
  output_code += ' */\n'
  output_code += 'class '+class_name+' {\n'
  output_code += 'public:\n'
  output_code += '  /*!\\brief Standard constrctor for initializing DNN\n'
  output_code += '   */\n'
  output_code += '  '+class_name+'();\n'
  output_code += '\n'
  output_code += '  enum class ActivationType{relu, elu, sigmoid};\n'
  output_code += '\n'
  output_code += '  /*!\\brief evaluates DNN score at a point input\n'
  output_code += '   */\n'
  output_code += '  float evaluate(std::vector<float> input) const;\n'
  output_code += '\n'
  output_code += 'private:\n'
  output_code += '  std::vector<float> scale(std::vector<float> input) const;\n'
  output_code += '  float dot_product(std::vector<float> x, std::vector<float> y) const;\n'
  output_code += '  float relu(float x) const;\n'
  output_code += '  float elu(float x) const;\n'
  output_code += '  float sigmoid(float x) const;\n'
  output_code += '\n'
  output_code += '  unsigned n_layer;\n'
  output_code += '  unsigned n_input;\n'
  output_code += '  std::vector<float> scale_lamb;\n'
  output_code += '  std::vector<float> scale_mean;\n'
  output_code += '  std::vector<float> scale_stdv;\n'
  output_code += '  std::vector<unsigned> n_unit;\n'
  output_code += '  std::vector<'+class_name+'::ActivationType> activation_type;\n'
  output_code += '  std::vector<std::vector<std::vector<float>>> weight;\n'
  output_code += '  std::vector<std::vector<float>> bias;\n'
  output_code += '};\n'
  output_code += '\n'
  output_code += ''+class_name+'::'+class_name+'() :\n'
  output_code += initialize_code
  output_code += '\n'
  output_code += 'float '+class_name+'::evaluate(std::vector<float> input) const {\n'
  output_code += '  std::vector<float> layer_input = scale(input);\n'
  output_code += '  for (unsigned ilayer = 0; ilayer < n_layer; ilayer++) {\n'
  output_code += '    std::vector<float> layer_output;\n'
  output_code += '    layer_output.resize(n_unit[ilayer]);\n'
  output_code += '    for (unsigned iunit = 0; iunit < n_unit[ilayer]; iunit++) {\n'
  output_code += '      float unit_input = dot_product(layer_input,weight[ilayer][iunit])+bias[ilayer][iunit];\n'
  output_code += '      if (activation_type[ilayer]=='+class_name+'::ActivationType::relu)\n'
  output_code += '        layer_output[iunit] = relu(unit_input);\n'
  output_code += '      else if (activation_type[ilayer]=='+class_name+'::ActivationType::sigmoid)\n'
  output_code += '        layer_output[iunit] = sigmoid(unit_input);\n'
  output_code += '      else if (activation_type[ilayer]=='+class_name+'::ActivationType::elu)\n'
  output_code += '        layer_output[iunit] = elu(unit_input);\n'
  output_code += '    }\n'
  output_code += '    layer_input = layer_output;\n'
  output_code += '  }\n'
  output_code += '  return layer_input[0];\n'
  output_code += '}\n'
  output_code += '\n'
  output_code += 'float '+class_name+'::relu(float x) const {\n'
  output_code += '  if (x < 0) return 0;\n'
  output_code += '  return x;\n'
  output_code += '}\n'
  output_code += '\n'
  output_code += 'float '+class_name+'::elu(float x) const {\n'
  output_code += '  if (x < 0) return (exp(x)-1.0);\n'
  output_code += '  return x;\n'
  output_code += '}\n'
  output_code += '\n'
  output_code += 'float '+class_name+'::sigmoid(float x) const {\n'
  output_code += '  return 1.0/(1.0+exp(-1.0*x));\n'
  output_code += '}\n'
  output_code += '\n'
  output_code += 'float '+class_name+'::dot_product(std::vector<float> x, std::vector<float> y) const {\n'
  output_code += '  float output = 0;\n'
  output_code += '  for (unsigned i = 0; i < x.size(); i++) {\n'
  output_code += '    output += x[i]*y[i];\n'
  output_code += '  }\n'
  output_code += '  return output;\n'
  output_code += '}\n'
  output_code += '\n'
  output_code += 'std::vector<float> '+class_name+'::scale(std::vector<float> input) const {\n'
  output_code += '  std::vector<float> output;\n'
  output_code += '  output.resize(n_input);\n'
  output_code += '  for (unsigned iin = 0; iin < n_input; iin++) {\n'
  if args.scaler_type=='normscale':
    output_code += '    if (scale_lamb[iin] > 0 || scale_lamb[iin] < 0)\n'
    output_code += '      output[iin] = (pow(input[iin], scale_lamb[iin])-1)/scale_lamb[iin];\n'
    output_code += '    else\n'
    output_code += '      output[iin] = log(input[iin]);\n'
    output_code += '    output[iin] = (output[iin]-scale_mean[iin])/scale_stdv[iin];\n'
  else:
    output_code += '    output[iin] = (input[iin]-scale_mean[iin])/scale_stdv[iin];\n'
  output_code += '  }\n'
  output_code += '  return output;\n'
  output_code += '}\n'
  output_code += '\n'
  output_code += '#endif'
  output_file = open(args.output_filename, 'w')
  output_file.write(output_code)
  output_file.close()
