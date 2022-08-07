#!/usr/bin/env python3
"""@package docstring
Small utility to generate a simple datacard that can be fed into Higgs combine

Usage: make_simple_datacard.py > datacard.txt

Currently, the signal_yields and background_yields variables in this file must
be edited. The script then generates a combine datacard with the appropriate 
expected yields in each bin. The observed yield is the sum of signal and
background, rounded to the nearest integer. 

One can then get the expected significance for this scenario by running

combine datacard.txt -M Significance
"""

#variables that specify the signal and background yields in each bin
#for now, manually edit these to set datacard data
#TODO: read input from table or similar
signal_yields = [1.0]
background_yields = [1.0]

if __name__=='__main__':
  #preamble
  print('imax '+str(len(signal_yields))+'  number of channels')
  print('jmax 1  number of backgrounds')
  print('kmax *  number of nuisance parameters')
  print('shapes * * FAKE')
  print('')
  
  #print observed yields
  print('                        bin',end='')
  for ibin in range(len(signal_yields)):
    print(('bin'+str(ibin)).rjust(26,' '),end='')
  print('')
  print('                Observation',end='')
  for ibin in range(len(signal_yields)):
    print(str(int(round(signal_yields[ibin]+background_yields[ibin]))).rjust(26,' '),end='')
  print('')
  print('')
  
  #print mc yields
  print('                        bin',end='')
  for ibin in range(len(signal_yields)):
    print(('bin'+str(ibin)).rjust(26,' '),end='')
    print(('bin'+str(ibin)).rjust(26,' '),end='')
  print('')
  print('                    process',end='')
  for ibin in range(len(signal_yields)):
    print(('sig').rjust(26,' '),end='')
    print(('bkg').rjust(26,' '),end='')
  print('')
  print('                    process',end='')
  for ibin in range(len(signal_yields)):
    print(('0').rjust(26,' '),end='')
    print(('1').rjust(26,' '),end='')
  print('')
  print('                       rate',end='')
  for ibin in range(len(signal_yields)):
    print(str(round(signal_yields[ibin],1)).rjust(26,' '),end='')
    print(str(round(background_yields[ibin],1)).rjust(26,' '),end='')
  print('')
  print('')
