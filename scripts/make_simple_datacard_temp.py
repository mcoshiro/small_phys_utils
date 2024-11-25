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
#signal_yields = [1.0]
#background_yields = [1.0]

#this block until the next empty line is custom, do not commit
#run 2 version
#signal_yields = [59.5,102.2,42.7,10.7]
#background_yields = [79740.2,24899.9,4133.5,278.8]
#photon pt in BDT
#signal_yields = [91.8,92.6,24.4,6.4]
#background_yields = [91974.1,15534.2,1426.49,117.6]
#signal_yields = [62.8,82.2,34.2,8.5]
#background_yields = [51047.7,11792.8,1693.6,97.7]
#photon pt split
#signal_yields = [75.9,97.4,22.6,19.3]
#background_yields = [83955.1,18638.6,5150.1,678.6]
#photon pt split, more bins
signal_yields = [46.0,73.0,54.3,11.0,21.6,9.3]
background_yields = [72206.3,22589.5,7797.9,4130.2,1480.2,218.4]

#summed
#signal_yields = [signal_yields[0]+signal_yields[1]+signal_yields[2]+signal_yields[3]]
#background_yields = [background_yields[0]+background_yields[1]+background_yields[2]+background_yields[3]]

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
