#!/usr/bin/env python3
"""@package docstring
Small utility to generate a simple datacard that can be fed into Higgs combine

"""

#TODO: read input from table or similar

def make_simple_datacard(filename,signal_yields,background_yields):
  '''Function that generates a simple combine datacard
  
  The script generates a combine datacard with the appropriate expected
  yields in each bin. The observed yield is the sum of signal and
  background, rounded to the nearest integer. 
  One can then get the expected significance for this scenario by running

  combine datacard.txt -M Significance

  @params
  filename - filename of datacard
  signal_yields - list of signal yields in each bin
  background_yields - list of background yields in each bin
  '''
  out_file = open(filename,'w')
  #preamble
  out_file.write('imax '+str(len(signal_yields))+'  number of channels\n')
  out_file.write('jmax 1  number of backgrounds\n')
  out_file.write('kmax *  number of nuisance parameters\n')
  out_file.write('shapes * * FAKE\n\n')
  
  #out_file.write observed yields
  out_file.write('                        bin')
  for ibin in range(len(signal_yields)):
    out_file.write(('bin'+str(ibin)).rjust(26,' '))
  out_file.write('\n')
  out_file.write('                Observation')
  for ibin in range(len(signal_yields)):
    out_file.write(str(int(round(signal_yields[ibin]+background_yields[ibin]))).rjust(26,' '))
  out_file.write('\n\n')
  
  #out_file.write mc yields
  out_file.write('                        bin')
  for ibin in range(len(signal_yields)):
    out_file.write(('bin'+str(ibin)).rjust(26,' '))
    out_file.write(('bin'+str(ibin)).rjust(26,' '))
  out_file.write('\n')
  out_file.write('                    process')
  for ibin in range(len(signal_yields)):
    out_file.write(('sig').rjust(26,' '))
    out_file.write(('bkg').rjust(26,' '))
  out_file.write('\n')
  out_file.write('                    process')
  for ibin in range(len(signal_yields)):
    out_file.write(('0').rjust(26,' '))
    out_file.write(('1').rjust(26,' '))
  out_file.write('\n')
  out_file.write('                       rate')
  for ibin in range(len(signal_yields)):
    out_file.write(str(round(signal_yields[ibin],1)).rjust(26,' '))
    out_file.write(str(round(background_yields[ibin],1)).rjust(26,' '))
  out_file.write('\n\n')
  out_file.close()
