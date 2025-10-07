#!/usr/bin/env python3
"""@package docstring
Converts electron reco weights from root to correctionlib JSON
"""

from argparse import ArgumentParser
from correctionlib import schemav2, CorrectionSet
import json

#constants

ETA_BINS = [-2.5, -2.0, -1.566, -1.444, -1.0, 0.0, 1.0, 1.444, 1.566, 2.0, 2.5]
PT_BINS = [10.0, 20.0, 45.0, 75.0, 100.0, 500.0]

def fix_correctionlib_json(json_texts):
  '''Fixes the format of correctionlib json created using corr.json, since 
  it is not properly formatted by default
  '''
  corr = []
  for json_text in json_texts:
    corr.append(json.loads(json_text))
  json_dict = {
    'schema_version' : 2,
    'description' : '',
    'corrections' : corr
  }
  return json.dumps(json_dict,indent=2)

def make_correction() -> schemav2.Correction:
  '''Generates correction object with fake photon corrections

  Returns:
    schemva2 correction object
  '''
  pt_bins = [15.0,20.0,30.0,500.0]
  abseta_bins = [0.0,0.8,1.5,2.0,2.5]
  run2_jetph_sfs = schemav2.MultiBinning(
      nodetype='multibinning',
      inputs=['pt','abseta'],
      edges=[pt_bins, abseta_bins],
      content=[0.87784536,0.86090972,0.83359961,1.00905351,
               1.20656524,1.15353438,1.28262424,1.15435585,
               0.52378337,0.87839728,0.63154299,1.13529129],
      flow='clamp',
      )
  run2_puph_sfs = schemav2.MultiBinning(
      nodetype='multibinning',
      inputs=['pt','abseta'],
      edges=[pt_bins, abseta_bins],
      content=[1.38054607,1.29143011,1.22872195,1.26953528,
               1.11197443,1.34174906,1.19494895,1.25048567,
               2.59609818,1.86582932,3.16422829,1.41518427],
      flow='clamp',
      )
  run3_jetph_sfs = schemav2.MultiBinning(
      nodetype='multibinning',
      inputs=['pt','abseta'],
      edges=[pt_bins, abseta_bins],
      content=[0.73009776,1.02800363,0.76673460,1.29045474,
               1.03924426,1.12586378,0.90955539,1.10089699,
               0.74353182,0.34798856,0.79232260,0.69352779],
      flow='clamp',
      )
  run3_puph_sfs = schemav2.MultiBinning(
      nodetype='multibinning',
      inputs=['pt','abseta'],
      edges=[pt_bins, abseta_bins],
      content=[1.26284615,1.17956303,1.28810042,1.39744805,
               1.45429152,1.37924462,1.47740148,1.47429645,
               2.11270871,3.20168632,1.77614373,1.49915009],
      flow='clamp',
      )

  return schemav2.Correction(
    name='fakephoton_corrections',
    version=1,
    inputs=[schemav2.Variable(name='run', type='string', description='run2 or run3'),
            schemav2.Variable(name='isjet', type='int', description='isjet'),
            schemav2.Variable(name='pt', type='real', description='pt'),
            schemav2.Variable(name='abseta', type='real', description='eta'),
            ],
    output=schemav2.Variable(name='sf', type='real', description='fake photon SF'),
    data=schemav2.Category(
      nodetype='category',
      input='run',
      content=[
        schemav2.CategoryItem(
          key='run2',
          value=schemav2.Category(
            nodetype='category',
            input='isjet',
            content=[
              schemav2.CategoryItem(
                key=1,
                value=run2_jetph_sfs
              ),
              schemav2.CategoryItem(
                key=0,
                value=run2_puph_sfs
              ),
            ]
          ),
        ),
        schemav2.CategoryItem(
          key='run3',
          value=schemav2.Category(
            nodetype='category',
            input='isjet',
            content=[
              schemav2.CategoryItem(
                key=1,
                value=run3_jetph_sfs
              ),
              schemav2.CategoryItem(
                key=0,
                value=run3_puph_sfs
              ),
            ]
          ),
        ),
        ]
      ),
    )

if __name__ == '__main__':

  #write correction
  corr = make_correction()
  with open('json/fakephoton.json','w') as output_file:
    output_file.write(fix_correctionlib_json([corr.json(exclude_unset=True)]))

  #check corrections
  ceval = CorrectionSet.from_file('json/fakephoton.json')
  print(ceval['fakephoton_corrections'].evaluate('run2',1,17.5,0.4))
  print(ceval['fakephoton_corrections'].evaluate('run2',1,17.5,1.0))
  print(ceval['fakephoton_corrections'].evaluate('run2',1,17.5,1.8))
  print(ceval['fakephoton_corrections'].evaluate('run2',1,17.5,2.25))
  print(ceval['fakephoton_corrections'].evaluate('run2',0,25.0,0.4))
  print(ceval['fakephoton_corrections'].evaluate('run2',0,25.0,1.0))
  print(ceval['fakephoton_corrections'].evaluate('run2',0,25.0,1.8))
  print(ceval['fakephoton_corrections'].evaluate('run2',0,25.0,2.25))
  print(ceval['fakephoton_corrections'].evaluate('run3',0,50.0,0.4))
  print(ceval['fakephoton_corrections'].evaluate('run3',0,50.0,1.0))
  print(ceval['fakephoton_corrections'].evaluate('run3',0,50.0,1.8))
  print(ceval['fakephoton_corrections'].evaluate('run3',0,50.0,2.25))
