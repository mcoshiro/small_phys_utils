#!/usr/bin/env python3
"""@package docstring
Make plots of b-tagging efficiencies for object review
"""

from root_plot_lib import RplPlot
import ROOT
from correctionlib import CorrectionSet
from argparse import ArgumentParser
from array import array

PT_BINS = [30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,
           350.0,400.0,500.0,600.0,800.0,1000.0]
ETA_BINS = [-2.4,-1.92,-1.44,-0.96,-0.48,0.0,0.48,0.96,1.44,1.92,2.4]
PT_BINS_ARRAY = array('d',PT_BINS)
ETA_BINS_ARRAY = array('d',ETA_BINS)
YEARS = ['2016APV','2016','2017','2018','2022','2022EE','2023','2023BPix']
PICO_FOLDERS = {'2016APV' : '2016preVFP_UL',
                '2016' : '2016postVFP_UL',
                '2017' : '2017_UL',
                '2018' : '2018_UL',
                '2022' : '2022',
                '2022EE' : '2022EE',
                '2023' : '2023',
                '2023BPix' : '2023BPix'}
BTAG_EFF_NAMES = ['Btag_b_WPmedium_MCeff', 'Btag_c_WPmedium_MCeff', 
                  'Btag_uds_WPmedium_MCeff']
FLAV_NAMES = ['b','c','l']
YEAR_LUMI = {'2016APV' : 19.51,
             '2016' : 16.80,
             '2017' : 41.48,
             '2018' : 59.83,
             '2022' : 8.17,
             '2022EE' : 27.01,
             '2023' : 17.61,
             '2023BPix' : 9.53}
YEAR_ENERGY = {'2016APV' : 13,
               '2016' : 13,
               '2017' : 13,
               '2018' : 13,
               '2022' : 13.6,
               '2022EE' : 13.6,
               '2023' : 13.6,
               '2023BPix' : 13.6}

def make_plot_year(year: str):
  """Makes b-tagging efficiency plot for year
  """
  in_file = ('/homes/oshiro/analysis/nano2pico/data/zgamma/'
             +f'{PICO_FOLDERS[year]}/btag_mceff.json')
  correction_set = CorrectionSet.from_file(in_file)
  hists = []
  for iflavor in range(3):
    hists.append(ROOT.TH2D(f'btageff_{FLAV_NAMES[iflavor]}_{year}',
        '; p_{T} [GeV]; #eta; '
        +f'{FLAV_NAMES[iflavor]} jets b-tag WP medium MC efficiency', 
        len(PT_BINS)-1, PT_BINS_ARRAY, len(ETA_BINS)-1, ETA_BINS_ARRAY))
    for ipt in range(len(PT_BINS)-1):
      pt_bin_mean = (PT_BINS[ipt]+PT_BINS[ipt+1])/2.0
      for ieta in range(len(ETA_BINS)-1):
        eta_bin_mean = (ETA_BINS[ieta]+ETA_BINS[ieta+1])/2.0
        eff = correction_set[BTAG_EFF_NAMES[iflavor]].evaluate('effmc',
            eta_bin_mean, pt_bin_mean)
        hists[iflavor].SetBinContent(ipt+1, ieta+1, eff)
    plot = RplPlot()
    plot.lumi_data = [(YEAR_LUMI[year],YEAR_ENERGY[year])]
    if iflavor == 2:
      plot.log_z = True
    plot.z_min = 0.0
    plot.z_max = 1.0
    plot.plot_colormap(hists[iflavor])
    plot.draw(f'plots/btageff_{FLAV_NAMES[iflavor]}_{year}.pdf')

if __name__ == '__main__':
  for year in YEARS:
    make_plot_year(year)
