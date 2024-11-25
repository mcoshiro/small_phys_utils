#!/usr/bin/env python3
"""@package docstring
Make plots related to photon preselection weights
"""

from root_plot_lib import RplPlot
import ROOT
from correctionlib import CorrectionSet
from argparse import ArgumentParser
from array import array

if __name__ == '__main__':

  argument_parser = ArgumentParser(prog='plot_photon_preselection_weights.py',
      description='Generates plots of photon preselection corrections used by H->Zgamma analysis.')
  argument_parser.add_argument('-j','--json_file')
  args = argument_parser.parse_args()

  pt_bins = [15,17.5,20,25,35,50,100]
  eta_bins = [0.0,0.8,1.5,2.0,2.5]
  x_vals = array('d')
  x_uncs = array('d')
  y_vals = [array('d') for ivar2 in range(len(eta_bins)-1)]
  y_unup = [array('d') for ivar2 in range(len(eta_bins)-1)]
  y_undn = [array('d') for ivar2 in range(len(eta_bins)-1)]
  correction_set = CorrectionSet.from_file(args.json_file)
  pt_bins_array = array('d')
  for ipt in range(len(pt_bins)-1):
    pt_bin_mean = (pt_bins[ipt]+pt_bins[ipt+1])/2.0
    pt_bin_width = (pt_bins[ipt+1]-pt_bins[ipt])/2.0
    x_vals.append(pt_bin_mean)
    x_uncs.append(pt_bin_width)
    for ieta in range(len(eta_bins)-1):
      eta_bin_mean = (eta_bins[ieta+1]+eta_bins[ieta])/2.0
      y_val = correction_set['photon_preselection_sf_pass_nom'].evaluate(
          pt_bin_mean, eta_bin_mean)
      y_up = correction_set['photon_preselection_sf_pass_up'].evaluate(
          pt_bin_mean, eta_bin_mean)
      y_dn = correction_set['photon_preselection_sf_pass_down'].evaluate(
          pt_bin_mean, eta_bin_mean)
      y_vals[ieta].append(y_val)
      y_unup[ieta].append(y_up-y_val)
      y_undn[ieta].append(y_val-y_dn)

  graphs = [ROOT.TGraphAsymmErrors(len(x_vals), x_vals, y_vals[ieta], x_uncs, 
                                   x_uncs, y_unup[ieta], y_undn[ieta])
                                   for ieta in range(len(eta_bins)-1)]

  for ieta in range(len(eta_bins)-1):
    graphs[ieta].SetTitle(str(eta_bins[ieta])+'<|#eta|<'+str(eta_bins[ieta+1]))
    graphs[ieta].GetXaxis().SetTitle('p_{T} [GeV]')

  plot = RplPlot()
  plot.lumi_data = [(60,13)]
  plot.y_title = 'Preselection Scale Factor'
  plot.x_min = 15.0
  plot.x_max = 100.0
  plot.log_x = True
  for graph in graphs:
    plot.plot_graph(graph)
  plot.draw('plots/photon_presel_weights.pdf')
