#!/usr/bin/env python3
"""@package docstring
Toys for stats studies
"""

from argparse import ArgumentParser
#from correctionlib import schemav2 
#from tnp_utils import do_tnp_fit, integrate_tgraph
from root_plot_lib import RplPlot
from math import sqrt, hypot
from statistics import fmean, stdev
from time import time_ns
import ROOT
import json
import subprocess
import ctypes

ROOT.gInterpreter.Declare("""

RooHist* cast_tobject_to_roohist(TObject* object)
{ return static_cast<RooHist*>(object); }

RooCurve* cast_tobject_to_roocurve(TObject* object)
{ return static_cast<RooCurve*>(object); }

""")

def integrate_tgraph_down(tgraph, low, high):
  '''
  Simple method to integrate a tgraph between two limits. 
  Assumes the points of the graph are in increasing x-order.
  '''
  xi_prev = 0.0
  yi_prev = 0.0
  xi = ctypes.c_double(0.0)
  yi = ctypes.c_double(0.0)
  integral = 0.0
  for i in range(tgraph.GetN()-4): #jank for this special case
    tgraph.GetPoint(i,xi,yi)
    if (xi.value < xi_prev):
      break
    if (xi.value > low and xi_prev < high): 
      integral += (yi.value+yi_prev)/2.0*(xi.value-xi_prev)
    xi_prev = xi.value
    yi_prev = yi.value
  return integral

def integrate_tgraph_up(tgraph, low, high):
  '''
  Simple method to integrate a tgraph between two limits. 
  Assumes the points of the graph are in increasing x-order.
  '''
  xi_prev = 0.0
  yi_prev = 0.0
  xi = ctypes.c_double(0.0)
  yi = ctypes.c_double(0.0)
  integral = 0.0
  for i in range(tgraph.GetN()-4): #jank
    tgraph.GetPoint(i,xi,yi)
    if (xi.value > xi_prev):
      yi_prev = yi.value
      xi_prev = xi.value
      continue
    #if (xi.value > low and xi_prev < high): 
    if (xi_prev > low and xi.value < high): 
      integral += (yi.value+yi_prev)/2.0*(xi_prev-xi.value)
    xi_prev = xi.value
    yi_prev = yi.value
  return integral

if __name__ == '__main__':

  ROOT.RooRandom.randomGenerator().SetSeed(time_ns())

  workspace = ROOT.RooWorkspace()
  mll = ROOT.RooRealVar('mll', 'm_{ll#gamma} [GeV]', 100.0, 160.0)
  mll.setRange('left', 100.0, 120.0)
  mll.setRange('right', 130.0, 160.0)

  gauss_mu = ROOT.RooRealVar('gauss_mu', 'Z peak Gaussian mean', 125.0, 125.0) 
  gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Z peak Gaussian width', 0.01, 15.0) 
  gauss_norm = ROOT.RooRealVar('sig_norm', 'Z peak normalization', 0.0, 1000000000.0) 
  erf_mu = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 90.0, 130.0)
  erf_sigma = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
  exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 100.0) 
  bak_norm = ROOT.RooRealVar('bak_norm', 'Nonresonant normalization', 0.0, 1000000000.0) 
  gauss_mu.setVal(125.0)
  gauss_sigma.setVal(3.0)
  #gauss_norm.setVal(0.5)
  gauss_norm.setVal(0.0)
  erf_mu.setVal(105.0)
  erf_sigma.setVal(0.15)
  exp_lambda.setVal(2.0)
  bak_norm.setVal(50.0)
  #getattr(workspace,'import')(gauss_mu)
  #getattr(workspace,'import')(gauss_sigma)
  #getattr(workspace,'import')(gauss_norm)
  #getattr(workspace,'import')(erf_mu)
  #getattr(workspace,'import')(erf_sigma)
  #getattr(workspace,'import')(exp_lambda)
  #getattr(workspace,'import')(bak_norm)
  pdf_s  = ROOT.RooGaussian('pdf_s','pdf_s', mll, gauss_mu, gauss_sigma)
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
      '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-100.0)/60.0)',
      ROOT.RooArgList(mll, erf_mu, erf_sigma, exp_lambda))

  erf_mu_ref = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 90.0, 130.0)
  erf_sigma_ref = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
  exp_lambda_ref = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 100.0) 
  erf_mu_ref.setVal(105.0)
  erf_sigma_ref.setVal(0.15)
  exp_lambda_ref.setVal(2.0)
  pdf_b_ref = ROOT.RooGenericPdf('pdf_b','pdf_b',
      '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-100.0)/60.0)',
      ROOT.RooArgList(mll, erf_mu_ref, erf_sigma_ref, exp_lambda_ref))

  #pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(gauss_norm, bak_norm))
  pdf_sb = ROOT.RooRealSumPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(gauss_norm, bak_norm), True)

  plot_truth = mll.frame()
  pdf_sb.plotOn(plot_truth,ROOT.RooFit.Name('truth_sb'))
  truth_canvas = ROOT.TCanvas()
  plot_truth.Draw()
  truth_canvas.SaveAs('plots/statstest_truth.pdf')

  plot_index = 0
  #fits_per_point = 4000
  fits_per_point = 40
  samples = [50]
  sb_uncs = [[] for i in range(len(samples))]
  bs_uncs = [[] for i in range(len(samples))]
  print_all = True
  print_index = (99,1)

  for i in range(len(samples)):
    for j in range(fits_per_point):
      gauss_mu.setVal(125.0)
      gauss_sigma.setVal(3.0)
      gauss_norm.setVal(0.0)
      erf_mu.setVal(105.0)
      erf_sigma.setVal(0.15)
      exp_lambda.setVal(2.0)
      bak_norm.setVal(50.0)

      test_dataset_unbinned = pdf_b.generate(ROOT.RooArgSet(mll),samples[i])
      ref_hist = ROOT.TH1D('ref_hist','',30,100,160)
      test_dataset = ROOT.RooDataHist('test_dataset','Toys',ROOT.RooArgList(mll),ref_hist)
      test_dataset.add(test_dataset_unbinned)

      #blinded_dataset = getattr(test_dataset,'reduce')('mll<120||mll>130')

      fit_result = pdf_b.fitTo(test_dataset, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True))
      #print(fit_result)
      if (fit_result.covQual() != 3) or (fit_result.status() != 0):
        continue

      dummy = ''

      if print_all:
        #print('Signal yield: ',end='')
        #print(gauss_norm.getValV(),end='')
        #print('+-')
        #print(gauss_norm.getError())
        #print('Background yield: ',end='')
        #print(bak_norm.getValV(),end='')
        #print('+-')
        #print(bak_norm.getError())
        plot_postfit_b = mll.frame()
        test_dataset.plotOn(plot_postfit_b)
        pdf_b_ref.plotOn(plot_postfit_b,ROOT.RooFit.Name('truth'),ROOT.RooFit.LineColor(ROOT.TColor.GetColor('#5790fc')))
        pdf_b.plotOn(plot_postfit_b,ROOT.RooFit.Name('postfit_b'),ROOT.RooFit.VisualizeError(fit_result,1,False),ROOT.RooFit.FillColor(ROOT.TColor.GetColor('#f89c20')))
        postfit_b_canvas = ROOT.TCanvas()
        plot_postfit_b.Draw()
        postfit_b_canvas.SaveAs('plots/statstest_bfit.pdf')
        truth_curve = ROOT.cast_tobject_to_roocurve(postfit_b_canvas.FindObject('truth'))
        fit_curve = ROOT.cast_tobject_to_roocurve(postfit_b_canvas.FindObject('postfit_b'))
        dn_yield = integrate_tgraph_down(fit_curve, 122.5, 127.5)
        up_yield = integrate_tgraph_up(fit_curve, 122.5, 127.5)
        true_yield = integrate_tgraph_down(truth_curve, 122.5, 127.5)
        nom_yield = (up_yield+dn_yield)/2.0
        print('True yield: '+str(true_yield))
        print('Fit yield: '+str(nom_yield))
        print('Difference (%): '+str((true_yield-nom_yield)/nom_yield))
        print('Stat uncertainty (%): '+str(sqrt(nom_yield)/nom_yield))
        print('Syst uncertainty (%): '+str((up_yield-nom_yield)/nom_yield))
        dummy = input()

      if dummy != 's':
        sb_uncs[i].append(gauss_norm.getError()*gauss_norm.getValV())
      #if sb_uncs[i][-1] > 10.0:
      #  sb_uncs[i].pop()

      #gauss_mu.setVal(125.0)
      #gauss_sigma.setVal(3.0)
      #gauss_norm.setVal(0.0)
      #erf_mu.setVal(105.0)
      #erf_sigma.setVal(0.15)
      #exp_lambda.setVal(2.0)
      #bak_norm.setVal(50.0)
      #gauss_norm.setConstant(True)
      #gauss_sigma.setConstant(True)
      #gauss_mu.setConstant(True)

      #fit_result = pdf_b.fitTo(blinded_dataset, ROOT.RooFit.Extended(True), ROOT.RooFit.SumCoefRange('left,right'), ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save(True))
      ##fit_result2 = pdf_sb.fitTo(blinded_dataset, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True), ROOT.RooFit.Range('left,right'))
      ##fit_result2 = pdf_sb.fitTo(blinded_dataset, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True))
      ##print(fit_result2)
      ##if (fit_result2.covQual() != 3): #or (fit_result.status() != 0):
      ##  continue
      #if (fit_result.covQual() != 3) or (fit_result.status() != 0):
      #  continue

      #gauss_norm.setVal(10.0)
      #erf_mu.setConstant(True)
      #erf_sigma.setConstant(True)
      #exp_lambda.setConstant(True)
      #bak_norm.setConstant(True)
      #gauss_norm.setConstant(False)
      #gauss_sigma.setConstant(False)
      #gauss_mu.setConstant(False)

      #fit_result = pdf_sb.fitTo(test_dataset, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True))
      ##print(fit_result)
      #if (fit_result.covQual() != 3) or (fit_result.status() != 0):
      #  continue

      #if print_all:
      #  print('Signal yield: ',end='')
      #  print(gauss_norm.getValV(),end='')
      #  print('+-')
      #  print(gauss_norm.getError())
      #  print('Background yield: ',end='')
      #  print(bak_norm.getValV(),end='')
      #  print('+-')
      #  print(bak_norm.getError())
      #  plot_prefit_sb = mll.frame()
      #  test_dataset.plotOn(plot_prefit_sb)
      #  pdf_sb.plotOn(plot_prefit_sb,ROOT.RooFit.Name('prefit_sb'))
      #  prefit_sb_canvas = ROOT.TCanvas()
      #  plot_prefit_sb.Draw()
      #  prefit_sb_canvas.SaveAs('plots/statstest_bsfit.pdf')
      #  dummy = input()

      #if dummy != 's':
      #  bs_uncs[i].append(gauss_norm.getError()*gauss_norm.getValV())
      ##if bs_uncs[i][-1] > 10.0:
      ##  bs_uncs[i].pop()

      #if i==print_index[0] and j==print_index[1]:
      #  print('Uncs:')
      #  print(sb_uncs[i][j])
      #  print(bs_uncs[i][j])
      #  print('S and B norms:')
      #  print(gauss_norm.getValV())
      #  print(bak_norm.getValV())
      #  plot_postfit_sb = mll.frame()
      #  test_dataset.plotOn(plot_postfit_sb)
      #  pdf_sb.plotOn(plot_postfit_sb,ROOT.RooFit.Name('postfit_sb'))
      #  postfit_sb_canvas = ROOT.TCanvas()
      #  plot_postfit_sb.Draw()
      #  postfit_sb_canvas.SaveAs('plots/statstest_postfit_sb.pdf')
      #  dummy = input()

  #for i in range(len(samples)):
  #  print('N samples: ',end='')
  #  print(samples[i])
  #  print('Simultaneous S+B: ',end='')
  #  print(fmean(sb_uncs[i]),end='')
  #  print('+-',end='')
  #  print(stdev(sb_uncs[i])/sqrt(len(sb_uncs[i])))
  #  print('Fixed B, then S: ',end='')
  #  print(fmean(bs_uncs[i]),end='')
  #  print('+-',end='')
  #  print(stdev(bs_uncs[i])/sqrt(len(bs_uncs[i])))

  #print(sb_uncs)
