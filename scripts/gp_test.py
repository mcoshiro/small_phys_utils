#script to test Gaussian Process regression

from array import array
import math
import numpy
import random
import ROOT
from root_plot_lib import RplPlot
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from time import time_ns

#ROOT.gInterpreter.Declare("""
#
#""")

def get_points_from_hist(hist, x_condition):
  """Return lists of x and y points from ROOT TH1D. 
  Condition allows one to provide a lambda, which can be used to blind part of 
  the histogram."""
  x_list = []
  y_list = []
  for ibin in range(1,hist.GetNbinsX()+1):
    x_lo = hist.GetBinLowEdge(ibin)
    x_wd = hist.GetBinWidth(ibin)
    x = x_lo+x_wd/2.0
    y = hist.GetBinContent(ibin)
    if (x_condition(x)):
      x_list.append(x)
      y_list.append(y)
  return x_list, y_list

def list_to_1darray(lst):
  """Convert list into 1D array (list of lists)"""
  return [[item] for item in lst]

def make_tgrapherrors(x_list, y_list, x_err_list, y_err_list):
  """Make TGraphErrors from python lists"""
  vx = array('d', x_list)
  vy = array('d', y_list)
  vex = array('d', x_err_list)
  vey = array('d', y_err_list)
  g = ROOT.TGraphErrors(len(x_list),vx,vy,vex,vey)
  return g

def fill_th1d(h, content):
  """Fill TH1D h with content from python list"""
  for ibin in range(len(content)):
    h.SetBinContent(ibin+1, content[ibin])

def double_bins(h, nbins_original, lo, hi):
  """Returns TH1D with twice as many bins and content filled by 
  interpolation"""
  h2 = ROOT.TH1D('','',nbins_original*2, lo, hi)
  for ibin in range(1,nbins_original+1):
    content1 = h.GetBinContent(ibin)
    content2 = h.GetBinContent(ibin+1)
    if (ibin == nbins_original):
      content2 = content1
    h2.SetBinContent(ibin*2-1, content1)
    h2.SetBinContent(ibin*2, (content1+content2)/2.0)
  return h2

if __name__=='__main__':

  ROOT.RooRandom.randomGenerator().SetSeed(time_ns())

  #make PDF model
  mllg = ROOT.RooRealVar('mllg', 'm_{ll#gamma} [GeV]', 100.0, 160.0)
  erf_mu = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 90.0, 130.0)
  erf_sigma = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
  exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 100.0) 
  erf_mu.setVal(105.0)
  erf_sigma.setVal(0.15)
  exp_lambda.setVal(2.0)
  s_mu = ROOT.RooRealVar('s_mu', 'Signal mean', 125.0, 125.0)
  s_sigma = ROOT.RooRealVar('s_sigma', 'Signal width', 1.8, 1.8) 
  shape_s = ROOT.RooGaussian('shape_s','',mllg,s_mu,s_sigma)
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
      '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-100.0)/60.0)',
      ROOT.RooArgList(mllg, erf_mu, erf_sigma, exp_lambda))
  norm_bgen = ROOT.RooRealVar('norm_bgen','Background norm',0.0,1.0e9)
  norm_sgen = ROOT.RooRealVar('norm_sgen','Signal norm',0.0,1.0) 
  norm_bgen.setVal(0.98)
  norm_sgen.setVal(0.02)
  shape_sb = ROOT.RooAddPdf('pdf_sb','',ROOT.RooArgList(shape_s,pdf_b),ROOT.RooArgList(norm_sgen))

  #make test dataset
  n_samples = 10000
  test_dataset_unbinned = shape_sb.generate(ROOT.RooArgSet(mllg),n_samples)
  test_data = test_dataset_unbinned.createHistogram('mllg',mllg,ROOT.RooFit.Binning(60,100,160))
  x_train, y_train = get_points_from_hist(test_data, lambda mll : (mll < 120 or mll > 130))
  #u_train = numpy.array([2.0*math.sqrt(y) for y in y_train])
  u_train = numpy.array([y for y in y_train]) #variance, not stddev
  x_test, y_test = get_points_from_hist(test_data, lambda x : True)
  x_train_gpr = list_to_1darray(x_train)
  x_test_gpr = list_to_1darray(x_test)

  ##generate dataset
  #n_samples = 10000
  #x_lo = [100.0,105.0,110.0,115.0,120.0]
  #x_hi = [130.0,135.0,140.0,145.0,150.0,155.0,160.0]
  #pdf_lo = [0.004,0.012,0.017,0.017,0.015]
  #pdf_hi = [0.0105,0.009,0.0075,0.0065,0.0055,0.0045,0.004]
  #pdf_norm = sum(pdf_lo)+sum(pdf_hi)
  #y_lo = [random.gauss(ipdf*n_samples/pdf_norm, math.sqrt(ipdf*n_samples/pdf_norm)) for ipdf in pdf_lo]
  #y_hi = [random.gauss(ipdf*n_samples/pdf_norm, math.sqrt(ipdf*n_samples/pdf_norm)) for ipdf in pdf_hi]
  #target = x_lo + [122.5, 125.0, 127.5] + x_hi
  #x_np = numpy.array([[ix] for ix in x_lo+x_hi])
  #y_np = numpy.array(y_lo+y_hi)
  #target_np = numpy.array([[itarget] for itarget in target])

  #do GPR
  #kernel = 1 * RBF(length_scale=0.1, length_scale_bounds=(1e-6, 1e2)) + WhiteKernel()
  kernel = 1 * RBF(length_scale=0.1, length_scale_bounds=(1e-6, 1e2))
  gaussian_process = GaussianProcessRegressor(kernel=kernel, alpha=u_train, n_restarts_optimizer=30)
  gaussian_process.fit(x_train_gpr, y_train)
  fit_test, unc_test = gaussian_process.predict(x_test_gpr,return_std=True)

  #make plot
  #x = x_test
  #y = [fit_test[i] for i in range(len(x_test))]
  #ux = [0.0 for i in x_test]
  #uy = [unc_test[i] for i in range(len(x_test))]
  #x = target
  #y = [y_target[i] for i in range(len(target))]
  #ux = [0.0 for i in target]
  #uy = [target_unc[i] for i in range(len(target))]
  ##x = x_lo
  ##y = y_lo
  ##ux = [0.0 for i in x_lo]
  ##uy = [0.0 for i in y_lo]
  ##x += target
  ##y += [y_target[i] for i in range(len(target))]
  ##ux += [0.0 for i in target]
  ##uy += [target_unc[i] for i in range(len(target))]
  ##x += x_hi
  ##y += y_hi
  ##ux += [0.0 for i in x_hi]
  ##uy += [0.0 for i in y_hi]
  #g = make_tgrapherrors(x, y, ux, uy)
  #c = ROOT.TCanvas()
  #g.SetFillColorAlpha(ROOT.TColor.GetColor('#5790fc'),0.5)
  #g.SetLineColor(ROOT.TColor.GetColor('#5790fc'))
  #g.SetLineWidth(2)
  #g.Draw('AL3')
  #test_data.Draw('SAME')
  #c.SaveAs('plots/gp_test.pdf')

  #do s+b fit
  x_template = [[100.125+0.25*i] for i in range(240)]
  template_th1d = ROOT.TH1D('template_th1d','',240,100.0,160.0)
  y_template = gaussian_process.predict(x_template, return_std=False)
  fill_th1d(template_th1d, y_template)
  template_th1d_2 = double_bins(template_th1d,240,100.0,160.0)
  template_roodatahist = ROOT.RooDataHist('','',ROOT.RooArgList(mllg),template_th1d_2)
  shape_b = ROOT.RooHistPdf('shape_b','',ROOT.RooArgSet(mllg),template_roodatahist)
  norm_b = ROOT.RooRealVar('norm_b','Background norm',0.0,1.0e9)
  norm_s = ROOT.RooRealVar('norm_s','Signal norm',0.0,1.0) 
  norm_b.setVal(0.98)
  norm_s.setVal(0.02)
  shape_sb = ROOT.RooAddPdf('pdf_sb','',ROOT.RooArgList(shape_s,shape_b),ROOT.RooArgList(norm_s))

  shape_sb.fitTo(test_dataset_unbinned)
  plot_postfit = mllg.frame()
  test_dataset_unbinned.plotOn(plot_postfit)
  shape_sb.plotOn(plot_postfit,ROOT.RooFit.Name('postfit_sb'))
  sig_norm_temp = norm_s.getValV()
  bak_norm_temp = 1.0-sig_norm_temp
  norm_s.setVal(0.0)
  shape_sb.plotOn(plot_postfit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
  norm_s.setVal(sig_norm_temp)
  postfit_canvas = ROOT.TCanvas()
  plot_postfit.Draw()
  postfit_canvas.SaveAs('plots/gptest_postfit.pdf')
  print(sig_norm_temp)

