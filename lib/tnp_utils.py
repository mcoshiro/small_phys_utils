#!/usr/bin/env python3
"""@package docstring
Package containing utilities useful for tag-and-probe analyses
"""

import ROOT
import math
import ctypes

ROOT.gInterpreter.Declare("""

RooHist* cast_tobject_to_roohist(TObject* object)
{ return static_cast<RooHist*>(object); }

RooCurve* cast_tobject_to_roocurve(TObject* object)
{ return static_cast<RooCurve*>(object); }

""")

rng = ROOT.TRandom3()

def integrate_tgraph(tgraph, low, high):
  '''
  Simple method to integrate a tgraph between two limits. 
  Assumes the points of the graph are in increasing x-order.
  '''
  xi_prev = 0.0
  yi_prev = 0.0
  xi = ctypes.c_double(0.0)
  yi = ctypes.c_double(0.0)
  integral = 0.0
  for i in range(tgraph.GetN()):
    tgraph.GetPoint(i,xi,yi)
    if (xi.value < xi_prev):
      raise ValueError('integrate method requires monotonic TGraph')
    if (xi.value > low and xi_prev < high): 
      integral += (yi.value+yi_prev)/2.0*(xi.value-xi_prev)
    xi_prev = xi.value
    yi_prev = yi.value
  return integral

def get_maximum_bin(hist, x_min, x_max):
  '''
  Find the bin with maximum content in a subrange
  '''
  x_min_bin = hist.FindBin(x_min)
  x_max_bin = hist.FindBin(x_max)
  max_bin = 0
  max_content = 0
  for ibin in range(x_min_bin,x_max_bin+1):
    bin_content = hist.GetBinContent(ibin)
    if bin_content > max_content:
      max_bin = ibin
      max_content = bin_content
  return max_bin

def get_value_tgraph(graph, x):
  '''Estimates the value of a TGraph at a point x via linear interpolation 
  between the nearest points

  graph  TGraph or similar (ex. RooCurve)
  x      point to get value at
  '''
  npoints = graph.GetN();
  xi = ctypes.c_double(0.0);
  yi = ctypes.c_double(0.0);
  graph.GetPoint(0,xi,yi)
  if (x < xi.value):
    raise ValueError('Requested x '+str(x)+' is less than graph minimum '+str(xi))
  graph.GetPoint(npoints-1,xi,yi)
  if (x > xi.value):
    raise ValueError('Requested x '+str(x)+' is greater than graph maximum '+str(xi))
  if (x == xi.value):
    return (yi.value+0.0)
  idx = 0
  graph.GetPoint(idx,xi,yi)
  while (xi.value <= x):
    graph.GetPoint(idx,xi,yi)
    idx += 1
  xi_prev = ctypes.c_double(0.0)
  yi_prev = ctypes.c_double(0.0)
  graph.GetPoint(idx-2,xi_prev,yi_prev)
  return (x-xi_prev.value)/(xi.value-xi_prev.value)*(yi.value-yi_prev.value)+yi_prev.value

def fit_quality(canvas, hist_name, curve_name, nbins, x_min, x_max):
  '''
  Alternative to chi-squared, KS, etc. for evaluating goodness of fit
  where the metric is the integrated absolute difference between the normalized curves

  canvas      TCanvas, canvas on which the RooPlot is drawn
  hist_name   string, name of the RooHist to compare
  curve_name  string, name of the RooCurve to compare
  nbins       int, number of bins for numeric integration
  x_min       float, lower range to consider
  x_max       float, upper range to consider
  '''
  data_roohist = ROOT.cast_tobject_to_roohist(canvas.FindObject(hist_name))
  fit_roocurve = ROOT.cast_tobject_to_roocurve(canvas.FindObject(curve_name))
  data_integral = integrate_tgraph(data_roohist,x_min,x_max)
  fit_integral = integrate_tgraph(fit_roocurve,x_min,x_max)
  #left here
  prev_diff = 0
  integral_diff = 0
  step_size = (x_max-x_min)/float(nbins)
  for ibin in range(nbins):
    x = x_min +step_size*float(ibin)
    diff = abs(get_value_tgraph(data_roohist, x)/data_integral-get_value_tgraph(fit_roocurve, x)/fit_integral)
    integral_diff += (diff+prev_diff)/2.0*step_size
    prev_diff = diff
  return integral_diff

def fit_cmsshape_cb(hist, return_type='pdf', verbose=False, output_dir='plots/'):
  '''Fit a CMSshape (erf*exp) plus crystal ball distribution to the given TH1
  params
  hist         TH1D, histogram to fit
  return_type  string, 'pdf' or 'curve'
  verbose      bool, print/save extra information
  output_dir   string, location ot save output

  returns 1 if the fit fails, or the RooAbsPdf ('pdf') or a tuple of the S+B 
  and B fits RooCruves ('curve') if the fit succeeds
  '''
  #check arguments
  if not return_type in ['pdf','curve']:
    raise ValueError('Invalid argument to fit_cmsshape_cb, must be "pdf" or "curve".')


  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  mll = ROOT.RooRealVar('mll', 'm_{ll} [GeV]', 50.0, 130.0)
  data = ROOT.RooDataHist('data','Di"photon" invariant mass', ROOT.RooArgList(mll), hist)
  getattr(workspace,'import')(mll)
  getattr(workspace,'import')(data)
  x_name = 'mll'

  gauss_mu = ROOT.RooRealVar('gauss_mu', 'Z peak Gaussian mean', 85.0, 95.0) 
  gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Z peak Gaussian width', 0.01, 15.0) 
  gauss_norm = ROOT.RooRealVar('sig_norm', 'Z peak normalization', 0.0, 100000.0) 
  cb_alphal = ROOT.RooRealVar('cb_alphal', 'Z peak CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('cb_nl', 'Z peak CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('cb_alphar', 'Z peak CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('cb_nr', 'Z peak CB right power', 0.1, 10.0) 
  erf_mu = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 30.0, 85.0)
  erf_sigma = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
  exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 10.0) 
  bak_norm = ROOT.RooRealVar('bak_norm', 'Nonresonant normalization', 0.0, 100000.0) 
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(gauss_norm)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(erf_mu)
  getattr(workspace,'import')(erf_sigma)
  getattr(workspace,'import')(exp_lambda)
  getattr(workspace,'import')(bak_norm)
  param_names = ['gauss_mu','gauss_sigma','sig_norm','cb_alphal','cb_nl',
                 'cb_alphar','cb_nr','erf_mu','erf_sigma','exp_lambda',
                 'bak_norm']

  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', mll, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
      '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-50.0)/80.0)',
         ROOT.RooArgList(mll, erf_mu, erf_sigma, exp_lambda))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(gauss_norm, bak_norm))
  getattr(workspace,'import')(pdf_sb)
  
  #pre-guess parameters to coax a good fit
  #make this an input in order to do other resonances (and same for sigma/edges)od fit
  gauss_sigma_initial = 3.0
  cb_alpha_initial = 1.8
  cb_n_initial = 2.0
  range_x = 80.0
  lowbound_x = 50.0
  gauss_mu_initial = hist.GetBinLowEdge(get_maximum_bin(hist,81.0,101.0))
  ref1_x = gauss_mu_initial-10.0
  ref2_x = gauss_mu_initial+10.0
  ref1_y = hist.GetBinContent(hist.FindBin(ref1_x))
  ref2_y = hist.GetBinContent(hist.FindBin(ref2_x))
  subpeak_value = (ref1_y+ref2_y)/2.0
  gauss_norm_initial = (hist.GetBinContent(hist.FindBin((ref2_x+ref1_x)/2.0))-subpeak_value)*math.sqrt(2.0*3.1416)*gauss_sigma_initial
  #determine pre-set values for exponential by solving system of equations
  #ref1_y = ref1_y-abs(ref1_y-ref2_y)*0.5
  ref1_x = 75.0
  ref1_y = hist.GetBinContent(hist.FindBin(ref1_x))
  exp_lambda_initial = (math.log(ref1_y)-math.log(ref2_y))/((ref2_x-ref1_x)/range_x)
  #temporary bak_norm, will re-set taking into automatic RooFit Normalization of pdf_b below
  bak_norm_initial = math.exp((math.log(ref1_y)*(ref2_x-lowbound_x)/range_x-math.log(ref2_y)*(ref1_x-lowbound_x)/range_x)/((ref2_x-ref1_x)/range_x))
  #solve another system of equations to get error function parameters
  background_peak_x = hist.GetBinLowEdge(get_maximum_bin(hist,50.0,75.0))
  ref1_x = 51.0
  ref2_x = 55.0
  if (background_peak_x > 55):
    ref1_x = background_peak_x-5
    ref2_x = background_peak_x-1
  ref1_y = hist.GetBinContent(hist.FindBin(ref1_x))
  ref2_y = hist.GetBinContent(hist.FindBin(ref2_x))
  ref1_expyield = bak_norm_initial*math.exp(-1.0*exp_lambda_initial*(ref1_x-lowbound_x)/range_x)
  ref2_expyield = bak_norm_initial*math.exp(-1.0*exp_lambda_initial*(ref2_x-lowbound_x)/range_x)
  bak_lowpedestal_initial = ref1_y/bak_norm_initial/10.0
  ref1_erfarg = ROOT.TMath.ErfInverse((ref1_y/ref1_expyield)*2.0-1.0)
  ref2_erfarg = ROOT.TMath.ErfInverse((ref2_y/ref2_expyield)*2.0-1.0)
  erf_mu_initial = (ref1_erfarg*ref2_x-ref2_erfarg*ref1_x)/(ref1_erfarg-ref2_erfarg+1.0)
  erf_sigma_initial = (ref1_erfarg-ref2_erfarg)/(ref1_x-ref2_x)
  #fix bak_norm
  bak_norm_initial  = hist.Integral()-gauss_norm_initial

  gauss_mu.setVal(gauss_mu_initial)
  gauss_sigma.setVal(gauss_sigma_initial)
  gauss_norm.setVal(gauss_norm_initial)
  cb_alphal.setVal(cb_alpha_initial)
  cb_nl.setVal(cb_n_initial)
  cb_alphar.setVal(cb_alpha_initial)
  cb_nr.setVal(cb_n_initial)
  erf_mu.setVal(erf_mu_initial)
  erf_sigma.setVal(erf_sigma_initial)
  exp_lambda.setVal(exp_lambda_initial)
  bak_norm.setVal(bak_norm_initial)

  workspace.saveSnapshot('initial_guess',','.join(param_names))

  max_tries = 10
  best_try = -1
  best_quality = 1.0
  jitter = 0.02
  quality_threshold = 0.06

  #try automatic fit several times
  #for itry in range(max_tries):

  #  workspace.loadSnapshot('initial_guess')
  #  quality = perform_automatic_fit(workspace, x_name, param_names, jitter, 
  #                                  True, 
  #                                  'plots/photon_corr/'+hist.GetName()+'_try'+str(itry))

  #  if (quality<best_quality):
  #    best_quality = quality
  #    best_try = itry

  #  if (quality>quality_threshold):
  #    print('Try '+str(itry)+' poor fit quality ('+str(quality)+'), refitting')
  #  else:
  #    plot_postfit = workspace.var(x_name).frame()
  #    data.plotOn(plot_postfit)
  #    pdf_sb.plotOn(plot_postfit,ROOT.RooFit.Name('postfit_sb'))
  #    sig_norm_temp = workspace.var('sig_norm').getValV()
  #    bak_norm_temp = workspace.var('bak_norm').getValV()
  #    workspace.var('sig_norm').setVal(0.0)
  #    pdf_sb.plotOn(plot_postfit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
  #    workspace.var('sig_norm').setVal(sig_norm_temp)
  #    postfit_canvas = ROOT.TCanvas()
  #    plot_postfit.Draw()

  #    postfit_canvas.SaveAs('plots/photon_corr/'+hist.GetName()+'_postfit.pdf')

  #    postfitsb = ROOT.cast_tobject_to_roocurve(postfit_canvas.FindObject('postfit_sb'))
  #    postfitb = ROOT.cast_tobject_to_roocurve(postfit_canvas.FindObject('postfit_b'))
  #    if (verbose): 
  #      print('fit quality: '+str(quality))
  #      print('Estimated Nsig: '+str(gauss_norm.getValV()))
  #      print('Estimated Nbak: '+str(bak_norm.getValV()))
  #      print('Entries in data: '+str(hist.Integral()))
  #    if return_type == 'curve':
  #      return (0,postfitsb, postfitb)
  #    else:
  #      return (0,pdf_sb)

  #try manual fit
  workspace.loadSnapshot('initial_guess')
  quality = perform_interactive_fit(workspace, x_name, param_names)
  if (quality<best_quality):
    best_quality = quality
    best_try = -2
  if (quality<quality_threshold):
    plot_postfit = workspace.var(x_name).frame()
    data.plotOn(plot_postfit)
    workspace.pdf('pdf_sb').plotOn(plot_postfit,ROOT.RooFit.Name('postfit_sb'))
    sig_norm_temp = workspace.var('sig_norm').getValV()
    bak_norm_temp = workspace.var('bak_norm').getValV()
    workspace.var('sig_norm').setVal(0.0)
    workspace.pdf('pdf_sb').plotOn(plot_postfit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
    workspace.var('sig_norm').setVal(sig_norm_temp)
    postfit_canvas = ROOT.TCanvas()
    plot_postfit.Draw()

    postfit_canvas.SaveAs('plots/photon_corr/'+hist.GetName()+'_postfit.pdf')

    postfitsb = ROOT.cast_tobject_to_roocurve(postfit_canvas.FindObject('postfit_sb'))
    postfitb = ROOT.cast_tobject_to_roocurve(postfit_canvas.FindObject('postfit_b'))
    if (verbose): 
      print('fit quality: '+str(quality))
      print('Estimated Nsig: '+str(workspace.var('sig_norm').getValV()))
      print('Estimated Nbak: '+str(workspace.var('bak_norm').getValV()))
      print('Entries in data: '+str(hist.Integral()))
    if return_type == 'curve':
      return (0,postfitsb, postfitb)
    else:
      return (0,pdf_sb)

  print('Fit failed, best automatic try: '+str(best_try)+', quality: '+str(best_quality))
  return (1,pdf_sb)

def do_tnp_fit(hist, function_type, output_name, return_type='pdf', verbose=False):
  '''Fit a CMSshape (erf*exp) plus crystal ball distribution to the given TH1
  params
  hist           TH1D, histogram to fit
  function_type  string, functions to use in fit - see below for supported
  return_type    string, 'pdf' or 'curve'
  verbose        bool, print/save extra information
  output_dir     string, location ot save output

  Supported functional forms: 'CB+CMS' 'CB+logisticexp' 'Gauss+CMS'

  returns 1 if the fit fails, or the RooAbsPdf ('pdf') or a tuple of the S+B 
  and B fits RooCruves ('curve') if the fit succeeds
  '''
  #check arguments
  if not return_type in ['pdf','curve']:
    raise ValueError('Invalid argument to fit_cmsshape_cb, must be "pdf" or "curve".')

  if not function_type in ['CB+CMS','CB+logisticexp','Gauss+CMS']:
    raise ValueError('Invalid function_type argument to do_tnp_fit.')

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  mll = ROOT.RooRealVar('mll', 'm_{ll} [GeV]', 50.0, 130.0)
  data = ROOT.RooDataHist('data','Di"photon" invariant mass', ROOT.RooArgList(mll), hist)
  getattr(workspace,'import')(mll)
  getattr(workspace,'import')(data)
  x_name = 'mll'
  param_names = []

  if function_type=='CB+CMS':
    gauss_mu = ROOT.RooRealVar('gauss_mu', 'Z peak Gaussian mean', 85.0, 95.0) 
    gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Z peak Gaussian width', 0.01, 15.0) 
    gauss_norm = ROOT.RooRealVar('sig_norm', 'Z peak normalization', 0.0, 1000000000.0) 
    cb_alphal = ROOT.RooRealVar('cb_alphal', 'Z peak CB left switchover', 0.1, 10.0) 
    cb_nl = ROOT.RooRealVar('cb_nl', 'Z peak CB left power', 0.1, 10.0) 
    cb_alphar = ROOT.RooRealVar('cb_alphar', 'Z peak CB right switchover', 0.1, 10.0) 
    cb_nr = ROOT.RooRealVar('cb_nr', 'Z peak CB right power', 0.1, 10.0) 
    erf_mu = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 30.0, 100.0)
    erf_sigma = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
    exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 100.0) 
    bak_pedestal = ROOT.RooRealVar('bak_pedestal', 'Background pedestal', -1.0, 1.0) 
    bak_norm = ROOT.RooRealVar('bak_norm', 'Nonresonant normalization', 0.0, 1000000000.0) 
    gauss_sigma.setVal(3.0)
    gauss_norm.setVal(50000.0)
    erf_mu.setVal(65.0)
    erf_sigma.setVal(0.2)
    exp_lambda.setVal(4.0)
    bak_norm.setVal(50000.0)
    getattr(workspace,'import')(gauss_mu)
    getattr(workspace,'import')(gauss_sigma)
    getattr(workspace,'import')(gauss_norm)
    getattr(workspace,'import')(cb_alphal)
    getattr(workspace,'import')(cb_nl)
    getattr(workspace,'import')(cb_alphar)
    getattr(workspace,'import')(cb_nr)
    getattr(workspace,'import')(erf_mu)
    getattr(workspace,'import')(erf_sigma)
    getattr(workspace,'import')(exp_lambda)
    getattr(workspace,'import')(bak_pedestal)
    getattr(workspace,'import')(bak_norm)
    param_names = ['gauss_mu','gauss_sigma','sig_norm','cb_alphal','cb_nl',
                   'cb_alphar','cb_nr','erf_mu','erf_sigma','exp_lambda',
                   'bak_pedestal','bak_norm']

    pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', mll, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
    pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
        '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-50.0)/80.0)+@4',
        ROOT.RooArgList(mll, erf_mu, erf_sigma, exp_lambda, bak_pedestal))
    pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(gauss_norm, bak_norm))
    getattr(workspace,'import')(pdf_sb)

  elif function_type=='Gauss+CMS':
    gauss_mu = ROOT.RooRealVar('gauss_mu', 'Z peak Gaussian mean', 85.0, 95.0) 
    gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Z peak Gaussian width', 0.01, 15.0) 
    gauss_norm = ROOT.RooRealVar('sig_norm', 'Z peak normalization', 0.0, 1000000000.0) 
    erf_mu = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 30.0, 100.0)
    erf_sigma = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
    exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 100.0) 
    bak_pedestal = ROOT.RooRealVar('bak_pedestal', 'Background pedestal', -1.0, 1.0) 
    bak_norm = ROOT.RooRealVar('bak_norm', 'Nonresonant normalization', 0.0, 1000000000.0) 
    gauss_sigma.setVal(3.0)
    gauss_norm.setVal(50000.0)
    erf_mu.setVal(65.0)
    erf_sigma.setVal(0.2)
    exp_lambda.setVal(4.0)
    bak_norm.setVal(50000.0)
    getattr(workspace,'import')(gauss_mu)
    getattr(workspace,'import')(gauss_sigma)
    getattr(workspace,'import')(gauss_norm)
    getattr(workspace,'import')(erf_mu)
    getattr(workspace,'import')(erf_sigma)
    getattr(workspace,'import')(exp_lambda)
    getattr(workspace,'import')(bak_pedestal)
    getattr(workspace,'import')(bak_norm)
    param_names = ['gauss_mu','gauss_sigma','sig_norm','erf_mu','erf_sigma',
                   'exp_lambda','bak_pedestal','bak_norm']

    pdf_s  = ROOT.RooGaussian('pdf_s','pdf_s', mll, gauss_mu, gauss_sigma)
    pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
        '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-50.0)/80.0)+@4',
           ROOT.RooArgList(mll, erf_mu, erf_sigma, exp_lambda, bak_pedestal))
    pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(gauss_norm, bak_norm))
    getattr(workspace,'import')(pdf_sb)

  elif function_type=='CB+logisticexp':
    gauss_mu = ROOT.RooRealVar('gauss_mu', 'Z peak Gaussian mean', 85.0, 95.0) 
    gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Z peak Gaussian width', 0.01, 15.0) 
    gauss_norm = ROOT.RooRealVar('sig_norm', 'Z peak normalization', 0.0, 1000000000.0) 
    cb_alphal = ROOT.RooRealVar('cb_alphal', 'Z peak CB left switchover', 0.1, 10.0) 
    cb_nl = ROOT.RooRealVar('cb_nl', 'Z peak CB left power', 0.1, 10.0) 
    cb_alphar = ROOT.RooRealVar('cb_alphar', 'Z peak CB right switchover', 0.1, 10.0) 
    cb_nr = ROOT.RooRealVar('cb_nr', 'Z peak CB right power', 0.1, 10.0) 
    logi_offset = ROOT.RooRealVar('logi_offset', 'Logistic curve offset', 30.0, 100.0) 
    logi_width = ROOT.RooRealVar('logi_width', 'Logistic curve width', 0.005, 30.0) 
    exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 100.0) 
    bak_pedestal = ROOT.RooRealVar('bak_pedestal', 'Background pedestal', -1.0, 1.0) 
    bak_norm = ROOT.RooRealVar('bak_norm', 'Nonresonant normalization', 0.0, 1000000000.0) 
    gauss_sigma.setVal(3.0)
    gauss_norm.setVal(50000.0)
    logi_offset.setVal(65.0)
    logi_width.setVal(0.2)
    exp_lambda.setVal(4.0)
    bak_norm.setVal(50000.0)
    getattr(workspace,'import')(gauss_mu)
    getattr(workspace,'import')(gauss_sigma)
    getattr(workspace,'import')(gauss_norm)
    getattr(workspace,'import')(cb_alphal)
    getattr(workspace,'import')(cb_nl)
    getattr(workspace,'import')(cb_alphar)
    getattr(workspace,'import')(cb_nr)
    getattr(workspace,'import')(logi_offset)
    getattr(workspace,'import')(logi_width)
    getattr(workspace,'import')(exp_lambda)
    getattr(workspace,'import')(bak_pedestal)
    getattr(workspace,'import')(bak_norm)
    param_names = ['gauss_mu','gauss_sigma','sig_norm','cb_alphal','cb_nl',
                   'cb_alphar','cb_nr','logi_offset','logi_width','exp_lambda',
                   'bak_pedestal','bak_norm']

    pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', mll, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
    pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
        'exp(-1.0*@3*(@0-50.0)/80.0)/(1.0+exp(-1.0*@2*(@0-@1)))+@4',
           ROOT.RooArgList(mll, logi_offset, logi_width, exp_lambda, bak_pedestal))
    pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(gauss_norm, bak_norm))
    getattr(workspace,'import')(pdf_sb)

  quality_threshold = 0.06

  #manual fit
  quality = perform_interactive_fit(workspace, x_name, param_names)
  #if (quality<quality_threshold):
  plot_postfit = workspace.var(x_name).frame()
  data.plotOn(plot_postfit)
  workspace.pdf('pdf_sb').plotOn(plot_postfit,ROOT.RooFit.Name('postfit_sb'))
  sig_norm_temp = workspace.var('sig_norm').getValV()
  bak_norm_temp = workspace.var('bak_norm').getValV()
  workspace.var('sig_norm').setVal(0.0)
  workspace.pdf('pdf_sb').plotOn(plot_postfit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
  workspace.var('sig_norm').setVal(sig_norm_temp)
  postfit_canvas = ROOT.TCanvas()
  plot_postfit.Draw()

  postfit_canvas.SaveAs(output_name)

  postfitsb = ROOT.cast_tobject_to_roocurve(postfit_canvas.FindObject('postfit_sb'))
  postfitb = ROOT.cast_tobject_to_roocurve(postfit_canvas.FindObject('postfit_b'))
  if (verbose): 
    print('fit quality: '+str(quality))
    print('Estimated Nsig: '+str(workspace.var('sig_norm').getValV()))
    print('Estimated Nbak: '+str(workspace.var('bak_norm').getValV()))
    print('Entries in data: '+str(hist.Integral()))
  if return_type == 'curve':
    return (0,postfitsb, postfitb)
  else:
    return (0,pdf_sb)

  #print('Fit failed, best automatic try: '+str(best_try)+', quality: '+str(best_quality))
  #return (1,pdf_sb)

def perform_automatic_fit(workspace, x_name, param_names, noise=0.0, verbose=False, output_name=''):
  '''Method to perform automated fit and quality check. Returns fit quality
  This method will overwrite workspace parameters, so save a snapshot beforehand
  if needed.

  workspace    RooWorkspace with data to fit and functional form of parameters
               must include a RooAbsHist data and RooAbsPdf pdf_sb with parameters 
               sig_norm and bak_norm
  x_name       name of fit variable
  param_names  list of strings of names of fit parameters
  noise        float how much random noise to add to initial parameters
  verbose      bool for debugging purposes
  output_name  name to use with verbose output plots
  '''
  initial_values = []
  for param_name in param_names:
    initial_values.append(workspace.var(param_name).getValV())
  if noise>0.0:
    for iparam in range(len(param_names)):
      workspace.var(param_names[iparam]).setVal(initial_values[iparam]*rng.Gaus(1.0,noise))

  data = workspace.data('data')
  pdf_sb = workspace.pdf('pdf_sb')

  #prefit
  plot_prefit = workspace.var(x_name).frame()
  data.plotOn(plot_prefit)
  pdf_sb.plotOn(plot_prefit,ROOT.RooFit.Name('prefit_sb'))
  sig_norm_temp = workspace.var('sig_norm').getValV()
  bak_norm_temp = workspace.var('bak_norm').getValV()
  workspace.var('sig_norm').setVal(0.0)
  pdf_sb.plotOn(plot_prefit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('prefit_b'))
  workspace.var('sig_norm').setVal(sig_norm_temp)
  prefit_canvas = ROOT.TCanvas()
  plot_prefit.Draw()

  pdf_sb.fitTo(data)

  #postfit
  plot_postfit = workspace.var(x_name).frame()
  data.plotOn(plot_postfit)
  pdf_sb.plotOn(plot_postfit,ROOT.RooFit.Name('postfit_sb'))
  sig_norm_temp = workspace.var('sig_norm').getValV()
  bak_norm_temp = workspace.var('bak_norm').getValV()
  workspace.var('sig_norm').setVal(0.0)
  pdf_sb.plotOn(plot_postfit,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
  workspace.var('sig_norm').setVal(sig_norm_temp)
  postfit_canvas = ROOT.TCanvas()
  plot_postfit.Draw()

  if (verbose): 
    prefit_canvas.SaveAs(output_name+'_prefit.pdf')
    postfit_canvas.SaveAs(output_name+'_postfit.pdf')

  #evaluate goodness of fit
  return fit_quality(postfit_canvas, 'h_data', 'postfit_sb', 78, 51.0, 129.0)


def perform_interactive_fit(workspace, x_name, param_names):
  '''Method that can be used to allow user to tune parameters via shell for fit
  This method will overwrite workspace parameters, so save a snapshot beforehand
  if needed.

  workspace    RooWorkspace with data to fit and functional form of parameters
               must include a RooAbsHist data and RooAbsPdf pdf_sb
  x_name       name of fit variable
  param_names  list of strings of names of fit parameters
  '''
  quit = False
  have_fit = False
  canvas = ROOT.TCanvas()
  data = workspace.data('data')
  pdf_sb = workspace.pdf('pdf_sb')
  quality = 1.0
  while not quit:
    user_input = input()
    user_input = user_input.split()
    if len(user_input)<1:
      continue
    if user_input[0]=='list' or user_input[0]=='l':
      for param_name in param_names:
        print(param_name+': '+str(workspace.var(param_name).getValV()))
      #workspace.Print('v') 
    if user_input[0]=='help' or user_input[0]=='h':
      print('This is an interactive fitting session. Commands include:')
      print('(f)it                attempt a fit')
      print('(l)ist               display values of variables')
      #print('(n)orm               automatically fix norms')
      print('(q)uit               exit interactive fitting session')
      print('(a)utoload <values>  set all variables at once (or output)')
      print('(r)evert             revert to prefit parameter values')
      print('(s)et <var> <value>  set variable <var> to <value>')
      print('(w)rite <fname>      write current canvas to a file')
    if user_input[0]=='set' or user_input[0]=='s':
      if len(user_input)<3:
        print('ERROR: (s)et takes two arguments: set <var> <value>')
        continue
      if not (user_input[1] in param_names):
        print('ERROR: unknown parameter '+user_input[1])
        continue
      try:
        float(user_input[2])
        workspace.var(user_input[1]).setVal(float(user_input[2]))
        plot = workspace.var(x_name).frame()
        data.plotOn(plot)
        pdf_sb.plotOn(plot,ROOT.RooFit.Name('prefit_sb'))
        sig_norm_temp = workspace.var('sig_norm').getValV()
        bak_norm_temp = workspace.var('bak_norm').getValV()
        workspace.var('sig_norm').setVal(0.0)
        pdf_sb.plotOn(plot,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('prefit_b'))
        workspace.var('sig_norm').setVal(sig_norm_temp)
        canvas.cd()
        plot.Draw()
        canvas.Update()
        quality = fit_quality(canvas, 'h_data', 'prefit_sb', 78, 51.0, 129.0)
        print('quality: '+str(quality))
      except ValueError:
        print('Unable to cast value, skipping')
    if user_input[0]=='autoload' or user_input[0]=='a':
      if len(user_input) < 2:
        #output autoload values
        first = True
        for param_name in param_names:
          if first:
            first = False
          else:
            print(',',end='')
          print(str(workspace.var(param_name).getValV()),end='')
        print('')
      else:
        #load autoload values
        param_values = []
        try:
          param_values = [float(value) for value in user_input[1].split(',')]
        except ValueError:
          print('unable to cast value, skipping')
        if len(param_values) != len(param_names):
          print('incorrect number of autoset values, skipping')
        else:
          for ipar in range(len(param_names)):
            workspace.var(param_names[ipar]).setVal(param_values[ipar])
          plot = workspace.var(x_name).frame()
          data.plotOn(plot)
          pdf_sb.plotOn(plot,ROOT.RooFit.Name('prefit_sb'))
          sig_norm_temp = workspace.var('sig_norm').getValV()
          bak_norm_temp = workspace.var('bak_norm').getValV()
          workspace.var('sig_norm').setVal(0.0)
          pdf_sb.plotOn(plot,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('prefit_b'))
          workspace.var('sig_norm').setVal(sig_norm_temp)
          canvas.cd()
          plot.Draw()
          canvas.Update()
          quality = fit_quality(canvas, 'h_data', 'prefit_sb', 78, 51.0, 129.0)
          print('quality: '+str(quality))
    if user_input[0]=='fit' or user_input[0]=='f':
      workspace.saveSnapshot('prefit',','.join(param_names))
      pdf_sb.fitTo(data)
      plot = workspace.var(x_name).frame()
      data.plotOn(plot)
      pdf_sb.plotOn(plot,ROOT.RooFit.Name('postfit_sb'))
      sig_norm_temp = workspace.var('sig_norm').getValV()
      bak_norm_temp = workspace.var('bak_norm').getValV()
      workspace.var('sig_norm').setVal(0.0)
      pdf_sb.plotOn(plot,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('postfit_b'))
      workspace.var('sig_norm').setVal(sig_norm_temp)
      canvas.cd()
      plot.Draw()
      canvas.Update()
      quality = fit_quality(canvas, 'h_data', 'postfit_sb', 78, 51.0, 129.0)
      print('quality: '+str(quality))
      have_fit = True
    if user_input[0]=='revert' or user_input[0]=='r':
      if not have_fit:
        print('ERROR: must have already performed a fit to load snapshot')
        continue
      workspace.loadSnapshot('prefit')
      plot = workspace.var(x_name).frame()
      data.plotOn(plot)
      pdf_sb.plotOn(plot,ROOT.RooFit.Name('prefit_sb'))
      sig_norm_temp = workspace.var('sig_norm').getValV()
      bak_norm_temp = workspace.var('bak_norm').getValV()
      workspace.var('sig_norm').setVal(0.0)
      pdf_sb.plotOn(plot,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('prefit_b'))
      workspace.var('sig_norm').setVal(sig_norm_temp)
      canvas.cd()
      plot.Draw()
      canvas.Update()
    #if user_input[0]=='norm' or user_input[0]=='n':
    #  norm_set = ROOT.RooArgSet()
    #  for param_name in param_names:
    #    norm_set.add(workspace.var(param_name))
    #  integral = pdf_sb.createIntegral(ROOT.RooArgSet(workspace.var(x_name)),norm_set).getVal()
    #  print(str(workspace.var('sig_norm').getValV()/integral))
    #  print(str(workspace.var('bak_norm').getValV()/integral))
    if user_input[0]=='write' or user_input[0]=='w':
      if len(user_input)<3:
        print('ERROR: (w)rite takes one argument: write <fname>')
        continue
      canvas.SaveAs(user_input[1])
    if user_input[0]=='quit' or user_input[0]=='q':
      quit = True
  return quality


