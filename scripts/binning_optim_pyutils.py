
from array import array
import ROOT
from root_plot_lib import RplPlot

ROOT.gInterpreter.ProcessLine('.L scripts/binning_optim_cpputils.cpp+')

HZG_COLOR_DICT = {'ggf' : ROOT.TColor.GetColor('#bd1f01'),
                  'vbf' : ROOT.TColor.GetColor('#94a4a2'),
                  'dyg' : ROOT.TColor.GetColor('#3f90da'),
                  'dyf' : ROOT.TColor.GetColor('#ffa90e')}

MIN_SIGNAL = 1.0
MIN_VALUE = 0.06143 #yield of single ZG event

def make_validation_plots(df):
  '''
  Makes simple plots to show signal and background

  df  RDataFrame
  '''
  bins = array('f', [i*2.0/50.0-1.0 for i in range(51)])
  for name, type_string in [('ggf','1'),('vbf','2'),('dyg','3'),('dyf','4')]:
    hist = df.Filter('type=={}'.format(type_string)).Histo2D(
        ('{}_hist'.format(name),'; VBF score; Inclusive score',
        50,-1.0,1.0,50,-1.0,1.0),
        'VBF','twoj','weight')
    ROOT.remove_zeros(hist.GetPtr(),MIN_VALUE)
    plot = RplPlot()
    plot.plot_colormap(hist.GetPtr())
    plot.draw('plots/{}_validate.pdf'.format(name))

def make_onedim_plots(ggf_hist, vbf_hist, dyg_hist, dyf_hist, name):
  '''Make plots of 1D discriminant

  ggf_hist  TH2D
  vbf_hist  TH2D
  dyg_hist  TH2D
  dyf_hist  TH2D
  name      string
  '''
  plot = RplPlot()
  plot.y_title = '% Events/bin'
  plot.plot_outline(ggf_hist,HZG_COLOR_DICT['ggf'])
  plot.plot_outline(vbf_hist,HZG_COLOR_DICT['vbf'])
  plot.plot_outline(dyg_hist,HZG_COLOR_DICT['dyg'])
  plot.plot_outline(dyf_hist,HZG_COLOR_DICT['dyf'])
  plot.draw('plots/binning_onedim{}.pdf'.format(name))

def make_cpp_vector_string(lst):
  '''Converts list of Python strings into C++ string vector

  lst  list of strings
  '''
  #pyROOT is supposed to be able to used Python lists, but something about the
  #way it infers types can cause issues, so we'll stick to C++ vectors for now
  cpp_vec = ROOT.std.vector('std::string')()
  for item in lst:
    cpp_vec.push_back(item)
  return cpp_vec

def analyze_2d_discriminant(df):
  '''Optimizes 2D binning

  df RDataFrame
  '''
  sig_hist_tst = df.Filter('type<=2&&(rdfentry_%2)==0').Histo2D(
      ('sig_hist','; VBF score; Inclusive score',50,-1.0,1.0,50,-1.0,1.0),
      'VBF','twoj','w')
  bkg_hist_tst = df.Filter('type>=3&&(rdfentry_%2)==0').Histo2D(
      ('bkg_hist','; VBF score; Inclusive score',50,-1.0,1.0,50,-1.0,1.0),
      'VBF','twoj','w')
  sig_hist_trn = df.Filter('type<=2&&(rdfentry_%2)==1').Histo2D(
      ('sig_hist','; VBF score; Inclusive score',50,-1.0,1.0,50,-1.0,1.0),
      'VBF','twoj','w')
  bkg_hist_trn = df.Filter('type>=3&&(rdfentry_%2)==1').Histo2D(
      ('bkg_hist','; VBF score; Inclusive score',50,-1.0,1.0,50,-1.0,1.0),
      'VBF','twoj','w')
  sig_hist_trn.Smooth()
  bkg_hist_trn.Smooth()
  #sig_hist_tst.Smooth()
  #bkg_hist_tst.Smooth()
  ROOT.optimize_2d_binning(sig_hist_trn.GetPtr(),bkg_hist_trn.GetPtr(),
                           sig_hist_tst.GetPtr(),bkg_hist_tst.GetPtr(),
                           MIN_SIGNAL,MIN_VALUE)

def get_likelihood_evaluator(df):
  '''Constructs likelihood discriminant evaluator. Returns 
  HistLikelihoodRatio pointer. Also does 2D binning optimization for
  reference

  df  RDataFrame
  '''
  #to avoid cherry picking bins with fluctuations for observation, use odd 
  #events
  sig_hist = df.Filter('type<=2&&(rdfentry_%2)==0').Histo2D(
      ('sig_hist','; VBF score; Inclusive score',24,-1.0,1.0,24,-1.0,1.0),
      'VBF','twoj','w')
  bkg_hist = df.Filter('type>=3&&(rdfentry_%2)==0').Histo2D(
      ('bkg_hist','; VBF score; Inclusive score',24,-1.0,1.0,24,-1.0,1.0),
      'VBF','twoj','w')
  return ROOT.HistLikelihoodRatio(sig_hist.GetPtr(), bkg_hist.GetPtr(), 
                                  MIN_VALUE)

def draw_likelihood_boundaries(likelihood_evaluator):
  '''Generates 2D histogram showing bin boundaries of 1D likelihood binning

  likelihood_evaluator  HistLikelihoodRatio pointer
  '''
  #show how 1D bins map onto 2D space
  disc_hist = ROOT.TH2D('disc_hist',
                        ';VBF Score; Inclusive Score; 1D Discriminant',
                        100,-1.0,1.0,100,-1.0,1.0)
  bin2d_hist = ROOT.TH2D('bin2d_hist',
                         ';VBF Score; Inclusive Score; Bin boundaries',
                         100,-1.0,1.0,100,-1.0,1.0)
  boundaries = [0.78,0.92,0.96]
  for x in range(100):
    x_val = x*0.02-1.0+0.01
    for y in range(100):
      y_val = y*0.02-1.0+0.01
      disc = likelihood_evaluator.evaluate(x_val, y_val)
      disc_hist.SetBinContent(x+1,y+1,disc)
      if (disc < boundaries[0]):
        bin2d_hist.SetBinContent(x+1,y+1,0.01)
      elif (disc < boundaries[1]):
        bin2d_hist.SetBinContent(x+1,y+1,1.0)
      elif (disc < boundaries[2]):
        bin2d_hist.SetBinContent(x+1,y+1,2.0)
      else:
        bin2d_hist.SetBinContent(x+1,y+1,3.0)
  boundary_plot = RplPlot()
  boundary_plot.plot_colormap(bin2d_hist)
  boundary_plot.draw('plots/likelihood_2dboundaries.pdf')
  ROOT.remove_zeros(disc_hist)
  disc_plot = RplPlot()
  disc_plot.plot_colormap(disc_hist)
  disc_plot.draw('plots/likelihood_2dvalue.pdf')

def analyze_1d_discriminant(df, column, nbins, lo_edge, hi_edge):
  '''Analyzes a 1D discriminant: produces 1D plot and optimizes binning

  df       dataframe
  column   discriminant column name
  nbins    number of histogram bins
  lo_edge  low edge of histograms
  hi_edge  high edge of histograms

  '''
  ggf_1d_hist_ptr = df.Filter('type==1&&(rdfentry_%2)==1').Histo1D((
      'ggf_1d_hist',
      'gg#rightarrow H#rightarrow Z#gamma; 1D score; % Events/bin',
      nbins, lo_edge, hi_edge),column,'w')
  vbf_1d_hist_ptr = df.Filter('type==2&&(rdfentry_%2)==1').Histo1D((
      'vbf_1d_hist',
      'VBF H#rightarrow Z#gamma; 1D score; % Events/bin',
      nbins, lo_edge, hi_edge),column,'w')
  dyg_1d_hist_ptr = df.Filter('type==3&&(rdfentry_%2)==1').Histo1D((
      'dyg_1d_hist','Z#gamma; 1D score; % Events/bin',
      nbins, lo_edge, hi_edge),column,'w')
  dyf_1d_hist_ptr = df.Filter('type==4&&(rdfentry_%2)==1').Histo1D((
      'dyf_1d_hist','Z+Fake photon; 1D score; % Events/bin',
      nbins, lo_edge, hi_edge),column,'w')
  sig_trn_1d_hist_ptr = df.Filter('type<=2&&(rdfentry_%2)==0').Histo1D((
      'sig_trn_1d_hist',
      'H#rightarrow Z#gamma; 1D score; % Events/bin',
      nbins, lo_edge, hi_edge),column,'w')
  bkg_trn_1d_hist_ptr = df.Filter('type>=3&&(rdfentry_%2)==0').Histo1D((
      'bkg_trn_1d_hist',
      'Background; 1D score; % Events/bin',
      nbins, lo_edge, hi_edge),column,'w')
  sig_trn_1d_hist = sig_trn_1d_hist_ptr.GetPtr()
  bkg_trn_1d_hist = bkg_trn_1d_hist_ptr.GetPtr()
  sig_trn_1d_hist.Smooth()
  bkg_trn_1d_hist.Smooth()
  ggf_1d_hist = ggf_1d_hist_ptr.GetPtr()
  vbf_1d_hist = vbf_1d_hist_ptr.GetPtr()
  dyg_1d_hist = dyg_1d_hist_ptr.GetPtr()
  dyf_1d_hist = dyf_1d_hist_ptr.GetPtr()
  sig_1d_hist = ggf_1d_hist.Clone()
  sig_1d_hist.Add(vbf_1d_hist)
  bkg_1d_hist = dyg_1d_hist.Clone()
  bkg_1d_hist.Add(dyf_1d_hist)
  ggf_1d_hist.Smooth()
  vbf_1d_hist.Smooth()
  dyg_1d_hist.Smooth()
  dyf_1d_hist.Smooth()
  ggf_1d_hist.Scale(1.0/ggf_1d_hist.Integral())
  vbf_1d_hist.Scale(1.0/vbf_1d_hist.Integral())
  dyg_1d_hist.Scale(1.0/dyg_1d_hist.Integral())
  dyf_1d_hist.Scale(1.0/dyf_1d_hist.Integral())
  make_onedim_plots(ggf_1d_hist, vbf_1d_hist, dyg_1d_hist, dyf_1d_hist, column)
  #sig_1d_hist.Smooth()
  #bkg_1d_hist.Smooth()
  ROOT.optimize_1d_binning(sig_trn_1d_hist,bkg_trn_1d_hist,
                           sig_1d_hist,bkg_1d_hist,MIN_SIGNAL,MIN_VALUE)

def main():
  '''Runs some studies related to 2D binning
  '''
  df = ROOT.RDataFrame('tree','ntuples/binning_tree.root')
  df = df.Define('w','weight*2')
  likelihood_evaluator = get_likelihood_evaluator(df)
  make_validation_plots(df)
  analyze_2d_discriminant(df)
  df = df.Define('disc',likelihood_evaluator.get_evaluator(),
                 make_cpp_vector_string(['VBF','twoj']))
  df = df.Define('simple_disc','VBF+twoj')
  analyze_1d_discriminant(df, 'disc', 50, 0.0, 1.0)
  analyze_1d_discriminant(df, 'simple_disc', 100, -2.0, 2.0)
  draw_likelihood_boundaries(likelihood_evaluator)

if __name__=='__main__':
  main()
