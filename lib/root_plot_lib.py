"""@package docstring
Python package to make nice-looking ROOT plots
"""

from array import array
#from os import path
import ROOT

def get_palette_hig21001():
  '''Returns list of colors used in HIG-21-001
  '''
  return [ROOT.TColor.GetColor('#de5a6a'), ROOT.TColor.GetColor('#ffcc66'), 
          ROOT.TColor.GetColor('#c0e684'), ROOT.TColor.GetColor('#64c0e8'), 
          ROOT.TColor.GetColor('#9999cc'), ROOT.TColor.GetColor('#ffccff')]

def get_palette_official(nplots):
  '''Returns list of colors to be used for plots

  @params
  nplots - number of different overlayed plots
  '''
  if (nplots <= 6):
    return [ROOT.TColor.GetColor('#5790fc'), ROOT.TColor.GetColor('#f89c20'), 
            ROOT.TColor.GetColor('#e42536'), ROOT.TColor.GetColor('#964a8b'), 
            ROOT.TColor.GetColor('#9c9ca1'), ROOT.TColor.GetColor('#7a21dd')]
  raise ValueError('More colors than allowed in official colors.')

def get_palette_lines(nplots):
  '''Returns list of colors to be used for plots

  @params
  nplots - number of different overlayed plots
  '''
  if (nplots <= 6):
    return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ffbd38'), 
            ROOT.TColor.GetColor('#aee82a'), ROOT.TColor.GetColor('#28aee8'), 
            ROOT.TColor.GetColor('#446bcc'), ROOT.TColor.GetColor('#7346cc')]
  if (nplots <= 8):
    return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ff9c38'), 
            ROOT.TColor.GetColor('#f6c24b'), ROOT.TColor.GetColor('#aee82a'), 
            ROOT.TColor.GetColor('#28aee8'), ROOT.TColor.GetColor('#446bcc'), 
            ROOT.TColor.GetColor('#7346cc'), ROOT.TColor.GetColor('#ee66ac')]
  return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ff9c38'), 
          ROOT.TColor.GetColor('#f6c24b'), ROOT.TColor.GetColor('#aee82a'), 
          ROOT.TColor.GetColor('#1b9e48'), ROOT.TColor.GetColor('#28aee8'), 
          ROOT.TColor.GetColor('#446bcc'), ROOT.TColor.GetColor('#50378f'),
          ROOT.TColor.GetColor('#7346cc'), ROOT.TColor.GetColor('#ee66ac')]

class RplPlot:
  '''Simple class to hold plot options and apply to ROOT plots
  '''
  def __init__(self):
    '''Initializer to set defaults
    '''
    self.plot_bottom = False
    self.y_title = 'Events/bin'
    self.y_title_lower = 'data/MC'
    self.x_title = 'x variable'
    self.x_min = -999.0
    self.x_max = -999.0
    self.y_min = -999.0
    self.y_max = -999.0
    self.y_min_lower = 0.75
    self.y_max_lower = 1.25
    self.log_x = False
    self.log_y = False
    self.log_y_bottom = False
    self.lumi_data = [(138,13)] #list of tuples in format (lumi [fbinv], energy [TeV])
    self.title_type = 'cms work in progress'
    self.legend_xlo = 0.19
    self.legend_xhi = 0.91
    self.legend_ylo = 0.78
    self.legend_yhi = 0.90
    self.legend_customsize = False
    self.legend_ncolumns = -1
    self.point_hists = []
    self.point_hist_color = []
    self.line_hists = []
    self.line_hist_color = []
    self.graphs = []
    self.graph_color = []
    self.bottom_plots = []
    self.bottom_is_ratio = True
    self.bottom_plot_color_index = []
    self.n_plots = 0

  def plot_points(self, hist, color=None):
    '''plots a TH1 as data points

    @params
    hist - root TH1 to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    '''
    self.point_hist_color.append(color)
    self.point_hists.append(hist)
    if self.n_plots == 0:
      xaxis_title = hist.GetXaxis().GetTitle()
      if self.x_title == 'x variable' and xaxis_title != '':
        self.x_title = xaxis_title
    self.n_plots += 1
    return self

  def plot_outline(self, hist, color=None):
    '''plots a TH1 as an outline

    @params
    hist - root TH1 to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    '''
    self.line_hist_color.append(color)
    self.line_hists.append(hist)
    if self.n_plots == 0:
      xaxis_title = hist.GetXaxis().GetTitle()
      if self.x_title == 'x variable' and xaxis_title != '':
        self.x_title = xaxis_title
    self.n_plots += 1
    return self

  def plot_graph(self, graph, color=None):
    '''plots a TGraph

    @params
    graph - root TGraph to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    '''
    self.graph_color.append(color)
    self.graphs.append(graph)
    if self.n_plots == 0:
      xaxis_title = hist.GetXaxis().GetTitle()
      if self.x_title == 'x variable' and xaxis_title != '':
        self.x_title = xaxis_title
    self.n_plots += 1
    return self

  def add_ratio(self, numerator_name, denominator_name):
    '''adds a ratio of two plot elements to the plot

    @params
    numerator_name - name of plot element to use as numerator
    denominator_name - name of plot element to use as denominator
    '''
    numerator_hist = None
    denominator_hist = None
    numerator_color_index = -1
    denominator_color_index = -1
    color_index = 0
    for hist in self.point_hists:
      if (hist.GetName()==numerator_name):
        numerator_hist = hist
        numerator_color_index = color_index
      if (hist.GetName()==denominator_name):
        denominator_hist = hist
        denominator_color_index = color_index
      color_index += 1
    for hist in self.line_hists:
      if (hist.GetName()==numerator_name):
        numerator_hist = hist
        numerator_color_index = color_index
      if (hist.GetName()==denominator_name):
        denominator_hist = hist
        denominator_color_index = color_index
      color_index += 1
    if (numerator_hist != None and denominator_hist != None):
      self.bottom_plots.append(numerator_hist.Clone())
      self.bottom_plots[-1].Divide(denominator_hist)
      self.bottom_plot_color_index.append(denominator_color_index)
      self.plot_bottom = True
    else:
      raise ValueError('Could not find numerator and denominator plots requested for ratio.')
    #TODO implement ratios of graphs
    return self

  def add_difference(self, minuend_name, subtrahend_name):
    '''adds a ratio of two plot elements to the plot

    @params
    minuend_name - name of plot element to use as minuend
    subtrahend_name - name of plot element to use as subtrahend
    '''
    self.bottom_is_ratio = False
    self.y_title_lower = 'data-MC'
    minuend_hist = None
    subtrahend_hist = None
    minuend_color_index = -1
    color_index = 0
    for hist in self.point_hists:
      if (hist.GetName()==minuend_name):
        minuend_hist = hist
        minuend_color_index = color_index
      if (hist.GetName()==subtrahend_name):
        subtrahend_hist = hist
      color_index += 1
    for hist in self.line_hists:
      if (hist.GetName()==minuend_name):
        minuend_hist = hist
        minuend_color_index = color_index
      if (hist.GetName()==subtrahend_name):
        subtrahend_hist = hist
      color_index += 1
    if (minuend_hist != None and subtrahend_hist != None):
      self.bottom_plots.append(minuend_hist.Clone())
      self.bottom_plots[-1].Add(subtrahend_hist, -1.0)
      self.bottom_plot_color_index.append(minuend_color_index)
      self.plot_bottom = True
    else:
      raise ValueError('Could not find minuend and subtrahend plots requested for difference.')
    #TODO implement differences of graphs
    return self

  def draw(self, filename='my_plot.pdf'):
    '''draws plot and saves to output file

    @params
    filename - name of plot to save file
    '''
    #remove stat box - not working??
    ROOT.gStyle.SetOptStat(0)

    #find maxima and minima
    if (self.y_max == -999.0 or self.y_min == -999.0):
      self.y_min = 0
      if self.log_y:
        self.y_min = 0.01
      self.y_max = 0 
      for point_hist in self.point_hists:
        if point_hist.GetMaximum()>self.y_max:
          self.y_max = point_hist.GetMaximum()
      for line_hist in self.line_hists:
        if line_hist.GetMaximum()>self.y_max:
          self.y_max = line_hist.GetMaximum()
      for graph in self.graphs:
        graph_max = ROOT.TMath.MaxElement(graph.GetN(), graph.GetY())
        if graph_max>self.y_max:
          self.y_max = graph_max
      if self.log_y:
        self.y_max = ((self.y_max/self.y_min)**1.45)*self.y_min
      else:
        self.y_max = self.y_max*1.45
    if (self.x_max == -999.0 or self.x_min == -999.0):
      self.x_min = 999999.0
      self.x_max = -999999.0
      #TODO implement something more robust, for now, just use first hist
      for point_hist in self.point_hists:
        lo_edge = point_hist.GetXaxis().GetBinLowEdge(1)
        hi_edge = point_hist.GetXaxis().GetBinUpEdge(point_hist.GetXaxis().GetNbins())
        if (lo_edge < self.x_min):
          self.x_min = lo_edge
        if (hi_edge > self.x_max):
          self.x_max = hi_edge
      for line_hist in self.line_hists:
        lo_edge = line_hist.GetXaxis().GetBinLowEdge(1)
        hi_edge = line_hist.GetXaxis().GetBinUpEdge(line_hist.GetXaxis().GetNbins())
        if (lo_edge < self.x_min):
          self.x_min = lo_edge
        if (hi_edge > self.x_max):
          self.x_max = hi_edge
      for graph in self.graphs:
        lo_edge = ROOT.TMath.MinElement(graph.GetN(), graph.GetX())
        hi_edge = ROOT.TMath.MaxElement(graph.GetN(), graph.GetX())
        if (lo_edge < self.x_min):
          self.x_min = lo_edge
        if (hi_edge > self.x_max):
          self.x_max = hi_edge

    #set up dummy histograms
    dummy_hist_upper = ROOT.TH1D('','',100,self.x_min,self.x_max)
    dummy_hist_upper.SetMinimum(self.y_min)
    dummy_hist_upper.SetMaximum(self.y_max)
    dummy_hist_upper.SetLabelSize(0.028,'y')
    dummy_hist_upper.SetTitleSize(0.032,'y')
    dummy_hist_upper.GetYaxis().SetTitle(self.y_title)
    dummy_hist_upper.GetYaxis().SetNdivisions(606)
    dummy_hist_upper.GetXaxis().SetNdivisions(606)
    if (not self.plot_bottom):
      dummy_hist_upper.GetXaxis().SetTitle(self.x_title)
      dummy_hist_upper.SetTitleSize(0.032,'x')
    else:
      dummy_hist_upper.SetLabelSize(0,'x')

    dummy_hist_lower = ROOT.TH1D('','',100,self.x_min,self.x_max)
    dummy_hist_lower.SetMinimum(self.y_min_lower)
    dummy_hist_lower.SetMaximum(self.y_max_lower)
    dummy_hist_lower.SetLabelSize(0.028,'x')
    dummy_hist_lower.SetLabelSize(0.028,'y')
    dummy_hist_lower.SetTitleOffset(1.5,'y')
    #dummy_hist_lower.SetTitleOffset(2.0,'y')
    #dummy_hist_lower.SetTitleSize(0.45/float(len(self.y_title_lower)),'y')
    dummy_hist_lower.SetTitleSize(0.032,'x')
    dummy_hist_lower.GetYaxis().SetTitle(self.y_title_lower)
    dummy_hist_lower.GetYaxis().SetNdivisions(606)
    dummy_hist_lower.GetYaxis().SetLimits(self.y_min_lower,self.y_max_lower)
    dummy_hist_lower.GetXaxis().SetTitle(self.x_title)
    dummy_hist_lower.GetXaxis().SetNdivisions(606)

    #set up pads
    canvas_name = filename[:filename.rfind('.')]
    if canvas_name.rfind('/') != -1:
      canvas_name = canvas_name[canvas_name.rfind('/')+1:]
    can = ROOT.TCanvas('c_'+canvas_name,'c',600,600)
    top_pad = ROOT.TPad('top_pad','',0.0,0.0,1.0,1.0)
    top_pad.SetTicks(1,1)
    if (self.plot_bottom):
      top_pad.SetMargin(0.14,0.06,0.31,0.05)
    else:
      top_pad.SetMargin(0.14,0.06,0.1,0.05)
    top_pad.SetFillStyle(4000)
    top_pad.SetLogy(self.log_y)
    top_pad.SetLogx(self.log_x)

    bot_pad = ROOT.TPad('bot_pad','',0.0,0.0,1.0,1.0)
    bot_pad.SetTicks(1,1)
    if (self.plot_bottom):
      bot_pad.SetMargin(0.14,0.06,0.1,0.71)
    else:
      bot_pad.SetMargin(0.14,0.06,0.0,1.0)
    bot_pad.SetFillStyle(4000)
    bot_pad.SetLogy(self.log_y_bottom)
    bot_pad.SetLogx(self.log_x)

    #draw upper plot
    top_pad.Draw()
    top_pad.cd()
    dummy_hist_upper.Draw()

    #draw plots and legend
    if (self.n_plots<0):
      raise ValueError('Must have positive number of plots.')
    if (self.n_plots>10):
      raise ValueError('Not enough colors for the number of plots.')
    palette = get_palette_official(self.n_plots)
    if (self.legend_ncolumns == -1):
      self.legend_ncolumns = self.n_plots//4+1
    if (not self.legend_customsize):
      if (self.n_plots < self.legend_ncolumns*4):
        self.legend_ylo = 0.9-0.03*(self.n_plots//self.legend_ncolumns+1)
    leg = ROOT.TLegend(self.legend_xlo,self.legend_ylo,self.legend_xhi,self.legend_yhi)
    leg.SetEntrySeparation(0)
    leg.SetTextSize(0.03)
    if (self.legend_ncolumns > -1):
      leg.SetNColumns(self.legend_ncolumns)
    else:
      n_columns = self.n_plots//4+1
      leg.SetNColumns(n_columns)
    color_index = 0
    custom_colors = self.point_hist_color+self.line_hist_color+self.graph_color
    for point_hist in self.point_hists:
      if (custom_colors[color_index]!=None):
        palette[color_index] = custom_colors[color_index]
      point_hist.SetLineWidth(3)
      point_hist.SetLineColor(palette[color_index])
      point_hist.SetMarkerColor(palette[color_index])
      point_hist.SetMarkerStyle(ROOT.kFullCircle)
      point_hist.Draw('same P')
      leg.AddEntry(point_hist, point_hist.GetTitle(), 'LP')
      color_index += 1
    for line_hist in self.line_hists:
      if (custom_colors[color_index]!=None):
        palette[color_index] = custom_colors[color_index]
      line_hist.SetLineWidth(3)
      line_hist.SetLineColor(palette[color_index])
      line_hist.Draw('same hist')
      leg.AddEntry(line_hist, line_hist.GetTitle(), 'F')
      color_index += 1
    for graph in self.graphs:
      if (custom_colors[color_index]!=None):
        palette[color_index] = custom_colors[color_index]
      graph.SetLineWidth(3)
      graph.SetLineColor(palette[color_index])
      graph.SetMarkerColor(palette[color_index])
      graph.SetMarkerStyle(ROOT.kFullCircle)
      graph.Draw('same P')
      leg.AddEntry(graph, graph.GetTitle(), 'LP')
      color_index += 1
    leg.SetBorderSize(0)
    leg.Draw('same')

    #draw CMS and lumi labels
    label = ROOT.TLatex()
    label.SetTextSize(0.032)
    label.SetNDC(ROOT.kTRUE)
    label.SetTextAlign(11)
    if (self.title_type == 'cms preliminary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}')
    elif (self.title_type == 'cms work in progress'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Work in Progress}}')
    elif (self.title_type == 'cms supplementary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Supplementary}}')
    elif (self.title_type == 'cms simulation'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}')
    elif (self.title_type == 'cms simulation supplementary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}')
    elif (self.title_type == 'cms simulation preliminary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Preliminary}}')
    elif (self.title_type == 'cms private work'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Private Work}}')
    elif (self.title_type == 'cms'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS}')
    label.SetTextAlign(31)
    label.SetTextSize(0.03)
    lumi_energy_string = '#font[42]{'
    first = True
    for lumi_datum in self.lumi_data:
      if first:
        first = False
      else:
        lumi_energy_string += ' + '
      lumi_energy_string += str(lumi_datum[0])+' fb^{-1} ('+str(lumi_datum[1])+' TeV)'
    lumi_energy_string += '}'
    label.DrawLatex(0.93,0.96,lumi_energy_string)
    top_pad.Modified()

    #draw lower plot
    if (self.plot_bottom):
      can.cd()
      bot_pad.Draw('same')
      bot_pad.cd()
      dummy_hist_lower.Draw()
      bottom_ref_value = 1.0
      if not self.bottom_is_ratio:
        bottom_ref_value = 0.0
      line = ROOT.TLine(self.x_min,bottom_ref_value,self.x_max,bottom_ref_value)
      line.SetNDC(ROOT.kFALSE)
      line.SetLineStyle(2)
      line.SetLineColor(ROOT.kBlack)
      line.SetLineWidth(2)
      line.Draw('SAME')
      bottom_index = 0
      for bottom_plot in self.bottom_plots:
        bottom_plot.SetLineWidth(3)
        bottom_plot.SetLineColor(palette[self.bottom_plot_color_index[bottom_index]])
        bottom_plot.SetMarkerColor(palette[self.bottom_plot_color_index[bottom_index]])
        bottom_plot.SetMarkerStyle(ROOT.kFullCircle)
        bottom_plot.Draw('same P')
        bottom_index += 1
      bot_pad.Modified()

    #draw everything and save
    can.Draw()
    can.SaveAs(filename)
    return self

