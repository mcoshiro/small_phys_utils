import ROOT, array

pen_PSGczDt = '#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}'
pen_PSGczDo = '#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcpmDo = '#lower[-0.12]{#tilde{#chi}}#lower[-0.1]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcmpDo = '#lower[-0.12]{#tilde{#chi}}#lower[-0.1]{#scale[0.85]{^{#mp}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcpDo = '#lower[-0.12]{#tilde{#chi}}#lower[0.1]{#scale[0.85]{^{+}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcmDo = '#lower[-0.12]{#tilde{#chi}}#lower[0.1]{#scale[0.85]{^{-}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
n1n2_string = pen_PSGczDo+pen_PSGczDt
cn1_string = pen_PSGcpmDo+pen_PSGczDo
cn2_string = pen_PSGcpmDo+pen_PSGczDt
cc_string = pen_PSGcpDo+pen_PSGcmDo

hist_left_edge = 3.0e1
hist_right_edge = 1.0e3
total_xs = 104.0e-3 #ATLAS 2022

class Point():

  def __init__(self,x_,y_,text_):
    self.x = x_
    self.y = y_
    self.text = text_
    self.manual_adjust_x = 1.0
    self.manual_adjust_y = 1.0

  def set_manual_adjust(self,adjust_x,adjust_y):
    self.manual_adjust_x = adjust_x
    self.manual_adjust_y = adjust_y
    return self

class LineGraph():

  def __init__(self, x_list, y_list, text_):
    if len(x_list) != len(y_list):
      print('ERROR: unequal length lists')
    x_array = array.array('d',x_list)
    y_array = array.array('d',y_list)
    self.first_x = x_list[0]
    self.first_y = y_list[0]
    self.manual_adjust_x = 1.0
    self.manual_adjust_y = 1.0
    self.graph = ROOT.TGraph(len(x_list),x_array,y_array)
    self.text = text_

  def set_color(self, color):
    self.graph.SetLineColor(color)
    return self

class Arrow():

  def __init__(self, y1_, text_):
    self.x1 = hist_left_edge*1.5
    self.y1 = y1_
    self.arrow = ROOT.TArrow(self.x1,y1_,hist_left_edge,y1_,0.008,'|>')
    self.manual_adjust_x = 1.0
    self.manual_adjust_y = 1.0
    self.text = text_

  def set_manual_adjust(self,adjust_x,adjust_y):
    self.manual_adjust_x = adjust_x
    self.manual_adjust_y = adjust_y
    return self

#Create the canvas and graph
ROOT.gStyle.SetOptStat(0)
mycanvas = ROOT.TCanvas('dm_canvas','dm_limits',600,400)
mycanvas.cd()
mycanvas.SetLogx()
mycanvas.SetLogy()
base_histo = ROOT.TH1D('base_histo','LHC process energy vs. cross section, #sqrt{s}=13 TeV; (Minimum) system energy/mass [GeV]; Cross section [b]',100,hist_left_edge,hist_right_edge)
#base_histo.GetXaxis().SetLimits(1.0e1, 1.0e4) #GeV
base_histo.GetXaxis().CenterTitle(True)
base_histo.GetXaxis().SetTitleOffset(1.35)
base_histo.GetXaxis().SetTickLength(-0.01)
base_histo.GetYaxis().CenterTitle(True)
base_histo.GetYaxis().SetTitleOffset(1.35)
base_histo.GetYaxis().SetTickLength(-0.01)
base_histo.SetMinimum(1.0e-15) #1 fb
base_histo.SetMaximum(1.0e0) #1 b
base_histo.Draw()

right_axis = ROOT.TGaxis(hist_right_edge,1.0e-15,hist_right_edge,1.0e0,1.0e-15/total_xs,1.0/total_xs,510,'GLS')
right_axis.SetTitle('Probability')
right_axis.SetTitleFont(42)
right_axis.SetTitleOffset(-1.35)
right_axis.SetTitleSize(0.035)
right_axis.CenterTitle(True)
right_axis.SetTickLength(0.01)
right_axis.SetLabelFont(42)
right_axis.SetLabelSize(0.035)
right_axis.SetLabelOffset(-0.01)
right_axis.Draw()

texts = ROOT.TLatex()
texts.SetTextSize(0.03)

#draw arrows
arrows_to_draw = [Arrow(1.0e-1,'Total cross section'),Arrow(70.0e-3,'Inelastic collisions').set_manual_adjust(1.0,0.5)]

for arrow in arrows_to_draw:
  arrow.arrow.Draw()
  texts.DrawLatex(arrow.x1*arrow.manual_adjust_x,arrow.y1*arrow.manual_adjust_y,' '+arrow.text)

#draw graphs
graphs_to_draw = [LineGraph([50.0,100.0,200.0,500.0],[246.0e-6,28.0e-6,1.71e-6,32.0e-9],'Jets'),
                  #LineGraph([1000.0,1500.0,2000.0,3000.0,4000.0,5000.0],[33.8e-12,2.77e-12,385.0e-15,15.7e-15,1.01e-15,76.80e-18],'#tilde{g}#tilde{g}'),
                  #LineGraph([300.0,600.0,1000.0,2000.0],[3.83e-12,285.0e-15,33.8e-15,968.0e-18],'#tilde{H}#tilde{H}')]
                  LineGraph([300.0,600.0,1000.0,2000.0],[715.0e-15,51.0e-15,5.9e-15,164.0e-18],n1n2_string)]

for graph in graphs_to_draw:
  graph.graph.Draw('C SAME')
  texts.DrawLatex(graph.first_x*graph.manual_adjust_x,graph.first_y*graph.manual_adjust_y,' '+graph.text)

#draw points
points_to_draw = [Point(81.0,20.5e-9*9.0,'W'),Point(91.0,6.07e-9*10.0,'Z'),Point(173.0*2.0,831.0e-12,'t#bar{t}'),
                  Point(173.0,113.0e-12+67.0e-12,'t'),Point(125.0,48.58e-12,'H'),
                  Point(81.0*2.0,118.7e-12,'WW').set_manual_adjust(0.83,1.0),Point(81.0+91.0,51.11e-12,'WZ'),Point(91.0*2.0,16.91e-12,'ZZ'),
                  Point(173.0*2.0+81.0,243.0e-15*3.0,'t#bar{t}W').set_manual_adjust(0.9,1.5),Point(173.0*2.0+91.0,540.0e-15,'t#bar{t}Z').set_manual_adjust(1.0,1.4),
                  Point(81.0+125.0,1.35e-12,'WH'),Point(91.0+125.0,884.0e-15,'ZH').set_manual_adjust(1.0,0.5),
                  Point(173.0*2.0+125.0,507.0e-15,'t#bar{t}H').set_manual_adjust(1.0,0.45),
                  Point(173.0*4.0,8.21e-15,'t#bar{t}t#bar{t}').set_manual_adjust(1.025,0.52),
                  Point(125.0*2.0,31.05e-15,'HH')]
markers = []

for point in points_to_draw:
  markers.append(ROOT.TMarker(point.x,point.y,20))
  markers[len(markers)-1].SetMarkerSize(0.5)
  markers[len(markers)-1].Draw()
  texts.DrawLatex(point.x*point.manual_adjust_x,point.y*point.manual_adjust_y,' '+point.text)

mycanvas.Draw()
mycanvas.SaveAs('lhc_xs.pdf')
