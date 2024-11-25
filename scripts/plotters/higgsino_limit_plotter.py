import ROOT, array

LOG_MASS_DIFFERENCE = True

pen_PSGczDt = '#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}'
pen_PSGczDo = '#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGczDot = '#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1/2}}}'
pen_PSGcpmDo = '#lower[-0.12]{#tilde{#chi}}#lower[-0.1]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcmpDo = '#lower[-0.12]{#tilde{#chi}}#lower[-0.1]{#scale[0.85]{^{#mp}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcpDo = '#lower[-0.12]{#tilde{#chi}}#lower[0.1]{#scale[0.85]{^{+}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
pen_PSGcmDo = '#lower[-0.12]{#tilde{#chi}}#lower[0.1]{#scale[0.85]{^{-}}}#kern[-1.3]{#scale[0.85]{_{1}}}'
n1n2_string = pen_PSGczDo+pen_PSGczDt
cn_string = pen_PSGcpmDo+pen_PSGczDot
cn1_string = pen_PSGcpmDo+pen_PSGczDo
cn2_string = pen_PSGcpmDo+pen_PSGczDt
cc_string = pen_PSGcpDo+pen_PSGcmDo

class LineGraph():

  def __init__(self, x_list, y_list, text_):
    if len(x_list) != len(y_list):
      print('ERROR: unequal length lists')
    y_list_new = []
    if LOG_MASS_DIFFERENCE:
      y_list_new = y_list
    else:
      for i in range(len(x_list)):
        y_list_new.append(x_list[i]-y_list[i])
    x_array = array.array('d',x_list)
    y_array = array.array('d',y_list_new)
    self.first_x = x_list[0]
    self.first_y = y_list_new[0]
    self.manual_adjust_x = 1.0
    self.manual_adjust_y = 1.0
    self.graph = ROOT.TGraph(len(x_list),x_array,y_array)
    self.text = text_

  def set_color(self, color):
    self.graph.SetLineColor(color)
    return self

  def set_manual_adjust(self,adjust_x,adjust_y):
    self.manual_adjust_x = adjust_x
    self.manual_adjust_y = adjust_y
    return self

#Create the canvas and graph
ROOT.gStyle.SetOptStat(0)
mycanvas = ROOT.TCanvas('higgsino_canvas','higgsino_limits',600,400)
mycanvas.cd()
#mycanvas.SetLogx()
hist_title = 'Limits on higgsino production; (Heavy) higgsino mass [GeV]; LSP mass [GeV]'
if LOG_MASS_DIFFERENCE:
  mycanvas.SetLogy()
  hist_title = 'Limits on higgsino production; (Heavy) higgsino mass [GeV]; (Heavy) higgsino mass - LSP mass [GeV]'
base_histo = ROOT.TH1D('base_histo',hist_title,100,0.0,1000.0)
#base_histo.GetXaxis().SetLimits(1.0e1, 1.0e4) #GeV
base_histo.GetXaxis().CenterTitle(True)
base_histo.GetXaxis().SetLabelOffset(0.014)
base_histo.GetXaxis().SetTitleOffset(1.35)
base_histo.GetXaxis().SetTickLength(-0.02)
base_histo.GetYaxis().CenterTitle(True)
base_histo.GetYaxis().SetLabelOffset(0.014)
base_histo.GetYaxis().SetTitleOffset(1.35)
base_histo.GetYaxis().SetTickLength(-0.01)
if LOG_MASS_DIFFERENCE:
  base_histo.SetMinimum(1.5e-1) #100 MeV
  base_histo.SetMaximum(1.0e3) #1 TeV
else:
  base_histo.SetMinimum(0)
  base_histo.SetMaximum(1.0e3) #1 TeV
base_histo.Draw()

texts = ROOT.TLatex()
texts.SetTextSize(0.03)

physical_boundary_y = [0.0,0.00,0.00,0.00,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
lep_limits_x = [102.9,103.6]
lep_limits_y = [0.100,103.6]
if LOG_MASS_DIFFERENCE:
  physical_boundary_y = [0.1,0.18,0.32,0.56,1.0,1.8,3.2,5.6,10,18,32,56,100,180,320,560,999]
  lep_limits_x = [102.9,92.8,94.9,94.9,94.9,102.5,102.5,103.6]
  lep_limits_y = [0.100,0.20,0.30,1.00,2.00,5.000,10.00,103.6]

graphs_to_draw = [
    LineGraph([0.1,0.18,0.32,0.56,1.0,1.8,3.2,5.6,10,18,32,56,100,180,320,560,999],
              physical_boundary_y,'Physical region boundary').set_color(ROOT.kBlack),
    LineGraph(lep_limits_x,
              lep_limits_y,'LEP limits').set_color(ROOT.kCyan),
    LineGraph([610.0,358.0,150.0,100.0],
              [0.150,0.200,0.300,0.330],'CMS disappearing track').set_color(ROOT.kRed),
    #LineGraph([93.0,200.0,400.0,600.0,710.0],
    #          [0.378,0.305,0.227,0.178,0.15],'ATLAS disappearing track').set_color(ROOT.kRed),
    LineGraph([175.4,200.0,292.7,400.0,500.0,700.0,750.0,800.0,808.0],
              [175.4,167.4,192.7,230.3,305.8,503.5,579.1,734.9,808.0],'EWK combo WH').set_color(ROOT.kViolet),
    LineGraph([266.0,282.0,300.0,306.6],
              [266.0,267.3,290.1,306.6],'HH 4b').set_color(ROOT.kGreen),
    LineGraph([155.0,188.0,205.0,150.0,125.0,110.0,125.0,150.0,187.5,175.0,150.0,175.0,200.0,300.0,400.0,500.0,575.0,600.0],
              [3.000,5.000,7.500,16.00,20.00,30.00,38.00,46.50,65.60,79.00,100.0,105.0,115.0,167.0,215.0,322.0,465.0,600.0],'EWK combo WZ').set_color(ROOT.kBlue)] #175 was above 50 GeV
    #LineGraph([100.0,132.0,184.0,196.2,180.0,153.0,149.0,154.0,169.0,186.0,199.0,206.0,212.0],
    #          [2.667,4.000,8.000,10.67,14.00,20.00,25.00,30.00,32.00,35.00,40.00,50.00,60.00],'ATLAS soft').set_color(ROOT.kOrange)]
    #LineGraph([500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,2750,3000],[33800,11300,4320,1810,812,385,71.5,15.7,3.85,1.01,0.275,0.0768,0.0216,0.00612],'#color[800]{#tilde{g}#tilde{g}}').set_color(ROOT.kOrange),
    #              #LineGraph([100,150,200,300,400,500,600,700,800,900,1000,1200,1400],[22670,5180,1807,387,121,46.4,20.1,9.51,4.76,2.50,1.34,0.416,0.131],'#color[432]{Wino '+cn1_string+'}').set_color(ROOT.kCyan),
    #              LineGraph([127,150,200,300,400,500,600,700,800,900,1000,1200,1400],[7602,3832,1336,285,88.7,33.8,14.7,6.90,3.46,1.81,0.969,0.299,0.078],'#color[632]{Higgsino '+n1n2_string+'+'+cn1_string+'+'+cn2_string+'+'+cc_string+'}').set_color(ROOT.kRed),
    #              LineGraph([127,150,200,300,400,500,600,700,800,900,1000,1200,1400],[1447,715,244,51.0,15.7,5.91,2.53,1.18,0.586,0.306,0.164,0.0516,0.0175],'#color[880]{Higgsino '+n1n2_string+'}').set_color(ROOT.kViolet).set_manual_adjust(0.7,1.04)] #.set_manual_adjust(1,0.05),
    #              #LineGraph([100,150,200,300,400,500,600,700,800,900,1000],[366,87.1,30.3,6.25,1.86,0.674,0.276,0.123,0.0586,0.0292,0.0150],'#color[600]{#tilde{l}_{L,R}#tilde{l}_{L,R}}').set_color(ROOT.kBlue).set_manual_adjust(1,0.15)]

for graph in graphs_to_draw:
  graph.graph.SetLineWidth(2)
  if (graph.text == 'Physical region boundary'):
    graph.graph.SetLineStyle(ROOT.kDashed)
  graph.graph.Draw('C SAME')
  #texts.DrawLatex(graph.first_x*graph.manual_adjust_x,graph.first_y*graph.manual_adjust_y,' '+graph.text)

if LOG_MASS_DIFFERENCE:
  texts.DrawLatex(250,50,'#color[632]{'+cn1_string+','+cn2_string+','+cc_string+' LLP} CMS disappearing track [2004.05153]')
  texts.DrawLatex(250,30,'#color[600]{'+cn2_string+'#rightarrow WZ+p_{T}^{miss}} CMS combination [CMS-PAS-SUS-21-008] reinterpretation')
  texts.DrawLatex(250,18,'#color[880]{'+cn1_string+','+cn2_string+'#rightarrow WH+p_{T}^{miss}} CMS combination [CMS-PAS-SUS-21-008]')
  texts.DrawLatex(250,10.8,'#color[416]{'+n1n2_string+'#rightarrow HH+p_{T}^{miss}} CMS 4b [2201.04206]')
  texts.DrawLatex(250,6.48,'#color[432]{'+cn1_string+','+cn2_string+','+cc_string+'} LEP combination [LEPSUSYWG/0.2-0.4.1,01-03.1]')
else:
  texts.DrawLatex(50,900,'#color[632]{'+cn1_string+','+cn2_string+','+cc_string+' LLP} CMS disappearing track [2004.05153]')
  texts.DrawLatex(50,850,'#color[600]{'+cn2_string+'#rightarrow WZ+p_{T}^{miss}} CMS combination [CMS-PAS-SUS-21-008] reinterpretation')
  texts.DrawLatex(50,800,'#color[880]{'+cn1_string+','+cn2_string+'#rightarrow WH+p_{T}^{miss}} CMS combination [CMS-PAS-SUS-21-008]')
  texts.DrawLatex(50,750,'#color[416]{'+n1n2_string+'#rightarrow HH+p_{T}^{miss}} CMS 4b [2201.04206]')
  texts.DrawLatex(50,700,'#color[432]{'+cn1_string+','+cn2_string+','+cc_string+'} LEP combination [LEPSUSYWG/0.2-0.4.1,01-03.1]')

mycanvas.Draw()
filename = 'higgsino_limits_lin.pdf'
if LOG_MASS_DIFFERENCE:
  filename = 'higgsino_limits.pdf'
mycanvas.SaveAs(filename)
