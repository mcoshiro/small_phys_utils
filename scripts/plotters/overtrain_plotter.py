import ROOT, array

#settings to change
x_list = []
y_list_a = []
y_list_b = []
plot_type_list = []
plot_xaxis_list = []
plot_xaxis_title = []
#nevents auc, csi, ks
x_list.append([5000,10000,25000,50000,100000,150000])
x_list.append([5000,10000,25000,50000,100000,150000])
x_list.append([5000,10000,25000,50000,100000,150000])
y_list_a.append([0.865, 0.848, 0.836, 0.823, 0.816, 0.814])
y_list_b.append([0.768, 0.782, 0.795, 0.802, 0.805, 0.805])
y_list_a.append([2.55, 2.14, 2.11, 2.07, 1.95, 1.94])
y_list_b.append([1.63, 1.67, 1.75, 1.83, 1.86, 1.86])
y_list_a.append([3.0e-17,1.1e-13,3.4e-7,0.0007,0.06,0.16])
y_list_b.append([0.007,0.009,0.0002,0.0007,0.09,0.32])
plot_type_list.append(0)
plot_type_list.append(1)
plot_type_list.append(2)
plot_xaxis_list.append('nevents')
plot_xaxis_list.append('nevents')
plot_xaxis_list.append('nevents')
plot_xaxis_title.append('Number of training events')
plot_xaxis_title.append('Number of training events')
plot_xaxis_title.append('Number of training events')
#ntrees auc, csi ks
x_list.append([10,25,50,100,200,400,850])
x_list.append([10,25,50,100,200,400,850])
x_list.append([10,25,50,100,200,400,850])
y_list_a.append([0.778,0.791,0.795,0.8,0.804,0.809,0.813])
y_list_b.append([0.778,0.79,0.794,0.798,0.801,0.804,0.805])
y_list_a.append([1.68,1.76,1.81,1.85,1.89,1.92,1.93])
y_list_b.append([1.69,1.75,1.80,1.84,1.86,1.88,1.85])
y_list_a.append([0.86,0.98,0.98,0.8,0.64,0.42,0.23])
y_list_b.append([0.73,0.83,0.81,0.86,0.78,0.64,0.12])
plot_type_list.append(0)
plot_type_list.append(1)
plot_type_list.append(2)
plot_xaxis_list.append('ntrees')
plot_xaxis_list.append('ntrees')
plot_xaxis_list.append('ntrees')
plot_xaxis_title.append('Number of trees')
plot_xaxis_title.append('Number of trees')
plot_xaxis_title.append('Number of trees')
#treedepth auc, csi ks
x_list.append([1.0,2.0,3.0,4.0])
x_list.append([1.0,2.0,3.0,4.0])
x_list.append([1.0,2.0,3.0,4.0])
y_list_a.append([0.788,0.807,0.813,0.818])
y_list_b.append([0.786,0.802,0.806,0.805])
y_list_a.append([1.76,1.89,1.93,1.98])
y_list_b.append([1.76,1.86,1.85,1.85])
y_list_a.append([0.83,0.25,0.24,0.009])
y_list_b.append([0.89,0.72,0.13,0.02])
plot_type_list.append(0)
plot_type_list.append(1)
plot_type_list.append(2)
plot_xaxis_list.append('treedepth')
plot_xaxis_list.append('treedepth')
plot_xaxis_list.append('treedepth')
plot_xaxis_title.append('Maximum tree depth')
plot_xaxis_title.append('Maximum tree depth')
plot_xaxis_title.append('Maximum tree depth')
#minnodesize auc, csi, ks
x_list.append([1.0,2.5,5.0,10.0,15.0])
x_list.append([1.0,2.5,5.0,10.0,15.0])
x_list.append([1.0,2.5,5.0,10.0,15.0])
y_list_a.append([0.813,0.813,0.810,0.805,0.798])
y_list_b.append([0.806,0.806,0.804,0.800,0.795])
y_list_a.append([1.93,1.93,1.92,1.86,1.81])
y_list_b.append([1.85,1.85,1.85,1.81,1.78])
y_list_a.append([0.24,0.23,0.26,0.48,0.66])
y_list_b.append([0.20,0.12,0.66,0.46,0.20])
plot_type_list.append(0)
plot_type_list.append(1)
plot_type_list.append(2)
plot_xaxis_list.append('minnodesize')
plot_xaxis_list.append('minnodesize')
plot_xaxis_list.append('minnodesize')
plot_xaxis_title.append('Minimum node size [%]')
plot_xaxis_title.append('Minimum node size [%]')
plot_xaxis_title.append('Minimum node size [%]')
#nvariables auc, csi, ks
x_list.append([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0])
x_list.append([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0])
x_list.append([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0])
y_list_a.append([0.727,0.755,0.788,0.791,0.797,0.805,0.810,0.811,0.812,0.812,0.813])
y_list_b.append([0.726,0.751,0.785,0.791,0.794,0.800,0.804,0.805,0.805,0.805,0.806])
y_list_a.append([1.29,1.38,1.76,1.82,1.86,1.88,1.90,1.94,1.93,1.92,1.93])
y_list_b.append([1.28,1.36,1.75,1.80,1.82,1.81,1.85,1.85,1.85,1.85,1.85])
y_list_a.append([0.45,0.34,0.36,0.71,0.79,0.63,0.52,0.35,0.55,0.27,0.24])
y_list_b.append([0.11,0.91,0.96,0.59,0.64,0.88,0.65,0.44,0.31,0.28,0.13])
plot_type_list.append(0)
plot_type_list.append(1)
plot_type_list.append(2)
plot_xaxis_list.append('nvariables')
plot_xaxis_list.append('nvariables')
plot_xaxis_list.append('nvariables')
plot_xaxis_title.append('Number of variables')
plot_xaxis_title.append('Number of variables')
plot_xaxis_title.append('Number of variables')

def make_plot(x_values, y_values_a, y_values_b, plot_type, plot_xaxis_short, plot_xaxis_title):
  is_auc = False
  is_csi = False
  is_ks = False
  if (plot_type==0):
    is_auc = True
  elif (plot_type==1):
    is_csi = True
  else:
    is_ks = True
  label_a = 'Training sample'
  label_b = 'Test sample'
  color_a = ROOT.kRed
  color_b = ROOT.kBlue
  plot_yaxis_title = 'AUC'
  log_y = False
  manual_y_axis = False
  y_axis_max = 1.0
  y_axis_min = 0.1*min(min(y_values_a),min(y_values_b))
  plot_type_str = 'auc'
  if is_csi:
    plot_yaxis_title = 'CSI'
    plot_type_str = 'csi'
  elif is_ks:
    label_a = 'Signal samples'
    label_b = 'Background samples'
    color_a = ROOT.kViolet
    color_b = ROOT.kCyan
    plot_yaxis_title = 'Kolmogorov-Smirnov p-value'
    log_y = True
    manual_y_axis = True
    plot_type_str = 'ks'
    if y_axis_min < 1.0e-8:
      y_axis_min = 1.0e-8
  plot_title = ''
  hist_title_str = plot_title+'; '+plot_xaxis_title+'; '+plot_yaxis_title
  
  #Create the canvas and graph
  ROOT.gStyle.SetOptStat(0)
  mycanvas = ROOT.TCanvas('canvas','',600,400)
  mycanvas.cd()
  if log_y:
    mycanvas.SetLogy()
  base_histo = ROOT.TH1D('base_histo',hist_title_str,100,min(x_values),max(x_values))
  base_histo.GetXaxis().CenterTitle(True)
  base_histo.GetXaxis().SetLabelOffset(0.014)
  base_histo.GetXaxis().SetTitleOffset(1.35)
  base_histo.GetXaxis().SetTickLength(-0.02)
  base_histo.GetYaxis().CenterTitle(True)
  base_histo.GetYaxis().SetLabelOffset(0.014)
  base_histo.GetYaxis().SetTitleOffset(1.35)
  base_histo.GetYaxis().SetTickLength(-0.01)
  if manual_y_axis:
    base_histo.SetMinimum(y_axis_min)
    base_histo.SetMaximum(y_axis_max)
  else:
    if log_y:
      base_histo.SetMinimum(0.1*min(min(y_values_a),min(y_values_b)))
      base_histo.SetMaximum(10.0*max(max(y_values_a),max(y_values_b)))
    else:
      base_histo.SetMinimum(0.9*min(min(y_values_a),min(y_values_b)))
      base_histo.SetMaximum(1.1*max(max(y_values_a),max(y_values_b)))
  base_histo.Draw()
  
  if len(x_values) != len(y_values_a) or len(x_values) != len(y_values_b):
    print('ERROR: unequal length lists')
  graph_a = ROOT.TGraph(len(x_values),array.array('d',x_values),array.array('d',y_values_a))
  graph_b = ROOT.TGraph(len(x_values),array.array('d',x_values),array.array('d',y_values_b))
  graph_a.SetLineColor(color_a)
  graph_b.SetLineColor(color_b)
  graph_a.SetLineWidth(2)
  graph_b.SetLineWidth(2)
  graph_a.Draw('SAME')
  graph_b.Draw('SAME')
  
  mylegend = ROOT.TLegend(0.6, 0.11, 0.89, 0.25)
  mylegend.AddEntry(graph_a,label_a,'l')
  mylegend.AddEntry(graph_b,label_b,'l')
  mylegend.SetTextSize(0.03)
  mylegend.SetBorderSize(0)
  mylegend.Draw('SAME')
  
  mycanvas.Draw()
  filename = 'overtrain_'+plot_type_str+'_'+plot_xaxis_short+'.pdf'
  mycanvas.SaveAs(filename)

for i in range(len(x_list)):
  make_plot(x_list[i],y_list_a[i],y_list_b[i],plot_type_list[i],plot_xaxis_list[i],plot_xaxis_title[i])
