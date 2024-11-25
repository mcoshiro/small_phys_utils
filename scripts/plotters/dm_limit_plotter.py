import ROOT, array

#Declare x values, y values, and uncertainties
#lz_mass = [9.0, 19.0, 29.0, 40.0, 55.0, 78.0, 110.0, 219.0, 308.0, 508.0, 716.0, 1008.0, 2000.0, 5000.0, 10000.0] #GeV
#lz_limit = [3.88e-46, 9.51e-48, 6.61e-48, 9.50e-48, 1.44e-47, 2.19e-47, 3.04e-47, 6.34e-47, 8.85e-47, 1.46e-46, 2.01e-46, 2.89e-46, 5.76e-46, 1.45e-45, 2.80e-45] #cm^2

interp_mass_points = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 23.0, 27.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 5000.0, 10000.0]

#lz 2022
lz_mass = [9.0, 11.0, 13.0, 16.0, 17.0, 19.0, 21.0, 23.0, 26.0, 29.0, 30.0, 32.0, 36.0, 40.0, 43.0, 46.0, 55.0, 65.0, 78.0, 91.0, 110.0, 129.0, 155.0, 182.0, 219.0, 256.0, 308.0, 361.0, 434.0, 508.0, 612.0, 716.0, 862.0, 1008.0, 1214.0, 1420.0, 1710.0, 2000.0, 5000.0, 10000.0]
lz_limit = [3.88141e-46, 1.05545e-46, 4.19377e-47, 1.52720e-47, 1.28035e-47, 9.50528e-48, 7.66097e-48, 6.94547e-48, 6.68266e-48, 6.61324e-48, 6.53944e-48, 6.91931e-48, 8.17366e-48, 9.49545e-48, 1.02165e-47, 1.13143e-47, 1.43786e-47, 1.76318e-47, 2.19292e-47, 2.53863e-47, 3.03790e-47, 3.64334e-47, 4.47954e-47, 5.15116e-47, 6.34341e-47, 7.36760e-47, 8.85464e-47, 1.01762e-46, 1.24217e-46, 1.46316e-46, 1.77592e-46, 2.00948e-46, 2.42561e-46, 2.89452e-46, 3.48215e-46, 4.09394e-46, 4.86699e-46, 5.75993e-46, 1.45172e-45, 2.79882e-45] 

#neutrino floor
#nf_mass = [0.2, 0.4, 1.0, 2.0, 4.0, 6.0, 7.0, 9.0, 10.0, 11.0, 13.0, 16.0, 20.0, 30.0, 50.0, 100.0, 1000.0, 10000.0]
#nf_limit = [6.0e-44, 5.0e-44, 1.0e-44, 2.8e-45, 3.0e-45, 4.0e-45, 2.0e-45, 7.0e-46, 2.5e-47, 4.5e-48, 3.5e-49, 1.5e-49, 1.0e-49, 8.0e-50, 1.0e-49, 2.5e-49, 2.0e-48, 2.0e-47]
nf_mass = [0.1, 0.5, 1.0, 2.0, 3.0, 6.0, 10.0, 30.0, 50.0, 1000.0, 10000.0]
nf_limit = [3.5e-45, 1.2e-44, 6.0e-45, 2.6e-45, 2.5e-45, 3.0e-45, 2.0e-48, 7.0e-50, 1.0e-49, 2.0e-48, 2.0e-47]

#pandax-4t (2021)
px_mass = [6.0, 10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 1000.0, 10000.0]
px_limit = [5.0e-45, 5.0e-46, 6.0e-47, 4.2e-47, 3.9e-47, 5.0e-47, 6.5e-47, 5.0e-46, 4.5e-45]

#pandax-4t 8B (2023)
pb_mass = [3.2, 5.0, 10.0]
pb_limit = [1.0e-42, 1.0e-44, 4.2e-46]

#pandax-4t low mass (2022)
ps_mass = [3.0, 5.0, 10.0]
ps_limit = [1.72e-43, 1.89e-44, 1.0e-44]

#xenon-1t (2018)
xn_mass = [6.0, 8.0, 10.0, 20.0, 30.0, 50.0, 100.0, 1000.0, 10000.0]
xn_limit = [2.5e-44, 2.2e-45, 5.5e-46, 6.0e-47, 4.1e-47, 5.1e-47, 9.0e-47, 9.0e-46, 8.0e-45]

#lz expected
lz_exp_mass = [10.0, 20.0, 30.0, 40.0, 100.0, 1000.0, 10000.0]
lz_exp_limit = [2.0e-46, 3.0e-48, 1.1e-48, 1.1e-48, 2.0e-48, 1.8e-47, 1.8e-46]

#xenon-1t S2 (2019)
xs_mass = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 18.0, 20.0, 30.0, 40.0, 50.0, 100.0]
xs_limit = [3.6e-41, 6.4e-43, 1.2e-43, 4.8e-44, 2.7e-44, 8.1e-45, 5.2e-45, 3.8e-45, 2.7e-45, 2.0e-45, 5.2e-46, 4.8e-46, 4.1e-46, 4.6e-46, 5.2e-46, 8.9e-46]

#xenon-1t 8B (2021)
xb_mass = [3.0, 5.0, 8.0, 12.0]
xb_limit = [3.2e-42, 1.5e-44, 1.5e-45, 3.1e-46]

#darkside-50 s2 (2022)
ds_mass = [1.2, 2.0, 6.0, 10.0]
ds_limit = [2.0e-41, 1.4e-42, 1.0e-43, 2.0e-43]

#darwin expected
dr_mass = [6.0, 10.0, 40.0, 100.0, 1000.0, 10000.0]
dr_limit = [6.0e-46, 7.0e-48, 2.1e-49, 3.2e-49, 2.0e-48, 2.0e-47]

interp_mass_array = array.array('d',interp_mass_points)

#Create the canvas and graph
mycanvas = ROOT.TCanvas('dm_canvas','dm_limits',600,400)
mycanvas.cd()
mycanvas.SetLogx()
mycanvas.SetLogy()
lz_mass_array = array.array('d',lz_mass)
lz_limit_array = array.array('d',lz_limit)
lz_graph = ROOT.TGraph(len(lz_mass), lz_mass_array, lz_limit_array)
lz_graph.SetTitle('; WIMP Mass [GeV/c^{2}]; SI WIMP-nucleon cross section [cm^{2}]')
lz_graph.SetLineWidth(2)
lz_graph.SetLineColor(ROOT.kOrange-3)
lz_graph.Draw('AC') 
lz_graph.GetXaxis().SetLimits(1.0e-1, 1.0e4) #GeV
lz_graph.GetXaxis().CenterTitle(True)
lz_graph.GetXaxis().SetTitleOffset(1.35)
lz_graph.GetXaxis().SetTickLength(-0.01)
lz_graph.GetYaxis().CenterTitle(True)
lz_graph.GetYaxis().SetTitleOffset(1.35)
lz_graph.GetYaxis().SetTickLength(-0.01)
lz_graph.GetHistogram().SetMinimum(1.0e-50) #cm^2
lz_graph.GetHistogram().SetMaximum(1.0e-41) #cm^2

nf_mass_array = array.array('d',nf_mass)
nf_limit_array = array.array('d',nf_limit)
nf_graph = ROOT.TGraph(len(nf_mass), nf_mass_array, nf_limit_array)
nf_graph.SetTitle('; WIMP Mass [GeV/c^{2}]; SI WIMP-nucleon cross section [cm^{2}]')
nf_graph.SetLineWidth(2)
nf_graph.SetLineColor(ROOT.kGray)
nf_graph.Draw('C SAME') 

def make_limit_graph(mass_points, limit_points, line_color, line_style):
  lim_mass_array = array.array('d',mass_points)
  lim_limit_array = array.array('d',limit_points)
  lim_graph = ROOT.TGraph(len(mass_points), lim_mass_array, lim_limit_array)
  lim_graph.SetTitle('; WIMP Mass [GeV/c^{2}]; SI WIMP-nucleon cross section [cm^{2}]')
  lim_graph.SetLineWidth(2)
  lim_graph.SetLineStyle(line_style)
  lim_graph.SetLineColor(line_color)
  return lim_graph

px_graph = make_limit_graph(px_mass, px_limit, ROOT.kBlue, ROOT.kSolid)
px_graph.Draw('C SAME')
xn_graph = make_limit_graph(xn_mass, xn_limit, ROOT.kRed, ROOT.kSolid)
xn_graph.Draw('C SAME')
lz_exp_graph = make_limit_graph(lz_exp_mass, lz_exp_limit, ROOT.kOrange-3, ROOT.kDashed)
lz_exp_graph.Draw('C SAME')
#xs_graph = make_limit_graph(xs_mass, xs_limit, ROOT.kRed, ROOT.kSolid)
#xs_graph.Draw('C SAME')
xb_graph = make_limit_graph(xb_mass, xb_limit, ROOT.kRed, ROOT.kSolid)
xb_graph.Draw('C SAME')
ds_graph = make_limit_graph(ds_mass, ds_limit, ROOT.kGreen+2, ROOT.kSolid)
ds_graph.Draw('C SAME')
dr_graph = make_limit_graph(dr_mass, dr_limit, ROOT.kRed+1, ROOT.kDashed)
dr_graph.Draw('C SAME')
pb_graph = make_limit_graph(pb_mass, pb_limit, ROOT.kBlue, ROOT.kSolid)
pb_graph.Draw('C SAME')
ps_graph = make_limit_graph(ps_mass, ps_limit, ROOT.kBlue, ROOT.kSolid)
ps_graph.Draw('C SAME')

#nf_graph_smoother = ROOT.TGraphSmooth('normal')
#nf_graph_smooth1 = nf_graph_smoother.Approx(nf_graph,'linear',len(interp_mass_points),interp_mass_array,6.0e-44,2.0e-47)
#nf_graph_smooth2 = nf_graph_smoother.SmoothKern(nf_graph_smooth1,'normal',2.0)
#nf_graph_smooth1.SetLineColor(ROOT.kGray)
#nf_graph_smooth1.Draw('L SAME') 
#nf_graph.GetXaxis().SetLimits(1.0e-1, 1.0e4) #GeV
#nf_graph.GetHistogram().SetMinimum(1.0e-50) #cm^2
#nf_graph.GetHistogram().SetMaximum(1.0e-37) #cm^2

mycanvas.Draw()
mycanvas.SaveAs('dm_limits.pdf')
