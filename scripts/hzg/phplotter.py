
from array import array
from root_plot_lib import RplPlot
import ROOT

eta_edges = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
ph_sfs = [1.0026715628389615, 0.9589399238464549, 0.9599955091347154, 0.9771583854366053, 0.9912104271076769, 0.9725093784891765, 0.9239176798268427, 0.965953120237465]
ph_uncs = [0.039219331721924225, 0.06475203962696309, 0.020404670112930956, 0.05235199407671402, 0.018441773709409325, 0.02260404821016293, 0.03135454521098947, 0.026469287058301615]
el_sfs = [0.9885971288901321,
         1.0161441664401092,
         0.9929995866856755,
         1.049892231087501,
         1.0285176897299955,
         1.013764959207935,
         0.8767926368416462,
         1.0153553529946038]
el_uncs = [0.02453792231415865,
           0.035565743458650774,
           0.050198762597549816,
           0.037788307666465584,
           0.057750780531404,
           0.01603256338717693,
           0.04373913284560969,
           0.021530143562328484]

if __name__=='__main__':
  eta_centers = []
  eta_widths = []
  nbins = len(eta_edges)-1
  for ieta in range(nbins):
    eta_centers.append((eta_edges[ieta+1]+eta_edges[ieta])/2.0)
    eta_widths.append((eta_edges[ieta+1]-eta_edges[ieta])/2.0)
  eta_array = array('d',eta_centers)
  eta_unc_array = array('d',eta_widths)
  el_array = array('d',el_sfs)
  el_unc_array = array('d',el_uncs)
  ph_array = array('d',ph_sfs)
  ph_unc_array = array('d',ph_uncs)
  el_graph = ROOT.TGraphErrors(nbins,eta_array,el_array,eta_unc_array,
                               el_unc_array)
  ph_graph = ROOT.TGraphErrors(nbins,eta_array,ph_array,eta_unc_array,
                               ph_unc_array)
  el_graph.SetTitle('Electrons')
  ph_graph.SetTitle('Photons')
  plot = RplPlot()
  plot.lumi_data = [(60.0,13)]
  plot.plot_graph(el_graph)
  plot.plot_graph(ph_graph)
  plot.x_title = 'Photon #eta'
  plot.y_title = 'Scale factor'
  plot.draw('plots/elphlowpt.pdf')
