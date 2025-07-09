
from array import array
from root_plot_lib import RplPlot
import ROOT

if __name__=='__main__':
  train_loss = [0.0638, 0.0634, 0.0633, 0.0633, 0.0633, 0.0632, 0.0632, 0.0632,
                0.0632, 0.0632, 0.0632, 0.0632, 0.0632, 0.0632, 0.0632, 0.0632,
                0.0632, 0.0632, 0.0632, 0.0632, 0.0632, 0.0632]
  valid_loss = [0.2132, 0.1854, 0.1255, 0.1023, 0.0854, 0.0788, 0.0702, 0.0684,
                0.0655, 0.0641, 0.0638, 0.0635, 0.0630, 0.0630, 0.0630, 0.0631,
                0.0630, 0.0631, 0.0630, 0.0631, 0.0631, 0.0631]
  nbins = len(train_loss)
  epoch = [i for i in range(nbins)]

  epoch_array = array('d',epoch)
  valid_array = array('d',valid_loss)
  train_array = array('d',train_loss)

  train_graph = ROOT.TGraph(nbins,epoch_array,train_array)
  valid_graph = ROOT.TGraph(nbins,epoch_array,valid_array)
  train_graph.SetTitle('Training loss; epoch')
  valid_graph.SetTitle('Validation loss; epoch')
  plot = RplPlot()
  plot.lumi_data = [(118.0,13)]
  plot.plot_graph(train_graph)
  plot.plot_graph(valid_graph)
  plot.y_title = 'BCE loss'
  plot.draw('plots/lossgraph.pdf')
