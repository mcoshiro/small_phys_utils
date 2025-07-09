from array import array
import ROOT

#setup colors
COLORS = [ROOT.TColor.GetColor('#3f90da'), ROOT.TColor.GetColor('#ffa90e'), 
          ROOT.TColor.GetColor('#bd1f01'), ROOT.TColor.GetColor('#832db6'), 
          ROOT.TColor.GetColor('#94a4a2'), ROOT.TColor.GetColor('#a96b59'),
          ROOT.TColor.GetColor('#e76300'), ROOT.TColor.GetColor('#b9ac70'),
          ROOT.TColor.GetColor('#717581'), ROOT.TColor.GetColor('#92dadd')]

def make_pie(title: str, vals: list[float], labels: list[str], output: str):
  """Makes a pie chart
  """
  c = ROOT.TCanvas('c','c',600,600)
  vals_array = array('d', vals)
  colors_array = array('i', COLORS[:len(vals)])
  pie_chart = ROOT.TPie('pie_chart', title, len(vals), vals_array, 
                        colors_array)
  for ientry in range(len(vals)):
    pie_chart.SetEntryLabel(ientry, labels[ientry])
  pie_chart.SetY(0.32)
  pie_chart.SetRadius(0.3);
  pie_chart.SetLabelsOffset(0.01);
  pie_chart.SetLabelFormat('%perc');
  pie_chart.Draw('nol <')
  legend = pie_chart.MakeLegend()
  legend.SetY1(0.6)
  legend.SetY2(0.9)
  legend.SetBorderSize(0)
  legend.Draw('same')
  c.SaveAs(output)

if __name__=='__main__':
  el_vals = [21.98, 5.2, 2.08, 0.39, 0.25]
  mu_vals = [29.48, 2.68, 0.28, 0.24, 0.08]
  labels = ['Pass dilepton trigger', '#geq 1 lepton not reconstructed',
            'Lepton fails ID/Iso', 'Lower leg fails p_{T} cut',
            'Upper leg fails p_{T} cut']
  make_pie('H#rightarrowZ#gamma#rightarrowee#gamma HLT decision', 
           el_vals, labels, 'plots/eltrigdecision.pdf')
  make_pie('H#rightarrowZ#gamma#rightarrow#mu#mu#gamma HLT decision', 
           mu_vals, labels, 'plots/mutrigdecision.pdf')

  
