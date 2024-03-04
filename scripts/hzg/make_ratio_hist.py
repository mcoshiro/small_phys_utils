'''
Script that makes ratio of histograms
'''

from argparse import ArgumentParser
from array import array
import ROOT

if __name__=='__main__':

  #parse arguments
  argument_parser = ArgumentParser(prog='make_ratio_hist',
      description='Generates a ratio of histograms and saves to a cpp file')
  argument_parser.add_argument('-n','--num_filename')
  argument_parser.add_argument('-d','--den_filename')
  argument_parser.add_argument('-o','--output_filename',default='pteta_hist.hpp')
  args = argument_parser.parse_args()

  pt_bins_list = [15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,60,65,70,75,100,125,150,400]
  abseta_bins_list = [d*0.25 for d in range(11)]

  pt_bins_centered = [pt_bins_list[0]]+[(pt_bins_list[ipt]+pt_bins_list[ipt+1])/2.0 for ipt in range(len(pt_bins_list)-1)]+[pt_bins_list[-1]]
  abseta_bins_centered = [abseta_bins_list[0]]+[(abseta_bins_list[iabseta]+abseta_bins_list[iabseta+1])/2.0 for iabseta in range(len(abseta_bins_list)-1)]+[abseta_bins_list[-1]]

  pt_bins = array('f',pt_bins_centered)
  abseta_bins = array('f',abseta_bins_centered)

  df_num = ROOT.RDataFrame('tree',args.num_filename).Define('ph_sc_abseta','TMath::Abs(ph_sc_eta)')
  df_den = ROOT.RDataFrame('tree',args.den_filename).Define('ph_sc_abseta','TMath::Abs(ph_sc_eta)')
  hist_num_ptr = df_num.Filter('pair_mass>81&&pair_mass<101').Histo2D(
      ('num_hist','',len(pt_bins)-1,pt_bins,len(abseta_bins)-1,abseta_bins),'ph_et','ph_sc_abseta')
  hist_den_ptr = df_den.Filter('pair_mass>81&&pair_mass<101').Histo2D(
      ('den_hist','',len(pt_bins)-1,pt_bins,len(abseta_bins)-1,abseta_bins),'ph_et','ph_sc_abseta')
  hist_num = hist_num_ptr.GetValue()
  hist_den = hist_den_ptr.GetValue()
  hist_num.Scale(1.0/hist_num.Integral())
  hist_den.Scale(1.0/hist_den.Integral())
  hist_num.Divide(hist_den)

  #save picture
  ROOT.gStyle.SetOptStat(0)
  can = ROOT.TCanvas('c','c',600,400)
  can.SetLogx()
  hist_num.Draw('colz')
  can.SaveAs('ratio_hist.pdf')

  #write c++ file
  output_code = '\n'
  output_code += 'const std::vector<float> pt_bins = {'
  is_first = True
  for pt_edge in pt_bins_list:
    if is_first:
      is_first = False
    else:
      output_code += ','
    output_code += str(pt_edge)
  output_code += '};\n'
  output_code += 'const std::vector<float> eta_bins = {'
  is_first = True
  for eta_edge in abseta_bins_list:
    if is_first:
      is_first = False
    else:
      output_code += ','
    output_code += str(eta_edge)
  output_code += '};\n'
  output_code += 'const std::vector<std::vector<float>> hist_content = {'
  for ipt in range(len(pt_bins_list)):
    if ipt != 0:
      output_code += ','
    output_code += '{'
    for ieta in range(len(abseta_bins_list)):
      if ieta != 0:
        output_code += ','
      output_code += str(hist_num.GetBinContent(ipt+1,ieta+1))
    output_code += '}'
  output_code += '};\n\n'
  output_code += 'float get_hist_interp(float pt, float abseta) {\n'
  output_code += '  if (pt<'+str(pt_bins_list[0])+' || pt>'+str(pt_bins_list[-1])
  output_code += ' || abseta<'+str(abseta_bins_list[0])+' || abseta>'+str(abseta_bins_list[-1])+') '
  output_code += 'return 1.0;\n'
  output_code += '  unsigned ipt(0), iabseta(0);\n'
  output_code += '  while (!(pt >= pt_bins[ipt] && pt < pt_bins[ipt+1]))\n'
  output_code += '    ipt++;\n'
  output_code += '  while !(abseta >= abseta_bins[iabseta] && abseta < abseta_bins[iabseta+1])\n'
  output_code += '    iabseta++;\n'
  output_code += '  float pt_interp_loeta = (pt-pt_bins[ipt])/(pt_bins[ipt+1]-pt_bins[ipt])*(hist_content[ipt+1][iabseta]-hist_content[ipt][iabseta])+hist_content[ipt][iabseta];\n'
  output_code += '  float pt_interp_hieta = (pt-pt_bins[ipt])/(pt_bins[ipt+1]-pt_bins[ipt])*(hist_content[ipt+1][iabseta+1]-hist_content[ipt][iabseta+1])+hist_content[ipt][iabseta+1];\n'
  output_code += '  return (abseta-abseta_bins[iabseta])/(abseta_bins[iabseta]+abseta_bins[iabseta+1])*(pt_interp_hieta-pt_interp_loeta)+pt_interp_loeta;\n'
  output_code += '}\n'
  output_file = open(args.output_filename, 'w')
  output_file.write(output_code)
  output_file.close()

