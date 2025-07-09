#!/usr/bin/env python3
"""@package docstring
Small utility for converting ROOT-generated PDFs to PNGs using pdftoppm 21.01.0

Usage: rootpdf_to_png.py PDF1
Usage: rootpdf_to_png.py PDF1 PDF2 ... [-f FIGURES_PER_ROW] [-o OUTPUT_NAME]

For whatever reason, running pdftoppm out of the box generates strange
whitespace when converting PDFs generated with CERN ROOT. This utility
circumvents this by first inserting the PDFs into a LaTeX document.

If multiple inputs are provided, one output png will be created with all the
included pdfs tiled. The f and o options can be used to specify the number
of figures per row and the output png name.

Binary Reqs: pdftoppm, lualatex
LaTeX Package Reqs: graphicx, float
"""

import sys
import subprocess

def strip_extension(filename):
  '''
  Function that removes extension from filename if an extension exists

  Param: filename - string of filename
  Returns: string that is filename with extension removed
  '''
  if '.' in filename:
    dot_pos = [i for i in range(len(filename)) if filename[i]=='.']
    return filename[0:dot_pos[-1]]
  else:
    return filename

if __name__=='__main__':
  if (len(sys.argv)<2):
    print('ERROR insufficient arguments.')
    print('usage: rootpdf_to_png.py <pdf1> <pdf2> ...')
  #parse arguments
  figures_per_row = 3
  output_name = 'combined_plots'
  input_filenames = []
  single_plot = False
  if (len(sys.argv)==2):
    figures_per_row = 1
    input_filenames.append(strip_extension(sys.argv[1]))
    output_name = input_filenames[0]
    single_plot = True
  else:
    row_number_flag = False
    output_name_flag = False
    for arg_idx in range(1,len(sys.argv)):
      if row_number_flag:
        figures_per_row = int(sys.argv[arg_idx])
        row_number_flag = False
      elif output_name_flag:
        output_name = strip_extension(sys.argv[arg_idx])
        output_name_flag = False
      elif sys.argv[arg_idx]=='-f':
        row_number_flag = True
      elif sys.argv[arg_idx]=='-o':
        output_name_flag = True
      else:
        #assume name of input
        input_filenames.append(strip_extension(sys.argv[arg_idx]))
  #calculate figure width
  figure_width = 1.0-0.01*figures_per_row
  figure_width = figure_width/figures_per_row
  figure_width_string = '{:.4f}'.format(figure_width)
  #write latex file
  latex_file = open('rootpdf_to_png_latexdoc.tex','w')
  #latex_file.write('\\documentclass[varwidth=true, border=0pt]{standalone}\n')
  latex_file.write('\\documentclass[10pt,oneside]{report}\n')
  latex_file.write('\\usepackage{graphicx,float}\n')
  latex_file.write('\\usepackage[active,tightpage]{preview}\n')
  latex_file.write('\\begin{document}\n')
  latex_file.write('\\begin{preview}\n')
  latex_file.write('\\begin{figure}[H]\n')
  for input_filename in input_filenames:
    latex_file.write('\\includegraphics[width='+figure_width_string)
    latex_file.write('\\textwidth]{'+input_filename+'.pdf}\n')
  latex_file.write('\\end{figure}')
  latex_file.write('\\end{preview}')
  latex_file.write('\\end{document}')
  latex_file.close()
  #compile latex document
  subprocess.run(['pdflatex','rootpdf_to_png_latexdoc.tex'])
  #generate and name png
  temp_filename = input_filenames[0]
  subprocess.run(['/data2/oshiro/linux/opt/xpdf-4.05/xpdf/pdftopng','-r',str(200*figures_per_row),'rootpdf_to_png_latexdoc.pdf',temp_filename])
  subprocess.run(['mv',temp_filename+'-000001.png',output_name+'.png'])
  #clean up
  subprocess.run(['rm','rootpdf_to_png_latexdoc.tex'])
  subprocess.run(['rm','rootpdf_to_png_latexdoc.aux'])
  subprocess.run(['rm','rootpdf_to_png_latexdoc.log'])
  if not single_plot:
    subprocess.run(['mv','rootpdf_to_png_latexdoc.pdf',output_name+'.pdf'])
  else:
    subprocess.run(['rm','rootpdf_to_png_latexdoc.pdf'])



