#!/usr/bin/env python3
"""@package docstring
Small utility to shuffle the order of events in a ROOT TTree

Usage: shuffle_tree.py my_file.root tree

The first argument is the filename, the second argument is the TTree name in the file
"""
import os
import ROOT
import root_utils
import spu_utils
import sys

def shuffle_events(filename, tree_name):
  '''
  Function to shuffle the TTree tree_name in file filename

  Parameters
  filename - ROOT file containing the tree of interest
  tree_name - name of TTree in file
  '''
  new_filename = spu_utils.get_path(filename)+'/shuffled_'+spu_utils.get_filename_in_directory(filename)
  columns = root_utils.get_tree_columns(filename, tree_name)
  out_file = open('shuffler_macro.cxx','w')
  out_file.write('void shuffler_macro() {\n')
  out_file.write('  TFile in_file("'+filename+'","READ");\n')
  out_file.write('  TTree* tree = static_cast<TTree*>(in_file.Get("'+tree_name+'"));\n')
  out_file.write('  TFile* out_file = new TFile("'+new_filename+'","RECREATE");\n')
  out_file.write('  TTree out_tree("'+tree_name+'","description");\n')
  out_file.write('\n')
  for column in columns:
      out_file.write('  '+root_utils.get_cpp_type(column[1])+' '+column[0]+' = 0;\n');
      if column[1]=='Float_t' or column[1]=='Double_t' or column[1]=='Int_t':
        out_file.write('  tree->SetBranchAddress("'+column[0]+'", &'+column[0]+');\n');
        out_file.write('  out_tree.Branch("'+column[0]+'", &'+column[0]+',"'+column[0]+'/'+root_utils.get_root_type(column[1])+'");\n');
  out_file.write('\n')
  out_file.write('  std::cout << "Constructing random sequence\\n";\n')
  out_file.write('  std::vector<long> event_order;\n')
  out_file.write('  for (long ievt = 0; ievt < tree->GetEntries(); ievt++) {\n')
  out_file.write('    event_order.push_back(ievt);\n')
  out_file.write('  }\n')
  out_file.write('  std::random_shuffle(event_order.begin(),event_order.end());\n')
  out_file.write('  std::cout << "Shuffling events\\n";\n')
  out_file.write('  tree->LoadBaskets();\n')
  out_file.write('  for (long ievt = 0; ievt < tree->GetEntries(); ievt++) {\n')
  out_file.write('    tree->GetEntry(event_order[ievt]);\n')
  out_file.write('    out_tree.Fill();\n')
  out_file.write('  }\n')
  out_file.write('  out_tree.Write();\n')
  out_file.write('  std::cout << "Wrote '+new_filename+'" << "\\n";\n')
  out_file.write('  out_file->Close();\n')
  out_file.write('}\n')
  out_file.close()
  os.system('root -l -q shuffler_macro.cxx')
  os.system('rm shuffler_macro.cxx')

if __name__=='__main__':
  if len(sys.argv)<3:
    print('Insufficient arguments')
    exit()
  shuffle_events(sys.argv[1],sys.argv[2])
