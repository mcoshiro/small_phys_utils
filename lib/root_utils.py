"""@package docstring
Small python utilities for small_phys_utils that depend on ROOT
"""
import ROOT

def get_tree_columns(filename, tree_name):
  '''Returns list of 2-tuples where first element is column(leaf) name and second element is C++ type
  
  Parameters:
  filename - filename to read
  tree_name - name of TTree in file
  '''
  in_file = ROOT.TFile(filename)
  tree = in_file.Get(tree_name)
  columns = [(key.GetName(),key.GetTypeName()) for key in tree.GetListOfLeaves()]
  in_file.Close()
  return columns

def get_cpp_type(type_string):
  '''Returns C++ type given the type string returned by TLeaf::GetTypeName
  
  Parameters:
  type_string - type as specified by TLeaf::GetTypeName
  '''
  if (type_string=='Float_t'):
    return 'Float_t'
  else:
    print('ERROR: unknown type')

def get_root_type(type_string):
  '''Returns ROOT type given the type string returned by TLeaf::GetTypeName
  
  Parameters:
  type_string - type as specified by TLeaf::GetTypeName
  '''
  if (type_string=='Float_t'):
    return 'F'
  else:
    print('ERROR: unknown type')

