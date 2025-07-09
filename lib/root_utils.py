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
  columns = [(key.GetName(),key.GetTypeName()) for key in 
              tree.GetListOfLeaves()]
  in_file.Close()
  return columns

def get_entries(filename, tree_name, cuts=''):
  '''Returns number of entries in TTree
  
  Parameters:
  filename - filename to read
  tree_name - name of TTree in file
  cuts - only count entries passing these criteria
  '''
  in_file = ROOT.TFile(filename)
  tree = in_file.Get(tree_name)
  if (cuts==''):
    return tree.GetEntries()
  return tree.GetEntries(cuts)

def get_column_types(column_names, filename, tree_name):
  '''Returns C++ types for given columns in TTree

  column_names - list of names of branches
  filename - file containing TTree
  tree_name - TTree name
  '''
  columns = get_tree_columns(filename, tree_name)
  column_types = []
  for column_name in column_names:
    for column in columns:
      if column[0]==column_name:
        column_types.append(column[1])
        continue
  return column_types

def get_cpp_type(type_string):
  '''Returns C++ type given the type string returned by TLeaf::GetTypeName
  
  Parameters:
  type_string - type as specified by TLeaf::GetTypeName
  '''
  if (type_string=='Float_t'):
    return 'Float_t'
  if (type_string=='Double_t'):
    return 'Double_t'
  if (type_string=='Int_t'):
    return 'Int_t'
  else:
    raise TypeError('unknown type string')

def get_root_type(type_string):
  '''Returns ROOT type given the type string returned by TLeaf::GetTypeName
  
  Parameters:
  type_string - type as specified by TLeaf::GetTypeName
  '''
  if (type_string=='Float_t'):
    return 'F'
  if (type_string=='Double_t'):
    return 'D'
  if (type_string=='Int_t'):
    return 'I'
  else:
    raise TypeError('unknown type string')

