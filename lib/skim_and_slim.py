'''
Function that generates skimmed/slimmed ROOT n-tuples, primarily for MVA training
'''
import ROOT

def write_ntuples(filenames, cuts, out_name, defines=[], tree_name='tree', branches=()):
  '''Generate ROOT n-tuple from existing n-tuple
  
  Parameters:
  filenames - list of filenames of ROOT n-tuples
  cuts - list of cuts expressed as strings in order they should be applied
  out_name - output filename of n-tuple
  defines - list of 2-tuples describing new branches to define in the format (name, expr)
            note that these must be defineable before cuts
  tree_name - name of tree in ROOT file
  branches - tuple of branches to save; if empty all branches are saved
  '''
  filenames_vec = ROOT.std.vector('string')()
  for filename in filenames:
    filenames_vec.push_back(filename)
  df = ROOT.RDataFrame('tree',filenames_vec)
  for define in defines:
    df = df.Define(define[0],define[1])
  for cut in cuts:
    df = df.Filter(cut)
  if (branches == ()):
    df.Snapshot(tree_name,'ntuples/'+out_name)
  else:
    df.Snapshot(tree_name,'ntuples/'+out_name,branches)
  print('Wrote ntuples/'+out_name)
