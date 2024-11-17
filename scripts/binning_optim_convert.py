import pandas as pd
import uproot

def convert_csv_to_ntuple():
  '''Converts data stored as csv into ROOT ntuple (TTree)
  '''
  df_ggf = pd.read_csv('ntuples/GGF.csv')
  df_ggf = df_ggf.drop('Unnamed: 0',axis=1)
  df_ggf['type'] = 1
  df_vbf = pd.read_csv('ntuples/VBF.csv')
  df_vbf = df_vbf.drop('Unnamed: 0',axis=1)
  df_vbf['type'] = 2
  df_dyg = pd.read_csv('ntuples/SMZg.csv')
  df_dyg = df_dyg.drop('Unnamed: 0',axis=1)
  df_dyg['type'] = 3
  df_dyf = pd.read_csv('ntuples/DY.csv')
  df_dyf = df_dyf.drop('Unnamed: 0',axis=1)
  df_dyf['type'] = 4
  
  df_full = pd.concat([df_ggf, df_vbf], axis=0)
  df_full = pd.concat([df_full, df_dyg], axis=0)
  df_full = pd.concat([df_full, df_dyf], axis=0)
  
  with uproot.recreate('ntuples/binning_tree.root') as output_file:
    output_file['tree'] = df_full

if __name__=='__main__':
  convert_csv_to_ntuple()
