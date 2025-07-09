import ROOT
import subprocess 
from typing import Tuple, List, Union

class RDataFrameSet:
  """Wrapper for collections of RDataFrames ex. when some columns need to
  be defined separately for each data set

  Attributes:
    dataframes: list of RDataframes
  """

  def __init__(self, dataframes: List[ROOT.RDataFrame]):
    """Constructor
    """
    self.dataframes = dataframes

  def RunMethod(self, method: str, *args):
    """Applies method to each dataframe and returns list of results
    """
    results = []
    for dataframe in self.dataframes:
      results.append(getattr(dataframe,method)(*args))
    return results

  def RunMethodEach(self, method: str, args):
    """Applies method to each dataframe and returns list of results
    """
    if len(defines) != len(self.dataframes):
      raise ValueError('RunMethodEach must have one arg set per dataframe')
    results = []
    for dataframe, arg in zip(self.dataframes,args):
      results.append(getattr(dataframe,method)(*arg))
    return results

  def Filter(self, *args) -> 'RDataFrameSet':
    """Applies filter to each dataframe and returns filtered dataframes
    """
    return RDataFrameSet(self.RunMethod('Filter',*args))

  def Define(self, *args) -> 'RDataFrameSet':
    """Applies define to each dataframe and returns dataframes
    """
    return RDataFrameSet(self.RunMethod('Define',*args))

  def Redefine(self, *args) -> 'RDataFrameSet':
    """Applies Redefine to each dataframe and returns dataframes
    """
    return RDataFrameSet(self.RunMethod('Redefine',*args))

  def DefineEach(self, defines: List) -> 'RDataFrameSet':
    """Applies define to each dataframe and returns dataframes
    """
    return RDataFrameSet(self.RunMethodEach('Define',defines))

  def Histo1D(self, *args) -> List[ROOT.RDF.RResultPtr]:
    """Calls Histo1D on each dataframe in the set
    """
    return self.RunMethod('Histo1D',*args)

  def Histo2D(self, *args) -> List[ROOT.RDF.RResultPtr]:
    """Calls Histo2D on each dataframe in the set
    """
    return self.RunMethod('Histo2D',*args)

  def Histo3D(self, *args) -> List[ROOT.RDF.RResultPtr]:
    """Calls Histo3D on each dataframe in the set
    """
    return self.RunMethod('Histo3D',*args)

  def Sum(self, *args) -> List[ROOT.RDF.RResultPtr]:
    """Calls Sum on each dataframe in the set
    """
    return self.RunMethod('Sum',*args)

def merge_hist_ptrs(hist_ptrs: List[ROOT.RDF.RResultPtr]) -> ROOT.TH1D:
  """Combines lists of histogram results pointers into a single histogram

  Args:
    hist_ptrs: list of RResultPtrs to histograms

  Returns:
    sum of histograms pointed to
  """
  sum_hist = hist_ptrs[0].GetPtr().Clone()
  for ihist in range(1, len(hist_ptrs)):
    sum_hist.Add(hist_ptrs[ihist].GetPtr())
  return sum_hist

def merge_float_ptrs(float_ptrs: List[ROOT.RDF.RResultPtr]) -> float:
  """Combines lists of RResultPtrs to float values into a single value

  Args:
    float_ptrs: list of RResultPtrs to floats

  Returns:
    sum of values
  """
  float_sum = 0.0
  for ptr in float_ptrs:
    value = ptr.GetValue()
    float_sum += value
  return float_sum

def write_dataframe_set(dataframe_set: RDataFrameSet, filename: str, 
                        columns: Tuple[str]=()):
  """Saves dataframe set to a single file

  Args:
    dataframe_set: set of dataframes
    filename: name of output file, should end in 'root'
    columns: columns to save; empty tuple saves all columns
  """
  filename_base = filename[:-5]
  if (columns==()):
    columns = ''
  idf = 0
  for df in dataframe_set.dataframes:
    df.Snapshot('tree', f'{filename_base}_{idf}.root', columns)
    idf += 1
  #merge outputs into a single file and clean
  combined_df = ROOT.RDataFrame('tree',f'{filename_base}_*.root')
  combined_df.Snapshot('tree',filename)
  for idf in range(len(dataframe_set.dataframes)):
    subprocess.run((f'rm {filename_base}_{idf}.root').split())

def merge_rdataframesets(dataframe_sets: List[RDataFrameSet]) -> RDataFrameSet:
  """Merges RDataFrameSets together

  Args:
    dataframe_sets: dataframesets to merge

  Returns:
    returns single merged RDataFrameSet
  """
  dataframe_list = []
  for dataframe_set in dataframe_sets:
    dataframe_list += dataframe_set.dataframes
  return RDataFrameSet(dataframe_list)
