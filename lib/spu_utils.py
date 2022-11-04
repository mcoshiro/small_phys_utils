"""@package docstring
Generic Python utilities for small_phys_utils
"""

def get_path(filename):
  '''Returns path to a file given filename
  
  Parameters:
  filename - filename possibly including path
  '''
  last_slash = filename.rfind('/')
  return filename[:last_slash]

def get_filename_in_directory(filename):
  '''Returns filename, removing any prefixed path
  
  Parameters:
  filename - filename possibly including path
  '''
  last_slash = filename.rfind('/')
  return filename[last_slash+1:]

