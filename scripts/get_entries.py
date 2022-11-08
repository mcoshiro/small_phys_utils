#!/usr/bin/env python3
"""@package docstring
Small utility to get entries in a TTree

Usage: shuffle_tree.py my_file.root tree [cuts]

The first argument is the filename, the second argument is the TTree name in the file, the third argument, if present is a string of cuts to apply
"""
import root_utils
import sys

if __name__=='__main__':
  if len(sys.argv)<3:
    print('Insufficient arguments')
    exit()
  if len(sys.argv)<4:
    print(root_utils.get_entries(sys.argv[1],sys.argv[2]))
  else:
    print(root_utils.get_entries(sys.argv[1],sys.argv[2],sys.argv[3]))
