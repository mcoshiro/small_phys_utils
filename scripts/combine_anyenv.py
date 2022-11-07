#!/usr/bin/env python3
"""@package docstring
Script to run combine regardless of currently set environment variables
"""
import subprocess
import sys

if __name__=='__main__':
  args = ' '.join(sys.argv[1:])
  subprocess.run(('env -i scripts/combine_anyenv.sh '+args).split())
