from argparse import ArgumentParser
from os.path import isdir
from subprocess import run

basedir = '/net/cms11/cms11r0/pico/'
nano_versions = ['NanoAODv9UCSB1/','NanoAODv9UCSB2/','NanoAODv9/',
                 'NanoAODv11/','NanoAODv11p9/','NanoAODv12/']
years = ['2016APV','2016','2017','2018','2022','2022EE','2023','2023BPix']
subfolders = ['data','mc']
subsubfolders = ['raw_pico','wgt_sums','corrections','unskimmed','skim_ll',
                 'skim_llg','merged_zgdata_ll','merged_zgmc_ll']

if __name__=='__main__':
  argument_parser = ArgumentParser(prog='clean_production',
      description='cleans nano2pico production')
  argument_parser.add_argument('production_name')
  argument_parser.add_argument('-r','--remove',action='store_true')
  argument_parser.add_argument('-f','--force',action='store_true')
  args = argument_parser.parse_args()

  for nano_version in nano_versions:
    prod_dir = basedir+nano_version+args.production_name+'/'
    if isdir(prod_dir):
      for year in years:
        year_dir = prod_dir+year+'/'
        if isdir(year_dir):
          for subfolder in subfolders:
            subfolder_dir = year_dir+subfolder+'/'
            if isdir(subfolder_dir):
              for subsubfolder in subsubfolders:
                subsubfolder_dir = subfolder_dir+subsubfolder
                if isdir(subsubfolder_dir):
                  flags = '-r'
                  if args.force:
                    flags = '-rf'
                  command = 'rm {} {}'.format(flags,subsubfolder_dir)
                  print(command)
                  if args.remove:
                    run(command.split())


              
      

  
