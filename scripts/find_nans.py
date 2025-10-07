#find picos with NaNs
#MO 2025-08-01

import ROOT
from glob import glob

pico_base_dirs = ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0/2016APV/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0/2016/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0/2017/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0/2018/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/']

pico_base_dirs = [
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-50to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-50to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/raw_pico_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/raw_pico_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root',
 '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_*.root']

pico_base_dirs = ['/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/raw_pico/raw_pico_DYGto2LG-1Jets_MLL-50_PTG-10to50_TuneCP5*.root']

pico_base_dirs = ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2016APV/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2016/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2017/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_redwood_v0_systsignal/2018/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2022/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2022EE/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2023/mc/raw_pico/',
                  '/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0_systsignal/2023BPix/mc/raw_pico/']

if __name__=='__main__':
  ibase = 1
  for base_dir in pico_base_dirs:
    print(f'In base_dir {base_dir} ({ibase}/{len(pico_base_dirs)})')
    slim_dir = base_dir.replace('raw_pico','merged_zgmc_llg')
    pico_list = glob(f'{slim_dir}*.root')
    #pico_list = glob(f'{base_dir}*.root')
    n_picos = len(pico_list)
    ipico = 1
    for pico_name in pico_list:
      print(f'Checking file {ipico}/{n_picos}')
      #print(f'{pico_name}')
      pico = ROOT.TFile(pico_name,'READ')
      nan_entries = pico.tree.GetEntries('isnan(sys_el[1])||isinf(sys_el[1])||isnan(sys_el[0])||isinf(sys_el[0])')
      if (nan_entries>0):
        print(f'Found nan in {pico_name}')
      pico.Close()
      ipico += 1
    ibase += 1


