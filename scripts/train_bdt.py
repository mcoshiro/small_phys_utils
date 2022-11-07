"""@package docstring
Script to automate training and testing of (AdaBoost) BDTs

"""
import bdt_utils

if __name__=='__main__':
  bdt_options = bdt_utils.BdtOptions()
  bdt_options.MaxDepth = 2
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity','photon_pt_mass']
  base_name = 'shuffled_kinbdt_masscut'
  #bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options)
  bdt_utils.evaluate_bdt(base_name,variables,'w_lumi')
  #bdt_utils.clean_bdt(base_name)
