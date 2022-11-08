"""@package docstring
Script to automate training and testing of (AdaBoost) BDTs

"""
import bdt_utils

if __name__=='__main__':
  bdt_options = bdt_utils.BdtOptions()
  #bdt_options.BoostType = 'Grad'
  #bdt_options.NTrees = 1000
  #bdt_options.MaxDepth = 2
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity']
  #in order
  #all_variables = ['photon_mva','min_dR','pt_mass','photon_res','max_dR','cosTheta','photon_rapidity','costheta','l1_rapidity','phi','l2_rapidity']
  base_name = 'shuffled_kinbdt_masscut'
  tag = '_phpt0'
  evaluations = []
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)
  tag = '_phpt1'
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity','photon_pt_mass']
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)
  base_name = 'shuffled_kinbdt_masscut_phptl35'
  tag = '_phpt2'
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity']
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)
  base_name = 'shuffled_kinbdt_masscut_phptg35'
  tag = '_phpt3'
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)
  

  #ntrains = [5000,10000,25000,50000,100000,150000]
  #tag = '_grad'
  #for ntrain in ntrains:
  #  tag = '_ntrain'+str(ntrain)+'_grad'
  #  bdt_options.nTrain_Signal = ntrain
  #  bdt_options.nTrain_Background = ntrain
  #  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  #  evaluation = bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  #  evaluations.append(evaluation)
  #  evaluation.print_summary()
  #  bdt_utils.clean_bdt(base_name,tag)
  #for trial in zip(ntrains,evaluations):
  #  print('For '+str(trial[0])+' training events: ')
  #  trial[1].print_summary()
