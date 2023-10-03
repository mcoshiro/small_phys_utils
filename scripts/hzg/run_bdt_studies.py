"""@package docstring
Script to automate training and testing of (AdaBoost) BDTs

"""
import bdt_utils

def final_photon_pt_studies():
  bdt_options = bdt_utils.BdtOptions()
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity']
  base_name = 'shuffled_kinbdt_masscut'
  tag = '_phpt0'
  evaluations = []
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  tag = '_phpt1'
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity','photon_pt_mass']
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  base_name = 'shuffled_kinbdt_masscut_phptl35'
  tag = '_phpt2'
  variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity']
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()
  base_name = 'shuffled_kinbdt_masscut_phptg35'
  tag = '_phpt3'
  bdt_utils.train_bdt(base_name,variables,'w_lumi',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi',bdt_options,tag)).print_summary()

def run_assocprod_studies():
  bdt_options = bdt_utils.BdtOptions()
  bdt_options.NTrees = 100
  #tth lep
  #variables = ['photon_mva','pt_mass','min_dR','min_lepton_quality','assoc_lepton_quality','max_lepton_miniso','assoc_lepton_miniso','min_lepton_pt','assoc_lepton_pt','assoc_lepton_pdgid','ht','met','njet_f','nbdfl_f','nbdfm_f','nbdft_f','nlep_f']
  #base_name = 'shuffled_tthlep_bdt'
  #tag = ''
  #tth had
  #variables = ['photon_mva','pt_mass','min_dR','ht','met','njet_f','nbdfl_f','nbdfm_f','nbdft_f','top_mass']
  #base_name = 'shuffled_tthhad_bdt'
  #tag = ''
  #vhlep
  #variables = ['photon_mva','pt_mass','min_dR',
  #    'min_lepton_quality','max_lepton_miniso','min_lepton_pt',
  #    'assoc_lepton_quality','assoc_lepton_miniso','assoc_lepton_pt',
  #    'assoc_lepton_pdgid','assoc_lepton_mt','ht','met','njet_f','nlep_f']
  #base_name = 'shuffled_vhlep_bdt'
  #tag = ''
  #vhmet
  variables = ['photon_mva','pt_mass','min_dR',
              'dphi_h_met','met','njet_f','nlep_f']
  base_name = 'shuffled_vhmet_bdt'
  tag = ''
  #bdt_utils.train_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)
  #(bdt_utils.evaluate_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)
  bdt_utils.make_bdt_friend_ttree(base_name, variables, 'w_lumi_year', bdt_options)

def mem_studies():
  bdt_options = bdt_utils.BdtOptions()
  #variables = ['min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_rapidity']
  variables = ['cosTheta','costheta','phi']
  base_name = 'shuffled_memtree'
  tag = '_run2'
  #bdt_utils.train_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)
  #(bdt_utils.evaluate_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)).print_summary()
  bdt_utils.clean_bdt(base_name,tag)

def phmva_studies():
  bdt_options = bdt_utils.BdtOptions()
  #train photon idmva
  variables = ['photon_ptransverse','photon_rapidity','ph_irel','ph_r9','ph_sieie','ph_hoe','photon_mva'] #add photon_res
  base_name = 'shuffled_phtree'
  tag = '_ph'
  #bdt_utils.train_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)
  #(bdt_utils.evaluate_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)
  variables = ['photon_newmva','pt_mass','min_dR','max_dR','cosTheta','costheta','phi']
  base_name = 'shuffled_phtree_bdt'
  tag = '_comb'
  bdt_utils.train_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)
  (bdt_utils.evaluate_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)).print_summary()
  #bdt_utils.clean_bdt(base_name,tag)

if __name__=='__main__':
  phmva_studies()
  #run_assocprod_studies()
  #bdt_options = bdt_utils.BdtOptions()
  #variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity']
  #bdt_utils.optimize_binning('shuffled_vbftree',variables,'w_lumi_year',bdt_options,'_run2','shuffled_kinbdt_masscut_idmvacut',3,1.0,340.0/138.0)
  #bdt_options = bdt_utils.BdtOptions()
  ##bdt_options.BoostType = 'Grad'
  ##bdt_options.MaxDepth = 2
  #bdt_options.NTrees = 50
  ##bdt_options.MinNodeSize = 4.0
  #variables = ['photon_mva','min_dR','max_dR','pt_mass','cosTheta','costheta','phi','photon_res','photon_rapidity','l1_rapidity','l2_rapidity']
  ###in order
  ###variables = ['photon_mva','min_dR','pt_mass','photon_res','max_dR','cosTheta','photon_rapidity','costheta','l1_rapidity','l2_rapidity','phi']
  ###variables = ['pt_mass', 'min_dR', 'photon_mva', 'cosTheta', 'costheta', 'photon_rapidity', 'l2_rapidity', 'max_dR', 'l1_rapidity', 'photon_res', 'phi']
  ###variables = ['pt_mass', 'min_dR', 'photon_mva', 'cosTheta', 'photon_rapidity']
  ##base_name = 'shuffled_kinbdt_masscut_idmvacut'
  ##tag = '_grad'
  #tag = '_run2'
  ##bdt_utils.train_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)
  ##(bdt_utils.evaluate_bdt(base_name,variables,'w_lumi_year',bdt_options,tag)).print_summary()
  #bdt_utils.optimize_binning('shuffled_kinbdt_masscut_idmvacut',variables,'w_lumi_year',bdt_options,tag,'shuffled_kinbdt_masscut_idmvacut',5,2.0,340.0/138.0)
  ##bdt_utils.clean_bdt(base_name,tag)
  ##evaluations = []
  ##for i in range(len(variables)+1,1,-1):
  ##  tag = '_missingvar'+str(i)
  ##  sub_variables = variables[:i]
  ##  bdt_utils.train_bdt(base_name,sub_variables,'w_lumi_year',bdt_options,tag)
  ##  evaluation = bdt_utils.evaluate_bdt(base_name,sub_variables,'w_lumi_year',bdt_options,tag)
  ##  evaluations.append(evaluation)
  ##  evaluation.print_summary()
  ##  bdt_utils.clean_bdt(base_name,tag)
  ##for i in zip(reversed(variables),evaluations):
  ##  print('For missing variable '+i[0])
  ##  i[1].print_summary()
