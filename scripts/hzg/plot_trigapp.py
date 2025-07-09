#!/usr/bin/env python3
"""@package docstring
Make plots related to photon weights
"""

from root_plot_lib import RplPlot
from array import array
import ROOT
import gc

ROOT.gInterpreter.Declare("""
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using std::vector;

//only for electrons
float get_costheta(RVec<float> el_pt, RVec<float> el_eta, RVec<float> el_phi, RVec<int> el_charge, RVec<float> ll_pt, RVec<float> ll_eta, RVec<float> ll_phi, RVec<float> ll_m, RVec<int> ll_i1, RVec<int> ll_i2) {
  TLorentzVector p_lm, p_negz;
  if (el_charge[ll_i1[0]]==-1) {
    p_lm.SetPtEtaPhiM(el_pt[ll_i1[0]],el_eta[ll_i1[0]],el_phi[ll_i1[0]],0.000511);
  }
  else {
    p_lm.SetPtEtaPhiM(el_pt[ll_i2[0]],el_eta[ll_i2[0]],el_phi[ll_i2[0]],0.000511);
  }
  p_negz.SetPtEtaPhiM(ll_pt[0],ll_eta[0],ll_phi[0],ll_m[0]);
  p_negz.SetPxPyPzE(-1.0*p_negz.Px(),-1.0*p_negz.Py(),-1.0*p_negz.Pz(),p_negz.E());
  TVector3 boost = p_negz.BoostVector();
  p_negz.Boost(boost);
  p_lm.Boost(boost);
  TVector3 p3_negz = p_negz.Vect();
  TVector3 p3_lm = p_lm.Vect();
  float ret = p3_negz.Dot(p3_lm)/p3_negz.Mag()/p3_lm.Mag();
  if (ret > 1.0) ret = 1.0;
  if (ret < -1.0) ret = -1.0;
  return ret;
}

""")

#MC_FILE = '/data2/oshiro/ntuples/triggervalidate/*.root'
MC_FILE = '/data2/oshiro/ntuples/triggervalidate/raw_pico_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL18NanoAODv9__106X_upgrade2018_realistic_v16_L1v1-v2__230000__1574B1FB-8C40-A24E-B059-59A80F397A0F.root'

MU_DATA_FILES = ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_SingleMuon*.root',
              '/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_DoubleMuon*.root']

EL_DATA_FILES = ['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/data/merged_zgdata_ll/merged_raw_pico_ll_EGamma*.root']

LEAD_BINS_MU = [0.0,20.0,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,47.5,
                50.0,52.5,55.0,57.5,60.0,65.0,70.0,90.0]
SUBL_BINS_MU = [0.0,10.0,15.0,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,
                42.5,45.0,50.0,60.0,90.0]
LEAD_BINS_EL = [0.0,20.0,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,47.5,
                50.0,52.5,55.0,57.5,60.0,65.0,70.0,90.0]
SUBL_BINS_EL = [10.0,15.0,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,
                45.0,50.0,60.0,90.0]
#ZPT_BINS = [0.0,2.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,16.0,18.0,
#            20.0,22.0,24.0,28.0,32.0,36.0,40.0,45.0,50.0,55.0,60.0,70.0,80.0,
#            90.0,100.0,500.0]
#ZABSETA_BINS = [0.0,0.5,1.0,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,4.0,4.5,5.0,
#                5.5,6.0,10.0]
#ZPT_BINS = [0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,25.0,30.0,35.0,
#            40.0,50.0,60.0,80.0,100.0,500.0]
#ZABSETA_BINS = [0.0,1.0,2.0,2.5,2.75,3.0,3.25,3.5,4.0,5.0,6.0,10.0]
ZPT_BINS = [0.0,2.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,
            24.0,26.0,28.0,30.0,35.0,40.0,50.0,60.0,80.0,100.0,500.0]
ZABSETA_BINS = [0.0,1.0,2.0,2.5,2.75,3.0,3.25,3.5,4.0,5.0,6.0,10.0]
ZABSCOSTHETA_BINS = [0.0,0.2,0.4,0.6,0.8,0.9,0.95,1.0]

if __name__ == '__main__':

  channel = 'el'
  data_files = MU_DATA_FILES
  lead_bins = LEAD_BINS_MU
  subl_bins = SUBL_BINS_MU
  if channel=='el':
    data_files = EL_DATA_FILES
    lead_bins = LEAD_BINS_EL
    subl_bins = SUBL_BINS_EL

  #generate histograms
  #ROOT.EnableImplicitMT()
  df_sim = ROOT.RDataFrame('tree', MC_FILE)
  df_dat = ROOT.RDataFrame('tree', data_files)

  mc_norm = df_sim.Sum('w_lumi')
  mc_norm_str = str(mc_norm.GetValue())
  print(mc_norm_str)

  df_sim = df_sim.Define('w_lumi_year','w_lumi*6077.22*60000.0/'+mc_norm_str)
  df_sim = df_sim.Define('w_lumi_year_lep','w_lumi_year*w_lep*w_pu*w_prefire')
  if channel=='mu':
    df_sim = df_sim.Filter('nll>=1&&nmu==2&&ll_m[0]>85&&ll_m[0]<95')
    df_dat = df_dat.Filter('nll>=1&&nmu==2&&ll_m[0]>85&&ll_m[0]<95')
    df_sim = df_sim.Define('lead_lep_pt','mu_pt[0]')
    df_sim = df_sim.Define('subl_lep_pt','mu_pt[1]')
    df_dat = df_dat.Define('lead_lep_pt','mu_pt[0]')
    df_dat = df_dat.Define('subl_lep_pt','mu_pt[1]')
    df_sim = df_sim.Define('w_lumi_year_singlelep','w_lumi_year_lep*w_trig_singlemu')
    df_sim = df_sim.Define('w_lumi_year_doublelep','w_lumi_year_lep*w_trig_doublemu')
  else:
    df_sim = df_sim.Filter('nll>=1&&nel==2&&ll_m[0]>85&&ll_m[0]<95')
    df_dat = df_dat.Filter('nll>=1&&nel==2&&ll_m[0]>85&&ll_m[0]<95') 
    df_sim = df_sim.Filter('nll>=1&&ll_lepid[0]==11&&nel==2&&ll_m[0]>85&&ll_m[0]<95')
    df_dat = df_dat.Filter('nll>=1&&ll_lepid[0]==11&&nel==2&&ll_m[0]>85&&ll_m[0]<95') 
    df_sim = df_sim.Define('lead_lep_pt','el_pt[0]')
    df_sim = df_sim.Define('subl_lep_pt','el_pt[1]')
    df_dat = df_dat.Define('lead_lep_pt','el_pt[0]')
    df_dat = df_dat.Define('subl_lep_pt','el_pt[1]')
    #df_sim = df_sim.Define('lead_lep_pt','el_charge[0]==1 ? el_pt[0]:el_pt[1]')
    #df_sim = df_sim.Define('subl_lep_pt','el_charge[0]==1 ? el_pt[1]:el_pt[0]')
    #df_dat = df_dat.Define('lead_lep_pt','el_charge[0]==1 ? el_pt[0]:el_pt[1]')
    #df_dat = df_dat.Define('subl_lep_pt','el_charge[0]==1 ? el_pt[1]:el_pt[0]')
    df_sim = df_sim.Define('lead_lep_eta','el_eta[0]')
    df_sim = df_sim.Define('subl_lep_eta','el_eta[1]')
    df_dat = df_dat.Define('lead_lep_eta','el_eta[0]')
    df_dat = df_dat.Define('subl_lep_eta','el_eta[1]')
    df_sim = df_sim.Define('w_lumi_year_singlelep','w_lumi_year_lep*w_trig_singleel')
    df_sim = df_sim.Define('w_lumi_year_doublelep','w_lumi_year_lep*w_trig_doubleel')
  #df_sim = df_sim.Define('w_trig_orel_fix','isnan(w_trig_orel) || isinf(w_trig_orel) ? 1.0 : w_trig_orel')
  df_sim = df_sim.Define('w_lumi_year_lep_trig','w_lumi_year_lep*w_trig')
  df_sim = df_sim.Define('z_pt','ll_pt[0]')
  df_dat = df_dat.Define('z_pt','ll_pt[0]')
  df_sim = df_sim.Define('z_eta','ll_eta')
  df_dat = df_dat.Define('z_eta','ll_eta')
  df_sim = df_sim.Define('z_abseta','fabs(ll_eta[0])')
  df_dat = df_dat.Define('z_abseta','fabs(ll_eta[0])')
  df_sim = df_sim.Define('z_costheta','get_costheta(el_pt,el_eta,el_phi,el_charge,ll_pt,ll_eta,ll_phi,ll_m,ll_i1,ll_i2)')
  df_dat = df_dat.Define('z_costheta','get_costheta(el_pt,el_eta,el_phi,el_charge,ll_pt,ll_eta,ll_phi,ll_m,ll_i1,ll_i2)')
  df_sim = df_sim.Define('z_abscostheta','fabs(z_costheta)')
  df_dat = df_dat.Define('z_abscostheta','fabs(z_costheta)')
  df_sim = df_sim.Filter('fabs(lead_lep_eta)<1.4442||fabs(lead_lep_eta)>1.566')
  df_sim = df_sim.Filter('fabs(subl_lep_eta)<1.4442||fabs(subl_lep_eta)>1.566')
  df_dat = df_dat.Filter('fabs(lead_lep_eta)<1.4442||fabs(lead_lep_eta)>1.566')
  df_dat = df_dat.Filter('fabs(subl_lep_eta)<1.4442||fabs(subl_lep_eta)>1.566')
  df_sim = df_sim.Filter('fabs(lead_lep_eta)<2.2&&fabs(subl_lep_eta)<2.2')
  df_dat = df_dat.Filter('fabs(lead_lep_eta)<2.2&&fabs(subl_lep_eta)<2.2')
  #df_sim = df_sim.Filter('el_pt[0]>30&&el_pt[0]<35')
  #disp = df_sim.Display(('z_pt','z_abseta','z_abscostheta'),30)
  #disp.Print()
  #exit()

  leadlep_dat_ptr = {}
  subllep_dat_ptr = {}
  leadlep_sim_ptrs = {}
  subllep_sim_ptrs = {}
  loopvars = [('trig_single_mu','silep'),('trig_double_mu','dilep'),
              ('trig_single_mu||trig_double_mu','orlep')]
  lep_name = 'Muon'
  if channel=='el':
    loopvars = [('trig_single_el','silep'),('trig_double_el','dilep'),
                ('trig_single_el||trig_double_el','orlep')]
    lep_name = 'Electron'

  #Get Z pT corrections
  print('Generating Z pT corrections')
  df_dat_zpt = df_dat.Filter(loopvars[2][0])
  df_sim_zpt = df_sim.Filter(loopvars[2][0])
  ##temp cuts
  #df_dat_zpt = df_dat.Filter('el_pt[0]<30.0&&el_pt[0]>20.0')
  #df_sim_zpt = df_sim.Filter('el_pt[0]<30.0&&el_pt[0]>20.0')
  sim_zpt_ptr = df_sim_zpt.Histo2D(('sim_zpt','MC;Z p_{T} [GeV]; Z |#eta|',
      len(ZPT_BINS)-1,array('d',ZPT_BINS),
      #len(ZABSETA_BINS)-1,array('d',ZABSETA_BINS),
      len(ZABSCOSTHETA_BINS)-1,array('d',ZABSCOSTHETA_BINS)),
      'z_pt',
      #'z_abseta',
      'z_abscostheta',
      'w_lumi_year_lep_trig')
  dat_zpt_ptr = df_dat_zpt.Histo2D(('dat_zpt','MC;Z p_{T} [GeV]; Z |#eta|',
      len(ZPT_BINS)-1,array('d',ZPT_BINS),
      #len(ZABSETA_BINS)-1,array('d',ZABSETA_BINS),
      len(ZABSCOSTHETA_BINS)-1,array('d',ZABSCOSTHETA_BINS)),
      'z_pt',
      #'z_abseta',
      'z_abscostheta')
  ROOT.gStyle.SetOptStat(0)
  #c1 = ROOT.TCanvas()
  #c1.SetLogx(True)
  sim_zpt_hist = sim_zpt_ptr.GetPtr()
  dat_zpt_hist = dat_zpt_ptr.GetPtr()
  #sim_zpt_hist.Draw('colz')
  #c1.SaveAs('plots/zpt_sim.pdf')
  #dat_zpt_hist.Draw('colz')
  #c1.SaveAs('plots/zpt_dat.pdf')
  sim_zpt_hist.Scale(1.0/sim_zpt_hist.Integral())
  dat_zpt_hist.Scale(1.0/dat_zpt_hist.Integral())
  dat_zpt_hist.Divide(sim_zpt_hist)
  #dat_zpt_hist.SetMinimum(0.5)
  #dat_zpt_hist.SetMaximum(1.5)
  #dat_zpt_hist.Draw('colz text')
  #c1.SaveAs('plots/zptrwgt.pdf')
  ROOT.gDirectory.Add(dat_zpt_hist)
  ROOT.gInterpreter.Declare("""
  TH2D* dat_zpt_hist = static_cast<TH2D*>(gDirectory->Get("dat_zpt"));
  """)
  df_sim = df_sim.Define('w_zpt','(dat_zpt_hist->GetBinContent(dat_zpt_hist->FindBin(z_pt,z_abscostheta)))')
  df_sim = df_sim.Define('w_zpt_notrig','w_lumi_year_lep*w_zpt')
  df_sim = df_sim.Define('w_zpt_sitrig','w_lumi_year_singlelep*w_zpt')
  df_sim = df_sim.Define('w_zpt_ditrig','w_lumi_year_doublelep*w_zpt')
  df_sim = df_sim.Define('w_zpt_ortrig','w_lumi_year_lep_trig*w_zpt')

  print('Making validation histograms')
  #check variables after reweighting
  df_dat_val = df_dat.Filter(loopvars[2][0])
  df_sim_val = df_sim.Filter(loopvars[2][0])
  #df_sim_zpt = df_sim_zpt.Filter('el_pt[0]>25.0&&el_pt[0]<35.0')
  #df_dat_zpt = df_dat_zpt.Filter('el_pt[0]>25.0&&el_pt[0]<35.0')
  validation_ptrs = []
  validation_ptrs.append((df_dat_val.Histo1D(('ctheta_dat','Data; cos #theta',30,-1.0,1.0),'z_costheta'),
                          df_sim_val.Histo1D(('ctheta_sim','MC; cos #theta',30,-1.0,1.0),'z_costheta','w_zpt_ortrig'),
                          'ctheta'))
  validation_ptrs.append((df_dat_val.Histo1D(('eta_dat','Data; #eta',40,-8.0,8.0),'z_eta'),
                          df_sim_val.Histo1D(('eta_sim','MC; #eta',40,-8.0,8.0),'z_eta','w_zpt_ortrig'),
                          'eta'))
  validation_ptrs.append((df_dat_val.Histo1D(('pt_dat','Data; p_{T} [GeV]',30,0.0,80.0),'z_pt'),
                          df_sim_val.Histo1D(('pt_sim','MC; p_{T} [GeV]',30,0.0,80.0),'z_pt','w_zpt_ortrig'),
                          'pt'))
  gc.collect()
  gc.disable()
  for dat_ptr, sim_ptr, desc in validation_ptrs:
    costheta_plot = RplPlot()
    costheta_plot.lumi_data = [(60,13)]
    costheta_dat_hist = dat_ptr.GetPtr()
    costheta_sim_hist = sim_ptr.GetPtr()
    costheta_sim_hist.Scale(costheta_dat_hist.Integral()
                            /costheta_sim_hist.Integral())
    costheta_plot.plot_points(costheta_dat_hist,ROOT.kBlack)
    costheta_plot.plot_outline(costheta_sim_hist)
    costheta_plot.add_ratio('{}_dat'.format(desc),'{}_sim'.format(desc))
    costheta_plot.draw('plots/zpt_{}.pdf'.format(desc))
  gc.enable()

  print('Making comparison histograms')
  #df_dat = df_dat.Filter('fabs(el_eta[0])<1.5')
  #df_sim = df_sim.Filter('fabs(el_eta[0])<1.5')
  for trig, trig_string in loopvars:
    leadlep_sim_ptrs[trig_string] = []
    subllep_sim_ptrs[trig_string] = []
    df_dat_trig = df_dat.Filter(trig)
    df_sim_trig = df_sim.Filter(trig)
    leadlep_dat_ptr[trig_string] = df_dat_trig.Histo1D(
        ('leadlep_hist_'+trig_string,f'Data;Lead {lep_name} p_{{T}} [GeV]',
        len(lead_bins)-1,
        array('d',lead_bins)),'lead_lep_pt')
    subllep_dat_ptr[trig_string] = df_dat_trig.Histo1D(
        ('subllep_hist_'+trig_string,f'Data;Sublead {lep_name} p_{{T}} [GeV]',
        len(subl_bins)-1,
        array('d',subl_bins)),'subl_lep_pt')
    #leadlep_dat_ptr[trig_string] = df_dat_trig.Histo1D(
    #    ('leadlep_hist_'+trig_string,f'Data;Lead {lep_name} #eta',
    #    30, -2.5, 2.5), 'lead_lep_eta')
    #subllep_dat_ptr[trig_string] = df_dat_trig.Histo1D(
    #    ('subllep_hist_'+trig_string,f'Data;Sublead {lep_name} #eta',
    #    30, -2.5, 2.5), 'subl_lep_eta')
    for weight, w_string, sim_name in [
        ('w_zpt_notrig','lumi','MC (nominal)'),
        ('w_zpt_sitrig','silep','MC (w_single_trig)'),
        ('w_zpt_ditrig','dilep','MC (w_double_trig)'),
        ('w_zpt_ortrig','orlep','MC (w_comb_trig)')]:
      leadlep_sim_ptrs[trig_string].append(df_sim_trig.Histo1D(
          ('leadlep_hist_{}_{}'.format(trig_string,w_string),
          sim_name+f';Lead {lep_name} p_{{T}}[GeV]',len(lead_bins)-1,
          array('d',lead_bins)),
          'lead_lep_pt',weight))
      subllep_sim_ptrs[trig_string].append(df_sim_trig.Histo1D(
          ('subllep_hist_{}_{}'.format(trig_string,w_string),
          sim_name+f';Sublead {lep_name} p_{{T}}[GeV]',len(subl_bins)-1,
          array('d',subl_bins)),
          'subl_lep_pt',weight))
      #leadlep_sim_ptrs[trig_string].append(df_sim_trig.Histo1D(
      #    ('leadlep_hist_{}_{}'.format(trig_string,w_string),
      #    sim_name+f';Lead {lep_name} #eta',
      #    30, -2.5, 2.5), 'lead_lep_eta',weight))
      #subllep_sim_ptrs[trig_string].append(df_sim_trig.Histo1D(
      #    ('subllep_hist_{}_{}'.format(trig_string,w_string),
      #    sim_name+f';Sublead {lep_name} #eta',
      #    30, -2.5, 2.5), 'subl_lep_eta',weight))

  #write histograms to file for ease of use
  out_file = ROOT.TFile('temp_trighists.root','RECREATE')
  for trig_string in ['silep','dilep','orlep']:
    leadlep_dat_ptr[trig_string].Write()
    subllep_dat_ptr[trig_string].Write()
    for iweight in range(4):
      leadlep_sim_ptrs[trig_string][iweight].Write()
      subllep_sim_ptrs[trig_string][iweight].Write()
  out_file.Close()

  #generate sumw2
  for dat_ptr, sim_ptrs in [(leadlep_dat_ptr, leadlep_sim_ptrs),
                            (subllep_dat_ptr, subllep_sim_ptrs)]:
    for trig_string in ['silep','dilep','orlep']:
      dat_ptr[trig_string].Sumw2()
      for ibin in range(4):
        sim_ptrs[trig_string][ibin].Sumw2()

  #make plots
  gc.collect()
  gc.disable()
  for dat_ptr, sim_ptrs, leg in [
          (leadlep_dat_ptr, leadlep_sim_ptrs,'leadlep'),
          (subllep_dat_ptr, subllep_sim_ptrs,'subllep')]:
    for trig_string in ['silep','dilep','orlep']:
      plot = RplPlot()
      plot.lumi_data = [(60,13)]
      plot.y_title = 'Events/bin'
      plot.y_title_lower = 'MC/Data'
      for imc in range(4):
        sim_hist = sim_ptrs[trig_string][imc].GetPtr()
        sim_hist.Scale(dat_ptr[trig_string].GetPtr().Integral()
                       /sim_hist.Integral())
      plot.plot_points(dat_ptr[trig_string].GetPtr(),ROOT.kBlack)
      plot.plot_outline(sim_ptrs[trig_string][0].GetPtr())
      plot.plot_outline(sim_ptrs[trig_string][1].GetPtr())
      plot.plot_outline(sim_ptrs[trig_string][2].GetPtr())
      plot.plot_outline(sim_ptrs[trig_string][3].GetPtr())
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_lumi'.format(leg,trig_string),False)
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_silep'.format(leg,trig_string),False)
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_dilep'.format(leg,trig_string),False)
      plot.add_ratio('{}_hist_{}'.format(leg,trig_string),
                     '{}_hist_{}_orlep'.format(leg,trig_string),False)
      plot.draw('plots/trigappcomp_{}_{}.pdf'.format(leg,trig_string))
    ##for trig_string, trig_desc in [('silep','Single lepon trigger'),
    ##                               ('dilep','Dilepon trigger')]:
    #for trig_string, trig_desc in [('dilep','Dilepon trigger')]:
    #  plot = RplPlot()
    #  dat_hist = dat_ptr['orlep'].Clone('dathist')
    #  dat_hist.Divide(dat_ptr[trig_string].GetPtr())
    #  sim_hist_lumi = sim_ptrs['orlep'][0].Clone('simhist_lumi')
    #  sim_hist_lumi.Divide(sim_ptrs[trig_string][2].GetPtr())
    #  sim_hist_orlep = sim_ptrs['orlep'][3].Clone('simhist_orlep')
    #  sim_hist_orlep.Divide(sim_ptrs[trig_string][2].GetPtr())
    #  plot.lumi_data = [(60,13)]
    #  plot.y_title = 'OR of triggers/'+trig_desc
    #  plot.plot_points(dat_hist,ROOT.kBlack)
    #  plot.plot_outline(sim_hist_lumi)
    #  plot.plot_outline(sim_hist_orlep)
    #  plot.add_ratio('dathist','simhist_lumi')
    #  plot.add_ratio('dathist','simhist_orlep')
    #  plot.draw('plots/trigappcomp_ratio_{}_{}.pdf'.format(leg,trig_string))
  gc.enable()
