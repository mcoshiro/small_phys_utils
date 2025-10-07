//this file contains a large number of utilities for inclusion in Python
//python/RDataFrame

#include <cmath>
#include <vector>

#include "ROOT/RVec.hxx"

#include "rw_mmp_good5.hpp"
#include "rw_mmp_r3_okay6.hpp"

const rw_mmp dnn_r2;
const rw_mmp_r3 dnn_r3;

template <class C>
using RVec = ROOT::VecOps::RVec<C>;

bool isgood_hornveto(bool isgood, float pt, float eta) {
  return isgood && (pt > 50 || fabs(eta)<2.5 || fabs(eta)>3.0);
}

float get_photon_sf_run2(std::vector<float> vars, int pflavor) {
  if (pflavor!=1) return 1.0;
  //input pt, abseta, idmva, res
  float dnn_output = dnn_r2.evaluate(vars);
  float sf = dnn_output/(1.0-dnn_output);
  //if (sf > 5.0) return 5.0;
  return sf;
}

float get_photon_sf_run3(std::vector<float> vars, int pflavor) {
  if (pflavor!=1) return 1.0;
  //input pt, abseta, idmva, res
  float dnn_output = dnn_r3.evaluate(vars);
  float sf = dnn_output/(1.0-dnn_output);
  //if (sf > 5.0) return 5.0;
  return sf;
}

double deltaPhi(double phi1, double phi2){
  const double PI = acos(-1.);
  double dphi = fmod(fabs(phi2-phi1), 2.*PI);
  return dphi>PI ? 2.*PI-dphi : dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
  return hypot(deltaPhi(phi1, phi2), eta2-eta1);
}

bool pass_trigs(bool trig_single_el, bool trig_single_mu, bool trig_double_el,
                bool trig_double_mu, int nel, int nmu, RVec<float> el_pt,
                RVec<float> mu_pt, float year) {
  if (trig_single_el && nel>=1) {
    if (year < 2016.99 && el_pt[0]>30) return true;
    if (year > 2016.99 && el_pt[0]>35) return true;
  }
  if (trig_double_el && nel>=2) {
    if (el_pt[0]>25 && el_pt[1]>15) return true;
  }
  if (trig_single_mu && nmu>=1) {
    if (year > 2016.99 && year < 2017.99) {
      if (mu_pt[0]>28) return true;
    }
    else {
      if (mu_pt[0]>25) return true;
    }
  }
  if (trig_double_mu && nmu>=2) {
    if (mu_pt[0]>20 && mu_pt[1]>10) return true;
  }
  return false;
}

bool photon_isjet(int nphoton, RVec<float> photon_eta, RVec<float> photon_phi,
                  RVec<float> mc_pt, RVec<float> mc_eta, RVec<float> mc_phi, 
                  RVec<int> mc_id, RVec<int> mc_statusflag) {
  if (nphoton < 1) return false;
  for (unsigned imc = 0; imc < mc_eta.size(); imc++) {
    //check pt>15 and fromHardProcess
    if (mc_pt[imc] > 15 && (mc_statusflag[imc] & 0x100) != 0) {
      //skip promptly decaying particles (W Z t H)
      int abs_mc_id = abs(mc_id[imc]);
      if (abs_mc_id == 6 || abs_mc_id == 23 || abs_mc_id == 24
          || abs_mc_id ==25) continue;
      //use large radius to capture ex. fragmentation photons in jets
      if (deltaR(photon_eta[0],photon_phi[0],mc_eta[imc],mc_phi[imc])<0.4) {
        return true;
      }
    }
  }
  return false;
}

float get_photon_mht_dphi(RVec<float> photon_pt, RVec<float> photon_phi, 
    RVec<bool> photon_sig, RVec<float> el_pt, RVec<float> el_phi, 
    RVec<bool> el_sig, RVec<float> mu_pt, RVec<float> mu_phi, 
    RVec<bool> mu_sig, RVec<float> jet_pt, RVec<float> jet_eta, 
    RVec<float> jet_phi, RVec<bool> jet_isgood) {
  float mht_x(0.0), mht_y(0.0);
  for (unsigned iel = 0; iel < el_sig.size(); iel++) {
    if (el_sig.at(iel)) {
      float pt = el_pt.at(iel);
      mht_x += pt*cos(el_phi.at(iel));
      mht_y += pt*sin(el_phi.at(iel));
    }
  }
  for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
    if (mu_sig.at(imu)) {
      float pt = mu_pt.at(imu);
      mht_x += pt*cos(mu_phi.at(imu));
      mht_y += pt*sin(mu_phi.at(imu));
    }
  }
  for (unsigned iph = 0; iph < photon_sig.size(); iph++) {
    if (photon_sig.at(iph)) {
      float pt = photon_pt.at(iph);
      mht_x += pt*cos(photon_phi.at(iph));
      mht_y += pt*sin(photon_phi.at(iph));
    }
  }
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      float pt = jet_pt.at(ijet);
      mht_x += pt*cos(jet_phi.at(ijet));
      mht_y += pt*sin(jet_phi.at(ijet));
    }
  }
  float mht_phi = atan2(-1.0*mht_y, -1.0*mht_x);
  return deltaPhi(photon_phi[0], mht_phi);
}

float get_mht(RVec<float> photon_pt, RVec<float> photon_phi, 
    RVec<bool> photon_sig, RVec<float> el_pt, RVec<float> el_phi, 
    RVec<bool> el_sig, RVec<float> mu_pt, RVec<float> mu_phi, 
    RVec<bool> mu_sig, RVec<float> jet_pt, RVec<float> jet_eta, 
    RVec<float> jet_phi, RVec<bool> jet_isgood) {
  float mht_x(0.0), mht_y(0.0);
  for (unsigned iel = 0; iel < el_sig.size(); iel++) {
    if (el_sig.at(iel)) {
      float pt = el_pt.at(iel);
      mht_x += pt*cos(el_phi.at(iel));
      mht_y += pt*sin(el_phi.at(iel));
    }
  }
  for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
    if (mu_sig.at(imu)) {
      float pt = mu_pt.at(imu);
      mht_x += pt*cos(mu_phi.at(imu));
      mht_y += pt*sin(mu_phi.at(imu));
    }
  }
  for (unsigned iph = 0; iph < photon_sig.size(); iph++) {
    if (photon_sig.at(iph)) {
      float pt = photon_pt.at(iph);
      mht_x += pt*cos(photon_phi.at(iph));
      mht_y += pt*sin(photon_phi.at(iph));
    }
  }
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      float pt = jet_pt.at(ijet);
      mht_x += pt*cos(jet_phi.at(ijet));
      mht_y += pt*sin(jet_phi.at(ijet));
    }
  }
  return sqrt(mht_x*mht_x+mht_y*mht_y);
}

float get_ht(RVec<float> photon_pt, RVec<float> photon_phi, 
    RVec<bool> photon_sig, RVec<float> el_pt, RVec<float> el_phi, 
    RVec<bool> el_sig, RVec<float> mu_pt, RVec<float> mu_phi, 
    RVec<bool> mu_sig, RVec<float> jet_pt, RVec<float> jet_phi,  
    RVec<float> jet_eta, RVec<bool> jet_isgood) {
  float ht = 0.0;
  for (unsigned iel = 0; iel < el_sig.size(); iel++) {
    if (el_sig.at(iel)) {
      ht += el_pt.at(iel);
    }
  }
  for (unsigned imu = 0; imu < mu_sig.size(); imu++) {
    if (mu_sig.at(imu)) {
      ht += mu_pt.at(imu);
    }
  }
  for (unsigned iph = 0; iph < photon_sig.size(); iph++) {
    if (photon_sig.at(iph)) {
      ht += photon_pt.at(iph);
    }
  }
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      ht += jet_pt.at(ijet);
    }
  }
  return ht;
}

float get_w_jet(float year, int type, int njet, bool photon_isjet) {
  //DY+jet
  if (type >= 6000 && type < 7000 && photon_isjet) {
    //run 2
    if (year < 2020) {
      if (njet==0) return 0.9111605383159621;
      else if (njet==1) return 1.1376003438448439; 
      else if (njet==2) return 1.713219935063525; 
      else return 3.1788706672250915;
    }
    //run 3
    else {
      if (njet==0) return 0.9728015266693347; 
      else if (njet==1) return 1.0002059092739153; 
      else if (njet==2) return 1.3489124392588017; 
      else return 2.9050241586119054;
    }
  }
  //DY+PU
  if (type >= 6000 && type < 7000 && photon_isjet) {
    //run 2
    if (year < 2020) {
      if (njet==0) return 1.0258646318505702;
      else if (njet==1) return 0.9269326556305092;
      else if (njet==2) return 0.7428198678341376;
      else return 0.5975578113064747;
    }
    //run 3
    else {
      if (njet==0) return 1.0026662278276175;
      else if (njet==1) return 0.9723535601861952;
      else if (njet==2) return 0.9755613569795534;
      else return 0.767196839119967;
    }
  }
  //ZG
  else if (type >= 17000 && type < 18000) {
    //run 2
    if (year < 2020) {
      if (njet==0) return 1.002995001426364;
      else if (njet==1) return 0.9112824642345808; 
      else if (njet==2) return 1.0720911090823475; 
      else return 1.7016023204619064;
    }
    //run 3
    else {
      if (njet==0) return 1.0265093217058043; 
      else if (njet==1) return 0.885799943406867;
      else if (njet==2) return 1.0100105581168564;
      else return 1.5602249486008501;
    }
  }
  return 1.0;
}

float get_lead_jet_pt(RVec<bool> jet_isgood, RVec<float> jet_pt, 
                      RVec<float> jet_eta) {
  float max_pt = -999.0;
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      if (jet_pt[ijet] > max_pt) {
        max_pt = jet_pt[ijet];
      }
    }
  }
  return max_pt;
}

float get_lead_jet_eta(RVec<bool> jet_isgood, RVec<float> jet_pt, 
                       RVec<float> jet_eta) {
  float max_pt = -999.0;
  float lead_eta = -999.0;
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      if (jet_pt[ijet] > max_pt) {
        max_pt = jet_pt[ijet];
        lead_eta = jet_eta[ijet];
      }
    }
  }
  return lead_eta;
}

float get_sublead_jet_pt(RVec<bool> jet_isgood, RVec<float> jet_pt, 
                         RVec<float> jet_eta) {
  float max_pt = -999.0;
  float sublead_pt = -999.0;
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      if (jet_pt[ijet] > max_pt) {
        sublead_pt = max_pt;
        max_pt = jet_pt[ijet];
      }
      else if (jet_pt[ijet] > sublead_pt) {
        sublead_pt = jet_pt[ijet];
      }
    }
  }
  return sublead_pt;
}

float get_sublead_jet_eta(RVec<bool> jet_isgood, RVec<float> jet_pt, 
                          RVec<float> jet_eta) {
  float max_pt = -999.0;
  float sublead_pt = -999.0;
  float lead_eta = -999.0;
  float sublead_eta = -999.0;
  for (unsigned ijet = 0; ijet < jet_isgood.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      if (jet_pt[ijet] > max_pt) {
        sublead_pt = max_pt;
        sublead_eta = lead_eta;
        max_pt = jet_pt[ijet];
        lead_eta = jet_eta[ijet];
      }
      else if (jet_pt[ijet] > sublead_pt) {
        sublead_pt = jet_pt[ijet];
        sublead_eta = jet_eta[ijet];
      }
    }
  }
  return sublead_eta;
}

//apply weights for photons with pT<20 in pinnacles_v0
float get_w_photon_lowpt(RVec<bool> photon_sig, RVec<int> photon_pflavor, 
                         RVec<float> photon_pt, RVec<float> photon_eta,
                         float year) {
  if (year>2022.9) return 1.0;
  float w = 1.0;
  for (unsigned iph = 0; iph < photon_pflavor.size(); iph++) {
    if (photon_pflavor[iph]==1 && photon_pt[iph]>15 && photon_pt[iph]<20) {
      float eta = photon_eta[iph];
      bool sig = photon_sig[iph];
      if (year > 2015.9 && year < 2016.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0637626399875102;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.980941745169525;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9713806263485965; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 1.016705800510692;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 1.013255390999418;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 1.004959318390466; 
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9028498618621597;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 1.0078660115088238;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.6458359529885689;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.0620333583712274;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1284296993607925;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 0.8791299133785585;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 0.9078126534514417;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 0.9790640106166821;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.3181105865876797;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 0.9540851991993182;
      }
      else if (year > 2016.4 && year < 2016.6) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 0.9285812754096276;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 1.0339806085618999;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9610075707366368; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9631877412261108;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 1.0046065422234682;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 1.0217023475344034;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9611024969289972;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 0.9499203568038671;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 1.3690933970619903;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 0.8979370999676594;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1627555236944658;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.248132746431878;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 0.9690851995285417;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 0.9137491261701439;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.1156632360862162;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 1.2645235937513208;
      }
      else if (year > 2016.9 && year < 2017.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0038600498776258;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.973441200861423;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9566560634147776; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9949802364843384;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9822339389814418;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 0.9735398543153528;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9589739931516128;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 1.0588456176736205;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.9850895146655023;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.0711514790710508;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1481185013060329;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.028288938158254;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.0985109885634243;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 1.0891240442031525;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.1097718747709093;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 0.7665019031151598;
      }
      else if (year > 2017.9 && year < 2018.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0026715628389615;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.9589399238464549;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9599955091347154; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9771583854366053;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9912104271076769;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 0.9725093784891765;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9239176798268427;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 0.965953120237465;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.9897845251239432;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.105854764771766;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1355540272325784;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.1233277647359423;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.0459648424357058;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 1.0901855371506401;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.1947301894556555;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 1.1346617949968554;
      }
      else if (year > 2021.9 && year < 2022.1) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 0.9952717636788198;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.9969438529681379;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9101314138722953; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9696422862626193;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9857683103368693;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 1.0480699104733633;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 1.0242583178501214;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 1.0899023852844136;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 1.0144283882263259;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.0106494907269516;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.2020556700335518;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.0949431857039922;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.0471284971827695;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 0.8866928706108905;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 0.916163665320336;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 0.6775232430131973;
      }
      else if (year > 2022.4 && year < 2022.6) {
        if (sig && eta > -2.5 && eta < -2.0)      w *= 1.0052678054133821;
        else if (sig && eta > -2.0 && eta < -1.5) w *= 0.962797814171453;
        else if (sig && eta > -1.5 && eta < -0.8) w *= 0.9263783918441049; 
        else if (sig && eta > -0.8 && eta < 0.0)  w *= 0.9374179718000482;
        else if (sig && eta > 0.0 && eta < 0.8)   w *= 0.9664231011166211;
        else if (sig && eta > 0.8 && eta < 1.5)   w *= 0.961355060779199;
        else if (sig && eta > 1.5 && eta < 2.0)   w *= 0.9799981025513707;
        else if (sig && eta > 2.0 && eta < 2.5)   w *= 0.9599848186660465;
        else if (!sig && eta > -2.5 && eta < -2.0) w *= 0.9825607111336123;
        else if (!sig && eta > -2.0 && eta < -1.5) w *= 1.1309203201239058;
        else if (!sig && eta > -1.5 && eta < -0.8) w *= 1.1664761752729338;
        else if (!sig && eta > -0.8 && eta < 0.0)  w *= 1.1990325160412474;
        else if (!sig && eta > 0.0 && eta < 0.8)   w *= 1.1100673514134627;
        else if (!sig && eta > 0.8 && eta < 1.5)   w *= 1.0907780198621746;
        else if (!sig && eta > 1.5 && eta < 2.0)   w *= 1.073268543888653;
        else if (!sig && eta > 2.0 && eta < 2.5)   w *= 1.1424960034725913;
      }
    }
  }
  return w;
}

float get_w_fakephoton(float year, bool this_photon_isjet, float ph_pt, 
                       float ph_abseta, int photon_pflavor) {
  static const std::vector<float> fakeph_ptbins = {15.0,20.0,500.0};
  static const std::vector<float> fakeph_absetabins = {0.0,1.5,2.0,2.5};

  const std::vector<float> fakeph_jet_sfs_run2 = {
    0.87784536,0.86090972,0.83359961,1.00905351,
    1.20656524,1.15353438,1.28262424,1.15435585,
    0.52378337,0.87839728,0.63154299,1.13529129};
  const std::vector<float> fakeph_pu_sfs_run2 = {
    1.38054607,1.29143011,1.22872195,1.26953528,
    1.11197443,1.34174906,1.19494895,1.25048567,
    2.59609818,1.86582932,3.16422829,1.41518427};
  const std::vector<float> fakeph_jet_sfs_run3 = {
    0.73009776,1.02800363,0.76673460,1.29045474,
    1.03924426,1.12586378,0.90955539,1.10089699,
    0.74353182,0.34798856,0.79232260,0.69352779};
  const std::vector<float> fakeph_pu_sfs_run3 = {
    1.26284615,1.17956303,1.28810042,1.39744805,
    1.45429152,1.37924462,1.47740148,1.47429645,
    2.11270871,3.20168632,1.77614373,1.49915009};

  if (photon_pflavor == 1) return 1.0;
  const std::vector<float> *fakeph_sfs;
  int run = 2;
  if (year > 2020.0) run = 3;
  if (run == 2 && this_photon_isjet) fakeph_sfs = &fakeph_jet_sfs_run2;
  else if (run == 2 && !this_photon_isjet) fakeph_sfs = &fakeph_pu_sfs_run2;
  else if (run == 3 && this_photon_isjet) fakeph_sfs = &fakeph_jet_sfs_run3;
  else if (run == 3 && !this_photon_isjet) fakeph_sfs = &fakeph_pu_sfs_run3;
  else return 1.0;
  for (unsigned ipt = 0; ipt < (fakeph_ptbins.size()-1); ipt++) {
    for (unsigned ieta = 0; ieta < (fakeph_absetabins.size()-1); ieta++) {
      if (ph_pt >= fakeph_ptbins[ipt] && ph_pt < fakeph_ptbins[ipt+1] 
          && ph_abseta >= fakeph_absetabins[ieta] 
          && ph_abseta < fakeph_absetabins[ieta+1]) {
        return fakeph_sfs->at(ipt*(fakeph_absetabins.size()-1)+ieta);
        break;
      }
    }
  }
  return 1.0;
}

//Gets fake photon weight
float get_w_fakephoton(float year, bool this_photon_isjet, float ph_pt,
                       float ph_abseta, float ph_idmva, float ph_res,
                       int photon_pflavor) {
  static const std::vector<float> fakeph_ptbins = {15.0,20.0,30.0,500.0};
  static const std::vector<float> fakeph_absetabins = {0.0,0.8,1.5,2.0,2.5};
  static const std::vector<float> fakeph_idbins = {0.14,0.5,0.8,0.9,1.0};
  static const std::vector<float> fakeph_resbins = {0.0,0.025,0.05,1.0};

  static const std::vector<float> fakeph_jet_sfs_run2 = {
    0.92932809,0.90764058,0.89241154,0.98757244,
    1.09096266,1.03151791,1.20695159,1.17096665,
    0.82937735,1.11916178,0.78897045,0.97441101};
  static const std::vector<float> fakeph_pu_sfs_run2 = {
    1.40036539,1.31645123,1.18902640,1.21409934,
    1.14731318,1.38274917,1.16330613,1.12162214,
    0.80839717,1.02609707,2.43566385,1.45087098};
  static const std::vector<float> fakeph_jet_sfs_run3 = {
    0.70322141,1.11836606,0.96501603,1.24315472,
    1.12892899,1.14611670,0.99618923,1.06107820,
    0.90955194,0.65634468,0.82285115,0.91302723};
  static const std::vector<float> fakeph_pu_sfs_run3 = {
    1.33178585,1.19331543,1.26473489,1.28245803,
    1.46886337,1.43250058,1.40016276,1.38490387,
    2.31379676,4.12526371,1.90097515,1.49419810};

  static const std::vector<float> fakeph_idres_sfs_run2 = {
    0.86703542,1.13279172,1.03335977,
    0.99237755,1.18276864,1.01707066,
    1.18351245,1.28974280,1.14538561,
    1.37581143,1.39472730,0.97234205};
  static const std::vector<float> fakeph_idres_sfs_run3 = {
    1.08025324,1.16182780,1.42111267,
    1.17597234,1.24471569,1.33403825,
    1.09437428,1.51553394,2.26460969,
    1.14807467,1.31767408,2.92946390};

  if (photon_pflavor == 1) return 1.0;
  const std::vector<float> *fakeph_sfs;
  const std::vector<float> *fakeph_idres_sfs;
  int run = 2;
  if (year > 2020.0) run = 3;
  if (run == 2 && this_photon_isjet) {
    fakeph_sfs = &fakeph_jet_sfs_run2;
    fakeph_idres_sfs = &fakeph_idres_sfs_run2;
  }
  else if (run == 2 && !this_photon_isjet) {
    fakeph_sfs = &fakeph_pu_sfs_run2;
    fakeph_idres_sfs = &fakeph_idres_sfs_run2;
  }
  else if (run == 3 && this_photon_isjet) {
    fakeph_sfs = &fakeph_jet_sfs_run3;
    fakeph_idres_sfs = &fakeph_idres_sfs_run3;
  }
  else if (run == 3 && !this_photon_isjet) {
    fakeph_sfs = &fakeph_pu_sfs_run3;
    fakeph_idres_sfs = &fakeph_idres_sfs_run3;
  }
  else return 1.0;

  float sf_pteta = 1.0;
  float sf_idres = 1.0;
  for (unsigned ipt = 0; ipt < (fakeph_ptbins.size()-1); ipt++) {
    for (unsigned ieta = 0; ieta < (fakeph_absetabins.size()-1); ieta++) {
      if (ph_pt >= fakeph_ptbins[ipt] && ph_pt < fakeph_ptbins[ipt+1]
          && ph_abseta >= fakeph_absetabins[ieta]
          && ph_abseta < fakeph_absetabins[ieta+1]) {
        sf_pteta = fakeph_sfs->at(ipt*(fakeph_absetabins.size()-1)+ieta);
        break;
      }
    }
  }
  for (unsigned iid = 0; iid < (fakeph_idbins.size()-1); iid++) {
    for (unsigned ires = 0; ires < (fakeph_resbins.size()-1); ires++) {
      if (ph_idmva >= fakeph_idbins[iid] && ph_idmva < fakeph_idbins[iid+1]
          && ph_res >= fakeph_resbins[ires]
          && ph_res < fakeph_resbins[ires+1]) {
        sf_idres = fakeph_idres_sfs->at(iid*(fakeph_resbins.size()-1)+ires);
        break;
      }
    }
  }
  return sf_pteta*sf_idres;
}

float get_w_llph_pt(float year, int type, bool this_photon_isjet,
                    float llphoton_pt) {
  static const std::vector<float> llph_ptbins = {
      0.0,5.0,10.0,15.0,20.0,30.0,40.0,60.0,80.0,100.0,500.0};
  static const std::vector<float> llph_zg_sfs_run2 = {
      0.80095164,0.97383617,1.07866292,1.10863927,1.05103401,1.01643762,
      1.00813146,1.02230660,1.15133346,0.99455812};
  static const std::vector<float> llph_jt_sfs_run2 = {
      0.43717314,0.94298997,0.98101131,0.89560594,1.04299537,1.05185543,
      1.04641453,1.21183438,2.03966944,2.19504543};
  static const std::vector<float> llph_pu_sfs_run2 = {
      1.12644807,0.88215348,0.96358980,0.93909889,1.02689214,1.00167445,
      1.11180094,1.06075877,0.50111905,0.61013672};
  static const std::vector<float> llph_zg_sfs_run3 = {
      0.79753245,0.99954447,1.08457808,1.12089606,1.11421871,0.95435079,
      0.96804314,1.02642743,0.99670823,1.01414788};
  static const std::vector<float> llph_jt_sfs_run3 = {
      0.38981804,0.81881434,0.82474122,1.14330376,1.11759692,1.07523328,
      1.24399480,1.12016585,1.19354538,1.11502050};
  static const std::vector<float> llph_pu_sfs_run3 = {
      1.21388660,0.96262377,0.95149128,0.90487912,0.95656613,1.08625187,
      1.01390092,1.04426156,1.32007871,1.27830503};
  const std::vector<float> *llph_sfs;
  int run = 2;
  if (year > 2020.0) run = 3;
  bool is_dyjet = (type >= 6000 && type < 7000 && this_photon_isjet);
  bool is_dypu = (type >= 6000 && type < 7000 && !this_photon_isjet);
  bool is_zg = (type >= 17000 && type < 18000);
  if (run == 2 && is_dyjet) llph_sfs = &llph_jt_sfs_run2;
  else if (run == 2 && is_dypu) llph_sfs = &llph_pu_sfs_run2;
  else if (run == 2 && is_zg) llph_sfs = &llph_zg_sfs_run2;
  else if (run == 3 && is_dyjet) llph_sfs = &llph_jt_sfs_run3;
  else if (run == 3 && is_dypu) llph_sfs = &llph_pu_sfs_run3;
  else if (run == 3 && is_zg) llph_sfs = &llph_zg_sfs_run3;
  else return 1.0;
  for (unsigned ipt = 0; ipt < (llph_ptbins.size()-1); ipt++) {
    if (llphoton_pt >= llph_ptbins[ipt] && llphoton_pt < llph_ptbins[ipt+1]) {
      return llph_sfs->at(ipt);
      break;
    }
  }
  return 1.0;
}

float get_llpjj_pt(RVec<float> llphoton_pt, RVec<float> llphoton_phi, 
                   RVec<float> jet_pt, RVec<float> jet_eta, 
                   RVec<float> jet_phi, RVec<bool> jet_isgood) {
  float pt_x = llphoton_pt[0]*cos(llphoton_phi[0]);
  float pt_y = llphoton_pt[0]*sin(llphoton_phi[0]);
  float lead_jet_pt = -999;
  float subl_jet_pt = -999;
  float lead_jet_phi = 0.0;
  float subl_jet_phi = 0.0;
  for (unsigned ijet = 0; ijet<jet_pt.size(); ijet++) {
    if (isgood_hornveto(jet_isgood.at(ijet), jet_pt.at(ijet), 
                        jet_eta.at(ijet))) {
      if (jet_pt[ijet] > lead_jet_pt) {
        subl_jet_pt = lead_jet_pt;
        subl_jet_phi = lead_jet_phi;
        lead_jet_pt = jet_pt[ijet];
        lead_jet_phi = jet_phi[ijet];
      }
      else if (jet_pt[ijet] > subl_jet_pt) {
        subl_jet_pt = jet_pt[ijet];
        subl_jet_phi = jet_phi[ijet];
      }
    }
  }
  pt_x += lead_jet_pt*cos(lead_jet_phi);
  pt_y += lead_jet_pt*sin(lead_jet_phi);
  pt_x += subl_jet_pt*cos(subl_jet_phi);
  pt_y += subl_jet_pt*sin(subl_jet_phi);
  return hypot(pt_x, pt_y);
}

//Gets llphotonjj pt weight
//TODO change name to balance
float get_w_llpjj_pt(float year, int type, RVec<float> llphoton_dijet_balance, 
                     int njet) {
  if (njet < 2) return 1.0;
  float llpjj_pt = llphoton_dijet_balance[0];
  static const std::vector<float> llpjj_ptbins = {
      0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.4,0.5,10.0};
  static const std::vector<float> llpjj_zg_sfs_run2 = {
      0.70346858,0.72718689,0.93645495,0.88680704,0.99462464,1.31871172,
      1.05642511,1.22181320,1.89596285,1.32755365,2.43267476};
  static const std::vector<float> llpjj_dy_sfs_run2 = {
      0.78097254,1.09496531,0.79448848,1.02530319,0.99289026,0.78107072,
      1.34585629,1.37819036,0.89429356,1.42922194,0.54031790};
  static const std::vector<float> llpjj_zg_sfs_run3 = {
      0.70576360,0.74688393,0.85359557,0.96456520,0.91700918,1.09012621,
      1.48679351,1.77907810,1.57360867,1.64248637,1.80862298};
  static const std::vector<float> llpjj_dy_sfs_run3 = {
      0.51050819,0.96290584,0.90547237,0.80490161,1.32526269,1.02370260,
      0.94556750,0.73765466,1.19208280,1.31317049,1.53405749};
  int run = 2;
  if (year > 2020) run = 3;
  const std::vector<float> *llpjj_sfs;
  bool is_dy = (type >= 6000 && type < 7000);
  bool is_zg = (type >= 17000 && type < 18000);
  if (run == 2 && is_dy) llpjj_sfs = &llpjj_dy_sfs_run2;
  else if (run == 2 && is_zg) llpjj_sfs = &llpjj_zg_sfs_run2;
  else if (run == 3 && is_dy) llpjj_sfs = &llpjj_dy_sfs_run3;
  else if (run == 3 && is_zg) llpjj_sfs = &llpjj_zg_sfs_run3;
  else return 1.0;
  for (unsigned ipt = 0; ipt < (llpjj_ptbins.size()-1); ipt++) {
    if (llpjj_pt >= llpjj_ptbins[ipt] && llpjj_pt < llpjj_ptbins[ipt+1]) {
      return llpjj_sfs->at(ipt);
      break;
    }
  }
  return 1.0;
}

float get_l1_rapidity(RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> mu_pt, RVec<float> mu_eta, RVec<int> ll_lepid, 
    RVec<int> ll_i1, RVec<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i1[0]] : 
           el_eta[ll_i2[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i1[0]] : 
         mu_eta[ll_i2[0]];
}

float get_l2_rapidity(RVec<float> el_pt, RVec<float> el_eta, 
    RVec<float> mu_pt, RVec<float> mu_eta, RVec<int> ll_lepid, 
    RVec<int> ll_i1, RVec<int> ll_i2) {
  if (ll_lepid[0]==11) {
    return (el_pt[ll_i1[0]] > el_pt[ll_i2[0]]) ? el_eta[ll_i2[0]] : 
           el_eta[ll_i1[0]];
  }
  return (mu_pt[ll_i1[0]] > mu_pt[ll_i2[0]]) ? mu_eta[ll_i2[0]] : 
         mu_eta[ll_i1[0]];
}

