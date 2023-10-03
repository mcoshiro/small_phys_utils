/**
 * @brief Functions for using the matrix element method (MEM) or matrix element likelihood analysis (MELA) for the H->Zgamma analysis
 */
#include "hzg_mem.hpp"

#include <cmath>

#include "TLorentzVector.h"
#include "TRotation.h"
#include "TVector3.h"

//constants used by the qqzg_likelihood function
const float mz = 91.2;
const float mh = 125.0;
const float sinsqtw = 0.23146;
const float cossqtw = 1.0-sinsqtw;
const float ulcz = 0.5*cossqtw-1.0/6.0*sinsqtw;
const float dlcz = -0.5*cossqtw-1.0/6.0*sinsqtw;
const float urcz = -2.0/3.0*sinsqtw;
const float drcz = 1.0/3.0*sinsqtw;
const float elcz = -0.5*cossqtw+1.0/2.0*sinsqtw;
const float ercz = 1.0*sinsqtw;
const float go_coef_u = (ercz+elcz)*(urcz+ulcz);
const float gt_coef_u = (ercz-elcz)*(urcz-ulcz);
const float go_coef_d = (ercz+elcz)*(drcz+dlcz);
const float gt_coef_d = (ercz-elcz)*(drcz-dlcz);
const float mf_p_ss = pow(mz,4)+pow(mh,4); 
const float ms_x_s = mz*mz*mh*mh;
const float m_x_sqs = mz*mh;
const float ms_p_s = mz*mz+mh*mh;

/**
 * @brief Returns the likelihood (up to a constant) for qq->Zgamma given decay angles
 *        specifying kinematics. 
 *
 *        Taken from 1112.1405 since I can't be bothered to take traces of 
 *        eight gamma matrices. Assumes mass already fixed to Z/H
 *
 * @param cosTheta          angle between Z/photon and initial parton in H frame
 * @param costheta          angle between Z/photon and lepton in H frame
 * @param phi               angle between photon+lepton plane and photon+parton 
 *                            plane in H frame
 */
float qqzg_likelihood(float const & cosTheta, float const & costheta, float const & phi) {
  //predefine some trig functions
  float costwotheta = 2.0*costheta*costheta-1.0;
  float sintwotheta = 2.0*costheta*sqrt(1.0-costheta*costheta);
  float sinsqtheta  = 1.0-costheta*costheta;
  float sintheta    = sqrt(sinsqtheta);
  float cscsqTheta  = 1.0/(1.0-cosTheta*cosTheta);
  float cosphi      = cos(phi);
  float cotTheta    = cosTheta/sqrt(1.0-costheta*costheta);

  //evaluate gs
  float gone = mf_p_ss*(3.0+costwotheta)*(4.0*cscsqTheta-2.0)
               +8.0*ms_x_s*sinsqtheta*(2.0+cos(2.0*phi))
               +8.0*m_x_sqs*ms_p_s*cotTheta*sintwotheta*cosphi;
  float gtwo = 16.0*sqrt(cscsqTheta)*(mf_p_ss*costheta*cotTheta
               +m_x_sqs*ms_p_s*sintheta*cosphi);

  //probably should integrate against PDFs for relative u/d fraction
  //for now, just assume 50/50? or check truth in nanos/picos?
  return (go_coef_u+go_coef_d)*gone+(gt_coef_u+gt_coef_d)*gtwo;
}

/**
 * @brief Returns the likelihood (up to a constant) for H->Zgamma given decay angles
 *        specifying kinematics. 
 *
 *
 * @param costheta          angle between Z/photon and lepton in H frame
 */
float hzg_likelihood(float const & costheta) {
  return (1.0+costheta*costheta);
}

/**
 * @brief Returns signal-background discriminant based on matrix element likelihood analysis (MELA)
 *
 * @param cosTheta          angle between Z/photon and initial parton in H frame
 * @param costheta          angle between Z/photon and lepton in H frame
 * @param phi               angle between photon+lepton plane and photon+parton 
 *                            plane in H frame
 */
float mela_discriminant(float const & cosTheta, float const & costheta, float const & phi) {
  float l_sig = hzg_likelihood(costheta);
  float l_bkg = qqzg_likelihood(cosTheta, costheta, phi);
  return 1.0/(1.0+l_bkg/l_sig/5000000); //arbitrary normalization
}

/**
 * @brief Returns signal-background discriminant based on matrix element likelihood analysis (MELA)
 *        This version assumes initial partons can be distinguished
 *
 * @param cosTheta          angle between Z/photon and initial parton in H frame
 * @param costheta          angle between Z/photon and lepton in H frame
 * @param phi               angle between photon+lepton plane and photon+parton 
 *                            plane in H frame
 */
float mela_discriminant_ideal(float const & cosTheta, float const & costheta, float const & phi) {
  float l_sig = hzg_likelihood(costheta);
  float l_bkg = qqzg_likelihood(cosTheta, costheta, phi);
  return 1.0/(1.0+l_bkg/l_sig); //currently, normalization is floating
}

/**
 * @brief Returns 4-momentum of initial parton going along +z axis
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
TLorentzVector get_q1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
  TLorentzVector llg_p4 = lep_minus + lep_plus + gamma;
  TVector3 htran = llg_p4.BoostVector();
  htran.SetZ(0);
  llg_p4.Boost(-1*htran);
  float pz, E;
  pz = llg_p4.Pz() + llg_p4.E();
  E  = llg_p4.E()  + llg_p4.Pz();
  TLorentzVector k1;
  k1.SetPxPyPzE(0,0,pz/2,E/2);
  k1.Boost(htran);
  return k1;
}

/**
 * @brief Returns cosine of the angle between the Z/photon and initial parton in 
 *        the H frame
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
float get_cosTheta_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
  TLorentzVector gamma_p4 = gamma;
  TLorentzVector l1_p4 = lep_minus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;
  TLorentzVector llg_p4 = ll_p4 + gamma_p4;
  TLorentzVector q1_p4 = get_q1(lep_minus, lep_plus, gamma);
  // Boost to llg frame
  TVector3 llgBoost = llg_p4.BoostVector();
  ll_p4.Boost(-1*llgBoost);
  q1_p4.Boost(-1*llgBoost);

  float cosTheta = cos(ll_p4.Angle(q1_p4.Vect()));
  return cosTheta;
}

/**
 * @brief Returns cosine of the angle between the Z/photon and lepton in 
 *        the H frame
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
float get_costheta_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
  TLorentzVector gamma_p4 = gamma;
  TLorentzVector l2_p4 = lep_plus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;
  TLorentzVector llg_p4 = ll_p4 + gamma_p4;
  // Boost to llg frame
  TVector3 llgBoost = llg_p4.BoostVector();
  ll_p4.Boost(-1*llgBoost);
  l2_p4.Boost(-1*llgBoost);
  // Boost to ll frame
  TVector3 llBoost = ll_p4.BoostVector();
  l2_p4.Boost(-1*llBoost);
  float costheta = cos(ll_p4.Angle(l2_p4.Vect()));
  return costheta;
}

/**
 * @brief Returns angle between photon+lepton plane and photon+parton plane
 *        the H frame
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
float get_phi_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma) {
  TLorentzVector gamma_p4 = gamma;
  TLorentzVector l1_p4 = lep_minus;
  TLorentzVector ll_p4 = lep_minus + lep_plus;
  TLorentzVector llg_p4 = ll_p4 + gamma_p4;
  TLorentzVector q1_p4 = get_q1(lep_minus, lep_plus, gamma);
  // Boost to llg frame
  TVector3 llgBoost = llg_p4.BoostVector();
  ll_p4.Boost(-1*llgBoost);
  q1_p4.Boost(-1*llgBoost);
  l1_p4.Boost(-1*llgBoost);

  TVector3 zAxis = ll_p4.Vect().Unit();
  TVector3 yAxis = q1_p4.Vect().Cross(zAxis.Unit()).Unit();
  TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();
  TRotation rotation;
  rotation = rotation.RotateAxes(xAxis,yAxis,zAxis).Inverse();
  l1_p4.Transform(rotation);
  float phi = l1_p4.Phi();
  return phi;
}
