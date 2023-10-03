/**
 * @brief Functions for using the matrix element method (MEM) or matrix element likelihood analysis (MELA) for the H->Zgamma analysis
 */
#ifndef H_HZG_MEM
#define H_HZG_MEM

#include "TLorentzVector.h"

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
float qqzg_likelihood(float const & cosTheta, float const & costheta, float const & phi);

/**
 * @brief Returns the likelihood (up to a constant) for H->Zgamma given decay angles
 *        specifying kinematics. 
 *
 *
 * @param costheta          angle between Z/photon and lepton in H frame
 */
float hzg_likelihood(float const & costheta);

/**
 * @brief Returns signal-background discriminant based on matrix element likelihood analysis (MELA)
 *
 * @param cosTheta          angle between Z/photon and initial parton in H frame
 * @param costheta          angle between Z/photon and lepton in H frame
 * @param phi               angle between photon+lepton plane and photon+parton 
 *                            plane in H frame
 */
float mela_discriminant(float const & cosTheta, float const & costheta, float const & phi);

/**
 * @brief Returns 4-momentum of initial parton going along +z axis
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
TLorentzVector get_q1(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);

/**
 * @brief Returns cosine of the angle between the Z/photon and initial parton in 
 *        the H frame
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
float get_cosTheta_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);

/**
 * @brief Returns cosine of the angle between the Z/photon and lepton in 
 *        the H frame
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
float get_costheta_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);

/**
 * @brief Returns angle between photon+lepton plane and photon+parton plane
 *        the H frame
 *        From github.com/jaebak/llg_angles
 *
 * @param lep_minus         four-momentum of negative lepton
 * @param lep_plus          four-momentum of positive lepton
 * @param gamma             four-momentum of photon
 */
float get_phi_alt2(TLorentzVector const & lep_minus, TLorentzVector const & lep_plus, TLorentzVector const & gamma);

#endif
