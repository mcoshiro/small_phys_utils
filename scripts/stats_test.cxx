#include <cmath>
#include <iostream>
#include <vector>

#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"

void stats_test() {
  TRandom3 rng;
  rng.SetSeed();
  int dummy = 0;
  int ntoys = 2000;
  float f_ntoys = static_cast<float>(ntoys);
  float sig = 4.0;
  float bkg = 10.0;
  float true_r = 0.0;
  float sb = true_r*sig+bkg;
  TCanvas c;

  ////calculate test stat as a function of R
  //std::vector<float> r_values = {-2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0};
  //std::vector<float> q_values;
  //for (float r : r_values) {
  //  q_values.push_back(-2.0*std::log(TMath::Poisson(nobs,r*sig+bkg)/TMath::Poisson(nobs,bkg)));
  //}
  //TGraph my_graph(r_values.size(),&r_values[0],&q_values[0]);
  //my_graph.Draw();
  //c.Print("qplot.pdf");
  
  //std::vector<float> r_values = {-1.0, -0.5, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
  std::vector<float> r_values;
  for (float r = -4; r < 8; r += 0.2) {
    if (r < 1.0e-3 && r > -1.0e-3) continue; //skip 0
    r_values.push_back(r);
  }

  float prev_clspb = 0;
  std::vector<float> ci_lower_bounds;
  std::vector<float> ci_upper_bounds;
  for (int iexp = 0; iexp < 100; iexp++) {
    //do a 'real' experiment
    bool found_lower = false;
    bool found_upper = false;
    float true_nobs = static_cast<float>(rng.Poisson(sb));
    std::cout << true_nobs << " ";

    prev_clspb = 0;
    for (float r : r_values) {
      //float best_fit_r = (true_nobs-bkg)/sig;
      //if (best_fit_r < 0) best_fit_r = 0;
      //float q_obs = -2.0*std::log(TMath::Poisson(true_nobs,r*sig+bkg)/TMath::Poisson(true_nobs,best_fit_r*sig+bkg));
      //float q_obs = -2.0*std::log(TMath::Poisson(true_nobs,r*sig+bkg)/TMath::Poisson(true_nobs,bkg));
      float q_obs = true_nobs;

      //throw toys
      //TH1D q_hist("q","q",100,-16,64);
      float nrej = 0;
      for (int i=0; i<f_ntoys; i++) {
        float nobs = static_cast<float>(rng.Poisson(r*sig+bkg));
        //best_fit_r = (nobs-bkg)/sig;
        //if (best_fit_r < 0) best_fit_r = 0;
        //float q = -2.0*std::log(TMath::Poisson(nobs,r*sig+bkg)/TMath::Poisson(nobs,best_fit_r*sig+bkg));
        //float q = -2.0*std::log(TMath::Poisson(nobs,r*sig+bkg)/TMath::Poisson(nobs,bkg));
        float q = nobs;
        if (q>q_obs) nrej += 1;
        //q_hist.Fill(q);
      }
      //c.cd();
      //q_hist.Draw();
      //c.Draw();
      //c.Print("qplot.pdf");
      //float clspb = q_hist.Integral(q_hist.FindBin(true_q),101)/q_hist.Integral();
      float clspb = nrej/f_ntoys;
      if (!found_lower && clspb > 0.05 && prev_clspb < 0.05) {
        ci_lower_bounds.push_back(r);
        found_lower = true;
      }
      if (!found_upper && clspb < 0.05 && prev_clspb > 0.05) {
        ci_upper_bounds.push_back(r);
        found_upper = true;
      }
      prev_clspb = clspb;
      //std::cout << "CL(s+b)(r = " << r << "): " << clspb << std::endl;
      //std::cin >> dummy;
    }
    //std::cout << "---------------------------------" << std::endl;
    if (found_upper && found_lower)
      std::cout << "CI(r) is (" << ci_lower_bounds.back() << "," << ci_upper_bounds.back() << ")" << std::endl;
    else if (!found_upper && found_lower)
      ci_lower_bounds.pop_back();
    else if (found_upper && !found_lower)
      ci_upper_bounds.pop_back();
  }

}
