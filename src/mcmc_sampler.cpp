/**
 * @brief Metropolis-Hastings-based MCMC sampling algorithm and MCMC uniform integrator
 */
#include "mcmc_sampler.hpp"

#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TTree.h"

namespace spu_mcmc {

/**
 * @brief Integrates a function using *uniform* Monte Carlo integration, may 
 *        have convergence issues for sharply peaked PDFs
 *
 * @param n_samples      number of MCMC samples to generate
 * @param variables      random variables to sample
 * @param pdf            function providing density function (up to a constant)
 */
float mc_integrate(unsigned int n_samples, std::vector<SampledVariable> variables, 
                   std::function<float (const std::vector<float>&)> pdf) {
  TRandom3 rng;
  std::vector<float> random_vars;
  float total_measure = 1.0;
  float sum_pdf = 0.0;
  for (unsigned irv = 0; irv < variables.size(); irv++) {
    random_vars.push_back(0.0);
    total_measure *= (variables[irv].max-variables[irv].min);
  }
  for (unsigned ievt = 0; ievt < n_samples; ievt++) {
    for (unsigned irv = 0; irv < variables.size(); irv++) {
      random_vars[irv] = rng.Uniform(variables[irv].min,variables[irv].max);
    }
    sum_pdf += pdf(random_vars);
  }
  return sum_pdf/static_cast<float>(n_samples)*total_measure;
}

/**
 * @brief Samples a given probability distribution using Metropolis-Hastings 
 *        algorithm 
 *
 * @param n_samples      number of MCMC samples to generate
 * @param variables      random variables to sample
 * @param pdf            function providing density function (up to a constant)
 * @param aux_variables  auxiliary variables to calculate from random variables
 * @param write_hists    whether to generate histograms
 * @param write_ntuple   whether to generate ntuple
 * @param n_burnin       number of burn-in samples
 */
void sample_mcmc(unsigned int n_samples, std::vector<SampledVariable> variables, 
                 std::function<float (const std::vector<float>&)> pdf, 
                 bool write_hists, bool write_ntuple, 
                 std::vector<AuxiliaryVariable> aux_variables,
                 unsigned int n_burnin) {
  //setup
  TTree out_tree("tree","tree");
  TRandom3 rng;
  unsigned int n_current = 0;
  std::vector<float> prev_random_vars;
  std::vector<float> next_random_vars;
  std::vector<float> aux_vars;
  std::vector<TH1D> hist;
  std::vector<TH1D> aux_hist;
  for (unsigned irv = 0; irv < variables.size(); irv++) {
    prev_random_vars.push_back((variables[irv].max+variables[irv].min)/2.0);
    next_random_vars.push_back(0);
    hist.push_back(TH1D(("hist_"+variables[irv].name).c_str(),
                   variables[irv].name.c_str(),100,
                   variables[irv].min,variables[irv].max));
    out_tree.Branch(variables[irv].name.c_str(), &prev_random_vars[irv], 
                    (variables[irv].name+"/F").c_str());
  }
  for (unsigned iaux = 0; iaux < aux_variables.size(); iaux++) {
    aux_vars.push_back(0);
    aux_hist.push_back(TH1D(("hist_"+aux_variables[iaux].name).c_str(),
                       aux_variables[iaux].name.c_str(),100,
                       aux_variables[iaux].min,aux_variables[iaux].max));
    out_tree.Branch(aux_variables[iaux].name.c_str(), &aux_vars[iaux], 
                    (aux_variables[iaux].name+"/F").c_str());
  }
  float prev_pdf = pdf(prev_random_vars);
  float next_pdf = 0.0;

  //main loop - sample from distribution
  while (n_current < (n_samples+n_burnin)) {

    //generate next random variables
    for (unsigned irv = 0; irv < variables.size(); irv++) {
      next_random_vars[irv] = rng.Gaus(prev_random_vars[irv],variables[irv].step_size);
      while (next_random_vars[irv] > variables[irv].max) 
        next_random_vars[irv] -= (variables[irv].max-variables[irv].min);
      while (next_random_vars[irv] < variables[irv].min) 
        next_random_vars[irv] += (variables[irv].max-variables[irv].min);
    }

    //check if event accepted
    next_pdf =pdf(next_random_vars);
    float rng_acc = rng.Uniform();
    if (rng_acc <= next_pdf/prev_pdf) {
      prev_pdf = next_pdf;
      for (unsigned irv = 0; irv < variables.size(); irv++) {
        prev_random_vars[irv] = next_random_vars[irv];
        if (n_current>n_burnin)
          hist[irv].Fill(prev_random_vars[irv]);
      }
      if (n_current>n_burnin) {
        for (unsigned iaux = 0; iaux < aux_variables.size(); iaux++) {
          aux_vars[iaux] = aux_variables[iaux].function(prev_random_vars);
          aux_hist[iaux].Fill(aux_vars[iaux]);
        }
        if (write_ntuple)
          out_tree.Fill();
      }
      if ((n_current % 20000) == 0) 
        std::cout << "Processing event " << n_current << std::endl;
      n_current++;
    }
  }
  
  //write output - histograms
  if (write_hists) {
    std::vector<TCanvas*> canvas;
    for (unsigned irv = 0; irv < variables.size(); irv++) {
      canvas.push_back(new TCanvas(("canvas_"+variables[irv].name).c_str(),"c",600,400));
      canvas.back()->cd();
      hist[irv].Draw();
      canvas.back()->SaveAs(("hist_"+variables[irv].name+".pdf").c_str());
    }
    for (unsigned iaux = 0; iaux < aux_variables.size(); iaux++) {
      canvas.push_back(new TCanvas(("canvas_"+aux_variables[iaux].name).c_str(),"c",600,400));
      canvas.back()->cd();
      aux_hist[iaux].Draw();
      canvas.back()->SaveAs(("hist_"+aux_variables[iaux].name+".pdf").c_str());
    }
    for (unsigned irv = 0; irv < canvas.size(); irv++) {
      delete canvas[irv];
    }
  }

  //write output - ntuples
  if (write_ntuple) {
    TFile out_file("sampled_tree.root","RECREATE");
    out_tree.Write();
    out_file.Close();
  }
}

}
