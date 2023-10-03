/**
 * @brief Metropolis-Hastings-based MCMC sampling algorithm and MCMC uniform integrator
 */
#include <functional>
#include <string>
#include <vector>

namespace spu_mcmc {

/**
 * @brief Struct containing information needed for variables sampled by MCMC
 *        algorithm including variable name, minimum value, maximum value, and
 *        characteristic size
 */
struct SampledVariable {
  std::string name;
  float min;
  float max;
  float step_size;
};

/**
 * @brief Struct containing information used for auxiliary variables computed
 *        from sampled variables including variable name, function that 
 *        computes it from sampled variables, minimum for plots and maximum
 *        for plots
 */
struct AuxiliaryVariable {
  std::string name;
  std::function<float (const std::vector<float>&)> function;
  float min;
  float max;
};

/**
 * @brief Integrates a function using *uniform* Monte Carlo integration, may 
 *        have convergence issues for sharply peaked PDFs
 *
 * @param n_samples      number of MCMC samples to generate
 * @param variables      random variables to sample
 * @param pdf            function providing density function (up to a constant)
 */
float mc_integrate(unsigned int n_samples, std::vector<SampledVariable> variables, 
                   std::function<float (const std::vector<float>&)> pdf);

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
                 bool write_hists = true, bool write_ntuple = false, 
                 std::vector<AuxiliaryVariable> aux_variables = {},
                 unsigned int n_burnin = 10000);

}
