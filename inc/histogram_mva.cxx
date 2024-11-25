#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"
#include "TTree.h"

/**
 * @brief Returns position of nth occurence of target in string
 *
 * @param search_string    string to search
 * @param target           target to find
 * @param n                occurence number
 */
size_t find_nth(std::string search_string, std::string target, int n);

/**
 * @brief Class used to evaluate histogram-based MVA discriminant
 */
class HistMVAEvaluator  {
public:

  /**
   * @brief Standard constructor
   *
   * @param filename   file name from which to load histograms
   * @param var_names  names of variables to use in evaluation
   */
  HistMVAEvaluator(std::string filename, std::vector<std::string> var_names);

  /**
   * @brief Evaluates MVA-based discriminant
   *
   * @param inputs   input variable values
   */
  float evaluate(std::vector<float> const & inputs) const;

private:
  std::vector<TH1F*> hist_1d_sig;
  std::vector<TH1F*> hist_1d_bak;
  std::vector<TH2F*> hist_2d_sig;
  std::vector<TH2F*> hist_2d_bak;
  std::vector<unsigned> oned_var_idx_sig;
  std::vector<unsigned> oned_var_idx_bak;
  std::vector<unsigned> twod_var1_idx_sig;
  std::vector<unsigned> twod_var1_idx_bak;
  std::vector<unsigned> twod_var2_idx_sig;
  std::vector<unsigned> twod_var2_idx_bak;
};

/**
 * @brief Struct containing information on input features including the name
 *        of the branch in the TTree, the minimum to use on histograms, the
 *        maximum, and the number of bins
 */
struct Feature {
  std::string name;
  float min;
  float max;
  unsigned nbins;
}

/**
 * @brief Trains an MVA discriminant based on approximating likelihoods with 
 *        histograms. Note that unlike TMVA, no shuffling is provided
 *
 * @param signal_tree       tree of signal events to train on
 * @param background_tree   tree of background events to train on
 * @param features          input variables organized as a vector of vectors, 
 *                          where variables in the same subvector are treated 
 *                          as correlated. Currently only float branches are
 *                          supported
 * @param output_filename   name of ROOT file that will store MVA parameters
 * @param weight_name       name of weight branch to be used
 * @param ntrain_signal     number of events to use for signal training. 0 uses half
 * @param ntrain_background number of events to use for background training. 0 uses half
 */
void train_histogram_mva(TTree* signal_tree, TTree* background_tree, 
    std::vector<std::vector<Feature>> features, std::string output_filename,
    std::string weight_name = "", unsigned ntrain_signal, 
    unsigned ntrain_background);
