/**
 * @brief Small script to optimize binning cuts given ROOT histograms (TH1D) for signal and background
 */
#include <cmath>
#include <iostream>
#include <vector>
#include "TH1.h"

/**
 * @brief Function to find optimized binning cuts given distribution as histograms
 *
 * @param signal     pointer to ROOT histogram of signal
 * @param background pointer to ROOT histogram of background
 * @param nbins      number of bins (ncuts + 1)
 */
void binning_optimizer(TH1D* signal, TH1D* background, int nbins) {
  if (nbins < 2) {
    float bin_s = signal->Integral();
    float bin_b = background->Integral();
    std::cout << "Optimal cuts(sig,bak): no cuts(" << bin_s << "," << bin_b << ")" << std::endl;
    std::cout << "Estimated significance: " << bin_s/sqrt(bin_b) << std::endl;
    return;
  }
  int hist_nbins = signal->GetNbinsX();
  float max_significance = 0;
  std::vector<int> max_sig_cuts;
  //initialize vector
  //convention: cuts are specified as lowest histogram bin in upper binning
  std::vector<int> cuts(nbins-1,0);
  for (int i = 0; i<(nbins-1); i++) {
    cuts[i] = i+2;
  }
  //loop over all possible combinations of cuts to find optimum 
  bool finished = false;
  while (!finished) {
    //calculate significance for current cuts
    std::vector<int> extended_cuts;
    extended_cuts.push_back(1);
    extended_cuts.insert(extended_cuts.end(),cuts.begin(),cuts.end());
    float significance = 0;
    for (int i = 0; i<(extended_cuts.size()-1); i++) {
      float bin_s = signal->Integral(extended_cuts[i],extended_cuts[i+1]-1);
      float bin_b = background->Integral(extended_cuts[i],extended_cuts[i+1]-1);
      if (bin_s <= 0) bin_s = 0;
      if (bin_b <= 0) bin_b = 0.01;
      significance += bin_s*bin_s/bin_b;
    }
    significance = sqrt(significance);
    if (significance > max_significance) {
      max_significance = significance;
      max_sig_cuts = cuts;
    }
    //iterate to next set of cuts
    int place = 0;
    bool iterated = false;
    while (!iterated) {
      if (cuts[nbins-2-place] != hist_nbins-1-place) {
        cuts[nbins-2-place]++;
        int place_offset = 1;
        for (int lower_place = place-1; lower_place >= 0; lower_place--) {
          cuts[nbins-2-lower_place] = cuts[nbins-2-place]+place_offset;
          place_offset++;
        }
        iterated = true;
      }
      else {
        if (place == nbins-2) {
          iterated = true;
          finished = true;
        }
        else {
          place += 1;
        }
      }
    }//end loop used to increment cuts
  }//end loop over all possible cuts
  //print optimal cuts
  std::cout << "Optimal cuts(sig,bak): ";
  std::cout << signal->GetBinLowEdge(1);
  std::cout << "(" << signal->Integral(1,max_sig_cuts[0]-1);
  std::cout << "," << background->Integral(1,max_sig_cuts[0]-1) << ")";
  bool first = false;
  for (int i = 0; i<(nbins-1); i++) {
    if (!first) std::cout << ", ";
    std::cout << signal->GetBinLowEdge(max_sig_cuts[i]);
    int upper_lim = hist_nbins;
    if (i != nbins-2) {
      upper_lim = max_sig_cuts[i+1]-1;
    }
    std::cout << "(" << signal->Integral(max_sig_cuts[i],upper_lim);
    std::cout << "," << background->Integral(max_sig_cuts[i],upper_lim) << ")";
  }
  std::cout << std::endl;
  std::cout << "Estimated significance: " << max_significance << std::endl;
}
