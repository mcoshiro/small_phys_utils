/**
 * @brief Collection of functions for analyzing distributions
 */
#include "distribution_analyzer.hpp"

#include <cmath>
#include <iostream>
#include <vector>

#include "TH1.h"

/**
 * @brief Returns area under ROC curve for given distributions
 *
 * @param signal            pointer to ROOT histogram of signal
 * @param background        pointer to ROOT histogram of background
 */
void get_roc_auc(TH1D* signal, TH1D* background) {

    //integrate where measure is %bkg accepted
    int hist_nbins = (signal->GetNbinsX())+2;
    float total_bkg = background->Integral();
    float total_sig = signal->Integral();
    float integrated_area = 0;
    float casi = 0;
    float prev_bkg = 1.0;
    float prev_sig = 1.0;
    float bkg_accepted = 0.0;
    float sig_accepted = 0.0;
    for (int ibin = 1; ibin < hist_nbins; ibin++) {
      bkg_accepted = (background->Integral(ibin,hist_nbins))/total_bkg;
      sig_accepted = (signal->Integral(ibin,hist_nbins))/total_sig;
      integrated_area += (sig_accepted+prev_sig)/2.0*(prev_bkg-bkg_accepted);
      if (prev_bkg>bkg_accepted)
        casi += (sig_accepted-prev_sig)*(sig_accepted-prev_sig)/(prev_bkg-bkg_accepted);
      prev_bkg = bkg_accepted;
      prev_sig = sig_accepted;
    }
    std::cout << "CSI = Continuous significance improvement (approximation to significance improvement)." << std::endl;
    std::cout << "ROC AUC is: " << integrated_area << std::endl;
    std::cout << "ROC CSI is: " << sqrt(casi) << std::endl;
}

/**
 * @brief Function to find optimized binning cuts given distribution as histograms
 *
 * @param signal            pointer to ROOT histogram of signal
 * @param background        pointer to ROOT histogram of background
 * @param max_nbins         maximum number of bins to consider (ncuts + 1)
 * @param min_signal_yield  minimum signal yield allowed in a bin
 * @param scale             scale factor for signal and background to change effective lumi
 * @param nbins_throw       number of bins to remove from significance calculation
 */
std::vector<int> binning_optimizer(TH1D* signal, TH1D* background, 
    int max_nbins, float min_signal_yield, float scale, unsigned nbins_throw) {

  //for (int nbins = 1; nbins <= max_nbins; nbins++) {
  int nbins = max_nbins;
    std::cout << "With " << nbins << " bins: " << std::endl;
    if (nbins < 2) {
      float bin_s = signal->Integral()*scale;
      float bin_b = background->Integral()*scale;
      std::cout << "Optimal cuts(sig,bak): -1.1(" << bin_s << "," << bin_b << ")" << std::endl;
      std::cout << "Estimated significance: " << bin_s/sqrt(bin_b) << std::endl;
      return {1,signal->GetNbinsX()};
    }
    int hist_nbins = signal->GetNbinsX();
    float max_significance = 0;
    std::vector<int> max_sig_cuts;
    //initialize vector
    //convention: cuts are specified as lowest histogram bin in upper binning
    std::vector<int> cuts(nbins-1,0);
    for (int i = 0; i<(nbins-1); i++) {
      cuts[i] = i+1;
    }
    //loop over all possible combinations of cuts to find optimum 
    bool finished = false;
    while (!finished) {
      //calculate significance for current cuts
      std::vector<int> extended_cuts;
      extended_cuts.push_back(0);
      extended_cuts.insert(extended_cuts.end(),cuts.begin(),cuts.end());
      extended_cuts.push_back(hist_nbins+2);
      float significance = 0;
      bool can_throw = false;
      if (extended_cuts.size() > nbins_throw+1) can_throw = true;
      for (unsigned i = 0; i<(extended_cuts.size()-1); i++) {
        float bin_s = signal->Integral(extended_cuts[i],extended_cuts[i+1]-1)*scale;
        float bin_b = background->Integral(extended_cuts[i],extended_cuts[i+1]-1)*scale;
        if (bin_s <= min_signal_yield) bin_s = 0;
        if (bin_b <= 0) bin_b = 0.01;
        if (!(can_throw && i<nbins_throw))
          significance += bin_s*bin_s/bin_b;
      }
      significance = sqrt(significance);
      if (significance > max_significance) {
        max_significance = significance;
        max_sig_cuts = cuts;
      }
      //debug print cuts
      //std::cout << "Cuts: "; for (int cut : cuts) std::cout << ", " << cut; std::cout << std::endl;
      //iterate to next set of cuts
      int place = 0;
      bool iterated = false;
      while (!iterated) {
        if (cuts[nbins-2-place] != hist_nbins+1-place) {
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
    std::cout << "(" << signal->Integral(1,max_sig_cuts[0]-1)*scale;
    std::cout << "," << background->Integral(1,max_sig_cuts[0]-1)*scale << ")";
    bool first = false;
    for (int i = 0; i<(nbins-1); i++) {
      if (!first) std::cout << ", ";
      std::cout << signal->GetBinLowEdge(max_sig_cuts[i]);
      int upper_lim = hist_nbins;
      if (i != nbins-2) {
        upper_lim = max_sig_cuts[i+1]-1;
      }
      std::cout << "(" << signal->Integral(max_sig_cuts[i],upper_lim)*scale;
      std::cout << "," << background->Integral(max_sig_cuts[i],upper_lim)*scale << ")";
    }
    std::cout << std::endl;
    std::cout << "Estimated significance: " << max_significance << std::endl;
  //}
  return max_sig_cuts;
}

/**
 * @brief Function to find optimized cut for a histogram
 *
 * @param signal            pointer to ROOT histogram of signal
 * @param background        pointer to ROOT histogram of background
 * @param min_signal_yield  minimum signal yield allowed after cut
 * @param cut_is_upper      if the cut is a lower bound or an upper bound
 * @param scale             scale factor for signal and background to change effective lumi
 */
void cut_optimizer(TH1D* signal, TH1D* background,
    float min_signal_yield, bool cut_is_upper, float scale) {
  int hist_nbins = signal->GetNbinsX();
  float max_significance = 0;
  int max_sig_cut = 0;
  //convention: cuts are specified as lowest histogram bin in upper binning
  //loop over all possible combinations of cuts to find optimum 
  if (!cut_is_upper) {
    for (int cut = 0; cut < hist_nbins; cut++) {
      float cut_s = signal->Integral(cut,hist_nbins+1)*scale;
      float cut_b = background->Integral(cut,hist_nbins+1)*scale;
      if (cut_s < min_signal_yield) continue;
      float significance = cut_s/sqrt(cut_b);
      if (significance > max_significance) {
        max_significance = significance;
        max_sig_cut = cut;
      }
    }
  }
  if (cut_is_upper) {
    for (int cut = hist_nbins; cut > 0; cut--) {
      float cut_s = signal->Integral(0,cut)*scale;
      float cut_b = background->Integral(0,cut)*scale;
      if (cut_s < min_signal_yield) continue;
      float significance = cut_s/sqrt(cut_b);
      if (significance > max_significance) {
        max_significance = significance;
        max_sig_cut = cut;
      }
    }
  }
  //print optimal cuts
  std::cout << "Optimal cut: " << signal->GetBinLowEdge(max_sig_cut) << std::endl;
  std::cout << "Estimated significance: " << max_significance << std::endl;
}
