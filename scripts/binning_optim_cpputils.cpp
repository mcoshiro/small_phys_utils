
#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

#include "TAxis.h"
#include "TH2D.h"

using std::cout;
using std::function;
using std::vector;

/**
 * Gets lowest positive yield across 1D histogram
 */
double get_lowest_positive_yield(TH1D* hist) {
  int nbins = hist->GetNbinsX();
  double lowest_positive_count = 9.9e9;
  for (int ix = 0; ix <= nbins+1; ix++) {
    double bin_content = hist->GetBinContent(ix);
    if (bin_content > 0 && bin_content < lowest_positive_count) {
      lowest_positive_count = bin_content;
    }
  }
  return lowest_positive_count;
}

/**
 * Gets lowest positive yield across 2D histogram
 */
double get_lowest_positive_yield(TH2D* hist) {
  int nbins_x = hist->GetNbinsX();
  int nbins_y = hist->GetNbinsY();
  double lowest_positive_count = 9.9e9;
  for (int ix = 0; ix <= nbins_x+1; ix++) {
    for (int iy = 0; iy <= nbins_y+1; iy++) {
      double bin_content = hist->GetBinContent(ix, iy);
      if (bin_content > 0 && bin_content < lowest_positive_count) {
        lowest_positive_count = bin_content;
      }
    }
  }
  return lowest_positive_count;
}

/**
 * modifies a given TH2D (in place) to change any bins with <=0 events to the 
 * smallest nonzero yield across the histogram
 */
void remove_zeros(TH2D* hist) {
  int nbins_x = hist->GetNbinsX();
  int nbins_y = hist->GetNbinsY();
  double lowest_positive_count = get_lowest_positive_yield(hist);
  for (int ix = 0; ix <= nbins_x+1; ix++) {
    for (int iy = 0; iy <= nbins_y+1; iy++) {
      if (hist->GetBinContent(ix, iy) <= 0) {
        hist->SetBinContent(ix, iy, lowest_positive_count);
      }
    }
  }
}

/**
 * modifies a given TH2D (in place) to change any bins with <=0 events to the 
 * given value
 */
void remove_zeros(TH2D* hist, float min_value) {
  int nbins_x = hist->GetNbinsX();
  int nbins_y = hist->GetNbinsY();
  for (int ix = 0; ix <= nbins_x+1; ix++) {
    for (int iy = 0; iy <= nbins_y+1; iy++) {
      if (hist->GetBinContent(ix, iy) <= 0) {
        hist->SetBinContent(ix, iy, min_value);
      }
    }
  }
}


/**
 * simple significance approximation
 */
vector<float> simple_significance(float s, float us, float b, float ub,
                                  float min_positive) {
  if (b <= 0) {
    b = min_positive;
    ub = min_positive;
  }
  float signif = s/sqrt(b);
  float sig_unc = us/sqrt(b);
  float bkg_unc = 0.5*ub/sqrt(b)/sqrt(b)/sqrt(b);
  float uncertainty = sqrt(sig_unc*sig_unc+bkg_unc*bkg_unc);
  return {signif, uncertainty};
}

/**
 * 4 bin simple significance approximation
 */
vector<float> simple_significance_4bins(float s[4], double us[4], float b[4], 
                                        double ub[4], float min_positive) {
  float bin_signifs[4] = {0.0, 0.0, 0.0, 0.0};
  float bin_uncs[4] = {0.0, 0.0, 0.0, 0.0};
  float total_signif = 0.0;
  float total_unc = 0.0;
  for (int i = 0; i < 4; i++) {
    vector<float> bin_results = simple_significance(s[i], us[i], b[i], ub[i],
                                                    min_positive);
    bin_signifs[i] = bin_results[0];
    bin_uncs[i] = bin_results[1];
    total_signif += bin_signifs[i]*bin_signifs[i];
  }
  total_signif = sqrt(total_signif);
  for (int i = 0; i < 4; i++) {
    bin_uncs[i] = bin_signifs[i]/total_signif*bin_uncs[i];
    total_unc += bin_uncs[i]*bin_uncs[i];
  }
  total_unc = sqrt(total_unc);
  return {total_signif, total_unc};
}

/**
 * does a simple s/sqrt(b) optimization of cuts in 2D
 */
void optimize_2d_binning(TH2D* signal_hist, TH2D* background_hist, 
                         TH2D* signal_test_hist, TH2D* background_test_hist, 
                         float min_signal, float min_positive) {
  int nbins_x = signal_hist->GetNbinsX();
  int nbins_y = signal_hist->GetNbinsY();
  float best_signif = 0.0;
  float test_best_signif = 0.0;
  float best_unc = 0.0;
  float test_best_unc = 0.0;
  float best_x_cut = 0.0;
  float best_y_cut = 0.0;
  float best_sig_yields[4] = {0.0, 0.0, 0.0, 0.0};
  float best_bkg_yields[4] = {0.0, 0.0, 0.0, 0.0};
  float best_sig_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
  float best_bkg_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
  //loop over all possible cuts, calculating significance by summing
  //bin significance in quadrature. Save best result
  for (int ix = 2; ix <= nbins_x+1; ix++) {
    for (int iy = 2; iy <= nbins_y+1; iy++) {
      int x_lo[4] = {1, 1, ix, ix};
      int x_hi[4] = {ix-1, ix-1, nbins_x, nbins_x};
      int y_lo[4] = {1, iy, 1, iy};
      int y_hi[4] = {iy-1, nbins_y, iy-1, nbins_y};
      float sig_yields[4] = {0.0, 0.0, 0.0, 0.0};
      float bkg_yields[4] = {0.0, 0.0, 0.0, 0.0};
      float sig_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
      float bkg_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
      double sig_uncs[4] = {0.0, 0.0, 0.0, 0.0};
      double bkg_uncs[4] = {0.0, 0.0, 0.0, 0.0};
      double sig_test_uncs[4] = {0.0, 0.0, 0.0, 0.0};
      double bkg_test_uncs[4] = {0.0, 0.0, 0.0, 0.0};
      bool signal_less_than_min = false;
      for (int i = 0; i < 4; i++) {
        sig_yields[i] = signal_hist->IntegralAndError(x_lo[i],x_hi[i],y_lo[i],
            y_hi[i],sig_uncs[i]);
        bkg_yields[i] = background_hist->IntegralAndError(x_lo[i],x_hi[i],
            y_lo[i],y_hi[i],bkg_uncs[i]);
        sig_test_yields[i] = signal_test_hist->IntegralAndError(x_lo[i],
            x_hi[i],y_lo[i],y_hi[i],sig_test_uncs[i]);
        bkg_test_yields[i] = background_test_hist->IntegralAndError(x_lo[i],
            x_hi[i],y_lo[i],y_hi[i],bkg_test_uncs[i]);
        if (sig_yields[i] < min_signal)
          signal_less_than_min = true;
      }
      if (!signal_less_than_min) {
        vector<float> signif_results = simple_significance_4bins(sig_yields,
            sig_uncs, bkg_yields, bkg_uncs, min_positive);
        vector<float> test_signif_results = simple_significance_4bins(
            sig_test_yields, sig_test_uncs, bkg_test_yields, bkg_test_uncs,
            min_positive);
        if (signif_results[0] > best_signif) {
          best_signif = signif_results[0];
          test_best_signif = test_signif_results[0];
          best_unc = signif_results[1];
          test_best_unc = test_signif_results[1];
          best_x_cut = signal_hist->GetXaxis()->GetBinLowEdge(ix);
          best_y_cut = signal_hist->GetYaxis()->GetBinLowEdge(iy);
          for (int i = 0; i < 4; i++) {
            best_sig_yields[i] = sig_yields[i];
            best_bkg_yields[i] = bkg_yields[i];
            best_sig_test_yields[i] = sig_test_yields[i];
            best_bkg_test_yields[i] = bkg_test_yields[i];
          }
        }
      }
    }
  }
  cout << "Optimized significance (train): " << best_signif << "+-" 
       << best_unc << "\n";
  cout << "Optimized significance (test): " << test_best_signif << "+-" 
       << test_best_unc << "\n";
  cout << "Optimized cuts: " << best_x_cut << ", " << best_y_cut << "\n";
  cout << "Optimized bin yields (train): ";
  for (int i = 0; i < 4; i++) {
    cout << "(" << best_sig_yields[i] << "," << best_bkg_yields[i] << ")";
  }
  cout << "\n";
  cout << "Optimized bin yields (test): ";
  for (int i = 0; i < 4; i++) {
    cout << "(" << best_sig_test_yields[i] << "," << best_bkg_test_yields[i] 
         << ")";
  }
  cout << "\n";
}

/**
 * does a simple s/sqrt(b) optimization of 4 bins in 1D
 */
void optimize_1d_binning(TH1D* signal_hist, TH1D* background_hist, 
                         TH1D* signal_test_hist, TH1D* background_test_hist,
                         float min_signal, float min_positive) {
  int nbins = signal_hist->GetNbinsX();
  float best_signif = 0.0;
  float test_best_signif = 0.0;
  float best_unc = 0.0;
  float test_best_unc = 0.0;
  float best_cut1 = 0.0;
  float best_cut2 = 0.0;
  float best_cut3 = 0.0;
  float best_sig_yields[4] = {0.0, 0.0, 0.0, 0.0};
  float best_bkg_yields[4] = {0.0, 0.0, 0.0, 0.0};
  float best_sig_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
  float best_bkg_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
  for (int i0 = 2; i0 < nbins-1; i0++) {
    for (int i1 = i0+1; i1 < nbins; i1++) {
      for (int i2 = i1+1; i2 < nbins+1; i2++) {
        int bnd[5] = {1,i0,i1,i2,nbins+1};
        float sig_yields[4] = {0.0, 0.0, 0.0, 0.0};
        float bkg_yields[4] = {0.0, 0.0, 0.0, 0.0};
        float sig_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
        float bkg_test_yields[4] = {0.0, 0.0, 0.0, 0.0};
        double sig_uncs[4] = {0.0, 0.0, 0.0, 0.0};
        double bkg_uncs[4] = {0.0, 0.0, 0.0, 0.0};
        double sig_test_uncs[4] = {0.0, 0.0, 0.0, 0.0};
        double bkg_test_uncs[4] = {0.0, 0.0, 0.0, 0.0};
        bool signal_less_than_min = false;
        for (int i = 0; i < 4; i++) {
          sig_yields[i] = signal_hist->IntegralAndError(bnd[i],bnd[i+1]-1,
              sig_uncs[i]);
          bkg_yields[i] = background_hist->IntegralAndError(bnd[i],bnd[i+1]-1,
              bkg_uncs[i]);
          sig_test_yields[i] = signal_test_hist->IntegralAndError(bnd[i],
              bnd[i+1]-1,sig_test_uncs[i]);
          bkg_test_yields[i] = background_test_hist->IntegralAndError(bnd[i],
              bnd[i+1]-1,bkg_test_uncs[i]);
          if (sig_yields[i] < min_signal)
            signal_less_than_min = true;
        }
        if (!signal_less_than_min) {
          vector<float> signif_results = simple_significance_4bins(sig_yields,
              sig_uncs, bkg_yields, bkg_uncs, min_positive);
          vector<float> test_signif_results = simple_significance_4bins(
              sig_test_yields, sig_test_uncs, bkg_test_yields, bkg_test_uncs,
              min_positive);
          if (signif_results[0] > best_signif) {
            best_signif = signif_results[0];
            test_best_signif = test_signif_results[0];
            best_unc = signif_results[1];
            test_best_unc = test_signif_results[1];
            best_cut1 = signal_hist->GetXaxis()->GetBinLowEdge(i0);
            best_cut2 = signal_hist->GetXaxis()->GetBinLowEdge(i1);
            best_cut3 = signal_hist->GetXaxis()->GetBinLowEdge(i2);
            for (int i = 0; i < 4; i++) {
              best_sig_yields[i] = sig_yields[i];
              best_bkg_yields[i] = bkg_yields[i];
              best_sig_test_yields[i] = sig_test_yields[i];
              best_bkg_test_yields[i] = bkg_test_yields[i];
            }
          }
        }
      }
    }
  }
  cout << "Optimized significance (train): " << best_signif << "+-" 
       << best_unc << "\n";
  cout << "Optimized significance (test): " << test_best_signif << "+-" 
       << test_best_unc << "\n";
  cout << "Optimized cuts: " << best_cut1 << ", " << best_cut2 << ", " 
       << best_cut3 << "\n";
  cout << "Optimized bin yields (train): ";
  for (int i = 0; i < 4; i++) {
    cout << "(" << best_sig_yields[i] << "," << best_bkg_yields[i] << ")";
  }
  cout << "\n";
  cout << "Optimized bin yields (test): ";
  for (int i = 0; i < 4; i++) {
    cout << "(" << best_sig_test_yields[i] << "," << best_bkg_test_yields[i] 
         << ")";
  }
  cout << "\n";
}

/**
 * class that implements a histogram-based 2D likelihood ratio esimator
 */
class HistLikelihoodRatio {

public: 

  /**
   * constructor from a signal and background histogram
   * also includes a minimum value to set 0-yield bins to
   */
  HistLikelihoodRatio(TH2D* signal_hist, TH2D* background_hist, float min_positive) {
    ratio_hist_ = static_cast<TH2D*>(signal_hist->Clone());
    ratio_hist_->Smooth();
    ratio_hist_->Scale(1.0/ratio_hist_->Integral());
    TH2D* background_hist_nonzero = static_cast<TH2D*>(background_hist->Clone());
    remove_zeros(background_hist_nonzero, min_positive);
    background_hist_nonzero->Smooth();
    background_hist_nonzero->Scale(1.0/background_hist_nonzero->Integral());
    ratio_hist_->Divide(background_hist_nonzero);
    delete background_hist_nonzero;
  }

  /**
   * Copy constructor
   */
  HistLikelihoodRatio(HistLikelihoodRatio& other) {
    ratio_hist_ = static_cast<TH2D*>(other.ratio_hist_->Clone());
  }

  /**
   * evaluator, which returns likelihood-based discriminant, specifically
   * 1/(1+ (likelihood_ratio)^-1), which ranges from 0 to 1
   */
  double evaluate(double var1, double var2) const {
    double likelihood_ratio = ratio_hist_->GetBinContent(
        ratio_hist_->FindBin(var1, var2));
    return 1.0/(1.0+1.0/likelihood_ratio);
  }

  /**
   * returns wrapper around evaluator for pyROOT usage
   */
  function<double(double,double)> get_evaluator() {
    return [this](double var1, double var2) {
      return this->evaluate(var1, var2);
    };
  }

  /**
   * destructor
   */
  ~HistLikelihoodRatio() {
    delete ratio_hist_;
  }

  TH2D* ratio_hist_;

private:

};
