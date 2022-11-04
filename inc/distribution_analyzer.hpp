/**
 * @brief Collection of functions for analyzing distributions
 */

#include "TH1D.h"

/**
 * @brief Returns area under ROC curve for given distributions
 *
 * @param signal            pointer to ROOT histogram of signal
 * @param background        pointer to ROOT histogram of background
 */
void get_roc_auc(TH1D* signal, TH1D* background);

/**
 * @brief Function to find optimized binning cuts given distribution as histograms
 *
 * @param signal            pointer to ROOT histogram of signal
 * @param background        pointer to ROOT histogram of background
 * @param max_nbins         maximum number of bins to consider (ncuts + 1)
 * @param min_signal_yield  minimum signal yield allowed in a bin
 * @param scale             scale factor for signal and background to change effective lumi
 */
void binning_optimizer(TH1D* signal, TH1D* background, int max_nbins=6, 
    float min_signal_yield=0.5, float scale=1.0);
