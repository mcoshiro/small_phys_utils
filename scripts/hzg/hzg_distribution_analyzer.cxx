/**
 * @brief Distribution analyzer driver script for H->Zgamma analysis
 * Example usage:  root -l -q 'scripts/hzg/hzg_distribution_analyzer.cxx("../draw_pico/plots/zgboost_decorr__decorrbdt.root")'
 */
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"

#include "distribution_analyzer.hpp"

/**
 * @brief Distribution analyzer driver for H->Zgamma analysis
 *
 * 1 argument; filename of ROOT file with histogram from draw_pico
 */
int hzg_distribution_analyzer(std::string filename) {
  gSystem->AddIncludePath("-Iinc");
  gSystem->AddLinkedLibs("-Llib -lSmallPhysUtils");

  //if (argc<2) {std::cout << "Insufficient arguments.\n"; return 1;}

  //extract histograms from files, knowing draw_pico naming convention
  TFile* file = TFile::Open(filename.c_str(),"read");
  TCanvas* canvas = static_cast<TCanvas*>(file->Get("canvas"));
  TH1D *bkg_1 = nullptr, *bkg_2 = nullptr, *signal = nullptr;
  bool found_bkg_1 = false;
  bool found_bkg_2 = false;
  bool found_sig_1 = false;
  for (int iplot = 0; iplot < 200; iplot++) {
    if (!found_bkg_1) {
      bkg_1 = static_cast<TH1D*>(canvas->FindObject(
            ("bkg_Z+Fake Photon_"+std::to_string(iplot)).c_str()));
      if (bkg_1 != nullptr) found_bkg_1 = true;
    }
    if (!found_bkg_2) {
      bkg_2 = static_cast<TH1D*>(canvas->FindObject(
            ("bkg_Z+#gamma_"+std::to_string(iplot)).c_str()));
      if (bkg_2 != nullptr) found_bkg_2 = true;
    }
    if (!found_sig_1) {
      signal = static_cast<TH1D*>(canvas->FindObject(
            ("sig_H#rightarrow Z#gamma_"+std::to_string(iplot)).c_str()));
      if (signal != nullptr) found_sig_1 = true;
    }
  }
  if (!found_bkg_1) {std::cout << "Bkg 1 not found\n"; return 1;}
  if (!found_bkg_1) {std::cout << "Bkg 2 not found\n"; return 1;}
  if (!found_sig_1) {std::cout << "Sig 1 not found\n"; return 1;}
  TH1D* background = static_cast<TH1D*>(bkg_1->Clone());
  background->Add(bkg_2);

  binning_optimizer(signal, background);
  get_roc_auc(signal, background);
  return 0;
}
