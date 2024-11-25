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
size_t find_nth(std::string search_string, std::string target, int n) {
  size_t pos = 0;
  for (int i = 0; i < n; i++) {
    pos = search_string.find(target, pos);
    if (pos==std::string::npos) return pos;
    pos++;
  }
  return pos;
}

/**
 * @brief Standard constructor
 *
 * @param filename   file name from which to load histograms
 * @param var_names  names of variables to use in evaluation (in order)
 */
HistMVAEvaluator::HistMVAEvaluator(std::string filename, std::vector<std::string> var_names) {
  TFile output_file(filename.c_str(),"READ");
  std::string th1_name = "TH1F";
  std::string th2_name = "TH2F";
  for (TObject *hist : output_file.GetList()) {
    if (hist->ClassName() == th1_name) {
      std::string hist_name = static_cast<TH1F*>(hist)->GetName();
      if (hist_name.substr(0,10) == "hist1d_sig") {
        hist_1d_sig.push_back(static_cast<TH1F*>(hist->Clone()));
        hist_1d_sig.back()->SetDirectory(0);
        std::string hist_var_name = hist_name.substr(11, hist_name.length()-11);
        bool found_var = false;
        for (unsigned ivar = 0; ivar < var_names.size(); ivar++) {
          if (var_names[ivar]==hist_var_name) {
            oned_var_idx_sig.push_back(ivar);
            found_var = true;
          }
        }
        if (!found_var) {
          std::cout << "ERROR: no variable " << hist_var_name << " found" << std::endl;
          exit(1);
        }
      }
      if (hist_name.substr(0,10) == "hist1d_bak") {
        hist_1d_bak.push_back(static_cast<TH1F*>(hist->Clone()));
        hist_1d_bak.back()->SetDirectory(0);
        std::string hist_var_name = hist_name.substr(11, hist_name.length()-11);
        bool found_var = false;
        for (unsigned ivar = 0; ivar < var_names.size(); ivar++) {
          if (var_names[ivar]==hist_var_name) {
            oned_var_idx_bak.push_back(ivar);
            found_var = true;
          }
        }
        if (!found_var) {
          std::cout << "ERROR: no variable " << hist_var_name << " found" << std::endl;
          exit(1);
        }
      }
    }
    if (hist->ClassName() == th2_name) {
      std::string hist_name = static_cast<TH2F*>(hist)->GetName();
      if (hist_name.substr(0,10) == "hist2d_sig") {
        hist_2d_sig.push_back(static_cast<TH2F*>(hist->Clone()));
        hist_2d_sig.back()->SetDirectory(0);
        size_t pos_var2 = find_nth(hist_name, "_", 3)+1;
        std::string hist_var1_name = hist_name.substr(11, pos_var2-1-11);
        std::string hist_var2_name = hist_name.substr(pos_var_2, hist_name.length()-pos_var_2);
        bool found_var1 = false;
        bool found_var2 = false;
        for (unsigned ivar = 0; ivar < var_names.size(); ivar++) {
          if (var_names[ivar]==hist_var1_name) {
            twod_var1_idx_sig.push_back(ivar);
            found_var1 = true;
          }
          if (var_names[ivar]==hist_var2_name) {
            twod_var2_idx_sig.push_back(ivar);
            found_var2 = true;
          }
        }
        if (!found_var1) {
          std::cout << "ERROR: no variable " << hist_var1_name << " found" << std::endl;
          exit(1);
        }
        if (!found_var2) {
          std::cout << "ERROR: no variable " << hist_var2_name << " found" << std::endl;
          exit(1);
        }
      }
      if (hist_name.substr(0,10) == "hist2d_bak") {
        hist_2d_bak.push_back(static_cast<TH2F*>(hist->Clone()));
        hist_2d_bak.back()->SetDirectory(0);
        size_t pos_var2 = find_nth(hist_name, "_", 3)+1;
        std::string hist_var1_name = hist_name.substr(11, pos_var2-1-11);
        std::string hist_var2_name = hist_name.substr(pos_var_2, hist_name.length()-pos_var_2);
        bool found_var1 = false;
        bool found_var2 = false;
        for (unsigned ivar = 0; ivar < var_names.size(); ivar++) {
          if (var_names[ivar]==hist_var1_name) {
            twod_var1_idx_bak.push_back(ivar);
            found_var1 = true;
          }
          if (var_names[ivar]==hist_var2_name) {
            twod_var2_idx_bak.push_back(ivar);
            found_var2 = true;
          }
        }
        if (!found_var1) {
          std::cout << "ERROR: no variable " << hist_var1_name << " found" << std::endl;
          exit(1);
        }
        if (!found_var2) {
          std::cout << "ERROR: no variable " << hist_var2_name << " found" << std::endl;
          exit(1);
        }
      }
    }
  }
  output_file.Close();
}

/**
 * @brief Evaluates MVA-based discriminant
 *
 * @param inputs   input variable values
 */
float HistMVAEvaluator::evaluate(std::vector<float> const & inputs) const {
  float prob_sig = 1.0;
  float prob_bak = 1.0;
  for (unsigned ihist = 0; ihist < hist_1d_sig.size(); ihist++) {
    prob_sig *= hist_1d_sig[ihist]->GetBinContent(hist_1d_sig[ihist]
                                  ->FindBin(inputs[oned_var_idx_sig[ihist]]));
  }
  for (unsigned ihist = 0; ihist < hist_1d_bak.size(); ihist++) {
    prob_bak *= hist_1d_bak[ihist]->GetBinContent(hist_1d_bak[ihist]
                                  ->FindBin(inputs[oned_var_idx_bak[ihist]]));
  }
  for (unsigned ihist = 0; ihist < hist_2d_sig.size(); ihist++) {
    prob_sig *= hist_2d_sig[ihist]->GetBinContent(hist_2d_sig[ihist]
                                  ->FindBin(inputs[twod_var1_idx_sig[ihist]],
                                            inputs[twod_var2_idx_sig[ihist]]));
  }
  for (unsigned ihist = 0; ihist < hist_2d_bak.size(); ihist++) {
    prob_bak *= hist_2d_bak[ihist]->GetBinContent(hist_2d_bak[ihist]
                                  ->FindBin(inputs[twod_var1_idx_bak[ihist]],
                                            inputs[twod_var2_idx_bak[ihist]]));
  }
  return 1.0/(1.0-prob_bak/prob_sig);
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
    unsigned ntrain_background) {

  //initialize
  unsigned n_inputs = 0;
  std::vector<std::string> unrolled_names;
  //for hists, trees, etc., 0 is signal, 1 is background
  std::vector<std::vector<TH1F>> hist_1d;
  std::vector<std::vector<TH2F>> hist_2d;
  hist_1d.push_back(std::vector<TH1F>());
  hist_2d.push_back(std::vector<TH2F>());
  for (unsigned icorr = 0; icorr < features.size(); icorr++) {
    n_inputs += features[icorr].size();
    if (features[icorr].size() == 1) {
      hist_1d[0].push_back(TH1D("hist1d_sig_"+features[icorr][0].name).c_str(),"",
                           features[icorr][0].nbins, features[icorr][0].min,
                           features[icorr][0].max);
      hist_1d[1].push_back(TH1D("hist1d_bak_"+features[icorr][0].name).c_str(),"",
                           features[icorr][0].nbins, features[icorr][0].min,
                           features[icorr][0].max);
      unrolled_names.push_back(features[icorr][0].name);
    }
    else if (features[icorr].size() == 2) {
      hist_2d[0].push_back(TH2D("hist2d_sig_"+features[icorr][0].name+"_"
                                +features[icorr][1].name).c_str(),"",
                           features[icorr][0].nbins, features[icorr][0].min,
                           features[icorr][0].max, features[icorr][1].nbins, 
                           features[icorr][1].min, features[icorr][1].max);
      hist_2d[1].push_back(TH2D("hist2d_bak_"+features[icorr][0].name+"_"
                                +features[icorr][1].name).c_str(),"",
                           features[icorr][0].nbins, features[icorr][0].min,
                           features[icorr][0].max, features[icorr][1].nbins, 
                           features[icorr][1].min, features[icorr][1].max);
      unrolled_names.push_back(features[icorr][0].name);
      unrolled_names.push_back(features[icorr][1].name);
    }
    else {
      std::cout << "ERROR: Sets of 3+ correlated variables not yet supported" 
                << std::endl;
      return;
    }
  }
  std::vector<float> inputs;
  float weight;
  inputs.resize(n_inputs);
  for (unsigned iin = 0; iin < n_inputs; iin++) {
    signal_tree->SetBranchAddress(unrolled_names[iin].c_str(), &inputs[iin]);
    background_tree->SetBranchAddress(unrolled_names[iin].c_str(), &inputs[iin]);
  }
  if (weight_name != "") {
    signal_tree->SetBranchAddress(weight_name.c_str(), &weight);
    background_tree->SetBranchAddress(weight_name.c_str(), &weight);
  }
  std::vector<unsigned> ntrain = {ntrain_signal, ntrain_background};
  std::vector<TTree*> tree = {signal_tree, background_tree};

  //fill histograms ("train")
  for (unsigned isample = 0; isample < 2; isample++) {
    if (ntrain[isample]==0) ntrain[isample] = tree[isample]->GetEntries()/2;
    for (long ievt = 0; ievt < static_cast<long>(ntrain[isample]); ievt++) {
      tree[isample]->GetEntry(ievt);
      unsigned hist_1d_idx = 0;
      unsigned hist_2d_idx = 0;
      unsigned feature_idx = 0;
      for (unsigned icorr = 0; icorr < features.size(); icorr++) {
        if (features[icorr].size() == 1) {
          hist_1d[isample][hist_1d_idx].Fill(inputs[feature_idx]);
          hist_1d_idx++;
          features_idx += 1;
        }
        else if (features[icorr].size() == 2) {
          hist_2d[isample][hist_2d_idx].Fill(inputs[feature_idx], 
                                             inputs[feature_idx+1]);
          hist_2d_idx++;
          features_idx += 2;
        }
      }
    }
  }

  //save output
  TFile output_file(output_filename.c_str(),"RECREATE");
  for (unsigned isample = 0; isample < 2; isample++) {
    for (TH1F & hist : hist_1d[isample]) {
      hist.Scale(1.0/hist.Integral());
      output_file.Write(&hist);
    }
    for (TH2F & hist : hist_2d[isample]) {
      hist.Scale(1.0/hist.Integral());
      output_file.Write(&hist);
    }
  }
  output_file.Close();
  
  //evaluate
  HistMVAEvaluator mva(output_filename.c_str(), unrolled_names);
  std::vector<TH1D> train_hist;
  train_hist.push_back(TH1D("sig_train_hist","",100,0,1));
  train_hist.push_back(TH1D("bak_train_hist","",100,0,1));
  std::vector<TH1D> tests_hist;
  tests_hist.push_back(TH1D("sig_tests_hist","",100,0,1));
  tests_hist.push_back(TH1D("bak_tests_hist","",100,0,1));
  for (unsigned isample = 0; isample < 2; isample++) {
    //training data loop
    for (long ievt = 0; ievt < static_cast<long>(ntrain[isample]); ievt++) {
      tree[isample]->GetEntry(ievt);
      train_hist[isample].Fill(mva.evaluate(inputs));
    }
    //testing data loop
    for (long ievt = static_cast<long>(ntrain[isample]); ievt < tree[isample].GetEntries(); ievt++) {
      tree[isample]->GetEntry(ievt);
      tests_hist[isample].Fill(mva.evaluate(inputs));
    }
  }

  std::cout << "Training data: " << std::endl;
  binning_optimizer(&train_hist[0], &train_hist[1]);
  get_roc_auc(&train_hist[0], &train_hist[1]);
  std::cout << "Testing data: " << std::endl;
  binning_optimizer(&tests_hist[0], &tests_hist[1]);
  get_roc_auc(&tests_hist[0], &tests_hist[1]);
}

