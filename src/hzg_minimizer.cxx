/**
 * Script to do simultaneous fit to determine scale factors for H->Zgamma 
 * background
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "TFile.h"
#include "TH2.h"

using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::vector;
using ROOT::Minuit2::FCNBase;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnUserParameters;

//evaluates Chebyshev polynomial (first kind)
double chebyshev(double *coefs, int order, double x) {
  double result = 0.0;
  if (order>0)
    result += coefs[0];
  if (order>1)
    result += coefs[1]*x;
  if (order>2)
    result += coefs[2]*(2.0*x*x-1.0);
  if (order>3)
    result += coefs[3]*(4.0*x*x*x-3.0*x);
  if (order>4)
    result += coefs[4]*(8.0*x*x*x*x-8.0*x*x+1.0);
  return result;
}

//3D 2nd-order Chebyshev (27 coefficients)
double cheby2_3d(double* coefs, double x, double y, double z) {
  double result = 0.0;
  for (unsigned ix = 0; ix < 3; ix++) {
    double x_term = 1.0;
    if (ix==1)
      x_term = x;
    else if (ix==2)
      x_term = 2*x*x-1.0;
    for (unsigned iy = 0; iy < 3; iy++) {
      double y_term = 1.0;
      if (iy==1)
        y_term = y;
      else if (iy==2)
        y_term = 2*y*y-1.0;
      for (unsigned iz = 0; iz < 3; iz++) {
        double z_term = 1.0;
        if (iz==1)
          z_term = z;
        else if (iz==2)
          z_term = 2*z*z-1.0;
        result += coefs[ix*9+iy*3+iz]*x_term*y_term*z_term;
      }
    }
  }
  return result;
}

/**
 * Class defining function that is minimized (~chi2)
 */
class HistFit : public FCNBase {

public:

  //ad-hoc regularization to prefer smaller weights
  //approximatley in units of chi2/DoF
  const double reg_scale = 0.1;

  /**
   * Weight for true Z+photon
   */
  double weight_zpho(double zph_pt, vector<double> &coefs) {
    return chebyshev(&coefs[0], 4, log(zph_pt));
  }

  /**
   * Weight for Z+jet
   */
  double weight_dyjt(double ph_pt, double ph_abseta, 
                     double ph_id, vector<double> &coefs) {
    return chebyshev(&coefs[4], 4, log(zph_pt))
           *cheby2_3d(&coefs[8], ph_pt, ph_abseta, ph_id);
  }

  /**
   * Weight for Z+PU
   */
  double weight_dypu(double ph_pt, double ph_abseta, 
                     double ph_id, vector<double> &coefs) {
    return chebyshev(&coefs[35], 4, log(zph_pt))
           *cheby2_3d(&coefs[39], ph_pt, ph_abseta, ph_id);
  }

  /**
   * Function that returns loss for fit to determine coefficients
   */
  double operator()(const vector<double>& args) const {
    double loss = 0.0;
    for (unsigned izph_pt = 0; izph_pt < zph_pt_bin_means.size(); izph_pt++) {
      double zph_pt_mean = zph_pt_bin_means[izph_pt];
      for (unsigned iph_pt = 0; iph_pt < ph_pt_bin_means.size(); iph_pt++) {
        double ph_pt_mean = ph_pt_bin_means[iph_pt];
        for (unsigned iph_eta = 0; iph_eta < ph_eta_bin_means.size(); 
             iph_eta++) {
          double ph_eta_mean = ph_eta_bin_means[iph_eta];
          for (unsigned iph_id = 0; iph_id < ph_id_bin_means.size(); 
               iph_id++) {
            double ph_id_mean = ph_id_bin_means[iph_id];
            double sf_zpho = weight_zpho(zph_pt_mean, coefs);
            double sf_dyjt = weight_dyjt(ph_pt_mean, ph_eta_mean, ph_id_mean, 
                                         coefs);
            double sf_dypu = weight_dypu(ph_pt_mean, ph_eta_mean, ph_id_mean, 
                                         coefs);
            double mc_yield = 
                (1.0+sf_zpho)*zpho_yields[izph_pt][iph_pt][iph_eta][iph_id]
                +sf_dyjt*dyjt_yields[izph_pt][iph_pt][iph_eta][iph_id]
                +sf_dypu*dypu_yields[izph_pt][iph_pt][iph_eta][iph_id];
            if (mc_yield <= 0.0)
              loss -= data_yield/1.4;
            else {
              double diff = mc_yield
                  -data_yields[izph_pt][iph_pt][iph_eta][iph_id];
              //chi2 loss
              loss -= diff*diff/mc_yield;
            }
            //add a small amount of regularization to promote small sfs
            loss -= (sf_zpho*sf_zpho+sf_dyjt*sf_dyjt+sf_dypu*sf_dypu)
                    *reg_scale;
          }
        }
      }
    }
  }

  /**
   * load hists and bin means
   */
  void load_hists_and_bin_means() {
    TFile hist_file("bkg_weight_hists.root", "READ");
    for (unsigned izph_pt = 0; izph_pt < zph_pt_bin_means.size(); izph_pt++) {
      data_yields.push_back(vector<vector<vector<double>>>());
      zpho_yields.push_back(vector<vector<vector<double>>>());
      dyjt_yields.push_back(vector<vector<vector<double>>>());
      dypu_yields.push_back(vector<vector<vector<double>>>());
      for (unsigned iph_pt = 0; iph_pt < ph_pt_bin_means.size(); iph_pt++) {
        string hist_name = "hist"+to_string(izph_pt)+"_"+to_string(iph_pt);
        TH2D* hist_data = hist_file.Get<TH2D>(("data2_"+hist_name).c_str());
        TH2D* hist_zpho = hist_file.Get<TH2D>(("mczg2_"+hist_name).c_str());
        TH2D* hist_dyjt = hist_file.Get<TH2D>(("dyjt2_"+hist_name).c_str());
        TH2D* hist_dypu = hist_file.Get<TH2D>(("dypu2_"+hist_name).c_str());
        data_yields[izph_pt].push_back(vector<vector<double>>());
        zpho_yields[izph_pt].push_back(vector<vector<double>>());
        dyjt_yields[izph_pt].push_back(vector<vector<double>>());
        dypu_yields[izph_pt].push_back(vector<vector<double>>());
        for (unsigned iph_eta = 0; iph_eta < ph_eta_bin_means.size(); 
             iph_eta++) {
          data_yields[izph_pt][iph_pt].push_back(vector<double>());
          zpho_yields[izph_pt][iph_pt].push_back(vector<double>());
          dyjt_yields[izph_pt][iph_pt].push_back(vector<double>());
          dypu_yields[izph_pt][iph_pt].push_back(vector<double>());
          for (unsigned iph_id = 0; iph_id < ph_id_bin_means.size(); 
               iph_id++) {
            data_yields[izph_pt][iph_pt][iph_eta].push_back(
                hist_data->GetBinContent(iph_eta+1, iph_id+1));
            zpho_yields[izph_pt][iph_pt][iph_eta].push_back(
                hist_zpho->GetBinContent(iph_eta+1, iph_id+1));
            dyjt_yields[izph_pt][iph_pt][iph_eta].push_back(
                hist_dyjt->GetBinContent(iph_eta+1, iph_id+1));
            dypu_yields[izph_pt][iph_pt][iph_eta].push_back(
                hist_dypu->GetBinContent(iph_eta+1, iph_id+1));
          }
        }
      }
    }
    hist_file.Close();
  }

  //constructor
  HistFit() {
    load_hists_and_bin_means();
  }

  vector<vector<vector<vector<double>>>> data_yields;
  vector<vector<vector<vector<double>>>> zpho_yields;
  vector<vector<vector<vector<double>>>> dyjt_yields;
  vector<vector<vector<vector<double>>>> dypu_yields;
  //try to minimize dimensionality
  //pT Zph: 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 30.0, 40.0, 80.0, 200.0
  //pT: 15.0, 17.5, 20.0, 30.0, 50.0
  //abseta: 0.0, 0.8, 1.5, 2.0, 2.5
  //ID: 0.14, 0.33, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0
  const vector<double> zph_pt_bin_means = {2.0,6.0,10.0,14.0,18.0,25.0,35.0,
                                           60.0,140.0};
  const vector<double> ph_pt_bin_means = {16.25,18.75,25.0,40.0};
  const vector<double> ph_eta_bin_means = {0.4,1.15,1.75,2.25};
  const vector<double> ph_id_bin_means = {0.235,0.415,0.6,0.75,0.85,0.925,
                                          0.975};
}

const int n_params = (4*3+27*2);
const vector<double> param_defaults(n_params, 0.0);
const vector<double> param_errors(n_params, 0.1);

/**
 * Driver that runs minimization
 */
int main() {
  HistFit hist_fit;
  MnUserParameters upar;
  //4*3+27*2 parameters
  for (unsigned ipar = 0; ipar<(4*3+27*2); ipar++)
    upar.add("par"+to_string(ipar), param_defaults[ipar], param_errors[ipar])
  MnMigrad migrad(hist_fit, upar);
  FunctionMinimum min = migrad();
  cout << "Fit result: " << min << endl;
}
