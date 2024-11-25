//script to do test of Gaussian Process regression

//stdlib
#include <cmath>
#include <vector>
//ROOT
#include "TGraph.h" 
#include "TRandom.h"

using std::vector;

void gp_test() {

  //generate toy data
  TRandom3 rng;
  vector<double> x = {100.0,105.0,110.0,115.0,120.0,130.0,135.0,140.0,145.0,150.0,155.0,160.0};
  vector<double> y = {0.004,0.012,0.017,0.017,0.015,0.0105,0.009,0.0075,0.0065,0.0055,0.0045,0.004};
  int n_samples = 1000;
  for (unsigned iy = 0; iy < y.size(); iy++) {
    y[iy] = rng.Gaus(n_samples*y[iy],sqrt(n_samples*y[iy]));
  }
  TGraph g(x.size(), &x[0], &y[0]);

  //do GPR
  //find hyperparameters

}
