#include <iomanip>
#include <iostream>
#include <random>

#include <armadillo>

#include "simulate.h"

using namespace std;

const int ITERATIONS = 100000;
const int N = 100;

int main(int argc, char *argv[]) {
  arma::Col<double> t{1,1,1,1,1,1,0,0,0,0,0,0};
  arma::Col<double> z{1,1,1,1,1,0,1,0,0,0,0,0}; 
  cout << testBalance(z, t) << endl; // expect 0.08
  // arma::Col<double> t{1,1,1,1,1,1,0,0,0,0,0,0};
  // arma::Col<double> z{1,1,1,1,1,0,1,1,0,0,0,0}; 
  // cout << testBalance(z, t) << endl; // expect 0.24
  // arma::Col<double> t{1,1,1,1,1,1,0,0,0,0,0,0};
  // arma::Col<double> z{1,1,1,1,1,1,1,1,0,0,0,0}; 
  // cout << testBalance(z, t) << endl; // expect 0.06
  // arma::Col<double> t{1,1,1,1,1,1,0,0,0,0,0,0};
  // arma::Col<double> z{1,1,1,1,1,0,0,0,0,0,0,0}; 
  // cout << testBalance(z, t) << endl; // expect 0.01515152
  // arma::Col<double> t{1,1,1,1,1,1,1,0,0,0,0,0,0,0};
  // arma::Col<double> z{1,1,1,1,1,1,1,1,0,0,0,0,0,0}; 
  // cout << testBalance(z, t) << endl; // expect 0.004662

  // bernoulli_distribution bernoulli(0.5);
  // // set up generators
  // random_device randomDevice;
  // mt19937_64 rng; rng.seed(randomDevice());
  // arma::Col<double> t(2*N, arma::fill::zeros);
  // for (int n = 0; n < N; ++n) t(n) = 1;
  // arma::Col<double> z(2*N);
  // for (int n = 0; n < 2*N; ++n) z(n) = bernoulli(rng);
  // int successes = 0;
  // for (int i = 0; i < ITERATIONS; ++i) {
  //   shuffle(t.begin(), t.end(), rng);
  //   double pValue = testBalance(z, t);
  //   if (pValue <= 0.05) ++successes;
  // }
  // cout << (double) successes/ITERATIONS << endl;
  return 0;
}

