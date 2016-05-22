#include "simulate.h"

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <tuple>

#include <armadillo>

using namespace std;

// parameters
const double sigma = 1;         // standard deviation

int main(int argc, char *argv[]) {
  clock_t startTime = clock();
  ios::sync_with_stdio(false); cin.tie(NULL);
  // parse arguments
  int ITERATIONS = stoi(argv[1]); // simulation runs
  int N = stoi(argv[2]);          // 1/2 sample size
  bool interceptTerm = argc >= 4 && strcmp(argv[3], "intercept") == 0 ? true : false;
  
  // set up generators
  random_device randomDevice;
  mt19937_64 rng; 
  // rng.seed(randomDevice());
  rng.seed(2016);
  normal_distribution<double> normalMean0(0, sigma);

  // store and print results
  vector<tuple<double, double, int, int, double, double>> results(ITERATIONS); // initial p-value, "cheated" beta, set size, subset size, balance p-value, "cheated" p-value
  cout << setprecision(12) << fixed;
  cout << "2N\tp.value\tbeta\tset.size\tsubset.size\tbalance.p.value\tfake.p.value" << endl;

  // randomly assign treatment
  vector<int> X(2*N); 
  for (int i = 0; i < N; ++i) X[i] = 1;  
  shuffle(X.begin() , X.end(), rng);  
  arma::mat Z(2*N, 2*N, arma::fill::zeros); // pre allocate Z matrix
  for (int i = 0; i < ITERATIONS; ++i) {
    // generate new response for each iteration
    arma::Col<double> Y(2*N);
    for (int j = 0; j < 2*N; ++j) Y(j) = normalMean0(rng);
    results[i] = simulate(Y, X, sigma, false, Z, rng, interceptTerm);
    cout << 2*N << '\t' 
         << get<0>(results[i]) << '\t' 
         << get<1>(results[i]) << '\t' 
         << get<2>(results[i]) << '\t' 
         << get<3>(results[i]) << '\t'
         << get<4>(results[i]) << '\t'
         << get<5>(results[i]) << endl;
  }

  double duration = (clock() - startTime) / (double) CLOCKS_PER_SEC;
  cerr << "Time taken (seconds): " << duration << endl; 
  return 0;
}
