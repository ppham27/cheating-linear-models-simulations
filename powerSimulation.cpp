#include "simulate.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <tuple>

#include <armadillo>
#include <boost/math/distributions/normal.hpp>

using namespace std;

// constants
const double ALPHA = 0.05;

int main(int argc, char *argv[]) {
  clock_t startTime = clock();
  ios::sync_with_stdio(false); cin.tie(NULL);
  // parse arguments
  int ITERATIONS = stoi(argv[1]); // simulation runs
  int N = stoi(argv[2]);          // 1/2 sample size
  double POWER = stod(argv[3]);    // power

  // calculate variance used for generators
  double variance = calculateVariance(N, ALPHA, POWER);
  double standardDeviation = sqrt(variance);
  
  // set up generators
  random_device randomDevice;
  mt19937_64 rng; 
  rng.seed(randomDevice());
  normal_distribution<double> normalMean0(0, standardDeviation);  
  normal_distribution<double> normalMean1(1, standardDeviation);  

  // store and print results
  vector<tuple<double, double, int, int, double, double>> results(ITERATIONS); // initial p-value, "cheated" beta, set size, subset size, balance p-value, "cheated" p-value
  cout << setprecision(12) << fixed;
  cout << "2N\tpower\tp.value\tbeta\tset.size\tsubset.size\tbalance.p.value\tfake.p.value" << endl;
  // randomly assign treatment
  vector<int> X(2*N); 
  for (int i = 0; i < N; ++i) X[i] = 1;  
  shuffle(X.begin() , X.end(), rng);  
  arma::mat Z(2*N, 2*N, arma::fill::zeros); // pre allocate Z matrix
  for (int i = 0; i < ITERATIONS; ++i) {
    // generate response
    arma::Col<double> Y(2*N);
    for (int j = 0; j < 2*N; ++j) {
      Y(j) = X[j] == 0 ? normalMean0(rng) : normalMean1(rng);
    }            
    results[i] = simulate(Y, X, standardDeviation, true, Z, rng);
    cout << 2*N << '\t' << POWER << '\t' 
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
