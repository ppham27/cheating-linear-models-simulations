#ifndef SIMULATE_H
#define SIMULATE_H

#include <random>
#include <vector>
#include <tuple>

#include <armadillo>

using namespace std;

// return variance needed to achieve power with level alpha test given 2N observations
double calculateVariance(int N, double alpha, double power);

// return initial p value, actual beta, set size, subset size
tuple<double, double, int, int> simulate(const arma::Col<double> &Y, const vector<int> X, 
                                         double sigma, bool varianceKnown,
                                         arma::mat &Z, mt19937_64 &rng);

#endif
