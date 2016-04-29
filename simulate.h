#ifndef SIMULATE_H
#define SIMULATE_H

#include <random>
#include <vector>
#include <tuple>

#include <armadillo>

using namespace std;

// return variance needed to achieve power with level alpha test given 2N observations
double calculateVariance(int N, double alpha, double power);

// return initial p value, final beta, set size, subset size, balance p-value
tuple<double, double, int, int, double> simulate(const arma::Col<double> &Y, const vector<int> X, 
                                                 double sigma, bool varianceKnown,
                                                 arma::mat &Z, mt19937_64 &rng);

// both must vectors of 0s and 1s
double testBalance(const arma::Col<double> &z, const arma::Col<double> &t);

#endif
