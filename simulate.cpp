#include "simulate.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <utility>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace std;

const boost::math::normal STANDARD_NORMAL(0, 1);

double calculateVariance(int N, double alpha, double power) {
  double z = boost::math::quantile(boost::math::complement(STANDARD_NORMAL, alpha/2)); // critical value
  double lower = 1;
  double upper = 100000000;
  double mid = (upper + lower)/2;
  while (upper - lower > 0.00000001) {
    double val = boost::math::cdf(STANDARD_NORMAL, -z + sqrt(N)/sqrt(mid)) + boost::math::cdf(STANDARD_NORMAL, -z - sqrt(N)/sqrt(mid));
    if (val > power) {
      lower = mid;
    } else {
      upper = mid;
    }
    mid = (upper + lower)/2;
  }
  return mid;
}

pair<double, double> calculateBetaPValue(const arma::mat &Z, const arma::Col<double> &Y, 
                                         double sigma, bool varianceKnown) {
  arma::mat ZZ = Z.t()*Z;
  arma::mat beta = arma::solve(ZZ, Z.t()*Y);
  arma::Col<double> e(ZZ.n_rows, arma::fill::zeros); e(0) = 1;
  arma::mat inverseFirstColumn = arma::solve(ZZ, e);
  if (varianceKnown) {
    double standardDeviation = sigma*sqrt(inverseFirstColumn(0));
    double t = abs(beta(0))/standardDeviation;  
    return make_pair(2*cdf(STANDARD_NORMAL, -t), beta(0));
  } else {    
    int df = Y.n_rows - beta.n_rows;            // degrees of freedom
    double s = arma::norm(Y - Z*beta)/sqrt(df); // sample standard deviation
    double standardDeviation = s*sqrt(inverseFirstColumn(0));
    double t = abs(beta(0))/standardDeviation;  
    boost::math::students_t tDist(df);
    return make_pair(2*cdf(tDist, -t), beta(0));
  }
}

tuple<double, double, int, int> simulate(const arma::Col<double> &Y, const vector<int> X, 
                                 double sigma, bool varianceKnown,
                                 arma::mat &Z, mt19937_64 &rng) {
  bernoulli_distribution bernoulli(0.5);
  int N = X.size();
  Z.fill(0);
  copy(X.begin(), X.end(), Z.begin_col(0));
  vector<arma::uvec> bestColumns{arma::uvec{0}}; bestColumns.reserve(N - 1);  
  vector<pair<double, double>> bestPValues; bestPValues.reserve(N - 1);
  bestPValues.push_back(calculateBetaPValue(Z.cols(bestColumns.front()), Y, sigma, varianceKnown));
  if (bestPValues.front().first <= 0.05) {
    return make_tuple(bestPValues.front().first, bestPValues.front().second, 0 ,0);
  } else {                    // need more covariates
    bool done = false;
    int smallestSubsetSize = INT_MAX;
    for (int j = 1; j < N - 1; ++j) { // add covariates one-by-one
      for (int k = 0; k < N; ++k) Z(k, j) = bernoulli(rng);
      while (arma::rank(Z) <= j) {
        for (int k = 0; k < N; ++k) Z(k, j) = bernoulli(rng);
      }        
      for (int k = j; k >= 1; --k) { // loop through subset sizes
        pair<double, double> newPValue;
        if (k == j) {
          bestColumns.emplace_back(j + 1);
          for (int l = 0; l <= k; ++l) bestColumns.back()(l) = l;
          newPValue = calculateBetaPValue(Z.cols(0, k), Y, sigma, varianceKnown);
          bestPValues.push_back(newPValue);
        } else {
          arma::uvec columnSubset(k + 1); 
          for (int l = 0; l < k; ++l) columnSubset(l) = bestColumns[k - 1](l); 
          columnSubset(k) = j;
          newPValue = calculateBetaPValue(Z.cols(columnSubset), Y, sigma, varianceKnown);
          if (bestPValues[k].first > newPValue.first) {
            bestPValues[k] = newPValue;
            bestColumns[k] = columnSubset;
          }
        }
        if (newPValue.first <= 0.05) {
          done = true;
          smallestSubsetSize = k;
        }
      }
      if (done) return make_tuple(bestPValues.front().first, bestPValues[smallestSubsetSize].second, j, smallestSubsetSize); 
    }    
  }  
  return make_tuple(bestPValues.front().first, bestPValues.front().second, -1, -1);
}
