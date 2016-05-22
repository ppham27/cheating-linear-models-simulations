#include "simulate.h"

#include <algorithm>
#include <exception>
#include <climits>
#include <cmath>
#include <utility>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace std;

const boost::math::normal STANDARD_NORMAL(0, 1);

double testBalance(const arma::Col<double> &z, const arma::Col<double> &t) {
  if (z.n_rows != t.n_rows) throw length_error("Column vectors must be the same length.");
  int N = t.n_rows;             // population size
  int n = 0;                    // draws
  int K = 0;                    // "success" events
  int k = 0;                    // observed successes
  for (int i = 0; i < N; ++i) {
    if (!(t(i) == 0 || t(i) == 1)) throw logic_error("Treatment must be a vector of 0s and 1s.");
    if (!(z(i) == 0 || z(i) == 1)) throw logic_error("Covariate must be a vector of 0s and 1s.");
    if (t(i) == 1) ++n;         // assigning to treatment is draw
    if (z(i) == 1) ++K;         // number of possible success events
    if (t(i) == 1 && z(i) == 1) ++k; // observed successes
  }    
  boost::math::hypergeometric_distribution<> X(K, n, N);
  // find Prob(|X - K/N*n| >= |k - K/N*n|) = Prob(X <= K/N*n - |k - K/N*n|) + Prob(X >= K/N*n + |k - K/N*n|)
  double p = ((double) K)/N;  
  double mean = p*n;
  double delta = abs(k - p*n);
  double pValue = 0;  
  int L = mean - delta;
  int U = mean + delta - 1; // complement is not inclusive
  return U < L ? 1 : boost::math::cdf(X, L) + boost::math::cdf(boost::math::complement(X, U));
}

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

tuple<double, double, int, int, double, double> simulate(const arma::Col<double> &Y, const vector<int> X, 
                                                         double sigma, bool varianceKnown,
                                                         arma::mat &Z, mt19937_64 &rng,
                                                         bool interceptTerm) {
  bernoulli_distribution bernoulli(0.5);
  int N = X.size();
  Z.fill(0);
  // bestColumns[k] keeps track of the k + 1 or k + 2 columns that produce the smallest p-value depending on interceptTerm
  vector<arma::uvec> bestColumns; bestColumns.reserve(N - 1); 
  if (interceptTerm) { // make intercept term last column of Z
    fill(Z.begin_col(N - 1), Z.end_col(N - 1), 1);
    copy(X.begin(), X.end(), Z.begin_col(0));
    bestColumns.push_back(arma::uvec{0, (unsigned long long) N - 1ULL}); 
  } else {
    copy(X.begin(), X.end(), Z.begin_col(0));
    bestColumns.push_back(arma::uvec{0});     
  }  
  // bestPValues[k] corresponds to p-value if the columns bestColumns[k] are used
  vector<pair<double, double>> bestPValues; bestPValues.reserve(N - 1);
  bestPValues.push_back(calculateBetaPValue(Z.cols(bestColumns.front()), Y, sigma, varianceKnown));
  if (bestPValues.front().first <= 0.05) {
    return make_tuple(bestPValues.front().first, bestPValues.front().second, 0, 0, -1, bestPValues.front().first);
  } else {                    // need more covariates
    bool done = false;
    int smallestSubsetSize = INT_MAX;
    /* add covariates one-by-one, we always include the treatment
     * if we're using the intercept two covariates are included by default
     */
    for (int j = 1; j < N - 2 || (j == N - 2 && !interceptTerm); ++j) { 
      for (int k = 0; k < N; ++k) Z(k, j) = bernoulli(rng);
      if (!interceptTerm) {
        while (arma::rank(Z) <= j) {
          for (int k = 0; k < N; ++k) Z(k, j) = bernoulli(rng);
        }        
      } else { // offset rank by 1 for intercept term
        while (arma::rank(Z) <= j + 1) {
          for (int k = 0; k < N; ++k) Z(k, j) = bernoulli(rng);
        }        
      }
      for (int k = j; k >= 1; --k) { // loop through subset sizes, k is the number of additional covariates
        pair<double, double> newPValue;
        if (k == j) {           // use all available covariates
          bestColumns.emplace_back(bestColumns.back().n_rows + 1); // add one more to biggest subset
          for (int l = 0; l < bestColumns.back().n_rows - 1; ++l) { 
            bestColumns.back()(l) = bestColumns[j - 1](l); // copy over from original subset
          }
          bestColumns.back()(bestColumns.back().n_rows - 1) = j; // add new covariate
          newPValue = calculateBetaPValue(Z.cols(bestColumns.back()), Y, sigma, varianceKnown);
          bestPValues.push_back(newPValue);
        } else {                // make a new subset of same size with new covariate
          arma::uvec columnSubset(bestColumns[k].n_rows); 
          for (int l = 0; l < columnSubset.n_rows - 1; ++l) 
            columnSubset(l) = bestColumns[k - 1](l); // copy over from smaller subset
          columnSubset(columnSubset.n_rows - 1) = j; // add new covariate
          newPValue = calculateBetaPValue(Z.cols(columnSubset), Y, sigma, varianceKnown);
          if (bestPValues[k].first > newPValue.first) { // if better subset replace
            bestPValues[k] = newPValue;
            bestColumns[k] = columnSubset;
          }
        }
        if (newPValue.first <= 0.05) { // stop when we reach significance
          done = true;
          smallestSubsetSize = k;
        }
      }
      if (done) {
        // compute balance p value in special case that only 1 covariate was needed
        double balancePValue = -1;
        if (smallestSubsetSize == 1 && !interceptTerm) {
          balancePValue = testBalance(Z.col(bestColumns[1](1)), Z.col(0));
        } else if (smallestSubsetSize == 1 && interceptTerm) {
          balancePValue = testBalance(Z.col(bestColumns[1](2)), Z.col(0));
        }
        return make_tuple(bestPValues.front().first, bestPValues[smallestSubsetSize].second, 
                          j, smallestSubsetSize, balancePValue, bestPValues[smallestSubsetSize].first); 
      }
    }    
  }  
  return make_tuple(bestPValues.front().first, bestPValues.front().second, -1, -1, -1, bestPValues.front().first);
}
