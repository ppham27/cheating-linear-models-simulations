#include "simulate.h"

#include <iomanip>
#include <iostream>

#include <armadillo>
#include <boost/math/distributions/normal.hpp>

using namespace std;

int main(int argc, char *argv[]) {
  int N = stoi(argv[1]);          // 1/2 sample size
  double ALPHA = stod(argv[2]);   // level of significance
  double POWER = stod(argv[3]); // power
  cout << calculateVariance(N, ALPHA, POWER) << endl;
  return 0;
}
