#!/usr/bin/env sh

./independentSimulation 1000 25 intercept > ../50_independent_beta_intercept.tsv
./independentSimulation 1000 50 intercept > ../100_independent_beta_intercept.tsv
./independentSimulation 1000 100 intercept > ../200_independent_beta_intercept.tsv
./independentSimulation 1000 200 intercept > ../400_independent_beta_intercept.tsv
# ./independentSimulation 1000 400 intercept > ../800_independent_beta_intercept.tsv
