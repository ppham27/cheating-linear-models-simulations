#!/usr/bin/env sh

OUT_DIR=analysis/output_data

mkdir -p $OUT_DIR

./independentSimulation 3000 25 intercept > $OUT_DIR/50_independent_beta_intercept.tsv
./independentSimulation 3000 50 intercept > $OUT_DIR/100_independent_beta_intercept.tsv
./independentSimulation 3000 100 intercept > $OUT_DIR/200_independent_beta_intercept.tsv
./independentSimulation 3000 200 intercept > $OUT_DIR/400_independent_beta_intercept.tsv
# ./independentSimulation 1000 400 intercept > $OUT_DIR/800_independent_beta_intercept.tsv
