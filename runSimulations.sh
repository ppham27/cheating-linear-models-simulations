#!/usr/bin/env sh

OUT_DIR=analysis/output_data

mkdir -p $OUT_DIR

./independentSimulation 1000 25 > $OUT_DIR/50_independent_beta.tsv
./independentSimulation 1000 50 > $OUT_DIR/100_independent_beta.tsv
./independentSimulation 1000 100 > $OUT_DIR/200_independent_beta.tsv
./independentSimulation 1000 200 > $OUT_DIR/400_independent_beta.tsv
# ./independentSimulation 1000 400 > $OUT_DIR/800_independent_beta.tsv

# ./powerSimulation 1000 25 0.1 > $OUT_DIR/50_01_R_beta.tsv
# ./powerSimulation 1000 50 0.1 > $OUT_DIR/100_01_R_beta.tsv
# ./powerSimulation 1000 100 0.1 > $OUT_DIR/200_01_R_beta.tsv
# ./powerSimulation 1000 200 0.1 > $OUT_DIR/400_01_R_beta.tsv
# ./powerSimulation 1000 400 0.1 > $OUT_DIR/800_01_R_beta.tsv

# ./powerSimulation 1000 25 0.3 > $OUT_DIR/50_03_R_beta.tsv
# ./powerSimulation 1000 50 0.3 > $OUT_DIR/100_03_R_beta.tsv
# ./powerSimulation 1000 100 0.3 > $OUT_DIR/200_03_R_beta.tsv
# ./powerSimulation 1000 200 0.3 > $OUT_DIR/400_03_R_beta.tsv
# ./powerSimulation 1000 400 0.3 > $OUT_DIR/800_03_R_beta.tsv

# ./powerSimulation 1000 25 0.5 > $OUT_DIR/50_05_R_beta.tsv
# ./powerSimulation 1000 50 0.5 > $OUT_DIR/100_05_R_beta.tsv
# ./powerSimulation 1000 100 0.5 > $OUT_DIR/200_05_R_beta.tsv
# ./powerSimulation 1000 200 0.5 > $OUT_DIR/400_05_R_beta.tsv
# ./powerSimulation 1000 400 0.5 > $OUT_DIR/800_05_R_beta.tsv

./powerSimulation 1000 25 0.7 > $OUT_DIR/50_07_R_beta.tsv
# ./powerSimulation 1000 50 0.7 > $OUT_DIR/100_07_R_beta.tsv
# ./powerSimulation 1000 100 0.7 > $OUT_DIR/200_07_R_beta.tsv
# ./powerSimulation 1000 200 0.7 > $OUT_DIR/400_07_R_beta.tsv
# ./powerSimulation 1000 400 0.7 > $OUT_DIR/800_07_R_beta.tsv

