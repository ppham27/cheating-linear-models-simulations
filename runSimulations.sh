#!/usr/bin/env sh

./independentSimulation 1000 25 > ../50_independent_beta.tsv
./independentSimulation 1000 50 > ../100_independent_beta.tsv
./independentSimulation 1000 100 > ../200_independent_beta.tsv
./independentSimulation 1000 200 > ../400_independent_beta.tsv
./independentSimulation 1000 400 > ../800_independent_beta.tsv

./powerSimulation 1000 25 0.1 > ../50_01_R_beta.tsv
./powerSimulation 1000 50 0.1 > ../100_01_R_beta.tsv
./powerSimulation 1000 100 0.1 > ../200_01_R_beta.tsv
./powerSimulation 1000 200 0.1 > ../400_01_R_beta.tsv
./powerSimulation 1000 400 0.1 > ../800_01_R_beta.tsv

./powerSimulation 1000 25 0.3 > ../50_03_R_beta.tsv
./powerSimulation 1000 50 0.3 > ../100_03_R_beta.tsv
./powerSimulation 1000 100 0.3 > ../200_03_R_beta.tsv
./powerSimulation 1000 200 0.3 > ../400_03_R_beta.tsv
./powerSimulation 1000 400 0.3 > ../800_03_R_beta.tsv

./powerSimulation 1000 25 0.5 > ../50_05_R_beta.tsv
./powerSimulation 1000 50 0.5 > ../100_05_R_beta.tsv
./powerSimulation 1000 100 0.5 > ../200_05_R_beta.tsv
./powerSimulation 1000 200 0.5 > ../400_05_R_beta.tsv
./powerSimulation 1000 400 0.5 > ../800_05_R_beta.tsv

./powerSimulation 1000 25 0.7 > ../50_07_R_beta.tsv
./powerSimulation 1000 50 0.7 > ../100_07_R_beta.tsv
./powerSimulation 1000 100 0.7 > ../200_07_R_beta.tsv
./powerSimulation 1000 200 0.7 > ../400_07_R_beta.tsv
./powerSimulation 1000 400 0.7 > ../800_07_R_beta.tsv

