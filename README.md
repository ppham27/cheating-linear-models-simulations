# Cheating Linear Models Simulations

In here is the code for my master's thesis. It consists of 
various statistical simulations related to linear models. 
[Boost](http://www.boost.org/) is used to for the statistical distributions. 
[Armadillo](http://arma.sourceforge.net/) is used for the linear algebra.

Each simulation generates `2N` observations. `N` observations are assigned to the treatment group, and the remaining `N` observations are the control group. The treatment `x` is a vector of 0s and 1s, where 1 denotes the treatment group. The response variable `y` for each obervation is normally distributed with mean and variance depending on user parameters. A simple linear regression is run with `y` as the dependent variable and `x` as the independent variable.

There are two possible modes:

1. **Independent mode**: In this case, `y` comes from a standard normal distribution and is independent of the treatment.
2. **Power mode**: In this case, the user specifies the statistical power `0.05 < 1 - beta <  1`. `y` from the control group has a normal distribution with mean 0. `y` from the treatment group has a normal distribution with mean 1. The variance is chosen such that probability of detecting the effect of the treatment is `1 - beta`, where we run a t-test at significance level 0.05.

Then, random vectors of 0s and 1s are generated, that is, each of the `2N` components of the vector follow a Bernoulli distribution with parameter 0.5 and are completely independent of the response and treatment. These vectors are generated one at a time to serve as additional predictor variables. After each vector is generated, a dynamic programming algorithm tries to choose an approximately minimal set of these vectors to make the treatment appear statistically significant. When a such a subset is found, the simulation stops.

Each simulation returns the following data:

- `2N`: 
- `p.value`:
- `beta`:
- `set.size`:
- `subset.size`:
- `balance.p.value`:
- `fake.p.value`:


## Compiling

This should be as simple as running the command

    make

## Generating the Data

`runSimulations.sh` and `runSimulationsWithIntercept.sh` provide examples of how to run the program. As noted in the name of script `runSimulationsWithIntercept.sh` specifies that we run the simulation an intercept term, which is just an additional vector of 1s as a predictor.

## Citations

- 2002. The Boost Graph Library: User Guide and Reference Manual. Addison-Wesley Longman Publishing Co., Inc., Boston, MA, USA.
- Conrad S. 2010. Armadillo: An open source C++ linear algebra library for fast prototyping and computationally intensive experiments. NICTA.

