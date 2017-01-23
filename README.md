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

- `2N`: `N` is the number of observations in the treatment and the control group, so `2N` is the number of total observations.
- `p.value`: A linear regression is done, `y = beta_x*x + epsilon` or `y = beta_x*x + beta_0 + epsilon` depending on if an intercept is used or not. Given a null hypothesis of `beta_x = 0`, this is the probability of getting an estimate of `beta_x` that is at least as large in absolute value as the estimate that we found.
- `beta`: The estimate for `beta_x` from doing the linear regression.
- `set.size`: The number of independent covariates generated before some subset of them could be used to fake statistical significance.
- `subset.size`: The number of covariates used from that set to fake statistical significance.
- `balance.p.value`: This is only generated in the case that that the subset size is 1. Call the one covariate in the subset `z`. We test if the covariate is evenly distributed among the treatment and control. This is modeled as a hypergeometric distribution. The total population size is `2N`. We make `N` draws from the population and assign them to the treatment. Assume `z` has `K` 1s. If there was perfect balance, both the treatment and control would have `K/2` 1s. Let `k` be the actual number of 1s in the treatment group, and let `k` be sampled from a random variable `X`. This number is the probability of assigning the treatment in a way that is just as unbalance or more unbalanced, that is, `Prob(X <= K/2 - |k - K/2|) + Prob(X >= K/2 + |k - K/2|)`.
- `fake.p.value`: The the probability of getting an estimate of `beta_x` that is at least as large in absolute value as the estimate that we found when we include the subset of independent covariates that results in statistical significance.
- `power`: This is only generated in power mode. This is the probability of being able to detect the effect of the treatment.


## Compiling

This should be as simple as running the command

    make

## Generating the Data

`runSimulations.sh` and `runSimulationsWithIntercept.sh` provide examples of how to run the program. As noted in the name of script `runSimulationsWithIntercept.sh` specifies that we run the simulation an intercept term, which is just an additional vector of 1s as a predictor.

The first argument is the total number of simulations to be run. The second argument is `N`. For independent simulations, there is an optional third argument. If you pass `intercept` as the third argument, the linear regression will be run with an intercept term. For power simulations, the third argument is required. It should be a number between 0.05 and 1, which specifies the probability of detecting the effect of the treatment.

## Citations

- 2002. The Boost Graph Library: User Guide and Reference Manual. Addison-Wesley Longman Publishing Co., Inc., Boston, MA, USA.
- Conrad S. 2010. Armadillo: An open source C++ linear algebra library for fast prototyping and computationally intensive experiments. NICTA.

