# 
# Project:      Estimating and simulating a random utility model
# Author:       F Bennhoff
# Description:  estimate a choice model on given data set
#

# clear all
rm(list = ls())

# load packages
# require tictoc to be installed
require(dplyr)
require(purrr)
require(tidyr)
require(utils)
require(EnvStats)
require(rlang)
require(tictoc)

# laod data
choices <- read.csv("./data/choices.csv")
choices <- choices[, 2:ncol(choices)]

# set parameters
params <- list(
  n_coef = 2,
  n_attr = 2
)


# Multinomial Probit: homogeneous coefficients
# Simulate choice probabilities
p_mprob <- function(j, coefs, attributes, cov0, nsim = 10000) {
  #
  # Description: Estimate the probability of j being chosen in a multinomial
  #              probit model of nc objects. Objects have na attributes and 
  #              individuals have (known) valuations for each attribute, given
  #              by a coefficient vector. Valuations are homogeneous.
  #
  # - coef (vector):        [1 x na] vector of slopes. Idividual specific, 
  #                         known, alternative-invariant.
  # - j (int):              alternative to be chosen 
  # - cov0 (matrix):        [nc x nc] matrix of error covariances
  # - attributes (matrix):  [nc x na] matrix of object specific characteristics
  # - nsim (int):           number of simulations to estimate probability
  #
  
  nc <- nrow(cov0) # number of choices in set (including j)
  
  # selector-matrix, S, used to create differences in noise:
  # define η_sj = ε_s - ε_j such that
  # η = Sε
  S = matrix(rep(0, nc^2), ncol = nc)
  S[, j] <- 1
  S <- diag(nc) - S
  S <- S[-j, ]
  
  # error variance of the η's
  cov_η = S %*% cov0 %*% t(S)
  
  # cholesky factorization: ΓΓ' = cov_η\
  # Γ is lower triangular
  Γ <- base::chol(cov_η) %>% t 
  
  V = attributes %*% coefs # observed value given parameters
  V = S %*% V # differences in observed between j and other choices
  bounds = V # these differences will be integration bounds! 
  
  # generate iid normally distributed variables using the probability transform
  # theorem, then simulate the probability that η = Γz <= bounds. Use 10e^3 reps.
  # this gives our choice probability given the parameter choice.
  #
  # Update: since we are assuming probit, the usual R routines are faster than
  # using the PIT.
  cp <- purrr::map_lgl(1:nsim, 
  function(x) return(all(Γ%*%rnorm(nc-1) <= bounds))
  ) %>% mean()
  
  return(cp)

}


## Not run:

# set some parameters and run the function for a choice of j ∈ {1, 2, 3}
attributes  <- matrix(c(-.2, 1 , .3, .2, -.1, -.6), ncol = 2)
coefs       <- c( 2  , 0.1)
cov0 <- matrix(c(  1,  .5,  .5,  
                  .5, 1  , 0.2, 
                  .5, 0.2,   1 ), ncol = 3)

j <- 2 # index of chosen alternative
tic("simulate a choice probability")
p_mprob(3, coefs, attributes, cov0, nsim = 1000000)
toc()

# set parameters as to simulate a fully arbitrary choice 
# and run the function for a choice of j ∈ {1, ..., 5}
# choice probability should be around .2
p_mprob(1, coefs = c(1), 
        attributes = matrix(rep(0, 5), ncol = 1), 
        cov0 = diag(5), 
        nsim = 1000000)


## End(Not run)