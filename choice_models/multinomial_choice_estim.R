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
require(mgcv)

# laod data
choices <- read.csv("./data/choices.csv")
choices <- choices[, 2:ncol(choices)]

# set parameters
params <- list(
  n_coef = 2,
  n_attr = 2
)

### ------------------------------------------------------------------------ ###
#                       Simulate choice probabilities
### ------------------------------------------------------------------------ ###
# Multinomial Probit: homogeneous coefficients ####
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

# Mixed Multinomial Logit Model: Heterogeneous coefficients ####

#   Probit coefficients
p_mmlogit_norm <- function(j, cmean, ccov, attributes, nsim = 10000, pc_id = NULL){
  # 
  # Description: Simulate the probability of objects j being chosen if coefficients
  #              on attributes (i.e. choice paramters) are normally distributed 
  #              in the population, with mean cmean and covariance ccov.
  #              Let na be the number of attributes of each alternative
  #              and nc be the number of alternatives to be chosen from.
  #
  # nsim (int):  number of samples drawn for simulation
  # ccov (mat):  [na x na] matrix
  # attributes (mat): [nc x na] matrix
  # j (int/vec): calculate choice probabilities for these objects. 
  #
  
  ζ <- rmvn(nsim, cmean, ccov) # sampled coefficients
  if (is.null(pc_id)) {
    sampled_probs <- calc_choice_prob_given_parameters_mmlogit(j, ζ, attributes) 
    if (is.vector(sampled_probs)) sampled_probs <- t(as.matrix(sampled_probs))
  }
  else {
    sampled_probs <- calc_choice_prob_given_parameters_mmlogit(0, ζ, attributes, by =pc_id) 
  }
  return(rowMeans(sampled_probs))
}

#   Latent class model
#       categorical coefficients: only for two attributes [deprecated ;)]
p_mmlogit_cat2  <- function(j, scat, fdist, attributes){

  # 
  # Description: Calculate the probability of objects j being chosen if coefficients
  #              on attributes (i.e. choice paramters) are drawn from scat (support) 
  #              categories in the population, with joint distribution fdist.
  #              Let na be the number of attributes of each alternative
  #              and nc be the number of alternatives to be chosen from.
  #
  # nsim (int):       number of samples drawn for simulation
  # fdist (mat):      [ncat x ncat] matrix, rows correspond to support of first attr.
  # j (int/vec):      calculate choice probabilities for these objects. 
  # attributes (mat): (observable) attributes of chooseable objects.
  #
  
  # expand possible coefficient combinations
  ζ <- tidyr::expand_grid(scat[,1], scat[,2]) %>% as.matrix 
  grid_probs <- calc_choice_prob_given_parameters_mmlogit(j, ζ, attributes)
  if (is.vector(grid_probs)) {
    grid_probs <- matrix(grid_probs, nrow = 1)
  }
  pweights <- utils::stack(data.frame(t(fdist)))$values
  sampled_probs <- purrr::modify(data.frame(t(grid_probs)), function(x) x * pweights)
  return(colSums(sampled_probs))
}

#       categorical coefficients: for any number of attributes
p_mmlogit_cat  <- function(j, scat_grid, attributes){
  
  # 
  # Description: Calculate the probability of objects j being chosen if coefficients
  #              on attributes (i.e. choice paramters) are drawn from scat (support) 
  #              categories in the population, with some joint distribution.
  #              Let na be the number of attributes of each alternative
  #              and nc be the number of alternatives to be chosen from.
  #
  # scat_grid (mat):  [ncat^na x (na + 1)] matrix, grid ranging over all combinations 
  #                   in the support of the joint distribution of coefficients. 
  #                   Rightmost column is the density evaluated at that point.
  # j (int/vec):      calculate choice probabilities for these objects. 
  # attributes (mat): (observable) attributes of chooseable objects.
  #
  
  #  possible coefficient combinations
  ζ <- scat2_grid[,-ncol(scat2_grid)]
  grid_probs <- calc_choice_prob_given_parameters_mmlogit(j, ζ, attributes) 
  choice_probs <- grid_probs %*% scat_grid[, ncol(scat_grid)]
  return(choice_probs)
}

# A helper function
calc_choice_prob_given_parameters_mmlogit <- function(j, coefs, attributes, by = NULL){
  # function is vectorized using df operations
  # function can be used to calculate choice probabilities in a single choice
  # set, in which case it requires input of vector j, the alternatives for which
  # to calculate such. It can also calculate probabilities across a range of
  # choice situations (hence one needs to supply a longer matrix of attributes), 
  # delineated by the "by" input variable.
  
  if (is.null(by)) {
    # j can be a single alternative or a vector of alternatives.
    # Two cases: coefficient vector is provided, or a matrix of coefficients with
    # each row corresponding to one sampled coefficient vector is provided.
    if (is.vector(coefs)) {
      vΩ <- exp(attributes %*% coefs)
      SvΩ <- vΩ %>% sum
      return(vΩ[j] / SvΩ )
    }
    else if (is.matrix(coefs)) {
      vΩ <- t(exp(attributes %*% t(coefs)))
      SvΩ <- vΩ %>% rowSums %>% as.matrix
      vΩ <- cbind(vΩ, SvΩ)
      return(t(vΩ[,j] / SvΩ[, ncol(SvΩ)]))
    }
    else {
      stop("Cannot work with input coefficients.")
    }
  }
  else {
      # if "by" is not null, j becomes redundant. Just output a matrix that calculates the 
      # choice probabilities in each choice set defined by "by" (rows), and across all 
      # draws of coefficients supplied (columns). This vectorized version of this function
      # is useful for faster computation of the likelihood of the entire sample, given a
      # (hyper-)parameter vector.
        
      # attributes <- as.matrix(data[, 5:6])
      # coefs <- ζ
      # by <- seq(1, 1500 + 2/3, 1/3) %>% floor()
    
      vΩ <- exp(attributes %*% t(coefs))
      SvΩ <- vΩ %>% rowsum(group = by) %>% as.matrix
      SvΩ <- SvΩ[rep(1:nrow(SvΩ), each=nrow(vΩ)/nrow(SvΩ)),]
      probs <- vΩ / SvΩ
      dimnames(probs)[[1]] <- paste0("cid_", by)
      dimnames(probs)[[2]] <- paste0("coefdraw_", 1:ncol(probs))
      
      return(probs)
  }

}



### ------------------------------------------------------------------------ ###
#                           Likelihood Functions
### ------------------------------------------------------------------------ ###
# note on identification:

LogLikelihoodFunction <- function(choices, p_id, pc_id, attr, model,...) {
  data <- data.frame(chosen = choices, 
                     p_id = p_id, 
                     pc_id = pc_id, 
                     attr)
  parameters <- list(...)
  # data <- data.frame(chosen = as.vector(choices[,"chosen"]), 
  #                    p_id = as.vector(choices[,"p_id"]), 
  #                    c_id = as.vector(choices[,"c_id"]), 
  #                    a_id = as.vector(choices[,"a_id"]), 
  #                    as.matrix(choices[,c("attr1", "attr2")]) )
  
  if (model == "p_mmlogit_norm") {
    choice_probabilities <- p_mmlogit_norm(0, 
                                           cmean = parameters$cmean, 
                                           ccov = parameters$ccov, 
                                           attributes = as.matrix(attr), 
                                           nsim = parameters$nsim, 
                                           pc_id = pc_id 
                                           )
  }
  
  # return simulated log likelihood
  (choices * log(choice_probabilities)) %>% sum %>% return
}

### ------------------------------ Sandbox --------------------------------- ###

## Not run:
# set some parameters and run the function for a choice of j ∈ {1, 2, 3}
# parameters for mprobit and mixed mlogit.
attributes  <- matrix(c(-.2, .12 , .3, .2, -.1, -.6), ncol = 2)
coefs       <- c( 2  , 0.1)
cov0        <- matrix(c(  1,  .5,  .5,  
                  .5, 1  , 0.2, 
                  .5, 0.2,   1 ), ncol = 3)
# parameters for mixed logit only.
cmean       <- c(0.2,0)
ccov        <- matrix(c(1, .5, .5, 1), ncol = 2)

# parameters/input for categorical model.
#   two dimensional attributes
scat <- matrix(c(seq(-.5, .5, by = 0.2), seq(-.5, .5, by = 0.2)), ncol = 2)
fdist <- matrix(purrr::rdunif(nrow(scat)^2, 10, 0), nrow = nrow(scat))
fdist <- fdist/sum(fdist)
#   lets also do this with four dimensions
attributes2 <- cbind(attributes, attributes + matrix(rnorm(n=6, sd = 0.2), ncol = 2))
scat2 <- cbind(scat, scat)
#   expand the support on a grid. Variables are called V1, ..., Vk, k = ncol(scat2)
scat2_grid <- do.call(tidyr::expand_grid, scat2 %>% as.data.frame %>% as.list)
scat2_grid$f <- runif(nrow(scat2_grid))
scat2_grid$f <- scat2_grid$f/sum(scat2_grid$f)
scat2_grid <- as.matrix(scat2_grid)

# 1. Try out the probability simulator for the multinomial mixed logit model.
cat("\n---- 1 ----", sep = "\n")
cat("P[ choice = j ] = \n")
tic("simulate all choice probabilities in the mixed multinomial logit model.\nTimed")
p_mmlogit_norm(1:nrow(attributes), cmean, ccov, attributes, nsim = 1000000) %>% cat(" ", sep = "\n")
toc()

# 2. Try out the probability simulator for the multinomial probit model.
j <- 2 # index of chosen alternative
cat("\n---- 2a ----", sep = "\n")
cat("P[ choice = j ] = \n")
tic("simulate a choice probability in the multinomial probit model.\nTimed")
p_mprob(3, coefs, attributes, cov0, nsim = 10000) %>% cat(" ", sep = "\n")
toc()

#   set parameters as to simulate a fully arbitrary choice 
#   and run the function for a choice of j ∈ {1, ..., 5}
#   choice probability should be around .2
cat("\n---- 2b ----", sep = "\n")
tic("simulate a choice probability in the multinomial probit model, \nequal likelihood case.\nTimed")
cat("P[ choice = j ] = \n")
p_mprob(1, coefs = c(1), 
        attributes = matrix(rep(0, 5), ncol = 1), 
        cov0 = diag(5), 
        nsim = 10000) %>% cat(" ", sep = "\n")
toc()

# 3. Try out calculating choice probabilities for the 
#    categorical mixed logit mulitnomial choice model

#   example of p_mmlogit_cat
cat("\n---- 3a ----", sep = "\n")
cat("P[ choice = j ] = \n")
tic("calculate a vector of choice probabiliites in the mutinomial \ncategorical mixed logit model, given a distribution of population coefficients. \nFour attributes case.\nTimed")
p_mmlogit_cat(1:nrow(attributes2), scat2_grid, attributes2) %>% cat(" ", sep = "\n")
toc()

#   example of p_mmlogit_cat2
cat("\n---- 3b ----", sep = "\n")
cat("P[ choice = j ] = \n")
tic("calculate a vector of choice probabiliites in the mutinomial \ncategorical mixed logit model, given a distribution of population coefficients. \nTwo attributes case.\nTimed")
p_mmlogit_cat2(1:nrow(attributes), scat = scat, fdist = fdist, attributes = attributes) %>% cat(" ", sep = "\n")
toc()

# 4. Calculate the simulated log likelihood for a given set of parameters in the
#    multinomial mixed logit model with normal mixing
cat("\n---- 4 ----", sep = "\n")
cat("L[ choice_1 = j_1, ...,   choice_1 = j_1 ; cmean, ccov] = \n")
tic("Calculate simulated log-likelihood for multinomial mixed logit model.")
LogLikelihoodFunction(
  choices = choices$chosen,
  p_id = choices$p_id,
  pc_id = as.matrix(choices$c_id)*10000+as.matrix(choices$p_id),
  attr = choices[,c(7,8)],
  model = "p_mmlogit_norm", 
  cmean = cmean, 
  ccov = ccov,
  nsim = 500) %>% cat(" ", sep = "\n")
toc()


purrr::map_dbl(1:50, function(x) LogLikelihoodFunction(
  choices = choices$chosen,
  p_id = choices$p_id,
  pc_id = as.matrix(choices$c_id)*10000+as.matrix(choices$p_id),
  attr = choices[,c(7,8)],
  model = "p_mmlogit_norm", 
  cmean = c(0.05,0.02), 
  ccov = matrix(c(1,.4271759,0.4271759,1), nrow = 2),
  nsim = 500) ) %>% sd

purrr::map_dbl(1:50, function(x) LogLikelihoodFunction(
  choices = choices$chosen,
  p_id = choices$p_id,
  pc_id = as.matrix(choices$c_id)*10000+as.matrix(choices$p_id),
  attr = choices[,c(7,8)],
  model = "p_mmlogit_norm", 
  cmean = c(0.05,0.02), 
  ccov = matrix(c(1,.4271759,0.4271759,1), nrow = 2),
  nsim = 10000) ) %>% sd

# To do:
# identification: make sure to normalize appropriately 
# (each i should e.g. have unit variance errors for choice 1)
# why does the likelihood vary so much across draws???

# 5. maximize the likelihood for a multinomial mixed logit with normal mixing
#    and two attributes

s11 <- 1
s22 <- seq(0.1, 4, 0.1)
s12 = list()
for (i in seq_along(s22)) {
  v <- seq(0, sqrt(s11*s22[i]), 0.1) %>% as.vector
  v <- c(-v[length(v):2], v)
  s12 <- list(i = v) %>% append(s12)
}


σ11 = s11
σ22 = s22
σ12 = σ21 <= sqrt(σ11*σ22)

## End(Not run)