# 
# Project:      Estimating and simulating a random utility model
# Author:       F Bennhoff
# Description:  Generate a simulated dataset of choices
#
# Status:       Can generate choice sets for normal mixing distribtion of prefe-
#               rence parameters.
# To be done:   - Extend to allow the normal mean/variance to depend on 
#                 individual characteristics
#               - Extend to generate data based on the latent class model.

# clear all
rm(list = ls())
set.seed(2209)

# load packages
require(dplyr)
require(purrr)
require(tidyr)
require(utils)
require(EnvStats)
require(rlang)

# setting parameters to generate dataset 
params <- list(
 alternatives_n = 10,
 choiceset_n = 3,
 choices_per_i = 1,
 players_n = 100,
 players_coef_var = matrix(c(1, 1/2, 1/2, 1), ncol = 2),
 players_coef_mean = c(1,-1)
)

gen_choice_situations   <- function(params) {
  #
  # Description: Generate a dataset of choice situations for each individual.
  #              Objects in each choice situations were randomly sampled.
  #
  
  # make a dummy dataframe to fill in
  data = data.frame(tidyr::expand_grid(1:params$players_n, 1:params$choices_per_i, 
                                       1:params$choiceset_n))
  names(data) <- c("p_id", "c_id", "a_id")
  
  # for each choice situation and individual, generate a set of alternatives
  samples <- purrr::map(1:(params$players_n*params$choices_per_i), 
             ~ sample(1:params$alternatives_n, params$choiceset_n)) 
  names(samples) <- 1:length(samples)
  samples <- utils::stack(samples) 
  samples %>% head()
  
  data$a_id <- samples$values
  data %>% head(n = 20)
  
  
  # for each alternative, generate two attributes. These are observable, hence
  # there is no need to impose any correlation between attributes.
  
  attributes <- matrix(rnorm(params$alternatives_n*2), ncol = 2, 
         dimnames = list(alternative = 1:params$alternatives_n, 
                         attributes = c("attr1", "attr2")))
  
  attributes <- cbind(a_id = 1:params$alternatives_n, attributes)
  attributes <- attributes %>% as.data.frame()
  
  data <- dplyr::left_join(data, attributes, by = c("a_id"))
  
  # return and exit
  return(data)
}
gen_person_coefs        <- function(params) {
  #
  # for each player, generate correlated taste parameters
  #
  cat("covariance matrix of parameters:\n")
  print(params$players_coef_var)
  
  len_rand <- nrow(params$players_coef_var)
  ?? <- base::chol(params$players_coef_var) %>% t
  coefs = purrr::map(1:params$players_n,   
                     ~ t(?? %*% rnorm(len_rand) + params$players_coef_mean
                     )) %>%
    do.call(what = rbind.data.frame)
  names(coefs) = paste0("coef_" , as.character(1:len_rand))
  coefs$p_id <- 1:params$players_n
  coefs %>% return
  
}
join_choices_and_coefs  <- function(data, coefs) {
  #
  # Description: Merge a frame of choice situations and personal coefficients
  #
  
  data %>% dplyr::left_join(coefs, by = c("p_id"))
}
gen_choice_utility      <- function(data, stable = TRUE, noise = "EV1", ...) {
  #
  # Description: Generate the noise in the personal utility of each choice.
  #              Noise can be interpreted as unobserved attributes and their 
  #              person specific effect.
  #
  
  args <- list(...)
  data_0 = data

  # we want to make the errors a-id-p_id specific, but repeated choice 
  # situations involving the same objects should have the same error (bc same
  # unobservables). This is meant by the option "stable".
  
  n <- nrow(data)
  if (noise == "EV1") {
    data$noise <- revd(n, location = args$location, scale = args$scale)
  }
  else {
    return(data_0)
    stop("No valid random number generator for noise provided.")
  }
  
  if (stable == TRUE) {
    df <- data %>% group_by(p_id, a_id)
    noise_by_p_a <- df %>% group_by(p_id, a_id) %>% 
      summarise(noise = first(noise)) %>%
      data.frame
    
    data <- data %>%
      dplyr::select(-c("noise")) %>% 
      dplyr::left_join(noise_by_p_a, by = c("p_id", "a_id"))
  }
  
  data_cols      <- colnames(data)
  attribute_cols <- data_cols[data_cols %>% grepl(pattern = "attr*")]
  coef_cols      <- data_cols[data_cols %>% grepl(pattern = "coef*")]
  
  data <- make_U(data, attribute_cols, coef_cols)
  return(data)
  
}
make_U                  <- function(data, coefficients, attributes) {
  #
  #  description: helper function to put attributes and coefficient columns 
  #               into a new utility column.
  #
  
  if (length(coefficients) != length(attributes)) {
    stop("number of attributes does not match number of coefficients")
  }
  data$U <- 0
  for (i in 1:length(coefficients)) {
    data["U"] <- data["U"] + data[coefficients[i]] * data[attributes[i]] 
  }
  data["U"] <- data["U"] + data$noise
  return(data)
}
gen_choices             <- function(data) {
  #
  # Description: Given data with choice utilities, find the objects chosen.
  #
  
  data <- data %>% group_by(p_id, c_id) %>% 
    summarise(maxU = max(U), a_id = a_id, chosen = U == max(U)) %>%
    right_join(data, by = c("p_id", "c_id", "a_id")) %>%
    as.data.frame() %>%
    select(p_id, c_id, a_id, chosen, U, maxU, starts_with("attr"), starts_with("coef")) 
  return(data)
}

main <- function(outname = "df") {
  #
  # Description: Generate a data set with choice data.
  #
  
  data <- gen_choice_situations(params)
  coefs <- gen_person_coefs(params) 
  data <- join_choices_and_coefs(data, coefs)
  data <- gen_choice_utility(data = data, noise = "EV1", location = 0, scale = 1) 
  data <- gen_choices(data)

  assign(paste0(outname), data, envir = global_env(),
       inherits = FALSE, immediate = TRUE)
    
}

main()
write.csv(df, file = "./data/choices.csv")




