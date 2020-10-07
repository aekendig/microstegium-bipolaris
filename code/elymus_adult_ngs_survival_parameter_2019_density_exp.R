##### info ####

# file: elymus_adult_ngs_survival_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/6/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_ngs_survival_model_2019_density_exp.rda")

# extract posterior distributions
evANgsSurvD1Samps <- posterior_samples(evANgsSurvD1Mod2)

# sample parameters
W_P_int <- sample(evANgsSurvD1Samps$b_Intercept, size = n_samps, replace = T)


#### survival function

W_P_fun <- function(iter) {
  
  # calculate survival
  W_P <- exp(W_P_int[iter])/(1 + exp(W_P_int[iter]))
  
  return(W_P)
}

