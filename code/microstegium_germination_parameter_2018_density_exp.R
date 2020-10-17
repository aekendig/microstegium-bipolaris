##### info ####

# file: microstegium_germination_2018_density_exp
# author: Amy Kendig
# date last edited: 10/16/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_germination_model_2018_density_exp.rda")

# extract posterior distribution
mvGermD1Samps <- posterior_samples(mvGermD1Mod2)

# sample parameters
g0_A_int_wat <- sample(mvGermD1Samps$b_Intercept, size = n_samps, replace = T)
g0_A_int_fun <- g0_A_int_wat + sample(mvGermD1Samps$b_fungicide, size = n_samps, replace = T)

# from Lili's experiment: The effect of litter on germination is weaker when Elymus is present, but converges with the estimate for when Elymus is absent when litter reaches ~2.5g. Using the model coefficients
beta_A_L <- -0.26 + 0.14


#### germination function ####

G_A_fun <- function(disease, litter, iter) {
  
  # maximum germination
  g_A_expr <- ifelse(disease == 1, g0_A_int_wat[iter] + beta_A_L * litter, g0_A_int_fun[iter] + beta_A_L * litter)
  
  G_A <- exp(g_A_expr) / (1 + exp(g_A_expr))
  
  return(G_A)
}