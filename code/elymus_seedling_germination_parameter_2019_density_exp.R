##### info ####

# file: elymus_seedling_germination_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/16/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_germination_model_2019_density_exp.rda")

# extract posterior distribution
evSGermD2Samps <- posterior_samples(evSGermD2Mod2)

# sample parameters
g0_S_int_wat <- sample(evSGermD2Samps$b_Intercept, size = n_samps, replace = T)
g0_S_int_fun <- g0_S_int_wat + sample(evSGermD2Samps$b_fungicide, size = n_samps, replace = T)

# from Lili's experiment: The presence of Microstegium litter reduced Elymus establishment by an average of 8.9%.
beta_SP_L <- 1 - 0.089


#### germination function ####

G_S_fun <- function(disease, litter_pres, iter) {
  
  # maximum germination
  g_S_expr <- ifelse(disease == 1, exp(g0_S_int_wat[iter])/(1 + exp(g0_S_int_wat[iter])), exp(g0_S_int_fun[iter])/(1 + exp(g0_S_int_fun[iter])))
  
  G_S <- ifelse(litter_pres == 1, g_S_expr * beta_SP_L, g_S_expr)
  
  return(G_S)
}