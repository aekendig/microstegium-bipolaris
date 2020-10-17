##### info ####

# file: elymus_adult_germination_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/16/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_germination_model_2019_density_exp.rda")

# extract posterior distribution
evAGermD2Samps <- posterior_samples(evAGermD2Mod2)

# sample parameters
g0_P_int_wat <- sample(evAGermD2Samps$b_Intercept, size = n_samps, replace = T)
g0_P_int_fun <- g0_P_int_wat + sample(evAGermD2Samps$b_fungicide, size = n_samps, replace = T)

# from Lili's experiment: The presence of Microstegium litter reduced Elymus establishment by an average of 8.9%.
beta_SP_L <- 1 - 0.089


#### germination function ####

G_P_fun <- function(disease, litter_pres, iter) {
  
  # maximum germination
  g_P_expr <- ifelse(disease == 1, exp(g0_P_int_wat[iter])/(1 + exp(g0_P_int_wat[iter])), exp(g0_P_int_fun[iter])/(1 + exp(g0_P_int_fun[iter])))
  
  G_P <- ifelse(litter_pres == 1, g_P_expr * beta_SP_L, g_P_expr)
  
  return(G_P)
}