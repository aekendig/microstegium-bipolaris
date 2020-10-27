##### info ####

# file: microstegium_germination_2018_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_germination_model_2018_density_exp.rda")

# extract posterior distribution
mvGermD1Samps <- posterior_samples(mvGermD1Mod2)

# sample parameters
g_A_df <- mvGermD1Samps[sample(nrow(mvGermD1Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)


#### germination function ####

g_A_fun <- function(disease, iter) {
  
  g_A_expr <- ifelse(disease == 1, g_A_df$int_wat[iter], g_A_df$int_fun[iter])
  
  g_A <- exp(g_A_expr) / (1 + exp(g_A_expr))
  
  return(g_A)
}
