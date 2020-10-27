##### info ####

# file: elymus_germination_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_germination_model_2019_density_exp.rda")

# extract posterior distribution
evGermD2Samps <- posterior_samples(evGermD2Mod2)

# sample parameters
g_S_df <- evGermD2Samps[sample(nrow(evGermD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)



#### germination function ####

g_S_fun <- function(disease, iter) {
  
  g_S_expr <- ifelse(disease == 1, g_S_df$int_wat[iter], g_S_df$int_fun[iter])
  
  g_S <- exp(g_S_expr) / (1 + exp(g_S_expr))
  
  return(g_S)
}
