##### info ####

# file: elymus_seedling_biomass_fung_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_biomass_fung_model_2019_density_exp.rda")

# extract posterior distributions
evSBioFuSamps <- posterior_samples(evSBioFuMod2)

# sample parameters
b_S_df <- evSBioFuSamps[sample(nrow(evSBioFuSamps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)
# remove the direct effect of fungicide by not including b_logv_treatmentfungicide 


#### biomass function ####

b_S_fung_fun <- function(disease, iter) {
  
  # calculate survival
  b_S_expr <- ifelse(disease == 1, 
                     b_S_df$int_wat[iter],
                     b_S_df$int_fun[iter])
  
  b_S <- exp(b_S_expr)
  
  return(b_S)
}

