##### info ####

# file: elymus_adult_biomass_fung_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_biomass_fung_model_2019_density_exp.rda")

# extract posterior distributions
evABioFuSamps <- posterior_samples(evABioFuMod2)

# sample parameters
b_P_df <- evABioFuSamps[sample(nrow(evABioFuSamps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)
# remove the direct effect of fungicide by not including b_logv_treatmentfungicide 


#### biomass function ####

b_P_fung_fun <- function(disease, iter) {
  
  # calculate survival
  b_P_expr <- ifelse(disease == 1, 
                     b_P_df$int_wat[iter],
                     b_P_df$int_fun[iter])
  
  b_P <- exp(b_P_expr)
  
  return(b_P)
}

