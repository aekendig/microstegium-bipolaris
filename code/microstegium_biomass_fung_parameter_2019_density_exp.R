##### info ####

# file: microstegium_biomass_fung_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_biomass_fung_model_2019_density_exp.rda")

# extract posterior distributions
mvBioFuD2Samps <- posterior_samples(mvBioFuD2Mod2)

# sample parameters
b_A_df <- mvBioFuD2Samps[sample(nrow(mvBioFuD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)
# remove the direct effect of fungicide by not including b_logv_treatmentfungicide 


#### biomass function ####

b_A_fung_fun <- function(disease, iter) {
  
  # calculate survival
  b_A_expr <- ifelse(disease == 1, 
                     b_A_df$int_wat[iter],
                     b_A_df$int_fun[iter])
  
  b_A <- exp(b_A_expr)
  
  return(b_A)
}

