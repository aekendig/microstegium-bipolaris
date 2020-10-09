##### info ####

# file: microstegium_seeds_per_biomass_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/8/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_seeds_per_biomass_model_2019_density_exp.rda")

# extract posterior distributions
mvSeedBioD2Samps <- posterior_samples(mvSeedBioD2Mod2)

# sample parameters
Y_A_int_wat <- sample(mvSeedBioD2Samps$b_Intercept, size = n_samps, replace = T)
Y_A_int_fun <- Y_A_int_wat + sample(mvSeedBioD2Samps$b_fungicide, size = n_samps, replace = T)

Y_A_bio_wat <- sample(mvSeedBioD2Samps$b_log_bio.g, size = n_samps, replace = T)
Y_A_bio_fun <- Y_A_bio_wat + sample(mvSeedBioD2Samps$"b_fungicide:log_bio.g", size = n_samps, replace = T)


#### seed function ####

Y_A_fun <- function(disease, biomass, iter) {
  
  # log-transform biomass
  log_bio <- log(biomass)
  
  # calculate seeds produced
  Y_A_lin_expr <- ifelse(disease == 1, 
                          Y_A_int_wat[iter] + Y_A_bio_wat[iter] * log_bio,
                          Y_A_int_fun[iter] + Y_A_bio_fun[iter] * log_bio)
  
  Y_A <- exp(Y_A_lin_expr )
  
  return(Y_A)
}

