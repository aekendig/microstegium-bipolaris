##### info ####

# file: elymus_seeds_per_biomass_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/8/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seeds_produced_model_2019_density_exp.rda")
load("output/elymus_seeds_per_biomass_model_2019_density_exp.rda")

# extract posterior distributions
evSeedD2Samps <- posterior_samples(evSeedD2Mod2)
evSeedBioD2Samps <- posterior_samples(evSeedBioD2Mod2)

# sample parameters
Y0_SP_int_wat <- sample(evSeedD2Samps$b_Intercept, size = n_samps, replace = T)
Y0_SP_int_fun <- Y0_SP_int_wat + sample(evSeedD2Samps$b_fungicide, size = n_samps, replace = T)

Y0_SP_bio_wat <- sample(evSeedD2Samps$b_log_bio.g, size = n_samps, replace = T)
Y0_SP_bio_fun <- Y0_SP_bio_wat + sample(evSeedD2Samps$"b_fungicide:log_bio.g", size = n_samps, replace = T)

Y_SP_int_wat <- sample(evSeedBioD2Samps$b_Intercept, size = n_samps, replace = T)
Y_SP_int_fun <- Y_SP_int_wat + sample(evSeedBioD2Samps$b_fungicide, size = n_samps, replace = T)

Y_SP_bio_wat <- sample(evSeedBioD2Samps$b_log_bio.g, size = n_samps, replace = T)
Y_SP_bio_fun <- Y_SP_bio_wat + sample(evSeedBioD2Samps$"b_fungicide:log_bio.g", size = n_samps, replace = T)


#### seed function ####

Y_SP_fun <- function(disease, biomass, iter) {
  
  # log-transform biomass
  log_bio <- ifelse(biomass > 0, log(biomass), 0)
  
  # calculate probability of producing seeds
  Y0_SP_lin_expr <- ifelse(disease == 1,
                           Y0_SP_int_wat[iter] + Y0_SP_bio_wat[iter] * log_bio,
                           Y0_SP_int_fun[iter] + Y0_SP_bio_fun[iter] * log_bio)
  
  Y0_SP <- rbinom(1, 1, exp(Y0_SP_lin_expr)/(1 + exp(Y0_SP_lin_expr)))
  
  # calculate seeds produced
  Y_SP_lin_expr <- ifelse(disease == 1, 
                          Y_SP_int_wat[iter] + Y_SP_bio_wat[iter] * log_bio,
                          Y_SP_int_fun[iter] + Y_SP_bio_fun[iter] * log_bio)
  
  Y_SP <- Y0_SP * exp(Y_SP_lin_expr )
  
  return(Y_SP)
}

