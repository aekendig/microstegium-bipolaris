##### info ####

# file: elymus_adult_biomass_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/7/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_biomass_model_2019_density_exp.rda")

# extract posterior distributions
evABioSamps <- posterior_samples(evABioMod2)

# sample parameters
V_P_int_wat <- sample(evABioSamps$b_logv_Intercept, size = n_samps, replace = T)
# remove the direct effect of fungicide by not including b_logv_treatmentfungicide in the calculation of V_P_int_fun
V_P_int_fun <- V_P_int_wat + sample(evABioSamps$b_logv_fungicide, size = n_samps, replace = T)

V_P_mv_dens_wat <- sample(evABioSamps$b_alphaA_treatmentcontrol, size = n_samps, replace = T)
V_P_mv_dens_fun <- sample(evABioSamps$b_alphaA_treatmentfungicide, size = n_samps, replace = T)

V_P_evS_dens_wat <- sample(evABioSamps$b_alphaS_treatmentcontrol, size = n_samps, replace = T)
V_P_evS_dens_fun <- sample(evABioSamps$b_alphaS_treatmentfungicide, size = n_samps, replace = T)

V_P_evA_dens_wat <- sample(evABioSamps$b_alphaP_treatmentcontrol, size = n_samps, replace = T)
V_P_evA_dens_fun <- sample(evABioSamps$b_alphaP_treatmentfungicide, size = n_samps, replace = T)


#### survival function

V_P_fun <- function(disease, A_dens, S_dens, P_dens, iter) {
  
  # calculate survival
  V_P_expr <- ifelse(disease == 1, 
                     V_P_int_wat[iter] - log(1 + V_P_mv_dens_wat[iter] * A_dens + V_P_evS_dens_wat[iter] * S_dens + V_P_evA_dens_wat[iter] * P_dens),
                     V_P_int_fun[iter] - log(1 + V_P_mv_dens_fun[iter] * A_dens + V_P_evS_dens_fun[iter] * S_dens + V_P_evA_dens_fun[iter] * P_dens))
  
  V_P <- exp(V_P_expr)
  
  return(V_P)
}

