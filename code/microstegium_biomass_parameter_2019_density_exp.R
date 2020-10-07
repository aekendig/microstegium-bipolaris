##### info ####

# file: microstegium_biomass_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/7/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_biomass_model_2019_density_exp.rda")

# extract posterior distributions
mvBioSamps <- posterior_samples(mvBioMod2)

# sample parameters
V_A_int_wat <- sample(mvBioSamps$b_logv_Intercept, size = n_samps, replace = T)
# remove the direct effect of fungicide by not including b_logv_treatmentfungicide in the calculation of V_A_int_fun
V_A_int_fun <- V_A_int_wat + sample(mvBioSamps$b_logv_fungicide, size = n_samps, replace = T)

V_A_mv_dens_wat <- sample(mvBioSamps$b_alphaA_treatmentcontrol, size = n_samps, replace = T)
V_A_mv_dens_fun <- sample(mvBioSamps$b_alphaA_treatmentfungicide, size = n_samps, replace = T)

V_A_evS_dens_wat <- sample(mvBioSamps$b_alphaS_treatmentcontrol, size = n_samps, replace = T)
V_A_evS_dens_fun <- sample(mvBioSamps$b_alphaS_treatmentfungicide, size = n_samps, replace = T)

V_A_evA_dens_wat <- sample(mvBioSamps$b_alphaP_treatmentcontrol, size = n_samps, replace = T)
V_A_evA_dens_fun <- sample(mvBioSamps$b_alphaP_treatmentfungicide, size = n_samps, replace = T)


#### survival function

V_A_fun <- function(disease, A_dens, S_dens, P_dens, iter) {
  
  # calculate survival
  V_A_expr <- ifelse(disease == 1, 
                     V_A_int_wat[iter] - log(1 + V_A_mv_dens_wat[iter] * A_dens + V_A_evS_dens_wat[iter] * S_dens + V_A_evA_dens_wat[iter] * P_dens),
                     V_A_int_fun[iter] - log(1 + V_A_mv_dens_fun[iter] * A_dens + V_A_evS_dens_fun[iter] * S_dens + V_A_evA_dens_fun[iter] * P_dens))
  
  V_A <- exp(V_A_expr)
  
  return(V_A)
}

