##### info ####

# file: microstegium_biomass_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_biomass_model_2019_density_exp.rda")

# extract posterior distributions
mvBioSamps <- posterior_samples(mvBioMod2)

# sample parameters
B_A_dens <- mvBioSamps[sample(nrow(mvBioSamps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_logv_Intercept,
         int_fun = int_wat + b_logv_fungicide,
         mv_dens_wat = b_alphaA_treatmentcontrol,
         mv_dens_fun = b_alphaA_treatmentfungicide,
         evS_dens_wat = b_alphaS_treatmentcontrol,
         evS_dens_fun = b_alphaS_treatmentfungicide,
         evA_dens_wat = b_alphaP_treatmentcontrol,
         evA_dens_fun = b_alphaP_treatmentfungicide)
# remove the direct effect of fungicide by not including b_logv_treatmentfungicide 


#### biomass function ####

B_A_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # calculate survival
  B_A_expr <- ifelse(disease == 1, 
                     B_A_dens$int_wat[iter] - log(1 + B_A_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + B_A_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + B_A_dens$evA_dens_wat[iter] * P_dens),
                     B_A_dens$int_fun[iter] - log(1 + B_A_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + B_A_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + B_A_dens$evA_dens_fun[iter] * P_dens))
  
  B_A <- exp(B_A_expr)
  
  return(B_A)
}

