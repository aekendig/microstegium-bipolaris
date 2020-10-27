##### info ####

# file: microstegium_seed_production_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_seed_model_2019_density_exp.rda")

# extract posterior distributions
mvSeedD2Samps <- posterior_samples(mvSeedD2Mod2)

# sample parameters
Y_A_dens <- mvSeedD2Samps[sample(nrow(mvSeedD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_maxS_Intercept,
         int_fun = int_wat + b_maxS_treatmentfungicide,
         mv_dens_wat = b_gammaA_treatmentcontrol,
         mv_dens_fun = b_gammaA_treatmentfungicide,
         evS_dens_wat = b_gammaS_treatmentcontrol,
         evS_dens_fun = b_gammaS_treatmentfungicide,
         evA_dens_wat = b_gammaP_treatmentcontrol,
         evA_dens_fun = b_gammaP_treatmentfungicide)


#### biomass function ####

Y_A_fun <- function(disease, A_dens, S_dens, P_dens, iter) {
  
  # calculate survival
  Y_A <- ifelse(disease == 1, 
                     Y_A_dens$int_wat[iter] / (1 + Y_A_dens$mv_dens_wat[iter] * A_dens + Y_A_dens$evS_dens_wat[iter] * S_dens + Y_A_dens$evA_dens_wat[iter] * P_dens),
                     Y_A_dens$int_fun[iter] / (1 + Y_A_dens$mv_dens_fun[iter] * A_dens + Y_A_dens$evS_dens_fun[iter] * S_dens + Y_A_dens$evA_dens_fun[iter] * P_dens))
  
  return(Y_A)
}

