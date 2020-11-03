##### info ####

# file: elymus_seedling_seed_production_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_seed_model_2019_density_exp.rda")

# extract posterior distributions
evSSeedD2Samps <- posterior_samples(evSSeedD2Mod2)

# sample parameters
Y_S_dens <- evSSeedD2Samps[sample(nrow(evSSeedD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_maxS_Intercept,
         int_fun = int_wat + b_maxS_treatmentfungicide,
         mv_dens_wat = b_gammaA_treatmentcontrol,
         mv_dens_fun = b_gammaA_treatmentfungicide,
         evS_dens_wat = b_gammaS_treatmentcontrol,
         evS_dens_fun = b_gammaS_treatmentfungicide,
         evA_dens_wat = b_gammaP_treatmentcontrol,
         evA_dens_fun = b_gammaP_treatmentfungicide)


#### biomass function ####

Y_S_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # calculate survival
  Y_S <- ifelse(disease == 1, 
                     Y_S_dens$int_wat[iter] / (1 + Y_S_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + Y_S_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + Y_S_dens$evA_dens_wat[iter] * P_dens),
                     Y_S_dens$int_fun[iter] / (1 + Y_S_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + Y_S_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + Y_S_dens$evA_dens_fun[iter] * P_dens))
  
  return(Y_S)
}

