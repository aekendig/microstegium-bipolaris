##### info ####

# file: microstegium_seed_production_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 11/12/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_seed_fung_model_2019_density_exp.rda")
load("output/microstegium_seed_model_2019_density_exp.rda")

# extract posterior distributions
mvSeedFuD2Samps <- posterior_samples(mvSeedFuD2Mod2)
mvSeedD2Samps <- posterior_samples(mvSeedD2Mod2)

# sample parameters
Y_A_fun_eff <- mvSeedFuD2Samps[sample(nrow(mvSeedFuD2Samps), size = n_samps, replace = T), ] %>%
  transmute(fun_eff = (b_Intercept + b_fungicide) / b_Intercept)

Y_A_dens <- mvSeedD2Samps[sample(nrow(mvSeedD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_logS_Intercept,
         mv_dens_wat = b_alphaA_treatmentcontrol,
         mv_dens_fun = b_alphaA_treatmentfungicide,
         evS_dens_wat = b_alphaS_treatmentcontrol,
         evS_dens_fun = b_alphaS_treatmentfungicide,
         evA_dens_wat = b_alphaP_treatmentcontrol,
         evA_dens_fun = b_alphaP_treatmentfungicide)


#### biomass function ####

Y_A_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # max seed production
  Y_A_max <- ifelse(disease == 1, 
                    Y_A_dens$int_wat[iter],
                    Y_A_dens$int_wat[iter] * Y_A_fun_eff$fun_eff[iter])
  
  # calculate seed production
  Y_A <- ifelse(disease == 1, 
                     Y_A_max - log(1 + Y_A_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + Y_A_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + Y_A_dens$evA_dens_wat[iter] * P_dens),
                     Y_A_max - log(1 + Y_A_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + Y_A_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + Y_A_dens$evA_dens_fun[iter] * P_dens))
  
  # correct negative values
  Y_A_out = ifelse(exp(Y_A) - 1 >= 0, exp(Y_A) - 1, 0)
  
  return(Y_A_out)
}

