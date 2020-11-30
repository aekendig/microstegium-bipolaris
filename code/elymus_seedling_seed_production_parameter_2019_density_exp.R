##### info ####

# file: elymus_seedling_seed_production_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 11/12/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_seed_fung_model_2019_density_exp.rda")
load("output/elymus_seedling_seed_model_2019_density_exp.rda")

# extract posterior distributions
evSSeedFuD2Samps <- posterior_samples(evSSeedFuD2Mod2)
evSSeedD2Samps <- posterior_samples(evSSeedD2Mod2)

# sample parameters
Y_S_fun_eff <- evSSeedFuD2Samps[sample(nrow(evSSeedFuD2Samps), size = n_samps, replace = T), ] %>%
  transmute(fun_eff = (b_Intercept + b_fungicide) / b_Intercept)

Y_S_dens <- evSSeedD2Samps[sample(nrow(evSSeedD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_logS_Intercept,
         mv_dens_wat = b_alphaA_treatmentcontrol,
         mv_dens_fun = b_alphaA_treatmentfungicide,
         evS_dens_wat = b_alphaS_treatmentcontrol,
         evS_dens_fun = b_alphaS_treatmentfungicide,
         evA_dens_wat = b_alphaP_treatmentcontrol,
         evA_dens_fun = b_alphaP_treatmentfungicide)


#### biomass function ####

Y_S_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # max seed production
  Y_S_max <- ifelse(disease == 1, 
                    Y_S_dens$int_wat[iter],
                    Y_S_dens$int_wat[iter] * Y_S_fun_eff$fun_eff[iter])
  
  # calculate seed production
  Y_S <- ifelse(disease == 1, 
                     Y_S_max - log(1 + Y_S_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + Y_S_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + Y_S_dens$evA_dens_wat[iter] * P_dens),
                     Y_S_max - log(1 + Y_S_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + Y_S_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + Y_S_dens$evA_dens_fun[iter] * P_dens))
  
  # correct negative values
  Y_S_out = ifelse(exp(Y_S) - 1 >= 0, exp(Y_S) - 1, 0)
  
  return(Y_S_out)
}

