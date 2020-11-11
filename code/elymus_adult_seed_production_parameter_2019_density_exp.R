##### info ####

# file: elymus_adult_seed_production_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_seed_fung_model_2019_density_exp.rda")
load("output/elymus_adult_seed_model_2019_density_exp.rda")

# extract posterior distributions
evASeedFuD2Samps <- posterior_samples(evASeedFuD2Mod2)
evASeedD2Samps <- posterior_samples(evASeedD2Mod2)

# sample parameters
Y_P_fun_eff <- evASeedFuD2Samps[sample(nrow(evASeedFuD2Samps), size = n_samps, replace = T), ] %>%
  transmute(fun_eff = exp(b_Intercept + b_fungicide) / exp(b_Intercept))

Y_P_dens <- evASeedD2Samps[sample(nrow(evASeedD2Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_maxS_Intercept,
         mv_dens_wat = b_gammaA_treatmentcontrol,
         mv_dens_fun = b_gammaA_treatmentfungicide,
         evS_dens_wat = b_gammaS_treatmentcontrol,
         evS_dens_fun = b_gammaS_treatmentfungicide,
         evA_dens_wat = b_gammaP_treatmentcontrol,
         evA_dens_fun = b_gammaP_treatmentfungicide)


#### biomass function ####

Y_P_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # max seed production
  Y_P_max <- ifelse(disease == 1, 
                    Y_P_dens$int_wat[iter],
                    Y_P_dens$int_wat[iter] * Y_P_fun_eff$fun_eff[iter])
  
  # calculate seed production
  Y_P <- ifelse(disease == 1, 
                     Y_P_max / (1 + Y_P_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + Y_P_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + Y_P_dens$evA_dens_wat[iter] * P_dens),
                     Y_P_max / (1 + Y_P_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + Y_P_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + Y_P_dens$evA_dens_fun[iter] * P_dens))
  
  return(Y_P)
}

