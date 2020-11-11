##### info ####

# file: microstegium_seed_production_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_seed_fung_model_2018_2019_density_exp.rda")
load("output/microstegium_seed_model_2018_2019_density_exp.rda")

# extract posterior distributions
mvSeedFuSamps <- posterior_samples(mvSeedFuMod2)
mvSeedSamps <- posterior_samples(mvSeedMod2)

# sample parameters
Y_A_fun_eff <- mvSeedFuSamps[sample(nrow(mvSeedFuSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_yearfyear2" = "b_fungicide:yearfyear2") %>%
  transmute(fun_eff_y1 = exp(b_Intercept + b_fungicide) / exp(b_Intercept),
            fun_eff_y2 = exp(b_Intercept + b_yearfyear2 + b_fungicide + b_fungicide_yearfyear2) / exp(b_Intercept + b_yearfyear2)) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")

# use water values from year 2 only (year 1 were not individual plants), and apply fungicide effect from year 1 to these
Y_A_dens1 <- mvSeedSamps[sample(nrow(mvSeedSamps), size = n_samps, replace = T), ] %>%
  transmute(int_wat_y2 = b_maxS_Intercept + b_maxS_yearfyear2,
            mv_dens_wat_y1 = b_gammaA_treatment_yearcontrol_year1,
            mv_dens_fun_y1 = b_gammaA_treatment_yearfungicide_year1,
            evS_dens_wat_y1 = b_gammaS_treatment_yearcontrol_year1,
            evS_dens_fun_y1 = b_gammaS_treatment_yearfungicide_year1,
            evA_dens_wat_y1 = b_gammaP_treatment_yearcontrol_year1,
            evA_dens_fun_y1 = b_gammaP_treatment_yearfungicide_year1,
            mv_dens_wat_y2 = b_gammaA_treatment_yearcontrol_year2,
            mv_dens_fun_y2 = b_gammaA_treatment_yearfungicide_year2,
            evS_dens_wat_y2 = b_gammaS_treatment_yearcontrol_year2,
            evS_dens_fun_y2 = b_gammaS_treatment_yearfungicide_year2,
            evA_dens_wat_y2 = b_gammaP_treatment_yearcontrol_year2,
            evA_dens_fun_y2 = b_gammaP_treatment_yearfungicide_year2)

# select the first half of water values to go with the fungicide effect from year 2 and the second half to go with the fungicide effect from year 1
Y_A_dens <- Y_A_dens1[1:(n_samps / 2), ] %>%
  mutate(int_wat_y1 = Y_A_dens1$int_wat_y2[(n_samps/2 + 1):n_samps]) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")



#### biomass function ####

Y_A_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # max seed production
  Y_A_max <- ifelse(disease == 1, 
                    Y_A_dens$int_wat[iter] * 100,
                    Y_A_dens$int_wat[iter] * 100 * Y_A_fun_eff$fun_eff[iter])
  
  # calculate seed production
  Y_A <- ifelse(disease == 1, 
                     Y_A_max / (1 + Y_A_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + Y_A_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + Y_A_dens$evA_dens_wat[iter] * P_dens),
                     Y_A_max / (1 + Y_A_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + Y_A_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + Y_A_dens$evA_dens_fun[iter] * P_dens))
  
  return(Y_A)
}

