##### info ####

# file: elymus_seedling_seed_production_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/12/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_seed_fung_model_2018_2019_density_exp.rda")
load("output/elymus_seedling_seed_model_2018_2019_density_exp.rda")

# extract posterior distributions
evSSeedFuSamps <- posterior_samples(evSSeedFuMod2)
evSSeedSamps <- posterior_samples(evSSeedMod2)

# sample parameters
Y_S_fun_eff <- evSSeedFuSamps[sample(nrow(evSSeedFuSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_yearfyear2" = "b_fungicide:yearfyear2") %>%
  transmute(fun_eff_y1 = (b_Intercept + b_fungicide) / b_Intercept,
            fun_eff_y2 = (b_Intercept + b_yearfyear2 + b_fungicide + b_fungicide_yearfyear2) / (b_Intercept + b_yearfyear2)) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")

Y_S_dens <- evSSeedSamps[sample(nrow(evSSeedSamps), size = n_samps / 2, replace = T), ] %>%
  transmute(int_wat_y1 = b_logS_Intercept,
            int_wat_y2 = int_wat_y1 + b_logS_yearfyear2,
            mv_dens_wat_y1 = b_alphaA_treatment_yearcontrol_year1,
            mv_dens_fun_y1 = b_alphaA_treatment_yearfungicide_year1,
            evS_dens_wat_y1 = b_alphaS_treatment_yearcontrol_year1,
            evS_dens_fun_y1 = b_alphaS_treatment_yearfungicide_year1,
            evA_dens_wat_y1 = b_alphaP_treatment_yearcontrol_year1,
            evA_dens_fun_y1 = b_alphaP_treatment_yearfungicide_year1,
            mv_dens_wat_y2 = b_alphaA_treatment_yearcontrol_year2,
            mv_dens_fun_y2 = b_alphaA_treatment_yearfungicide_year2,
            evS_dens_wat_y2 = b_alphaS_treatment_yearcontrol_year2,
            evS_dens_fun_y2 = b_alphaS_treatment_yearfungicide_year2,
            evA_dens_wat_y2 = b_alphaP_treatment_yearcontrol_year2,
            evA_dens_fun_y2 = b_alphaP_treatment_yearfungicide_year2) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")


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

