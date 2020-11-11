##### info ####

# file: elymus_adult_seed_production_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_seed_fung_model_2018_2019_density_exp.rda")
load("output/elymus_adult_seed_model_2018_2019_density_exp.rda")

# extract posterior distributions
evASeedFuSamps <- posterior_samples(evASeedFuMod2)
evASeedSamps <- posterior_samples(evASeedMod2)

# sample parameters
Y_P_fun_eff <- evASeedFuSamps[sample(nrow(evASeedFuSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_yearfyear2" = "b_fungicide:yearfyear2") %>%
  transmute(fun_eff_y1 = exp(b_Intercept + b_fungicide) / exp(b_Intercept),
            fun_eff_y2 = exp(b_Intercept + b_yearfyear2 + b_fungicide + b_fungicide_yearfyear2) / exp(b_Intercept + b_yearfyear2)) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")

Y_P_dens <- evASeedSamps[sample(nrow(evASeedSamps), size = n_samps / 2, replace = T), ] %>%
  transmute(int_wat_y1 = b_maxS_Intercept,
            int_wat_y2 = int_wat_y1 + b_maxS_yearfyear2,
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
            evA_dens_fun_y2 = b_gammaP_treatment_yearfungicide_year2) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")


#### biomass function ####

Y_P_fun <- function(disease, g.A, E.A, A_dens, g.S, E.S, S_dens, P_dens, iter) {
  
  # max seed production
  Y_P_max <- ifelse(disease == 1, 
                    Y_P_dens$int_wat[iter] * 10,
                    Y_P_dens$int_wat[iter] * 10 * Y_P_fun_eff$fun_eff[iter])
  
  # calculate seed production
  Y_P <- ifelse(disease == 1, 
                     Y_P_max / (1 + Y_P_dens$mv_dens_wat[iter] * g.A * E.A * A_dens + Y_P_dens$evS_dens_wat[iter] * g.S * E.S * S_dens + Y_P_dens$evA_dens_wat[iter] * P_dens),
                     Y_P_max / (1 + Y_P_dens$mv_dens_fun[iter] * g.A * E.A * A_dens + Y_P_dens$evS_dens_fun[iter] * g.S * E.S * S_dens + Y_P_dens$evA_dens_fun[iter] * P_dens))
  
  return(Y_P)
}

