##### info ####

# file: elymus_seedling_establishment_bh_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/4/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_gs_survival_bh_model_2018_2019_density_exp.rda")
load("output/elymus_litter_establishment_bh_model_2018_greenhouse_exp.rda")
load("output/microstegium_litter_establishment_bh_model_2018_litter_exp.rda")

# extract posterior distributions
evSGsSurvBhSamps <- posterior_samples(evSGsSurvBhMod2)
evLitEstBhSamps <- posterior_samples(evLitEstBhMod2)
evEstL1BhSamps <- posterior_samples(mvEstL1BhMod2)

# sample parameters
E_S_df <- evSGsSurvBhSamps[sample(nrow(evSGsSurvBhSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_year2" = "b_fungicide:yearfyear2") %>%
  transmute(int_wat_y1 = b_Intercept,
            int_fun_y1 = int_wat_y1 + b_fungicide,
            int_wat_y2 = int_wat_y1 + b_yearfyear2,
            int_fun_y2 = int_fun_y1 + b_yearfyear2 + b_fungicide_year2) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")

E_S_beta_wat <- sample(evLitEstBhSamps$b_betaL_Intercept, size = n_samps, replace = T)/10000

E_S_beta_fun <- evEstL1BhSamps[sample(nrow(evEstL1BhSamps), size = n_samps, replace = T), ] %>%
  mutate(sterilization_eff = b_betaL_sterilizedFsterilized/b_betaL_sterilizedFlive,
         live_litter = E_S_beta_wat * 10000,
         beta = (live_litter * sterilization_eff)/10000)


#### survival function

E_S_bh_fun <- function(disease, litter, iter) {
  
  # max survival
  E_S_lin_max = ifelse(disease == 1, 
                       E_S_df$int_wat[iter], 
                       E_S_df$int_fun[iter])
  E_S_max = exp(E_S_lin_max) / (1 + exp(E_S_lin_max))
  
  # add litter effect
  E_S_lit = ifelse(disease == 1, 
                   1 + E_S_beta_wat[iter] * litter, 
                   1 + E_S_beta_fun$beta[iter] * litter)
  
  E_S <- E_S_max / E_S_lit
  
  return(E_S)
}
