##### info ####

# file: microstegium_establishment_bh_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/4/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_gs_survival_bh_model_2018_2019_density_exp.rda")
load("output/microstegium_litter_establishment_bh_model_2018_greenhouse_exp.rda")
load("output/microstegium_litter_establishment_bh_model_2018_litter_exp.rda")

# extract posterior distributions
mvGsSurvBhSamps <- posterior_samples(mvGsSurvBhMod4)
mvLitEstBhSamps <- posterior_samples(mvLitEstBhMod2)
mvEstL1BhSamps <- posterior_samples(mvEstL1BhMod2)

# sample parameters
E_A_df <- mvGsSurvBhSamps[sample(nrow(mvGsSurvBhSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_year2" = "b_fungicide:yearfyear2") %>%
  transmute(int_wat_y1 = b_Intercept,
         int_fun_y1 = int_wat_y1 + b_fungicide,
         int_wat_y2 = int_wat_y1 + b_yearfyear2,
         int_fun_y2 = int_fun_y1 + b_yearfyear2 + b_fungicide_year2) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")

E_A_beta_wat <- sample(mvLitEstBhSamps$b_betaL_Intercept, size = n_samps, replace = T)/10000

E_A_beta_fun <- mvEstL1BhSamps[sample(nrow(mvEstL1BhSamps), size = n_samps, replace = T), ] %>%
  mutate(sterilization_eff = b_betaL_sterilizedFsterilized/b_betaL_sterilizedFlive,
         live_litter = E_A_beta_wat * 10000,
         beta = (live_litter * sterilization_eff)/10000)


#### survival function

E_A_bh_fun <- function(disease, litter, iter) {
  
  # max survival
  E_A_lin_max = ifelse(disease == 1, 
                       E_A_df$int_wat[iter], 
                       E_A_df$int_fun[iter])
  E_A_max = exp(E_A_lin_max) / (1 + exp(E_A_lin_max))
  
  # add litter effect
  E_A_lit = ifelse(disease == 1, 
                   1 + E_A_beta_wat[iter] * litter, 
                   1 + E_A_beta_fun$beta[iter] * litter)
  
  E_A <- E_A_max / E_A_lit
  
  return(E_A)
}

