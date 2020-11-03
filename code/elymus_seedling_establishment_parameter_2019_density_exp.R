##### info ####

# file: elymus_seedling_establishment_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_gs_survival_model_2019_density_exp.rda")
load("output/elymus_litter_establishment_model_2018_greenhouse_exp.rda")
load("output/microstegium_litter_establishment_model_2018_litter_exp.rda")

# extract posterior distributions
evSGsSurvD2Samps <- posterior_samples(evSGsSurvD2Mod2)
evLitEstSamps <- posterior_samples(evLitEstMod2)
evEstL1Samps <- posterior_samples(mvEstL1Mod2)

# sample parameters
E_S_dens <- evSGsSurvD2Samps[sample(nrow(evSGsSurvD2Samps), size = n_samps, replace = T), ] %>%
  rename("fungicide_mv_seedling_density" = "b_fungicide:mv_seedling_density",
         "fungicide_ev_seedling_density" = "b_fungicide:ev_seedling_density",
         "fungicide_ev_adult_density" = "b_fungicide:ev_adult_density") %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide,
         mv_dens_wat = b_mv_seedling_density,
         mv_dens_fun = mv_dens_wat + fungicide_mv_seedling_density,
         evS_dens_wat = b_ev_seedling_density,
         evS_dens_fun = evS_dens_wat + fungicide_ev_seedling_density,
         evA_dens_wat = b_ev_adult_density,
         evA_dens_fun = evA_dens_wat + fungicide_ev_adult_density)

E_S_lit_wat <- sample(evLitEstSamps$b_litter.g.cm2, size = n_samps, replace = T)/10000

E_S_lit_fun <- evEstL1Samps[sample(nrow(evEstL1Samps), size = n_samps, replace = T), ] %>%
  rename("interaction" = "b_litter.g.cm2:sterilized") %>%
  mutate(sterilization_eff = interaction/b_litter.g.cm2,
         live_litter = E_S_lit_wat,
         sterilized_litter = (live_litter + live_litter * sterilization_eff)/10000) %>%
  select(sterilized_litter)


#### survival function

E_S_fun <- function(disease, g.A, A_dens, g.S, S_dens, P_dens, litter, iter) {
  
  # calculate survival
  E_S_lin_expr <- ifelse(disease == 1, 
                         E_S_dens$int_wat[iter] + E_S_dens$mv_dens_wat[iter] * g.A * A_dens + E_S_dens$evS_dens_wat[iter] * g.S * S_dens + E_S_dens$evA_dens_wat[iter] * P_dens + E_S_lit_wat[iter] * litter,
                         E_S_dens$int_fun[iter] + E_S_dens$mv_dens_fun[iter] * g.A * A_dens + E_S_dens$evS_dens_fun[iter] * g.S * S_dens + E_S_dens$evA_dens_fun[iter] * P_dens + E_S_lit_fun$sterilized_litter[iter] * litter)
  
  E_S <- ifelse(exp(E_S_lin_expr) == Inf, 1, exp(E_S_lin_expr)/(1 + exp(E_S_lin_expr)))
  
  return(E_S)
}

