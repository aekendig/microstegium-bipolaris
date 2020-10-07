##### info ####

# file: elymus_adult_gs_survival_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/5/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_gs_survival_model_2019_density_exp.rda")

# extract posterior distributions
evAGsSurvD2Samps <- posterior_samples(evAGsSurvD2Mod2)

# sample parameters
H_P_int_wat <- sample(evAGsSurvD2Samps$b_Intercept, size = n_samps, replace = T)
H_P_int_fun <- H_P_int_wat + sample(evAGsSurvD2Samps$b_fungicide, size = n_samps, replace = T)

H_P_mv_dens_wat <- sample(evAGsSurvD2Samps$b_mv_seedling_density, size = n_samps, replace = T)
H_P_mv_dens_fun <- H_P_mv_dens_wat + sample(evAGsSurvD2Samps$"b_fungicide:mv_seedling_density", size = n_samps, replace = T)

H_P_evS_dens_wat <- sample(evAGsSurvD2Samps$b_ev_seedling_density, size = n_samps, replace = T)
H_P_evS_dens_fun <- H_P_evS_dens_wat + sample(evAGsSurvD2Samps$"b_fungicide:ev_seedling_density", size = n_samps, replace = T)

H_P_evA_dens_wat <- sample(evAGsSurvD2Samps$b_ev_adult_density, size = n_samps, replace = T)
H_P_evA_dens_fun <- H_P_evA_dens_wat + sample(evAGsSurvD2Samps$"b_fungicide:ev_adult_density", size = n_samps, replace = T)


#### survival function

H_P_fun <- function(disease, A_dens, S_dens, P_dens, iter) {
  
  # calculate survival
  H_P_lin_expr <- ifelse(disease == 1, 
                         H_P_int_wat[iter] + H_P_mv_dens_wat[iter] * A_dens + H_P_evS_dens_wat[iter] * S_dens + H_P_evA_dens_wat[iter] * P_dens,
                         H_P_int_fun[iter] + H_P_mv_dens_fun[iter] * A_dens + H_P_evS_dens_fun[iter] * S_dens + H_P_evA_dens_fun[iter] * P_dens)
  
  H_P <- exp(H_P_lin_expr)/(1 + exp(H_P_lin_expr))
  
  return(H_P)
}

