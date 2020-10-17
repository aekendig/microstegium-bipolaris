##### info ####

# file: elymus_seedling_ngs_survival_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/6/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_ngs_survival_model_2019_density_exp.rda")

# extract posterior distributions
evSNgsSurvD1Samps <- posterior_samples(evSNgsSurvD1Mod2)

# sample parameters
W_S_int_wat <- sample(evSNgsSurvD1Samps$b_Intercept, size = n_samps, replace = T)
W_S_int_fun <- W_S_int_wat + sample(evSNgsSurvD1Samps$b_fungicide, size = n_samps, replace = T)

W_S_mv_dens_wat <- sample(evSNgsSurvD1Samps$b_mv_seedling_density, size = n_samps, replace = T)
W_S_mv_dens_fun <- W_S_mv_dens_wat + sample(evSNgsSurvD1Samps$"b_fungicide:mv_seedling_density", size = n_samps, replace = T)

W_S_evS_dens_wat <- sample(evSNgsSurvD1Samps$b_ev_seedling_density, size = n_samps, replace = T)
W_S_evS_dens_fun <- W_S_evS_dens_wat + sample(evSNgsSurvD1Samps$"b_fungicide:ev_seedling_density", size = n_samps, replace = T)

W_S_evA_dens_wat <- sample(evSNgsSurvD1Samps$b_ev_adult_density, size = n_samps, replace = T)
W_S_evA_dens_fun <- W_S_evA_dens_wat + sample(evSNgsSurvD1Samps$"b_fungicide:ev_adult_density", size = n_samps, replace = T)


#### survival function

W_S_fun <- function(disease, A_dens, S_dens, P_dens, iter) {
  
  # calculate survival
  W_S_lin_expr <- ifelse(disease == 1, 
                         W_S_int_wat[iter] + W_S_mv_dens_wat[iter] * A_dens + W_S_evS_dens_wat[iter] * S_dens + W_S_evA_dens_wat[iter] * P_dens,
                         W_S_int_fun[iter] + W_S_mv_dens_fun[iter] * A_dens + W_S_evS_dens_fun[iter] * S_dens + W_S_evA_dens_fun[iter] * P_dens)
  
  W_S <- ifelse(exp(W_S_lin_expr) == Inf, 1, exp(W_S_lin_expr)/(1 + exp(W_S_lin_expr)))
  
  return(W_S)
}

