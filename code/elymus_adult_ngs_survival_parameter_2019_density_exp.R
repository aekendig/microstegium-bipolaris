##### info ####

# file: elymus_adult_ngs_survival_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_ngs_survival_model_2019_density_exp.rda")

# extract posterior distributions
evANgsSurvD1Samps <- posterior_samples(evANgsSurvD1Mod2)

# sample parameters
w_P_df <- evANgsSurvD1Samps[sample(nrow(evANgsSurvD1Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)


#### survival function ####

w_P_fun <- function(disease, iter) {
  
  # non-growing season survival
  w_P_lin_expr <- ifelse(disease == 1, 
                         w_P_df$int_wat[iter],
                         w_P_df$int_fun[iter])
  
  w_P <- exp(w_P_lin_expr)/(1 + exp(w_P_lin_expr))
  
  return(w_P)
}

