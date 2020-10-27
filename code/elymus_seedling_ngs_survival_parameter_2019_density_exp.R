##### info ####

# file: elymus_seedling_ngs_survival_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_seedling_ngs_survival_model_2019_density_exp.rda")

# extract posterior distributions
evSNgsSurvD1Samps <- posterior_samples(evSNgsSurvD1Mod2)

# sample parameters
w_S_df <- evSNgsSurvD1Samps[sample(nrow(evSNgsSurvD1Samps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)


#### survival function

w_S_fun <- function(disease, iter) {
  
  # calculate survival
  w_S_lin_expr <- ifelse(disease == 1, 
                         w_S_df$int_wat[iter],
                         w_S_df$int_fun[iter])
  
  w_S <- exp(w_S_lin_expr)/(1 + exp(w_S_lin_expr))
  
  return(w_S)
}

