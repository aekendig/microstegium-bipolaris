##### info ####

# file: elymus_adult_gs_survival_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/4/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_adult_gs_survival_fung_model_2018_2019_density_exp.rda")

# extract posterior distributions
evAGsSurvFuSamps <- posterior_samples(evAGsSurvFuMod2)

# sample parameters
u_P_df <- evAGsSurvFuSamps[sample(nrow(evAGsSurvFuSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_year2" = "b_fungicide:yearfyear2") %>%
  transmute(int_wat_y1 = b_Intercept,
            int_fun_y1 = int_wat_y1 + b_fungicide,
            int_wat_y2 = int_wat_y1 + b_yearfyear2,
            int_fun_y2 = int_fun_y1 + b_yearfyear2 + b_fungicide_year2) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")


#### survival function ####

u_P_fun <- function(disease, iter) {
  
  u_P_lin_expr <- ifelse(disease == 1, 
                         u_P_df$int_wat[iter],
                         u_P_df$int_fun[iter])
  
  u_P <- exp(u_P_lin_expr)/(1 + exp(u_P_lin_expr))
  
  return(u_P)
}
