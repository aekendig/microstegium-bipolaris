##### info ####

# file: elymus_germination_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/10/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/elymus_germination_model_2018_2019_density_exp.rda")

# extract posterior distribution
evGermSamps <- posterior_samples(evGermMod2)

# sample parameters
g_S_df <- evGermSamps[sample(nrow(evGermSamps), size = n_samps / 2, replace = T), ] %>%
  rename("b_fungicide_year2" = "b_fungicide:yearfyear2") %>%
  transmute(int_wat_y1 = b_Intercept,
            int_fun_y1 = int_wat_y1 + b_fungicide,
            int_wat_y2 = int_wat_y1 + b_yearfyear2,
            int_fun_y2 = int_fun_y1 + b_yearfyear2 + b_fungicide_year2) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")



#### germination function ####

g_S_fun <- function(disease, iter) {
  
  g_S_expr <- ifelse(disease == 1, g_S_df$int_wat[iter], g_S_df$int_fun[iter])
  
  g_S <- exp(g_S_expr) / (1 + exp(g_S_expr))
  
  return(g_S)
}
