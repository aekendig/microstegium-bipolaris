##### info ####

# file: elymus_adult_gs_survival_parameter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
# load("output/elymus_adult_gs_survival_model_2019_density_exp.rda")
load("output/elymus_adult_gs_survival_fung_model_2019_density_exp.rda")

# extract posterior distributions
# evAGsSurvD2Samps <- posterior_samples(evAGsSurvD2Mod2)
evAGsSurvD2FuSamps <- posterior_samples(evAGsSurvD2FuMod2)

# sample parameters
# U_P_dens <- evAGsSurvD2Samps[sample(nrow(evAGsSurvD2Samps), size = n_samps, replace = T), ] %>%
#   rename("fungicide_mv_seedling_density" = "b_fungicide:mv_seedling_density",
#          "fungicide_ev_seedling_density" = "b_fungicide:ev_seedling_density",
#          "fungicide_ev_adult_density" = "b_fungicide:ev_adult_density") %>%
#   mutate(int_wat = b_Intercept,
#          int_fun = int_wat + b_fungicide,
#          mv_dens_wat = b_mv_seedling_density,
#          mv_dens_fun = mv_dens_wat + fungicide_mv_seedling_density,
#          evS_dens_wat = b_ev_seedling_density,
#          evS_dens_fun = evS_dens_wat + fungicide_ev_seedling_density,
#          evA_dens_wat = b_ev_adult_density,
#          evA_dens_fun = evA_dens_wat + fungicide_ev_adult_density)

u_P_df <- evAGsSurvD2FuSamps[sample(nrow(evAGsSurvD2FuSamps), size = n_samps, replace = T), ] %>%
  mutate(int_wat = b_Intercept,
         int_fun = int_wat + b_fungicide)


#### survival function ####

# U_P_fun <- function(disease, A_dens, S_dens, P_dens, iter) {
#   
#   U_P_lin_expr <- ifelse(disease == 1, 
#                          U_P_dens$int_wat[iter] + U_P_dens$mv_dens_wat[iter] * A_dens + U_P_dens$evS_dens_wat[iter] * S_dens + U_P_dens$evA_dens_wat[iter] * P_dens,
#                          U_P_dens$int_fun[iter] + U_P_dens$mv_dens_fun[iter] * A_dens + U_P_dens$evS_dens_fun[iter] * S_dens + U_P_dens$evA_dens_fun[iter] * P_dens)
#   
#   U_P <- ifelse(exp(U_P_lin_expr) == Inf, 1, exp(U_P_lin_expr)/(1 + exp(U_P_lin_expr)))
#   
#   return(U_P)
# }

u_P_fun <- function(disease, iter) {
  
  u_P_lin_expr <- ifelse(disease == 1, 
                         u_P_df$int_wat[iter],
                         u_P_df$int_fun[iter])
  
  u_P <- exp(u_P_lin_expr)/(1 + exp(u_P_lin_expr))
  
  return(u_P)
}
