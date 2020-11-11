##### info ####

# file: microstegium_biomass_fung_parameter_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: sample from model coefficients to estimate survival


#### set-up ####

# tidyverse and brms packages must be loaded

# load model
load("output/microstegium_biomass_fung_model_2018_2019_density_exp.rda")

# extract posterior distributions
mvBioFuSamps <- posterior_samples(mvBioFuMod2)

# sample parameters
# use water values from year 2 only (year 1 were not individual plants), and apply fungicide effect from year 1 to these
b_A_df1 <- mvBioFuSamps[sample(nrow(mvBioFuSamps), size = n_samps, replace = T), ] %>%
  rename("b_fungicide_yearfyear2" = "b_fungicide:yearfyear2") %>%
  transmute(fun_eff_y1 = (b_Intercept + b_fungicide) / b_Intercept,
            int_wat_y2 = b_Intercept + b_yearfyear2,
            int_fun_y2 = b_Intercept + b_fungicide + b_yearfyear2 + b_fungicide_yearfyear2) %>%
  mutate(int_wat_y1 = int_wat_y2,
         int_fun_y1 = int_wat_y1 * fun_eff_y1) 

# select the first half of water values to go with the fungicide effect from year 1 and the second half to go with the fungicide effect from year 2
b_A_df <- tibble(int_wat_y1 = b_A_df1$int_wat_y1[1:(n_samps/2)],
                 int_fun_y1 = b_A_df1$int_fun_y1[1:(n_samps/2)],
                 int_wat_y2 = b_A_df1$int_wat_y2[(n_samps/2 + 1):n_samps],
                 int_fun_y2 = b_A_df1$int_fun_y2[(n_samps/2 + 1):n_samps]) %>%
  pivot_longer(everything(),
               names_to = c(".value", "year"),
               names_sep = "_y")
# remove the direct effect of fungicide by not including treatmentfungicide 


#### biomass function ####

b_A_fung_fun <- function(disease, iter) {
  
  # calculate survival
  b_A_expr <- ifelse(disease == 1, 
                     b_A_df$int_wat[iter],
                     b_A_df$int_fun[iter])
  
  b_A <- exp(b_A_expr)
  
  return(b_A)
}

