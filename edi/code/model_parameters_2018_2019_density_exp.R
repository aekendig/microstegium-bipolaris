##### outputs ####

# continuous_model_parameters_2018_2019_density_exp.csv
# discrete_model_parameters_2018_2019_density_exp.csv


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)


#### initial biomass ####

# divide minimum observed biomass by 10

# import data
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# edit data
init_bio_parms <- tibble(Parameter = c("b_A", "b_F", "b_P"),
                         Plant_group = c("annual",
                                         "fy perennial",
                                         "adult perennial"),
                         Estimate = c(min(mvBioD2Dat$biomass_weight.g, na.rm = T)/10,
                                      min(filter(evBioD2Dat, ID %in% c("1", "2", "3"))$weight, na.rm = T)/10,
                                      min(filter(evBioD2Dat, ID == "A")$weight, na.rm = T)/10),
                         Lower = rep(NA_real_, 3),
                         Upper = rep(NA_real_, 3),
                         Units = rep("g", 3),
                         Description = "initial biomass",
                         Source = rep("experiment", 3))


#### growth & interactions ####

# load model
load("output/focal_growth_biomass_model_2019_density_exp.rda")

# initial biomass values
b_A0 <- init_bio_parms %>% filter(Parameter == "b_A") %>% pull(Estimate)
b_F0 <- init_bio_parms %>% filter(Parameter == "b_F") %>% pull(Estimate)
b_P0 <- init_bio_parms %>% filter(Parameter == "b_P") %>% pull(Estimate)

# sample model
growth_mod_samps <- as_draws_df(growthD2Mod2)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(r_A_ctrl = Intercept - log(b_A0),
            r_F_ctrl = (Intercept + focs) - log(b_F0),
            r_P_ctrl = (Intercept + foca) - log(b_P0),
            r_A_fung = (Intercept + fungicide) - log(b_A0),
            r_F_fung = (Intercept + focs + fungicide + focs_fungicide) - log(b_F0),
            r_P_fung = (Intercept + foca + fungicide + foca_fungicide) - log(b_P0),
            alpha_AA_ctrl = plot_biomass,
            alpha_FA_ctrl = plot_biomass + focs_plot_biomass,
            alpha_PA_ctrl = plot_biomass + foca_plot_biomass,
            alpha_AF_ctrl = plot_biomass + plot_biomass_bgs,
            alpha_FF_ctrl = plot_biomass + focs_plot_biomass + plot_biomass_bgs + focs_plot_biomass_bgs,
            alpha_PF_ctrl = plot_biomass + foca_plot_biomass + plot_biomass_bgs + foca_plot_biomass_bgs,
            alpha_AP_ctrl = plot_biomass + plot_biomass_bga,
            alpha_FP_ctrl = plot_biomass + plot_biomass_bga + focs_plot_biomass + focs_plot_biomass_bga,
            alpha_PP_ctrl = plot_biomass + foca_plot_biomass + plot_biomass_bga + foca_plot_biomass_bga,
            alpha_AA_fung = plot_biomass + fungicide_plot_biomass,
            alpha_FA_fung = plot_biomass + fungicide_plot_biomass + focs_plot_biomass + focs_fungicide_plot_biomass,
            alpha_PA_fung = plot_biomass + fungicide_plot_biomass + foca_plot_biomass + foca_fungicide_plot_biomass,
            alpha_AF_fung = plot_biomass + plot_biomass_bgs + fungicide_plot_biomass + fungicide_plot_biomass_bgs,
            alpha_FF_fung = plot_biomass + focs_plot_biomass + plot_biomass_bgs + focs_plot_biomass_bgs + fungicide_plot_biomass + focs_fungicide_plot_biomass + fungicide_plot_biomass_bgs + focs_fungicide_plot_biomass_bgs,
            alpha_PF_fung = plot_biomass + foca_plot_biomass + plot_biomass_bgs + foca_plot_biomass_bgs + fungicide_plot_biomass + foca_fungicide_plot_biomass + fungicide_plot_biomass_bgs + foca_fungicide_plot_biomass_bgs,
            alpha_AP_fung = plot_biomass + plot_biomass_bga + fungicide_plot_biomass + fungicide_plot_biomass_bga,
            alpha_FP_fung = plot_biomass + plot_biomass_bga + focs_plot_biomass + focs_plot_biomass_bga + fungicide_plot_biomass + fungicide_plot_biomass_bga + focs_fungicide_plot_biomass + focs_fungicide_plot_biomass_bga,
            alpha_PP_fung = plot_biomass + foca_plot_biomass + plot_biomass_bga + foca_plot_biomass_bga + fungicide_plot_biomass + foca_fungicide_plot_biomass + fungicide_plot_biomass_bga + foca_fungicide_plot_biomass_bga) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate")

# growth rates and interaction coefficients
growth_mod_parms <- growth_mod_samps %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  ungroup() %>%
  transmute(Parameter = str_replace(parameter, "_fung", ""),
            Parameter = str_replace(Parameter, "_ctrl", ""),
            Plant_group = case_when(str_detect(parameter, "_AA") == T ~ "annual/annual",
                                    str_detect(parameter, "_FA") == T ~ "annual/fy perennial",
                                    str_detect(parameter, "_PA") == T ~ "annual/adult perennial",
                                    str_detect(parameter, "_AF") == T ~ "fy perennial/annual",
                                    str_detect(parameter, "_FF") == T ~ "fy perennial/fy perennial",
                                    str_detect(parameter, "_PF") == T ~ "fy perennial/adult perennial",
                                    str_detect(parameter, "_AP") == T ~ "adult perennial/annual",
                                    str_detect(parameter, "_FP") == T ~ "adult perennial/fy perennial",
                                    str_detect(parameter, "_PP") == T ~ "adult perennial/adult perennial",
                                    str_detect(parameter, "r_A") == T ~ "annual",
                                    str_detect(parameter, "r_F") == T ~ "perennial fy",
                                    str_detect(parameter, "r_P") == T ~ "perennial adult"),
            Description = case_when(str_detect(parameter, "alpha") == T ~ "interaction", 
                                    str_detect(parameter, "r") == T ~ "growth rate"),
            Treatment = if_else(str_detect(parameter, "_fung") == T, "fungicide", "control"),
            Estimate = case_when(str_detect(parameter, "r_") == T ~ estimate/160, # r may not be sig, but we need an estimate
                                 str_detect(parameter, "alpha_") == T & estimate > 0 ~ 0,
                                 TRUE ~ estimate * -1),
            Lower = case_when(str_detect(parameter, "r_") == T ~ .lower/160,
                              str_detect(parameter, "alpha_") == T & estimate > 0 ~ 0,
                              TRUE ~ .upper * -1), # switch order because of switched sign
            Upper = case_when(str_detect(parameter, "r_") == T ~ .upper/160,
                              str_detect(parameter, "alpha_") == T & estimate > 0 ~ 0,
                              TRUE ~ .lower * -1),
            Units = case_when(str_detect(parameter, "alpha") == T ~ "m^2^ g^-1^", 
                              str_detect(parameter, "r") == T ~ "day^-1^"),
            Source = "experiment")


#### alternative Mv comp coef ####

# because intraspecific competition coefficient is higher in fungicide than control

load("output/mv_plot_biomass_density_model_2019_dens_exp.rda")

# ratio of control:fungicide from plot biomass model
mv_alpha_ratio <- as_draws_df(mvBioDensMod)  %>%
  transmute(control = b_alpha_treatmentwater,
            fungicide = b_alpha_treatmentfungicide,
            ctrl_fung_ratio = control / fungicide) %>%
  mean_hdi(ctrl_fung_ratio) %>%
  ungroup() %>%
  pull(ctrl_fung_ratio)

# re-calculate control comp. coeff using fungicide 
mv_alpha_samps <- as_draws_df(growthD2Mod2)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(alpha_AA = (plot_biomass + fungicide_plot_biomass) * mv_alpha_ratio * -1) %>%
  mean_hdi(alpha_AA) %>%
  rename(estimate = "alpha_AA", lower = ".lower", upper = ".upper")

# replace comp. coeff
growth_mod_parms2 <- growth_mod_parms %>%
  mutate(Estimate = if_else(Parameter == "alpha_AA" & Treatment == "control",
                            mv_alpha_samps$estimate, Estimate),
         Lower = if_else(Parameter == "alpha_AA" & Treatment == "control",
                            mv_alpha_samps$lower, Lower),
         Upper = if_else(Parameter == "alpha_AA" & Treatment == "control",
                            mv_alpha_samps$upper, Upper))


#### alternative Ev comp coefs ####

# because intraspecific adult Ev competition coefficient is non-existent with fungicide

# use other intraspecific coefficients
ev_PP_fung <- growth_mod_parms2 %>%
  filter(Parameter == "alpha_FP" & Treatment == "fungicide")

# replace comp. coeff
growth_mod_parms3 <- growth_mod_parms2 %>%
  mutate(Estimate = case_when(Parameter == "alpha_PP" & Treatment == "fungicide" ~ ev_PP_fung$Estimate,
                              TRUE ~ Estimate),
         Lower = case_when(Parameter == "alpha_PP" & Treatment == "fungicide" ~ ev_PP_fung$Lower,
                              TRUE ~ Lower),
         Upper = case_when(Parameter == "alpha_PP" & Treatment == "fungicide" ~ ev_PP_fung$Upper,
                              TRUE ~ Upper))


#### transmission ####

# load model
trans_mod_coef <- read_csv("output/focal_severity_model_coefficients_2019_dens_exp.csv")

# make negative values zero
# multiply all by constant
trans_mod_parms <- trans_mod_coef %>%
  mutate(source = fct_recode(source,
                             "A" = "Surrounding invader",
                             "A" = "Invader (Mv)",
                             "F" = "1st yr comp. (Ev)",
                             "P" = "Adult comp. (Ev)"),
         plant_group = fct_recode(plant_group,
                                  "beta_A" = "Invader (Mv)",
                                  "beta_F" = "1st yr comp. (Ev)",
                                  "beta_P" = "Adult comp. (Ev)"),
         Parameter = paste0(plant_group, source),
         Treatment = if_else(fungicide == 1, "fungicide", "control"),
         Estimate = if_else(trend < 0, 0, trend * 5e-4)) %>%
  group_by(Parameter, Treatment) %>%
  summarize(Estimate = max(Estimate)) %>%
  ungroup() %>%
  mutate(Plant_group = case_when(str_detect(Parameter, "_AA") == T ~ "annual/annual",
                                 str_detect(Parameter, "_FA") == T ~ "annual/fy perennial",
                                 str_detect(Parameter, "_PA") == T ~ "annual/adult perennial",
                                 str_detect(Parameter, "_AF") == T ~ "fy perennial/annual",
                                 str_detect(Parameter, "_FF") == T ~ "fy perennial/fy perennial",
                                 str_detect(Parameter, "_PF") == T ~ "fy perennial/adult perennial",
                                 str_detect(Parameter, "_AP") == T ~ "adult perennial/annual",
                                 str_detect(Parameter, "_FP") == T ~ "adult perennial/fy perennial",
                                 str_detect(Parameter, "_PP") == T ~ "adult perennial/adult perennial"),
         Description = "transmission",
         Units = "m^2^ g^-1^ day^-1^",
         Source = "experiment")


#### germination ####

load("output/mv_germination_fungicide_model_2018_density_exp.rda")

mv_germ_parms <- as_draws_df(mvGermD1Mod2)  %>%
  transmute(g_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            g_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = "g_A",
            Plant_group = "annual",
            Description = "germination fraction",
            Treatment = if_else(parameter == "g_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "NA",
            Source = "experiment")


load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

ev_germ_parms <- as_draws_df(evGermMod)  %>%
  transmute(g_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            g_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>%
  transmute(Parameter = "g_P",
            Plant_group = "perennial",
            Description = "germination fraction",
            Treatment = if_else(parameter == "g_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "NA",
            Source = "experiment")


#### seeds per biomass ####

# import models
load("output/mv_seeds_per_biomass_untransformed_model_2019_density_exp.rda")
load("output/evS_seeds_per_biomass_untransformed_model_2019_density_exp.rda")
load("output/evA_seeds_per_biomass_untransformed_model_2019_density_exp.rda")

mv_seed_parms <- as_draws_df(mvSeedsBioD2Mod2)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(c_W = b_biomass_weight.g,
            c_F = b_biomass_weight.g + b_biomass_weight.g_fungicide) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = "c_A",
            Plant_group = "annual",
            Description = "seed conversion",
            Treatment = if_else(parameter == "c_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "seeds g^-1^",
            Source = "experiment")

evS_seed_parms <- as_draws_df(evSSeedsBioD2Mod2)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(c_W = b_biomass_weight.g,
            c_F = b_biomass_weight.g + b_biomass_weight.g_fungicide) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = "c_F",
            Plant_group = "fy perennial",
            Description = "seed conversion",
            Treatment = if_else(parameter == "c_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "seeds g^-1^",
            Source = "experiment")

evA_seed_parms <- as_draws_df(evASeedsBioD2Mod2)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(c_W = b_biomass_weight.g,
            c_F = b_biomass_weight.g + b_biomass_weight.g_fungicide) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = "c_P",
            Plant_group = "adult perennial",
            Description = "seed conversion",
            Treatment = if_else(parameter == "c_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "seeds g^-1^",
            Source = "experiment")


#### establishment ####

# import model
load("output/survival_fungicide_model_2019_density_exp.rda")

est_mod_parms <- as_draws_df(survFungD2Mod)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(e_A_W = exp(b_Intercept)/(1 + exp(b_Intercept)),
            e_A_F = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide)),
            e_P_W = exp(b_Intercept + b_focs)/(1 + exp(b_Intercept + b_focs)),
            e_P_F = exp(b_Intercept + b_focs + b_fungicide)/(1 + exp(b_Intercept + b_focs + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>%
  transmute(Parameter = str_sub(parameter, 1, 3),
            Plant_group = rep(c("annual", "fy perennial"), each = 2),
            Description = "establishment fraction",
            Treatment = if_else(str_detect(parameter, "W") == T, "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "NA",
            Source = "experiment")


#### adult survival ####

# import model
load("output/ev_adult_survival_fungicide_model_2018_2019_density_exp.rda")

evA_surv_parms <- as_draws_df(adultSurvD1Mod2)  %>%
  transmute(l_W = exp(b_Intercept)/(1 + exp(b_Intercept)),
            l_F = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>% # two separate hdi's - continuous is very similar
  transmute(Parameter = "l_P",
            Plant_group = "adult perennial",
            Description = "survival fraction",
            Treatment = if_else(parameter == "l_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "NA",
            Source = "experiment")


#### biomass mortality ####

# import data
ev_fung_gh <- read_csv("data/ev_biomass_dec_2019_fungicide_exp.csv")

# only have measurement for Ev
# apply to Mv as well
bio_mort_parms <- ev_fung_gh %>%
  select(-notes) %>%
  pivot_wider(names_from = biomass_type,
              values_from = weight.g) %>%
  mutate(prop_dead = dead / (dead + live)) %>%
  summarise(Prop_dead = mean(prop_dead, na.rm = T)) %>%
  expand_grid(tibble(Parameter = c("m_A", "m_F", "m_P"),
                     Plant_group = c("annual",
                                     "fy perennial",
                                     "adult perennial"),
                     Description = "biomass mortality")) %>%
  mutate(Estimate = Prop_dead/292, # experiment lasted 292 days
         Lower = NA,
         Upper = NA,
         Units = "NA",
         Source = "experiment") %>%
  select(-Prop_dead)


#### adult new biomass ####

# use greenhouse experiment and survival estimate
evA_new_bio_parms <- ev_fung_gh %>%
  select(-notes) %>%
  pivot_wider(names_from = biomass_type,
              values_from = weight.g) %>%
  mutate(prop_live = live / (dead + live)) %>%
  summarise(Prop_live = mean(prop_live, na.rm = T)) %>%
  expand_grid(evA_surv_parms %>%
                select(Plant_group, Treatment, Estimate)) %>%
  mutate(Parameter = "n_P",
         Description = "surviving biomass",
         Estimate = Prop_live * Estimate,
         Lower = NA,
         Upper = NA,
         Units = "NA",
         Source = "experiment") %>%
  select(-Prop_live)
# the proportion of live biomass at the end of the growing season (FL greenhouse)
# multiplied by the proportion of adults that survive


#### litter sensitivity ####

# import model
load("output/litter_reu_mv_establishment_model.rda")

mv_lit_parms <- as_draws_df(mv_bh_mod)  %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(gamma_A = beta_Intercept) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = parameter,
            Plant_group = "annual",
            Description = "sensitivity to litter",
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "m^2^ g^-1^",
            Source = "Benitez et al. 2021")

load("output/litter_reu_ev_establishment_model.rda")

ev_lit_parms <- as_draws_df(ev_bh_mod)  %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(gamma_P = beta_Intercept) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = parameter,
            Plant_group = "perennial",
            Description = "sensitivity to litter",
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "m^2^ g^-1^",
            Source = "Benitez et al. 2021")


#### literature/estimated parameters ####

lit_parms <- tibble(Parameter = c("h", 
                                  "d", 
                                  "s_P", "s_A",
                                  "beta_AC", "beta_FC", "beta_PC", 
                                  "v_A", "v_F", "v_P", 
                                  "a",
                                  "dis_thresh"),
                    Plant_group = c(NA_character_,
                                    NA_character_,
                                    "perennial",
                                    "annual",
                                    rep(c("annual", "fy perennial", "adult perennial"), 2),
                                    NA_character_,
                                    NA_character_),
                    Description = c("inoculum addition to litter",
                                    "litter decomposition fraction",
                                    "surviving seed fraction",
                                    "surviving seed fraction",
                                    "litter transmission",
                                    "litter transmission",
                                    "litter transmission",
                                    "infected tissue loss",
                                    "infected tissue loss",
                                    "infected tissue loss",
                                    "inoculum loss from litter",
                                    "disease threshold"),
                    Treatment = rep(NA_character_, 12),
                    Estimate = c(5500, 0.59, 0.05, 0.15, rep(1e-7, 3), rep(1e-4, 3), 0.95, 0.15),
                    Source = c("Benitez et al. 2021", "DeMeester and Richter 2010", "Garrison and Stier 2010", "Redwood et al. 2018", rep("NA", 8))) %>%
  mutate(Lower = NA_real_,
         Upper = NA_real_,
         Units = c("g^-1^", rep("NA", 3), rep ("day^-1^ g^-1^", 3), rep("day^-1^", 3), "day^-1^", "NA"))


#### combine ####
parms <- growth_mod_parms3 %>%
  full_join(trans_mod_parms) %>%
  full_join(mv_germ_parms) %>%
  full_join(ev_germ_parms) %>%
  full_join(mv_seed_parms) %>%
  full_join(evS_seed_parms) %>%
  full_join(evA_seed_parms) %>%
  full_join(est_mod_parms) %>%
  full_join(evA_surv_parms) %>%
  full_join(evA_new_bio_parms) %>%
  full_join(init_bio_parms) %>%
  full_join(bio_mort_parms) %>%
  full_join(mv_lit_parms) %>%
  full_join(ev_lit_parms) %>%
  full_join(lit_parms)

# continuous parameters
cont_parms <- growth_mod_parms3 %>%
  full_join(trans_mod_parms) %>%
  full_join(bio_mort_parms) %>%
  full_join(lit_parms %>%
              filter(Parameter %in% c("beta_AC", "beta_FC", "beta_PC",
                                      "v_A", "v_F", "v_P",
                                      "a", "h"))) %>%
  select(Parameter, Treatment, Estimate)

# discrete parameters
disc_parms <- mv_germ_parms %>%
  full_join(ev_germ_parms) %>%
  full_join(mv_seed_parms) %>%
  full_join(evS_seed_parms) %>%
  full_join(evA_seed_parms) %>%
  full_join(est_mod_parms) %>%
  full_join(evA_surv_parms) %>%
  full_join(evA_new_bio_parms) %>%
  full_join(init_bio_parms) %>%
  full_join(mv_lit_parms) %>%
  full_join(ev_lit_parms) %>%
  full_join(lit_parms %>%
              filter(Parameter %in% c("s_P", "s_A",
                                       "d", "h", "dis_thresh"))) %>%
  select(Parameter, Treatment, Estimate)


#### export ####
write_csv(parms, "output/model_parameters_2018_2019_density_exp.csv")
write_csv(cont_parms, "output/continuous_model_parameters_2018_2019_density_exp.csv")
write_csv(disc_parms, "output/discrete_model_parameters_2018_2019_density_exp.csv")
