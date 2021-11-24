##### info ####

# file: model_parameters_2018_2019_density_exp.R
# author: Amy Kendig
# date last edited: 11/22/21
# goal: estimate parameters for each treatment


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
  transmute(r_A = Intercept - log(b_A0),
            r_F = (Intercept + focs) - log(b_F0),
            r_P = (Intercept + foca) - log(b_P0),
            r_A_fung = (Intercept + fungicide) - log(b_A0),
            r_F_fung = (Intercept + focs + fungicide + focs_fungicide) - log(b_F0),
            r_P_fung = (Intercept + foca + fungicide + foca_fungicide) - log(b_P0),
            alpha_AA = plot_biomass,
            alpha_FA = plot_biomass + plot_biomass_focs,
            alpha_PA = plot_biomass + plot_biomass_foca,
            alpha_AF = plot_biomass + plot_biomass_bgs,
            alpha_FF = plot_biomass + plot_biomass_focs + plot_biomass_bgs + plot_biomass_focs_bgs,
            alpha_PF = plot_biomass + plot_biomass_foca + plot_biomass_bgs + plot_biomass_foca_bgs,
            alpha_AP = plot_biomass + plot_biomass_bga,
            alpha_FP = plot_biomass + plot_biomass_bga + plot_biomass_focs + plot_biomass_focs_bga,
            alpha_PP = plot_biomass + plot_biomass_foca + plot_biomass_bga + plot_biomass_foca_bga,
            alpha_AA_fung = plot_biomass + plot_biomass_fungicide,
            alpha_FA_fung = plot_biomass + plot_biomass_fungicide + plot_biomass_focs + plot_biomass_focs_fungicide,
            alpha_PA_fung = plot_biomass + plot_biomass_fungicide + plot_biomass_foca + plot_biomass_foca_fungicide,
            alpha_AF_fung = plot_biomass + plot_biomass_bgs + plot_biomass_fungicide + plot_biomass_bgs_fungicide,
            alpha_FF_fung = plot_biomass + plot_biomass_focs + plot_biomass_bgs + plot_biomass_focs_bgs + plot_biomass_fungicide + plot_biomass_focs_fungicide + plot_biomass_bgs_fungicide + plot_biomass_focs_bgs_fungicide,
            alpha_PF_fung = plot_biomass + plot_biomass_foca + plot_biomass_bgs + plot_biomass_foca_bgs + plot_biomass_fungicide + plot_biomass_foca_fungicide + plot_biomass_bgs_fungicide + plot_biomass_foca_bgs_fungicide,
            alpha_AP_fung = plot_biomass + plot_biomass_bga + plot_biomass_fungicide + plot_biomass_bga_fungicide,
            alpha_FP_fung = plot_biomass + plot_biomass_bga + plot_biomass_focs + plot_biomass_focs_bga + plot_biomass_fungicide + plot_biomass_bga_fungicide + plot_biomass_focs_fungicide + plot_biomass_focs_bga_fungicide,
            alpha_PP_fung = plot_biomass + plot_biomass_foca + plot_biomass_bga + plot_biomass_foca_bga + plot_biomass_fungicide + plot_biomass_foca_fungicide + plot_biomass_bga_fungicide + plot_biomass_foca_bga_fungicide) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate")

# growth rates and interaction alphaficients
growth_mod_parms <- growth_mod_samps %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  ungroup() %>%
  transmute(Parameter = parameter,
            Plant_group = rep(c("annual/annual",
                                "annual/fy perennial",
                                "annual/adult perennial",
                                "fy perennial/annual",
                                "fy perennial/fy perennial",
                                "fy perennial/adult perennial",
                                "adult perennial/annual",
                                "adult perennial/fy perennial",
                                "adult perennial/adult perennial",
                                "annual",
                                "perennial fy",
                                "perennial adult"), each = 2),
            Description = c(rep("interaction", 18), rep("growth rate", 6)),
            Treatment = rep(c("control (water)", "fungicide"), 12),
            Estimate = case_when(str_detect(parameter, "r_") == T ~ estimate/160, # r may not be sig, but we need an estimate
                                 str_detect(parameter, "alpha_") == T & estimate > 0 ~ 0,
                                 TRUE ~ estimate * -1),
            Lower = case_when(str_detect(parameter, "r_") == T ~ .lower/160,
                              str_detect(parameter, "alpha_") == T & estimate > 0 ~ 0,
                              TRUE ~ .upper * -1), # switch order because of switched sign
            Upper = case_when(str_detect(parameter, "r_") == T ~ .upper/160,
                              str_detect(parameter, "alpha_") == T & estimate > 0 ~ 0,
                              TRUE ~ .lower * -1),
            Units = c(rep("g^-1^", 18), rep("day^-1^", 6)),
            Source = "experiment")


#### transmission ####

# load model
load("output/plot_transmission_model_2019_density_exp.rda")

# sample model
trans_mod_samps <- as_draws_df(sevD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(beta_AA = bg_severity,
            beta_AA_fung = bg_severity + bg_severity_fungicide,
            beta_FA = bg_severity + bg_severity_focs,
            beta_FA_fung = bg_severity + bg_severity_focs + bg_severity_fungicide + bg_severity_focs_fungicide,
            beta_PA = bg_severity + bg_severity_foca,
            beta_PA_fung = bg_severity + bg_severity_foca + bg_severity_fungicide + bg_severity_foca_fungicide,
            beta_AF = bg_severity + bg_severity_bgs,
            beta_AF_fung = bg_severity + bg_severity_bgs + bg_severity_fungicide + bg_severity_bgs_fungicide,
            beta_FF = bg_severity + bg_severity_focs + bg_severity_bgs + bg_severity_focs_bgs,
            beta_FF_fung = bg_severity + bg_severity_focs + bg_severity_bgs + bg_severity_focs_bgs + bg_severity_fungicide + bg_severity_bgs_fungicide + bg_severity_focs_fungicide + bg_severity_focs_bgs_fungicide,
            beta_PF = bg_severity + bg_severity_foca + bg_severity_bgs + bg_severity_foca_bgs,
            beta_PF_fung = bg_severity + bg_severity_foca + bg_severity_bgs + bg_severity_foca_bgs + bg_severity_fungicide + bg_severity_bgs_fungicide + bg_severity_foca_fungicide + bg_severity_foca_bgs_fungicide,
            beta_AP = bg_severity + bg_severity_bga,
            beta_AP_fung = bg_severity + bg_severity_bga + bg_severity_fungicide + bg_severity_bga_fungicide,
            beta_FP = bg_severity + bg_severity_bga + bg_severity_focs + bg_severity_focs_bga,
            beta_FP_fung = bg_severity + bg_severity_bga + bg_severity_focs + bg_severity_focs_bga + bg_severity_fungicide + bg_severity_bga_fungicide + bg_severity_focs_fungicide + bg_severity_focs_bga_fungicide,
            beta_PP = bg_severity + bg_severity_foca + bg_severity_bga + bg_severity_foca_bga,
            beta_PP_fung = bg_severity + bg_severity_foca + bg_severity_bga + bg_severity_foca_bga + bg_severity_fungicide + bg_severity_bga_fungicide + bg_severity_foca_fungicide + bg_severity_foca_bga_fungicide) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate")

# transmission rates
trans_mod_parms <- trans_mod_samps %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  ungroup() %>%
  transmute(Parameter = parameter,
            Plant_group = rep(c("annual/annual",
                                "annual/fy perennial",
                                "annual/adult perennial",
                                "fy perennial/annual",
                                "fy perennial/fy perennial",
                                "fy perennial/adult perennial",
                                "adult perennial/annual",
                                "adult perennial/fy perennial",
                                "adult perennial/adult perennial"), each = 2),
            Description = "transmission",
            Treatment = rep(c("control (water)", "fungicide"), 9),
            Estimate = case_when(estimate < 0 ~ 0,
                                 TRUE ~ estimate/30), # days between censuses
            Lower = case_when(estimate < 0 ~ 0,
                              TRUE ~ .lower/30),
            Upper = case_when(estimate < 0 ~ 0,
                              TRUE ~ .upper/30),
            Units = "day^-1^ g^-1^",
            Source = "experiment")


#### germination ####

load("output/mv_germination_fungicide_model_2018_density_exp.rda")

mv_germ_parms <- as_draws_df(mvGermD1Mod3)  %>%
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
            Treatment = if_else(parameter == "g_W", "control (water)", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper,
            Units = "NA",
            Source = "experiment")


load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

ev_germ_parms <- as_draws_df(evGermMod2)  %>%
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
            Treatment = if_else(parameter == "g_W", "control (water)", "fungicide"),
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
            Treatment = if_else(parameter == "c_W", "control (water)", "fungicide"),
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
            Treatment = if_else(parameter == "c_W", "control (water)", "fungicide"),
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
            Treatment = if_else(parameter == "c_W", "control (water)", "fungicide"),
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
            Treatment = if_else(str_detect(parameter, "W") == T, "control (water)", "fungicide"),
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
            Treatment = if_else(parameter == "l_W", "control (water)", "fungicide"),
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
            Units = "g^-1^",
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
            Units = "g^-1^",
            Source = "Benitez et al. 2021")


#### literature/estimated parameters ####

lit_parms <- tibble(Parameter = c("h", 
                                  "d", 
                                  "s_P", "s_A", "s_A",
                                  "beta_AC", "beta_FC", "beta_PC", 
                                  "v_A", "v_F", "v_P", 
                                  "b"),
                    Plant_group = c(NA_character_,
                                    NA_character_,
                                    "perennial",
                                    rep("annual", 2),
                                    rep(c("annual", "fy perennial", "adult perennial"), 2),
                                    NA_character_),
                    Description = c("inoculum addition to litter",
                                    "litter decomposition fraction",
                                    "surviving seed fraction",
                                    "surviving seed fraction",
                                    "surviving seed fraction",
                                    "litter transmission",
                                    "litter transmission",
                                    "litter transmission",
                                    "infected tissue loss",
                                    "infected tissue loss",
                                    "infected tissue loss",
                                    "inoculum loss from litter"),
                    Treatment = c(rep(NA_character_, 3), "fungicide", "control (water)", rep(NA_character_, 7)),
                    Estimate = c(5500, 0.59, 0.05, 0.15, 0, rep(1e-3, 3), rep(1e-3, 3), 0.5),
                    Source = c("Benitez et al. 2021", "DeMeester and Richter 2010", "Garrison and Stier 2010", "Redwood et al. 2018", "experiment", rep("NA", 7))) %>%
  mutate(Lower = NA_real_,
         Upper = NA_real_,
         Units = c("g^-1^", rep("NA", 4), rep ("day^-1^ g^-1^", 3), rep("day^-1^", 3), "day^-1^"))


#### combine ####
parms <- growth_mod_parms %>%
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


#### export ####
write_csv(parms, "output/model_parameters_2018_2019_density_exp.csv")

# growing season days
gs_time <- 160