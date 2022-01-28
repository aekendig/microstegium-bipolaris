##### info ####

# file: plot_severity_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 9/29/21
# goal: effects of within and outside of plot severity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)
library(gridExtra)
library(car)
library(janitor)

# import data
d1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
d2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") 
# covariate_data_processing_2018_density_exp
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") 
# temp_humidity_data_processing_2019_density_exp
bioD2Dat <- read_csv("intermediate-data/plot_biomass_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")

# logit adjustment
logit_adjust <- 0.001


#### edit data ####

# # format edge severity
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(!(month %in% c("sep"))) %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity)

# one plot in one month is missing, use nearby plots
edgeSevD2Dat3 <- tibble(site = "D1", plot = 6, treatment = "fungicide", month = "early_aug") %>%
  mutate(edge_severity = edgeSevD2Dat2 %>%
           filter(site == "D1" & plot %in% c(4, 5) & treatment == "water" & month == "early_aug") %>%
           summarise(sev = mean(edge_severity)) %>%
           pull(sev)) %>%
  full_join(edgeSevD2Dat2) %>%
  filter(!is.na(edge_severity)) # use to check for missing values

# bg severity
bgSevD1Dat <- d1Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") | # select background measurements (includes focal)
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_tot = area_tot,
         bg_sp = sp,
         bg_age = age) %>%
  select(-c(severity, healthy))

bgSevD2Dat <- d2Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") |
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_tot = biomass_tot,
         bg_sp = sp,
         bg_age = age) %>%
  select(-c(severity, healthy))

# focal dataset (includes background when background is same species)
focNextSevD1Dat <- d1Dat %>%
  filter(month != "jul") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_next_healthy = healthy) %>%
  mutate(month = dplyr::recode(month, 
                               "late_aug" = "jul", # match prior month
                               "sep" = "late_aug")) %>%
  select(-c(severity, lesions, area_tot))

focNextSevD2Dat <- d2Dat %>%
  filter(month != "may") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_next_healthy = healthy) %>%
  mutate(month = dplyr::recode(month,  # match prior month
                               "jun" = "may",
                               "jul" = "jun",
                               "early_aug" = "jul",
                               "late_aug" = "early_aug")) %>%
  select(-c(severity, lesions, biomass_tot))

# prior severity
focSevD1Dat <- d1Dat %>%
  filter(month != "sep") %>% # remove last month
  rename(foc_sp = sp,
         foc_age = age,
         foc_healthy = healthy,
         foc_lesions = lesions,
         foc_tot = area_tot) %>%
  select(-severity)

focSevD2Dat <- d2Dat %>%
  filter(month != "late_aug") %>% # remove last month
  rename(foc_sp = sp,
         foc_age = age,
         foc_healthy = healthy,
         foc_lesions = lesions,
         foc_tot = biomass_tot) %>%
  select(-severity)

# 2018 data
sevD1Dat <- focNextSevD1Dat %>%
  left_join(focSevD1Dat) %>%
  left_join(bgSevD1Dat) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         bg_sp = case_when(plot == 1 ~ foc_sp, # make the focal species the bg in 1 plots
                           TRUE ~ bg_sp),
         bg_age = case_when(plot == 1 ~ foc_age,
                            TRUE ~ bg_age),
         focal = paste(foc_sp, foc_age) %>% fct_recode("Mv" = "Mv seedling") %>% as.character(),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = paste(bg_sp, bg_age) %>% fct_recode("Mv" = "Mv seedling") %>% as.character(),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         bg_lesions = case_when(focal == background & plot == 1 ~ foc_lesions, # use focal lesions as bg in 1 plots
                                TRUE ~ bg_lesions),
         bg_tot = case_when(focal == background & plot == 1 ~ foc_tot, # use focal total as bg in 1 plots
                                TRUE ~ bg_tot),
         bg_severity = bg_lesions / bg_tot,
         foc_healthy_change = log(foc_healthy / foc_next_healthy),
         foc_healthy_change = case_when(foc_healthy == 0 & foc_next_healthy == 0 ~ log(1), # 4 plots had no lesions, no change
                                        TRUE ~ foc_healthy_change),
         foc_healthy_change = case_when(month == "jul" ~ foc_healthy_change / 2, # account for longer lapse between samples 
                                        TRUE ~ foc_healthy_change)) %>%
  group_by(focal, background) %>%
  mutate(mean_bg_sev = mean(bg_severity, na.rm = T),
         sd_bg_sev = sd(bg_severity, na.rm = T)) %>%
  ungroup() %>%
  mutate(bg_sev_s = (bg_severity - mean_bg_sev) / sd_bg_sev,
         edge_severity = 0) %>%
  filter(!is.na(foc_next_healthy) & !is.na(foc_healthy) & !is.na(bg_severity))

# check for missing data
filter(sevD1Dat, is.na(foc_next_healthy)) # 0
filter(sevD1Dat, is.na(bg_severity)) # 30 missing
filter(sevD1Dat, is.na(foc_healthy)) # 7 missing
filter(sevD1Dat, is.na(foc_healthy_change) & !is.na(foc_next_healthy) & !is.na(foc_healthy))  
filter(sevD1Dat, foc_healthy_change %in% c(Inf, -Inf))

# 2019 data
sevD2Dat <- focNextSevD2Dat %>%
  left_join(focSevD2Dat) %>%
  left_join(bgSevD2Dat) %>%
  left_join(edgeSevD2Dat3) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         bg_sp = case_when(plot == 1 ~ foc_sp, # make the focal species the bg in 1 plots
                           TRUE ~ bg_sp),
         bg_age = case_when(plot == 1 ~ foc_age,
                            TRUE ~ bg_age),
         focal = paste(foc_sp, foc_age) %>%
           fct_recode("Mv" = "Mv seedling"),
         focal = paste(foc_sp, foc_age) %>% fct_recode("Mv" = "Mv seedling") %>% as.character(),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = paste(bg_sp, bg_age) %>% fct_recode("Mv" = "Mv seedling") %>% as.character(),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         bg_lesions = case_when(focal == background & plot == 1 ~ foc_lesions, # use focal lesions as bg in 1 plots
                                TRUE ~ bg_lesions),
         bg_tot = case_when(focal == background & plot == 1 ~ foc_tot, # use focal total as bg in 1 plots
                            TRUE ~ bg_tot),
         bg_severity = bg_lesions / bg_tot,
         foc_healthy_change = log(foc_healthy / foc_next_healthy),
         foc_healthy_change = case_when(foc_healthy == 0 & foc_next_healthy == 0 ~ log(1), # 4 plots had no lesions, no change
                                        TRUE ~ foc_healthy_change)) %>%
  group_by(focal, background) %>%
  mutate(mean_bg_sev = mean(bg_severity, na.rm = T),
         sd_bg_sev = sd(bg_severity, na.rm = T)) %>%
  ungroup() %>%
  mutate(bg_sev_s = (bg_severity - mean_bg_sev) / sd_bg_sev) %>%
  filter(!is.na(foc_next_healthy) & !is.na(foc_healthy) & !is.na(bg_severity))

# check for missing data
filter(sevD2Dat, is.na(foc_next_healthy))
filter(sevD2Dat, is.na(bg_severity)) # 189 missing
filter(sevD2Dat, is.na(foc_healthy)) # 185 missing
filter(sevD2Dat, is.na(foc_healthy_change) & !is.na(foc_next_healthy) & !is.na(foc_healthy)) 
filter(sevD2Dat, foc_healthy_change %in% c(Inf, -Inf))

# environmental data
# make month the final month instead of initial
envD1Dat2 <- envD1Dat %>%
  inner_join(sevD1Dat %>%
               mutate(month = fct_recode(month, sep = "late_aug", late_aug = "jul")))

envD2Dat2 <- envD2Dat %>%
  filter(month %in% c("early_aug", "late_aug")) %>%
  mutate(dew_c = (dew_intensity2 - mean(dew_intensity2, na.rm = T)) / sd(dew_intensity2, na.rm = T)) %>%
  inner_join(sevD2Dat %>%
               mutate(month = fct_recode(month, 
                                         late_aug = "early_aug", early_aug = "jul", jul = "jun", jun  = "may")))

envBioD2Dat <- envD2Dat  %>%
  filter(month %in% c("early_aug", "late_aug")) %>%
  mutate(dew_c = (dew_intensity2 - mean(dew_intensity2, na.rm = T)) / sd(dew_intensity2, na.rm = T)) %>%
  left_join(bioD2Dat %>%
              mutate(biomass_tot = biomass_bg + biomass_foc_mv + biomass_foc_evS + biomass_foc_evA,
                     biomass_c = (biomass_tot - mean(biomass_tot, na.rm = T)) / sd(biomass_tot, na.rm = T)) %>%
              select(site, treatment, plot, biomass_tot, biomass_c)) %>%
  mutate(plotf = paste0(site, plot, substr(treatment, 1, 1)))


#### visualizations ####

# healthy change distributions
ggplot(sevD1Dat, aes(foc_healthy_change)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(foc_healthy_change)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

# bg severity distributions
ggplot(sevD1Dat, aes(bg_severity)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD1Dat, aes(bg_sev_s)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(bg_severity)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(bg_sev_s)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

# pairwise combinations
ggplot(sevD1Dat, aes(bg_severity, foc_healthy_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(focal ~ background, scales = "free")

ggplot(sevD2Dat, aes(bg_severity, foc_healthy_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(focal ~ background, scales = "free")
# looked at high Ev adult plot (bg_les_s > 7) - no Ev scans for late August

# edge
ggplot(sevD2Dat, aes(edge_severity)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(edge_severity, foc_healthy_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(~ focal, scales = "free")

# edge and background
ggplot(edgeSevD2Dat2, aes(x = treatment, y = edge_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ month)

ggplot(filter(bgSevD2Dat, bg_sp == "Mv"), aes(x = treatment, y = bg_lesions/bg_tot)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ month)


#### fit models ####

# Mv
sevD1Mod <- brm(foc_healthy_change ~ bg_severity * foc * bg * fungicide + (1|plotf),
                   data = sevD1Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99))
# use plot instead of site for random effects because of multiple temporal sampling
mod_check_fun(sevD1Mod)

sevD2Mod <- brm(foc_healthy_change ~ bg_severity * foc * bg * fungicide + edge_severity + edge_severity:foc + edge_severity:fungicide + edge_severity:foc:fungicide + (1|plotf),
                data = sevD2Dat, family = gaussian,
                prior <- c(prior(normal(0, 1), class = "Intercept"),
                           prior(normal(0, 1), class = "b")), # use default for sigma
                iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                control = list(adapt_delta = 0.99))
mod_check_fun(sevD2Mod)

# save models
save(sevD1Mod, file = "output/plot_transmission_model_2018_density_exp.rda")
save(sevD2Mod, file = "output/plot_transmission_model_2019_density_exp.rda")


#### intraspecific vs. interspecific ####

# common terms on both sides of = were deleted
# inter listed first in name
# intra on left side of =

# Mv 
evS_mv_ctrl_hyp = "-bg_severity:focs = 0"
evS_mv_fung_hyp = "-bg_severity:focs - bg_severity:focs:fungicide = 0"
evA_mv_ctrl_hyp = "-bg_severity:foca = 0"
evA_mv_fung_hyp = "-bg_severity:foca - bg_severity:foca:fungicide = 0"

# EvS 
mv_evS_ctrl_hyp = "0.5 * (bg_severity:focs + bg_severity:focs:bgs + bg_severity:foca + bg_severity:foca:bgs) = 0"
mv_evS_fung_hyp = "0.5 * (bg_severity:focs + bg_severity:focs:bgs + bg_severity:focs:fungicide + bg_severity:focs:bgs:fungicide + bg_severity:foca + bg_severity:foca:bgs + bg_severity:foca:fungicide + bg_severity:foca:bgs:fungicide) = 0"

# EvA 
mv_evA_ctrl_hyp = "0.5 * (bg_severity:focs + bg_severity:focs:bga + bg_severity:foca + bg_severity:foca:bga) = 0"
mv_evA_fung_hyp = "0.5 * (bg_severity:focs + bg_severity:focs:bga + bg_severity:focs:fungicide + bg_severity:focs:bga:fungicide + bg_severity:foca + bg_severity:foca:bga + bg_severity:foca:fungicide + bg_severity:foca:bga:fungicide) = 0"

sevD1hyps <- hypothesis(sevD1Mod, 
                        c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
                          mv_evS_ctrl_hyp, mv_evS_fung_hyp, mv_evA_ctrl_hyp, mv_evA_fung_hyp))

sevD2hyps <- hypothesis(sevD2Mod, 
                        c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
                          mv_evS_ctrl_hyp, mv_evS_fung_hyp, mv_evA_ctrl_hyp, mv_evA_fung_hyp))
# none are significantly different from zero

write_csv(sevD1hyps[[1]], "output/plot_transmission_intra_vs_inter_2018_density_exp.csv")
write_csv(sevD2hyps[[1]], "output/plot_transmission_intra_vs_inter_2019_density_exp.csv")


#### baseline infection ####

mv_ctrl_base <- "Intercept = 0"
evS_ctrl_base <- "Intercept + focs = 0"
evA_ctrl_base <- "Intercept + foca = 0"

mv_fung_base <- "Intercept + fungicide = 0"
evS_fung_base <- "Intercept + focs + fungicide + focs:fungicide = 0"
evA_fung_base <- "Intercept + foca + fungicide + foca:fungicide = 0"

sevD1base <- hypothesis(sevD1Mod,
                        c(mv_ctrl_base, mv_fung_base, evS_ctrl_base, evS_fung_base, evA_ctrl_base, evA_fung_base))[[1]] %>%
  as_tibble() %>%
  mutate(treatment = rep(c("control", "fungicide"), 3),
         focal = rep(c("Mv", "Ev first-year", "Ev adult"), each = 2),
         background = "background") %>%
  select(focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper) %>%
  arrange(focal, treatment)

sevD2base <- hypothesis(sevD2Mod,
                        c(mv_ctrl_base, mv_fung_base, evS_ctrl_base, evS_fung_base, evA_ctrl_base, evA_fung_base))[[1]] %>%
  as_tibble() %>%
  mutate(treatment = rep(c("control", "fungicide"), 3),
         focal = rep(c("Mv", "Ev first-year", "Ev adult"), each = 2),
         background = "background") %>%
  select(focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper) %>%
  arrange(focal, treatment)


#### transmission coefficients (betas) ####

# are betas different than 0? (_beta)
# does fungicide treatment affect betas? (_trt_eff, fung_beta - ctrl_beta)

# Mv background
mv_mv_ctrl_beta = "bg_severity = 0"
mv_mv_fung_beta = "bg_severity + bg_severity:fungicide = 0"
mv_mv_trt_eff = "bg_severity:fungicide = 0"
evS_mv_ctrl_beta = "bg_severity + bg_severity:focs = 0"
evS_mv_fung_beta = "bg_severity + bg_severity:fungicide + bg_severity:focs + bg_severity:focs:fungicide = 0"
evS_mv_trt_eff = "bg_severity:fungicide + bg_severity:focs:fungicide = 0"
evA_mv_ctrl_beta = "bg_severity + bg_severity:foca = 0"
evA_mv_fung_beta = "bg_severity + bg_severity:fungicide + bg_severity:foca + bg_severity:foca:fungicide = 0"
evA_mv_trt_eff = "bg_severity:fungicide + bg_severity:foca:fungicide = 0"

# EvS background
evS_evS_ctrl_beta = "bg_severity + bg_severity:focs + bg_severity:bgs + bg_severity:focs:bgs = 0"
evS_evS_fung_beta = "bg_severity + bg_severity:focs + bg_severity:bgs + bg_severity:focs:bgs + bg_severity:fungicide + bg_severity:focs:fungicide + bg_severity:bgs:fungicide + bg_severity:focs:bgs:fungicide = 0"
evS_evS_trt_eff = "bg_severity:fungicide + bg_severity:focs:fungicide + bg_severity:bgs:fungicide + bg_severity:focs:bgs:fungicide = 0"
mv_evS_ctrl_beta = "bg_severity +  bg_severity:bgs = 0"
mv_evS_fung_beta = "bg_severity +  bg_severity:bgs + bg_severity:fungicide + bg_severity:bgs:fungicide = 0"
mv_evS_trt_eff = "bg_severity:fungicide + bg_severity:bgs:fungicide = 0"
evA_evS_ctrl_beta = "bg_severity +  bg_severity:bgs + bg_severity:foca + bg_severity:foca:bgs = 0"
evA_evS_fung_beta = "bg_severity +  bg_severity:bgs + bg_severity:foca + bg_severity:foca:bgs + bg_severity:fungicide +  bg_severity:bgs:fungicide + bg_severity:foca:fungicide + bg_severity:foca:bgs:fungicide = 0"
evA_evS_trt_eff = "bg_severity:fungicide +  bg_severity:bgs:fungicide + bg_severity:foca:fungicide + bg_severity:foca:bgs:fungicide = 0"

# EvA background
evA_evA_ctrl_beta = "bg_severity + bg_severity:foca + bg_severity:bga + bg_severity:foca:bga = 0"
evA_evA_fung_beta = "bg_severity + bg_severity:foca + bg_severity:bga + bg_severity:foca:bga + bg_severity:fungicide + bg_severity:foca:fungicide + bg_severity:bga:fungicide + bg_severity:foca:bga:fungicide = 0"
evA_evA_trt_eff = "bg_severity:fungicide + bg_severity:foca:fungicide + bg_severity:bga:fungicide + bg_severity:foca:bga:fungicide = 0"
mv_evA_ctrl_beta = "bg_severity +  bg_severity:bga = 0"
mv_evA_fung_beta = "bg_severity +  bg_severity:bga + bg_severity:fungicide + bg_severity:bga:fungicide = 0"
mv_evA_trt_eff = "bg_severity:fungicide + bg_severity:bga:fungicide = 0"
evS_evA_ctrl_beta = "bg_severity + bg_severity:bga + bg_severity:focs + bg_severity:focs:bga = 0"
evS_evA_fung_beta = "bg_severity + bg_severity:bga + bg_severity:focs + bg_severity:focs:bga + bg_severity:fungicide + bg_severity:bga:fungicide + bg_severity:focs:fungicide + bg_severity:focs:bga:fungicide = 0"
evS_evA_trt_eff = "bg_severity:fungicide + bg_severity:bga:fungicide + bg_severity:focs:fungicide + bg_severity:focs:bga:fungicide = 0"

sevD1betas <- hypothesis(sevD1Mod, 
                         c(mv_mv_ctrl_beta, mv_mv_fung_beta, 
                           evS_mv_ctrl_beta, evS_mv_fung_beta, 
                           evA_mv_ctrl_beta, evA_mv_fung_beta,
                           evS_evS_ctrl_beta, evS_evS_fung_beta,
                           mv_evS_ctrl_beta, mv_evS_fung_beta, 
                           evA_evS_ctrl_beta, evA_evS_fung_beta, 
                           evA_evA_ctrl_beta, evA_evA_fung_beta,
                           mv_evA_ctrl_beta, mv_evA_fung_beta, 
                           evS_evA_ctrl_beta, evS_evA_fung_beta))

sevD1TrtEff <- hypothesis(sevD1Mod, 
                         c(mv_mv_trt_eff, 
                           evS_mv_trt_eff, 
                           evA_mv_trt_eff,
                           evS_evS_trt_eff,
                           mv_evS_trt_eff, 
                           evA_evS_trt_eff, 
                           evA_evA_trt_eff,
                           mv_evA_trt_eff, 
                           evS_evA_trt_eff)) # none sig

sevD2betas <- hypothesis(sevD2Mod, 
                         c(mv_mv_ctrl_beta, mv_mv_fung_beta, 
                           evS_mv_ctrl_beta, evS_mv_fung_beta, 
                           evA_mv_ctrl_beta, evA_mv_fung_beta,
                           evS_evS_ctrl_beta, evS_evS_fung_beta,
                           mv_evS_ctrl_beta, mv_evS_fung_beta, 
                           evA_evS_ctrl_beta, evA_evS_fung_beta, 
                           evA_evA_ctrl_beta, evA_evA_fung_beta,
                           mv_evA_ctrl_beta, mv_evA_fung_beta, 
                           evS_evA_ctrl_beta, evS_evA_fung_beta))

sevD2TrtEff <- hypothesis(sevD2Mod, 
                          c(mv_mv_trt_eff, 
                            evS_mv_trt_eff, 
                            evA_mv_trt_eff,
                            evS_evS_trt_eff, # fung sig increase
                            mv_evS_trt_eff, 
                            evA_evS_trt_eff, 
                            evA_evA_trt_eff,
                            mv_evA_trt_eff, 
                            evS_evA_trt_eff))

# combine betas
betaDat <- sevD1betas[[1]] %>%
  mutate(year = "2018",
         foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                        "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                        "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                        "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                        "s_a_ctrl", "s_a_fung")) %>%
  full_join(sevD2betas[[1]] %>%
              mutate(year = "2019",
                     foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                                    "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                                    "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                                    "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                                    "s_a_ctrl", "s_a_fung"))) %>%
  select(-Hypothesis) %>%
  rowwise() %>%
  mutate(foc = str_split(foc_bg_trt, "_")[[1]][1],
         bg = str_split(foc_bg_trt, "_")[[1]][2],
         trt = str_split(foc_bg_trt, "_")[[1]][3]) %>%
  ungroup() %>%
  mutate(treatment = fct_recode(trt, "control (water)" = "ctrl",
                                "fungicide" = "fung"),
         sig = case_when((CI.Lower < 0 & CI.Upper < 0) | (CI.Lower > 0 & CI.Upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"))

# combine treatment effects
trtEffDat <- sevD1TrtEff[[1]] %>%
  mutate(year = "2018",
         foc_bg = c("m_m", "s_m", "a_m", 
                    "s_s", "m_s", "a_s",
                    "a_a", "m_a","s_a")) %>%
  full_join(sevD2TrtEff[[1]] %>%
              mutate(year = "2019",
                     foc_bg = c("m_m", "s_m", "a_m", 
                                "s_s", "m_s", "a_s",
                                "a_a", "m_a","s_a"))) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob, Star)) %>%
  rowwise() %>%
  mutate(foc = str_split(foc_bg, "_")[[1]][1],
         bg = str_split(foc_bg, "_")[[1]][2]) %>%
  ungroup()

# edit to save
betaD1DatSave <- bind_rows(sevD1base, betaDat %>%
                             filter(year == 2018) %>%
  left_join(sevD2Dat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(treatment = fct_recode(treatment, "control" = "control (water)")) %>%
  select(focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper) %>%
  arrange(background, focal, treatment))

betaD2DatSave <- bind_rows(sevD2base, betaDat %>%
                             filter(year == 2019) %>%
                             left_join(sevD2Dat %>%
                                         select(foc, focal, bg, background) %>%
                                         unique()) %>% 
                             mutate(treatment = fct_recode(treatment, "control" = "control (water)")) %>%
                             select(focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper) %>%
                             arrange(background, focal, treatment))

trtEffDatSave <- trtEffDat %>%
  left_join(sevD2Dat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  select(year, focal, background, Estimate, Est.Error, CI.Lower, CI.Upper) %>%
  arrange(year, background, focal)

# save
write_csv(betaD1DatSave, "output/plot_transmission_coefficients_2018_density_exp.csv")
write_csv(trtEffDatSave, "output/plot_transmission_fungicide_effects_2018_2019_density_exp.csv")


#### edge effects ####

mv_ctrl_edge <- "edge_severity = 0"
mv_fung_edge <- "edge_severity + fungicide:edge_severity = 0"
evS_ctrl_edge <- "edge_severity + focs:edge_severity = 0"
evS_fung_edge <- "edge_severity + focs:edge_severity + fungicide:edge_severity + focs:fungicide:edge_severity = 0"
evA_ctrl_edge <- "edge_severity + foca:edge_severity = 0"
evA_fung_edge <- "edge_severity + foca:edge_severity + fungicide:edge_severity + foca:fungicide:edge_severity = 0"

edgeD2hyps <- hypothesis(sevD2Mod, 
                        c(mv_ctrl_edge, mv_fung_edge,
                          evS_ctrl_edge, evS_fung_edge,
                          evA_ctrl_edge, evA_fung_edge))
# all are significantly different than zero

edgeD2hyps2 <- edgeD2hyps[[1]] %>%
  as_tibble() %>%
  mutate(treatment = rep(c("control (water)", "fungicide"), 3),
         foc = rep(c("m", "s", "a"), each = 2),
         sig = case_when((CI.Lower < 0 & CI.Upper < 0) | (CI.Lower > 0 & CI.Upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"))

# add to betas to save
betaD2DatSave2 <- rbind(betaD2DatSave, edgeD2hyps2 %>%
                          left_join(sevD2Dat %>%
                                      select(foc, focal) %>%
                                      unique()) %>%
                          mutate(treatment = fct_recode(treatment, "control" = "control (water)"),
                                 background = "Mv edge") %>%
                          select(focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper) %>%
                          arrange(focal, treatment))
  
# save
write_csv(betaD2DatSave2, "output/plot_transmission_coefficients_2019_density_exp.csv")


#### figure ####

# bg_severity function
bg_sev_fun <- function(f, b, trt, yr){
  
  if(yr == "2018"){
    dat <- sevD1Dat %>% filter(foc == f & bg == b & treatment == trt)
  }else{
    dat <- sevD2Dat %>% filter(foc == f & bg == b & treatment == trt)
  }
  
  bg_severity = seq(min(dat$bg_severity), max(dat$bg_severity), length.out = 100)
  
  return(bg_severity)
}

# predicted data
predD1Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2018",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(bg_severity = pmap(list(foc, bg, treatment, year), bg_sev_fun)) %>%
  unnest(bg_severity) %>%
  mutate(foc_healthy_change = fitted(sevD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(sevD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(sevD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predD2Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2019",
         plotf = "A",
         edge_severity = 0,
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(bg_severity = pmap(list(foc, bg, treatment, year), bg_sev_fun)) %>%
  unnest(bg_severity) %>%
  mutate(foc_healthy_change = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predDat <- predD1Dat %>%
  full_join(predD2Dat) %>%
  mutate(focal = fct_recode(foc, "Invader (Mv)" = "m", "1st yr competitor (Ev)" = "s", "Adult competitor (Ev)" = "a") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         background = fct_recode(bg, "Invader (Mv)" = "m", "1st yr competitor (Ev)" = "s", "Adult competitor (Ev)" = "a") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev(),
         bg_severity = bg_severity * 100)

# edge severity
sevD2Dat %>%
  group_by(treatment) %>%
  summarise(minEdge = min(edge_severity),
            maxEdge = max(edge_severity))

predD2EdgeDat <- tibble(foc = c("a", "s", "m"),
                        focal = c("Adult competitor (Ev)", "1st yr competitor (Ev)", "Invader (Mv)")) %>%
  expand_grid(tibble(treatment = rep(c("control (water)", "fungicide"), each = 100),
                     edge_severity = c(seq(0, 0.0856, length.out = 100),
                                       seq(0, 0.0639, length.out = 100)))) %>%
  mutate(year = "2019",
         plotf = "A",
         bg = "m",
         background = "Mv",
         bg_severity = 0,
         fungicide = case_when(treatment == "control (water)" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(foc_healthy_change = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"],
         edge_severity = edge_severity * 100,
         treatment = fct_rev(treatment),
         focal = fct_relevel(focal, "Invader (Mv)", "1st yr competitor (Ev)"))

# combine with preddat
predDat2 <- predDat %>%
  left_join(betaDat) %>%
  mutate(intra = case_when(foc == bg ~ "yes",
                           foc %in% c("a", "s") & bg %in% c("a", "s") ~ "yes",
                           TRUE ~ "no"))

# raw data
figDat <- sevD1Dat %>%
  select(month, site, plot, treatment, focal, background, bg_severity, edge_severity, foc_healthy_change) %>%
  mutate(year = "2018") %>%
  full_join(sevD2Dat %>%
              select(month, site, plot, treatment, focal, background, bg_severity, edge_severity, foc_healthy_change) %>%
              mutate(year = "2019")) %>%
  left_join(plotsD %>%
              select(plot, density_level) %>%
              unique()) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"),
         focal = fct_recode(focal, "Invader (Mv)" = "Mv",
                            "Adult competitor (Ev)" = "Ev adult",
                            "1st yr competitor (Ev)" = "Ev seedling") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         background = fct_recode(background, "Invader (Mv)" = "Mv",
                                 "Adult competitor (Ev)" = "Ev adult",
                                 "1st yr competitor (Ev)" = "Ev seedling") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         bg_severity = bg_severity * 100,
         edge_severity = edge_severity * 100,
         density_level = fct_relevel(density_level, "low", "medium"),
         intra = if_else(focal == background, "yes", "no"))

# betas for figure
# sig beta values
betaDat2 <- betaDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, bg, background) %>%
              summarise(bg_severity = max(bg_severity))) %>%
  left_join(predDat2 %>%
              select(foc, focal) %>% # add focal names
              unique()) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(upper = max(upper)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(dat = max(foc_healthy_change))) %>%
              rowwise() %>%
              mutate(foc_healthy_change = max(c(dat, upper))) %>%
              ungroup() %>%
              select(-c(dat, upper))) %>%
  rowwise() %>%
  mutate(parm = as.character(as.expression(substitute(beta~"="~est,
                                                      list(est = format(round_half_up(Estimate, 1), nsmall = 1)))))) %>%
  ungroup()

# add sig to pred data
predD2EdgeDat2 <- predD2EdgeDat %>%
  left_join(edgeD2hyps2 %>%
              select(foc, treatment, sig)) %>%
  mutate(intra = if_else(foc == "m", "yes", "no"))

# edge coefficients for figure
edgeD2hyps3 <- edgeD2hyps2 %>%
  filter(sig == "omits 0") %>%
  left_join(predD2EdgeDat %>%
              group_by(foc, focal) %>%
              summarise(edge_severity = max(edge_severity),
                        upper = max(upper),
                        lower = min(lower)) %>%
              ungroup()) %>%
  left_join(sevD2Dat %>%
              group_by(foc) %>%
              summarise(dat_max = max(foc_healthy_change)) %>%
              ungroup()) %>%
  rowwise() %>%
  mutate(foc_healthy_change = max(c(dat_max, upper))) %>%
  ungroup() %>%
  select(-c(dat_max, upper, lower)) %>%
  rowwise() %>%
  mutate(parm = as.character(as.expression(substitute(beta~"="~est,
                                                      list(est = format(round_half_up(Estimate, 1), nsmall = 1)))))) %>%
  ungroup() %>%
  mutate(foc_healthy_change = case_when(treatment == "fungicide" & foc == "m" ~ foc_healthy_change - 0.25,
                                        treatment == "fungicide" & foc == "a" ~ foc_healthy_change - 0.75,
                                        treatment == "fungicide" & foc == "s" ~ foc_healthy_change - 0.4,
                                        TRUE ~ foc_healthy_change),
         vjust = 1,
         hjust = 1,
         background = "Mv")
         
# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

textSize = 2.5

box_shade = "gray92"

# yearText <- tibble(year = c("2018", "2019"),
#                    bg_severity = c(0, 0),
#                    foc_healthy_change = c(2, 3),
#                    background = "Ev adult",
#                    focal = "Ev adult",
#                    treatment = "fungicide")



# legend
legFig <- ggplot(predDat2, aes(x = bg_severity, y = foc_healthy_change)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = figDat, aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.5) +
  facet_grid(rows = vars(focal),
             cols = vars(background)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density treatment") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical") +
  guides(shape = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

leg <- get_legend(legFig)

# 2018 figure
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = bg_severity, y = foc_healthy_change)) +
  geom_rect(aes(fill = intra), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, size = 0.15) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  # geom_text(data = filter(betaDat2, year == "2018"), 
  #           aes(label = paste("beta", " == ", parm, sep = ""), 
  #               color = treatment), parse = T, hjust = 1, vjust = 0, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c(col_pal, "white", box_shade)) +
  xlab("Initial disease severity (%)") +
  ylab("Change in infected tissue") +
  fig_theme

pdf("output/plot_transmission_pairwise_figure_2018_density_exp.pdf", width = 3.94, height = 4.33)
plot_grid(pairD1Fig, leg,
          nrow = 2, 
          rel_heights = c(1, 0.15))
dev.off()

# 2019 figure
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = bg_severity, y = foc_healthy_change)) +
  geom_rect(aes(fill = intra), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, size = 0.15) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(betaDat2, year == "2019"),
            aes(label = parm,
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c(col_pal, "white", box_shade)) +
  xlab("Initial disease severity (%)") +
  ylab("Change in infected tissue") +
  fig_theme

edgeD2Fig <- ggplot(predD2EdgeDat2, aes(x = edge_severity, y = foc_healthy_change)) +
  geom_rect(aes(fill = intra), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, size = 0.15) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = filter(figDat, year == "2019") %>% 
               mutate(background = "Mv"), 
             aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = edgeD2hyps3, 
            aes(label = parm, 
                color = treatment,
                hjust = hjust, vjust = vjust), 
            parse = T, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             #cols = vars(background),
             scales = "free",
             switch = "both") +
  # scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c(col_pal, "white", "gray85")) +
  xlab("Disease severity\nof invader (Mv)\nsurrounding plots (%)") +
  fig_theme +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_text(size = 7, margin = margin(t = 7, r = 0, b = 0, l = 0)))

# combine plots
combFig <- plot_grid(pairD2Fig, 
                     edgeD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     label_size = 10,
                     rel_widths = c(1, 0.33),
                     hjust = c(-0.5, 0.3))

# combine
pdf("output/plot_transmission_pairwise_figure_2019_density_exp.pdf", width = 5.12, height = 4.33)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(1, 0.15))
dev.off()


#### environmental models ####

# 2018
envSevD1Mod <- brm(foc_healthy_change ~ foc * fungicide * (soil_moisture_jun.prop + canopy_cover.prop) + (1|plotf),
                   data = envD1Dat2, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(envSevD1Mod)

save(envSevD1Mod, file = "output/environmental_infection_change_model_2018_density_exp.rda")

# 2019
ggplot(envD2Dat2, aes(x = dew_intensity2)) +
  geom_histogram()
ggplot(envD2Dat2, aes(x = dew_c)) +
  geom_histogram()

envSevD2Mod <- brm(foc_healthy_change ~ foc * dew_c + (1|plotf),
                   data = envD2Dat2, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(envSevD2Mod)

save(envSevD2Mod, file = "output/environmental_infection_change_model_2019_density_exp.rda")

# coefficients
mv_ctrl_soil_hyp <- "soil_moisture_jun.prop = 0"
mv_fung_soil_hyp <- "soil_moisture_jun.prop + fungicide:soil_moisture_jun.prop = 0"
evS_ctrl_soil_hyp <- "soil_moisture_jun.prop + focs:soil_moisture_jun.prop = 0"
evS_fung_soil_hyp <- "soil_moisture_jun.prop + fungicide:soil_moisture_jun.prop + focs:soil_moisture_jun.prop + focs:fungicide:soil_moisture_jun.prop = 0"
evA_ctrl_soil_hyp <- "soil_moisture_jun.prop + foca:soil_moisture_jun.prop = 0"
evA_fung_soil_hyp <- "soil_moisture_jun.prop + fungicide:soil_moisture_jun.prop + foca:soil_moisture_jun.prop + foca:fungicide:soil_moisture_jun.prop = 0"

mv_ctrl_canopy_hyp <- "canopy_cover.prop = 0"
mv_fung_canopy_hyp <- "canopy_cover.prop + fungicide:canopy_cover.prop = 0"
evS_ctrl_canopy_hyp <- "canopy_cover.prop + focs:canopy_cover.prop = 0"
evS_fung_canopy_hyp <- "canopy_cover.prop + fungicide:canopy_cover.prop + focs:canopy_cover.prop + focs:fungicide:canopy_cover.prop = 0"
evA_ctrl_canopy_hyp <- "canopy_cover.prop + foca:canopy_cover.prop = 0"
evA_fung_canopy_hyp <- "canopy_cover.prop + fungicide:canopy_cover.prop + foca:canopy_cover.prop + foca:fungicide:canopy_cover.prop = 0"

envD1Hyps <- hypothesis(envSevD1Mod, c(mv_ctrl_soil_hyp, mv_fung_soil_hyp, mv_ctrl_canopy_hyp, mv_fung_canopy_hyp,
                          evS_ctrl_soil_hyp, evS_fung_soil_hyp, evS_ctrl_canopy_hyp, evS_fung_canopy_hyp,
                          evA_ctrl_soil_hyp, evA_fung_soil_hyp, evA_ctrl_canopy_hyp, evA_fung_canopy_hyp))[[1]] %>%
  as_tibble() %>%
  mutate(treatment = rep(c("control", "fungicide"), 6),
         foc = rep(c("Mv", "Ev first-year", "Ev adult"), each = 4),
         env_var = rep(rep(c("soil moisture", "canopy cover"), each = 2), 3)) %>%
  select(foc, env_var, treatment, Estimate:CI.Upper) %>%
  arrange(env_var, foc, treatment)
# no significant effects

write_csv(envD1Hyps, "output/environmental_infection_change_coefficients_2018_density_exp.csv")

mv_ctrl_dew_hyp <- "dew_c = 0"
evS_ctrl_dew_hyp <- "dew_c + focs:dew_c = 0"
evA_ctrl_dew_hyp <- "dew_c + foca:dew_c = 0"

envD2Hyps <- hypothesis(envSevD2Mod, c(mv_ctrl_dew_hyp, evS_ctrl_dew_hyp, evA_ctrl_dew_hyp))[[1]] %>%
  as_tibble() %>%
  mutate(foc = c("Mv", "Ev first-year", "Ev adult")) %>%
  select(foc, Estimate:CI.Upper)

# Ev adult and dew intensity
ggplot(filter(envD2Dat2, foc == "a"), aes(dew_c, foc_healthy_change)) +
  geom_point()

# remove high Ev adult value and refit model
envD2Dat3 <- envD2Dat2 %>%
  filter(!(foc == "a" & foc_healthy_change > 3))

envSevD2Mod2 <- brm(foc_healthy_change ~ foc * dew_c + (1|plotf),
                   data = envD2Dat3, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)

save(envSevD2Mod2, file = "output/environmental_infection_change_model2_2019_density_exp.rda")

envD2Hyps2 <- envD2Hyps %>%
  full_join(hypothesis(envSevD2Mod2, evA_ctrl_dew_hyp)[[1]] %>%
              as_tibble() %>%
              mutate(foc = "Ev adult without highest value") %>%
              select(-c(Hypothesis, Evid.Ratio, Post.Prob, Star))) %>%
  arrange(foc)
# sig effect is lost

write_csv(envD2Hyps2, "output/environmental_infection_change_coefficients_2019_density_exp.csv")

# biomass effect on dew intensity
ggplot(envBioD2Dat, aes(x = biomass_tot)) +
  geom_histogram()

ggplot(envBioD2Dat, aes(x = biomass_c)) +
  geom_histogram()

envBioD2Mod <- brm(dew_c ~ biomass_c * month + (1|plotf),
                   data = envBioD2Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)

mod_check_fun(envBioD2Mod)

save(envBioD2Mod, file = "output/dew_biomass_model_2019_density_exp.rda")

hypothesis(envBioD2Mod, c("biomass_c = 0", "biomass_c + biomass_c:monthlate_aug = 0"))
# no sig effect


#### water (control) effect - compare to edge ####

# combine datasets
ctrlD2Dat <- edgeSevD2Dat2 %>%
  filter(plot %in% c(2:4) & month != "may") %>%
  mutate(plant_type = "edge (untreated)") %>%
  rename(severity = "edge_severity") %>%
  full_join(bgSevD2Dat %>%
              filter(bg_sp == "Mv") %>%
              mutate(severity = bg_lesions/bg_tot,
                     plant_type = "experimental")) %>%
  mutate(month = fct_recode(month, 
                            "early August" = "early_aug",
                            "July" = "jul",
                            "June" = "jun",
                            "late August" = "late_aug") %>%
           fct_relevel("June", "July", "early August"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         severity2 = case_when(severity == 0 ~ 1e-6, # lowest value is 2.4e-6
                               TRUE ~ severity))

# sample sizes
ctrlD2Dat %>%
  group_by(plant_type, site, treatment) %>%
  count() %>%
  data.frame()

ctrlD2Dat %>%
  group_by(plotf, month) %>%
  count() %>%
  filter(n != 2)

# model
ctrlD2Mod <- brm(severity2 ~ plant_type*treatment*month + (1|plotf),
                 data = ctrlD2Dat, family = Beta,
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3)

mod_check_fun(ctrlD2Mod)

# save model
save(ctrlD2Mod, file = "output/edge_experiment_mv_model_2019_density_exp.rda")

# logit to prob
logit2prob <- function(x){
  exp(x)/(1 + exp(x))
}

# do experimental plants have different severity than edge?
jun_ctrl_comp = "logit2prob(Intercept + plant_typeexperimental)*100 = logit2prob(Intercept)*100"
jun_fung_comp = "logit2prob(Intercept + plant_typeexperimental + treatmentfungicide + plant_typeexperimental:treatmentfungicide)*100 = logit2prob(Intercept + treatmentfungicide)*100"
jul_ctrl_comp = "logit2prob(Intercept + plant_typeexperimental + monthJuly + plant_typeexperimental:monthJuly)*100 = logit2prob(Intercept + monthJuly)*100"
jul_fung_comp = "logit2prob(Intercept + plant_typeexperimental + treatmentfungicide + monthJuly + plant_typeexperimental:treatmentfungicide + plant_typeexperimental:monthJuly + treatmentfungicide:monthJuly + plant_typeexperimental:treatmentfungicide:monthJuly)*100 = logit2prob(Intercept + treatmentfungicide + monthJuly + treatmentfungicide:monthJuly)*100"
eau_ctrl_comp = "logit2prob(Intercept + plant_typeexperimental + monthearlyAugust + plant_typeexperimental:monthearlyAugust)*100 = logit2prob(Intercept + monthearlyAugust)*100"
eau_fung_comp = "logit2prob(Intercept + plant_typeexperimental + treatmentfungicide + monthearlyAugust + plant_typeexperimental:treatmentfungicide + plant_typeexperimental:monthearlyAugust + treatmentfungicide:monthearlyAugust + plant_typeexperimental:treatmentfungicide:monthearlyAugust)*100 = logit2prob(Intercept + treatmentfungicide + monthearlyAugust + treatmentfungicide:monthearlyAugust)*100"
lau_ctrl_comp = "logit2prob(Intercept + plant_typeexperimental + monthlateAugust + plant_typeexperimental:monthlateAugust)*100 = logit2prob(Intercept + monthlateAugust)*100"
lau_fung_comp = "logit2prob(Intercept + plant_typeexperimental + treatmentfungicide + monthlateAugust + plant_typeexperimental:treatmentfungicide + plant_typeexperimental:monthlateAugust + treatmentfungicide:monthlateAugust + plant_typeexperimental:treatmentfungicide:monthlateAugust)*100 = logit2prob(Intercept + treatmentfungicide + monthlateAugust + treatmentfungicide:monthlateAugust)*100"

# estimates
ctrlD2hyps <- hypothesis(ctrlD2Mod, c(jun_ctrl_comp, jun_fung_comp, jul_ctrl_comp, jul_fung_comp,
                                      eau_ctrl_comp, eau_fung_comp, lau_ctrl_comp, lau_fung_comp))

ctrlD2hypsOut <- ctrlD2hyps[[1]] %>%
  as_tibble() %>%
  mutate(Month = rep(c("June", "July", "early August", "late August"), each = 2),
         treatment = rep(c("control (water)", "fungicide"), 4))

write_csv(ctrlD2hypsOut, "output/edge_experiment_mv_2019_density_exp.csv")

# figure
pdf("output/edge_experiment_mv_figure_2019_density_exp.pdf", width = 2.5, height = 8)
ggplot(ctrlD2Dat, aes(x = treatment, y = severity2*100, color = plant_type)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 2.5, position = position_dodge(0.2)) +
  facet_wrap(~ month, scales = "free_y", ncol = 1) +
  labs(x = "Plot treatment", y = "Disease severity (%)") +
  scale_color_manual(values = c("black", "#20A387FF"), name = "Plant type") +
  fig_theme +
  theme(legend.position = c(0.7, 0.94))
dev.off()


#### treatment effect on severity ####

# initial visualization
ggplot(d2Dat, aes(x = month, y = severity, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ paste0(sp, age))
# most disease in late august
# fungicide effect strongest for Ev in early August

# edit dataset
trtD2Dat <- d2Dat %>%
  filter((month == "late_aug" & sp == "Mv") |
           (month == "early_aug" & sp == "Ev")) %>%
  mutate(treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control"),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         plant_type = paste(sp, age, sep = "_"),
         fungicide = if_else(treatment == "fungicide", 1, 0),
         severity2 = severity + 1e-6) # lowest value is 2e-4

# min severity
min(filter(trtD2Dat, severity > 0)$severity)
max(trtD2Dat$severity)

# min severity
min(filter(d2Dat, severity > 0)$severity)

# model
trtD2Mod <- brm(severity2 ~ plant_type*fungicide + (1|plotf),
                 data = trtD2Dat, family = Beta,
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3)

mod_check_fun(trtD2Mod)

# save model
save(trtD2Mod, file = "output/severity_fungicide_model_2019_density_exp.rda")

# treatment effects
hypothesis(trtD2Mod, c("fungicide = 0",
                       "fungicide + plant_typeEv_seedling:fungicide = 0",
                       "fungicide + plant_typeMv_seedling:fungicide = 0"))
# effect of fungicide on Mv is significantly different from zero
