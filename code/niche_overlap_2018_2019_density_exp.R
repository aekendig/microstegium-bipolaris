##### info ####

# file:niche_overlap_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 5/19/21
# goal: fungicide effects on alphas and niche overlap


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# load models
load("output/mv_tillers_mv_density_model_2018_density_exp.rda")
load("output/evS_tillers_mv_density_model_2018_density_exp.rda")
load("output/evA_tillers_mv_density_model_2018_density_exp.rda")
load("output/mv_tillers_evS_density_model_2018_density_exp.rda")
load("output/evS_tillers_evS_density_model_2018_density_exp.rda")
load("output/evA_tillers_evS_density_model_2018_density_exp.rda")
load("output/mv_tillers_evA_density_model_2018_density_exp.rda")
load("output/evS_tillers_evA_density_model_2018_density_exp.rda")
load("output/evA_tillers_evA_density_model_2018_density_exp.rda")
load("output/mv_biomass_mv_density_model_2019_density_exp.rda")
load("output/evS_biomass_mv_density_model_2019_density_exp.rda")
load("output/evA_biomass_mv_density_model_2019_density_exp.rda")
load("output/mv_biomass_evS_density_model_2019_density_exp.rda")
load("output/evS_biomass_evS_density_model_2019_density_exp.rda")
load("output/evA_biomass_evS_density_model_2019_density_exp.rda")
load("output/mv_biomass_evA_density_model_2019_density_exp.rda")
load("output/evS_biomass_evA_density_model_2019_density_exp.rda")
load("output/evA_biomass_evA_density_model_2019_density_exp.rda")

load("output/mv_biomass_Ricker_model_2019_density_exp.rda")
load("output/evS_biomass_Ricker_model_2019_density_exp.rda")
load("output/evA_biomass_Ricker_model_2019_density_exp.rda")
load("output/mv_tillers_Ricker_model_2018_density_exp.rda")
load("output/evS_tillers_Ricker_model_2018_density_exp.rda")
load("output/evA_tillers_Ricker_model_2018_density_exp.rda")

load("output/evA_growing_season_survival_model_2018_density_exp.rda")
load("output/evA_growing_season_survival_model_2019_density_exp.rda")
load("output/evA_winter_survival_model_2018_density_exp.rda")


#### adult survival ####

(survD1 <- exp(summary(evASurvD1Mod)$fixed[1, 1]) / (1 + exp(summary(evASurvD1Mod)$fixed[1, 1])))
(survD2 <- exp(summary(evASurvD2Mod)$fixed[1, 1]) / (1 + exp(summary(evASurvD2Mod)$fixed[1, 1])))
(survW <- exp(summary(evAWinSurvD1Mod)$fixed[1, 1]) / (1 + exp(summary(evAWinSurvD1Mod)$fixed[1, 1])))

(survP <- round(mean(c(survD1, survD2)) * survW, 2))
(survS <- 1 - survP)


#### 2018 values ####

# check models
summary(evSEvSD1Mod)
summary(evSEvAD1Mod)
summary(evSMvD1Mod)

summary(evAEvSD1Mod)
summary(evAEvAD1Mod)
summary(evAMvD1Mod)

summary(mvEvSD1Mod)
summary(mvEvAD1Mod)
summary(mvMvD1Mod)

# combine alphas
alphasD1 <- posterior_samples(evSEvAD1Mod) %>%
  mutate(subscr = "sp") %>%
  full_join(posterior_samples(evSEvSD1Mod) %>%
              mutate(subscr = "ss")) %>%
  full_join(posterior_samples(evSMvD1Mod) %>%
            mutate(subscr = "sa")) %>%
  full_join(posterior_samples(evAEvAD1Mod) %>%
              mutate(subscr = "pp")) %>%
  full_join(posterior_samples(evAEvSD1Mod) %>%
              mutate(subscr = "ps")) %>%
  full_join(posterior_samples(evAMvD1Mod) %>%
              mutate(subscr = "pa")) %>%
  full_join(posterior_samples(mvEvAD1Mod) %>%
              mutate(subscr = "ap")) %>%
  full_join(posterior_samples(mvEvSD1Mod) %>%
              mutate(subscr = "as")) %>%
  full_join(posterior_samples(mvMvD1Mod) %>%
              mutate(subscr = "aa")) %>%
  select(subscr, b_alpha_treatmentfungicide, b_alpha_treatmentwater) %>%
  as_tibble() %>%
  pivot_longer(cols = -subscr,
               names_to = "treatment",
               values_to = "alpha",
               names_prefix = "b_alpha_treatment") %>%
  mutate(iter = rep(rep(1:15000, each = 2), 9)) %>%
  pivot_wider(names_from = subscr,
              values_from = alpha) %>%
  mutate(ap = ap * survP,
         as = as * survS,
         pa = pa * survP,
         sa = sa * survS,
         pp = pp * survP,
         ss = ss * survS,
         rho = sqrt((ap + as) * (pa + sa) / (aa * (ss + pp)))) 

rhoD1 <- alphasD1 %>%
  group_by(treatment) %>%
  mean_hdi(rho)

diffD1 <- alphasD1 %>%
  pivot_wider(names_from = treatment,
              values_from = c(sp, ss, sa, pp, ps, pa, ap, as, aa, rho),
              names_glue = "{.value}_{treatment}") %>%
  transmute(sp = (sp_fungicide - sp_water) / sp_water,
            ss = (ss_fungicide - ss_water) / ss_water,
            sa = (sa_fungicide - sa_water) / sa_water,
            pp = (pp_fungicide - pp_water) / pp_water,
            ps = (ps_fungicide - ps_water) / ps_water,
            pa = (pa_fungicide - pa_water) / pa_water,
            ap = (ap_fungicide - ap_water) / ap_water,
            as = (as_fungicide - as_water) / as_water,
            aa = (aa_fungicide - aa_water) / aa_water,
            rho = (rho_fungicide - rho_water) / rho_water) %>%
  pivot_longer(cols = everything(),
               names_to = "alpha",
               values_to = "diff") %>%
  group_by(alpha) %>%
  mean_hdi(diff)


#### 2019 values ####

# check models
summary(evSEvSD2Mod)
summary(evSEvAD2Mod)
summary(evSMvD2Mod)

summary(evAEvSD2Mod)
summary(evAEvAD2Mod)
summary(evAMvD2Mod)

summary(mvEvSD2Mod)
summary(mvEvAD2Mod)
summary(mvMvD2Mod)

# combine alphas
alphasD2 <- posterior_samples(evSEvAD2Mod) %>%
  mutate(subscr = "sp") %>%
  full_join(posterior_samples(evSEvSD2Mod) %>%
              mutate(subscr = "ss")) %>%
  full_join(posterior_samples(evSMvD2Mod) %>%
              mutate(subscr = "sa")) %>%
  full_join(posterior_samples(evAEvAD2Mod) %>%
              mutate(subscr = "pp")) %>%
  full_join(posterior_samples(evAEvSD2Mod) %>%
              mutate(subscr = "ps")) %>%
  full_join(posterior_samples(evAMvD2Mod) %>%
              mutate(subscr = "pa")) %>%
  full_join(posterior_samples(mvEvAD2Mod) %>%
              mutate(subscr = "ap")) %>%
  full_join(posterior_samples(mvEvSD2Mod) %>%
              mutate(subscr = "as")) %>%
  full_join(posterior_samples(mvMvD2Mod) %>%
              mutate(subscr = "aa")) %>%
  select(subscr, b_alpha_treatmentfungicide, b_alpha_treatmentwater) %>%
  as_tibble() %>%
  pivot_longer(cols = -subscr,
               names_to = "treatment",
               values_to = "alpha",
               names_prefix = "b_alpha_treatment") %>%
  mutate(iter = rep(rep(1:15000, each = 2), 9)) %>%
  pivot_wider(names_from = subscr,
              values_from = alpha) %>%
  mutate(ap = ap * survP,
         as = as * survS,
         pa = pa * survP,
         sa = sa * survS,
         pp = pp * survP,
         ss = ss * survS,
         rho = sqrt((ap + as) * (pa + sa) / (aa * (ss + pp))))

rhoD2 <- alphasD2 %>%
  group_by(treatment) %>%
  mean_hdi(rho)

diffD2 <- alphasD2 %>%
  pivot_wider(names_from = treatment,
              values_from = c(sp, ss, sa, pp, ps, pa, ap, as, aa, rho),
              names_glue = "{.value}_{treatment}") %>%
  transmute(sp = (sp_fungicide - sp_water) / sp_water,
            ss = (ss_fungicide - ss_water) / ss_water,
            sa = (sa_fungicide - sa_water) / sa_water,
            pp = (pp_fungicide - pp_water) / pp_water,
            ps = (ps_fungicide - ps_water) / ps_water,
            pa = (pa_fungicide - pa_water) / pa_water,
            ap = (ap_fungicide - ap_water) / ap_water,
            as = (as_fungicide - as_water) / as_water,
            aa = (aa_fungicide - aa_water) / aa_water,
            rho = (rho_fungicide - rho_water) / rho_water) %>%
  pivot_longer(cols = everything(),
               names_to = "alpha",
               values_to = "diff") %>%
  group_by(alpha) %>%
  mean_hdi(diff)


#### Ricker models ####

# combine alphas
alphasRick <- posterior_samples(mvRickD1Mod) %>%
  mutate(focal = "a", year = "2018") %>%
  full_join(posterior_samples(evSRickD1Mod) %>%
              mutate(focal = "s", year = "2018")) %>%
  full_join(posterior_samples(evARickD1Mod) %>%
              mutate(focal = "p", year = "2018")) %>%
  full_join(posterior_samples(mvRickD2Mod) %>%
              mutate(focal = "a", year = "2019")) %>%
  full_join(posterior_samples(evSRickD2Mod) %>%
              mutate(focal = "s", year = "2019")) %>%
  full_join(posterior_samples(evARickD2Mod) %>%
              mutate(focal = "p", year = "2019")) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  mutate(a_water = b_density,
         s_water = b_density + b_density_backgroundEv_seedling,
         p_water = b_density + b_density_backgroundEv_adult,
         a_fungicide = b_density + b_fungicide_density,
         s_fungicide = s_water + b_fungicide_density + b_fungicide_backgroundEv_seedling,
         p_fungicide = p_water + b_fungicide_density + b_fungicide_backgroundEv_adult) %>%
  select(year, focal, a_water, s_water, p_water, a_fungicide, s_fungicide, p_fungicide) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(year, focal),
               names_to = c(".value", "treatment"),
               names_pattern = "(.)_(.*)") %>%
  mutate(iter = rep(rep(1:15000, each = 2), 6)) %>%
  pivot_wider(names_from = focal,
              names_glue = "{focal}{.value}",
              values_from = c(a, s, p)) %>%
  mutate(aa = -1 * aa,
         sa = -1 * sa * survS,
         pa = -1 * pa * survP,
         as = -1 * as * survS,
         ss = -1 * ss * survS,
         ps = -1 * ps,
         ap = -1 * ap * survP,
         pp = -1 * pp * survP,
         rho = sqrt((ap + as) * (pa + sa) / (aa * (ss + pp))))

sum(is.na(alphasRick$rho))/60000
# almost half are lost because of negative (facilitative) values

rhoRick <- alphasRick %>%
  pivot_longer(cols = c(aa:rho),
               names_to = "param",
               values_to = "est") %>%
  group_by(year, treatment, param) %>%
  mean_hdi(est, na.rm = T) %>%
  mutate(sig = case_when(.lower < 0 & .upper < 0 ~ T,
                         .upper > 0 & .lower > 0 ~ T,
                         TRUE ~ F))
data.frame(rhoRick)

diffRick <- alphasRick %>%
  pivot_wider(names_from = treatment,
              values_from = c(sp, ss, sa, pp, ps, pa, ap, as, aa, rho),
              names_glue = "{.value}_{treatment}") %>%
  transmute(sp = (sp_fungicide - sp_water) / sp_water,
            ss = (ss_fungicide - ss_water) / ss_water,
            sa = (sa_fungicide - sa_water) / sa_water,
            pp = (pp_fungicide - pp_water) / pp_water,
            ps = (ps_fungicide - ps_water) / ps_water,
            pa = (pa_fungicide - pa_water) / pa_water,
            ap = (ap_fungicide - ap_water) / ap_water,
            as = (as_fungicide - as_water) / as_water,
            aa = (aa_fungicide - aa_water) / aa_water,
            rho = (rho_fungicide - rho_water) / rho_water) %>%
  pivot_longer(cols = everything(),
               names_to = "alpha",
               values_to = "diff") %>%
  group_by(alpha) %>%
  mean_hdi(diff, na.rm = T)
# no sig effect of fungicide