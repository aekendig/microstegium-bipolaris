##### info ####

# file: plot_severity_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/5/21
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

# import data
d1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
d2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R

# model functions
source("code/brms_model_fitting_functions.R")


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
         bg_sp = sp,
         bg_age = age) %>%
  select(-severity)

bgSevD2Dat <- d2Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") |
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-severity)

# focal dataset (includes background when background is same species)
focNextSevD1Dat <- d1Dat %>%
  filter(month != "jul") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_next_lesions = lesions) %>%
  mutate(month = dplyr::recode(month, 
                               "late_aug" = "jul", # match prior month
                               "sep" = "late_aug")) %>%
  select(-severity)

focNextSevD2Dat <- d2Dat %>%
  filter(month != "may") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_next_lesions = lesions) %>%
  mutate(month = dplyr::recode(month,  # match prior month
                               "jun" = "may",
                               "jul" = "jun",
                               "early_aug" = "jul",
                               "late_aug" = "early_aug")) %>%
  select(-severity)

# prior severity
focSevD1Dat <- d1Dat %>%
  filter(month != "sep") %>% # remove last month
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  select(-severity)

focSevD2Dat <- d2Dat %>%
  filter(month != "late_aug") %>% # remove last month
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
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
         foc_lesions_change = log(foc_next_lesions / foc_lesions),
         foc_lesions_change = case_when(foc_lesions == 0 & foc_next_lesions == 0 ~ log(1), # 4 plots had no lesions, no change
                                        TRUE ~ foc_lesions_change),
         foc_lesions_change = case_when(month == "jul" ~ foc_lesions_change / 2, # account for longer lapse between samples 
                                        TRUE ~ foc_lesions_change)) %>%
  group_by(focal, background) %>%
  mutate(mean_bg_les = mean(bg_lesions, na.rm = T),
         sd_bg_les = sd(bg_lesions, na.rm = T)) %>%
  ungroup() %>%
  mutate(bg_les_s = (bg_lesions - mean_bg_les) / sd_bg_les,
         edge_severity = 0) %>%
  filter(!is.na(foc_next_lesions) & !is.na(foc_lesions) & !is.na(bg_lesions))

# check for missing data
filter(sevD1Dat, is.na(foc_next_lesions))
filter(sevD1Dat, is.na(bg_lesions)) # 18 missing
filter(sevD1Dat, is.na(foc_lesions)) # 7 missing
filter(sevD1Dat, is.na(foc_lesions_change) & !is.na(foc_next_lesions) & !is.na(foc_lesions))  
filter(sevD1Dat, foc_lesions_change %in% c(Inf, -Inf))

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
         foc_lesions_change = log(foc_next_lesions / foc_lesions),
         foc_lesions_change = case_when(foc_lesions == 0 & foc_next_lesions == 0 ~ log(1), # 4 plots had no lesions, no change
                                        TRUE ~ foc_lesions_change),
         foc_lesions_change = case_when(foc_lesions == 0 & foc_next_lesions > 0 ~ log(foc_next_lesions / 1e-6),
                                        TRUE ~ foc_lesions_change),  # used close to smallest lesion amount (6.4e-6) as denominator
         foc_lesions_change = case_when(foc_lesions > 0 & foc_next_lesions == 0 ~ log(1e-6 / foc_lesions),
                                        TRUE ~ foc_lesions_change)) %>%  # used close to smallest lesion amount (6.4e-6) as numerator
  group_by(focal, background) %>%
  mutate(mean_bg_les = mean(bg_lesions, na.rm = T),
         sd_bg_les = sd(bg_lesions, na.rm = T)) %>%
  ungroup() %>%
  mutate(bg_les_s = (bg_lesions - mean_bg_les) / sd_bg_les) %>%
  filter(!is.na(foc_next_lesions) & !is.na(foc_lesions) & !is.na(bg_lesions))

# check for missing data
filter(sevD2Dat, is.na(foc_next_lesions))
filter(sevD2Dat, is.na(bg_lesions)) # 189 missing
filter(sevD2Dat, is.na(foc_lesions)) # 185 missing
filter(sevD2Dat, is.na(foc_lesions_change) & !is.na(foc_next_lesions) & !is.na(foc_lesions)) 
filter(sevD2Dat, foc_lesions_change %in% c(Inf, -Inf)) # 21 cases


#### visualizations ####

# lesions change distributions
ggplot(sevD1Dat, aes(foc_lesions_change)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(foc_lesions_change)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

# bg lesions distributions
ggplot(sevD1Dat, aes(bg_lesions)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD1Dat, aes(bg_les_s)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(bg_lesions)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(bg_les_s)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

# pairwise combinations
ggplot(sevD1Dat, aes(bg_les_s, foc_lesions_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(focal ~ background, scales = "free")

ggplot(sevD2Dat, aes(bg_les_s, foc_lesions_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(focal ~ background, scales = "free")
# looked at high Ev adult plot (bg_les_s > 7) - no Ev scans for late August

# edge
ggplot(sevD2Dat, aes(edge_severity)) +
  geom_density() +
  facet_wrap(~ focal, scales = "free")

ggplot(sevD2Dat, aes(edge_severity, foc_lesions_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(~ focal, scales = "free")


#### fit models ####

# Mv
sevD1Mod <- brm(foc_lesions_change ~ bg_les_s * foc * bg * fungicide + (1|plotf),
                   data = sevD1Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99))
# use plot instead of site for random effects because of multiple temporal sampling
mod_check_fun(sevD1Mod)

sevD2Mod <- brm(foc_lesions_change ~ bg_les_s * foc * bg * fungicide + edge_severity + edge_severity:foc + edge_severity:fungicide + edge_severity:foc:fungicide + (1|plotf),
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
evS_mv_ctrl_hyp = "-bg_les_s:focs = 0"
evS_mv_fung_hyp = "-bg_les_s:focs - bg_les_s:focs:fungicide = 0"
evA_mv_ctrl_hyp = "-bg_les_s:foca = 0"
evA_mv_fung_hyp = "-bg_les_s:foca - bg_les_s:foca:fungicide = 0"

# EvS 
mv_evS_ctrl_hyp = "0.5 * (bg_les_s:focs + bg_les_s:focs:bgs + bg_les_s:foca + bg_les_s:foca:bgs) = 0"
mv_evS_fung_hyp = "0.5 * (bg_les_s:focs + bg_les_s:focs:bgs + bg_les_s:focs:fungicide + bg_les_s:focs:bgs:fungicide + bg_les_s:foca + bg_les_s:foca:bgs + bg_les_s:foca:fungicide + bg_les_s:foca:bgs:fungicide) = 0"

# EvA 
mv_evA_ctrl_hyp = "0.5 * (bg_les_s:focs + bg_les_s:focs:bga + bg_les_s:foca + bg_les_s:foca:bga) = 0"
mv_evA_fung_hyp = "0.5 * (bg_les_s:focs + bg_les_s:focs:bga + bg_les_s:focs:fungicide + bg_les_s:focs:bga:fungicide + bg_les_s:foca + bg_les_s:foca:bga + bg_les_s:foca:fungicide + bg_les_s:foca:bga:fungicide) = 0"

sevD1hyps <- hypothesis(sevD1Mod, 
                        c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
                          mv_evS_ctrl_hyp, mv_evS_fung_hyp, mv_evA_ctrl_hyp, mv_evA_fung_hyp))

sevD2hyps <- hypothesis(sevD2Mod, 
                        c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
                          mv_evS_ctrl_hyp, mv_evS_fung_hyp, mv_evA_ctrl_hyp, mv_evA_fung_hyp))
# none are significantly different from zero

write_csv(sevD1hyps[[1]], "output/plot_transmission_intra_vs_intra_2018_density_exp.csv")
write_csv(sevD2hyps[[1]], "output/plot_transmission_intra_vs_intra_2019_density_exp.csv")


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
# none are significantly different than zero

write_csv(edgeD2hyps[[1]], "output/edge_transmission_2019_density_exp.csv")


#### figure ####

# bg_les_s function
bg_les_fun <- function(f, b, trt, yr){
  
  if(yr == "2018"){
    dat <- sevD1Dat %>% filter(foc == f & bg == b & treatment == trt)
  }else{
    dat <- sevD2Dat %>% filter(foc == f & bg == b & treatment == trt)
  }
  
  bg_les_s = seq(min(dat$bg_les_s), max(dat$bg_les_s), length.out = 100)
  
  return(bg_les_s)
}

# predicted data
predD1Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2018",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(bg_les_s = pmap(list(foc, bg, treatment, year), bg_les_fun)) %>%
  unnest(bg_les_s) %>%
  mutate(foc_lesions_change = fitted(sevD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
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
  mutate(bg_les_s = pmap(list(foc, bg, treatment, year), bg_les_fun)) %>%
  unnest(bg_les_s) %>%
  mutate(foc_lesions_change = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(sevD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predDat <- predD1Dat %>%
  full_join(predD2Dat) %>%
  mutate(focal = fct_recode(foc, Mv = "m", "Ev seedling" = "s", "Ev adult" = "a"),
         background = fct_recode(bg, Mv = "m", "Ev seedling" = "s", "Ev adult" = "a"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

# Mv background
mv_mv_ctrl_beta = "bg_les_s = 0"
mv_mv_fung_beta = "bg_les_s + bg_les_s:fungicide = 0"
evS_mv_ctrl_beta = "bg_les_s + bg_les_s:focs = 0"
evS_mv_fung_beta = "bg_les_s + bg_les_s:fungicide + bg_les_s:focs + bg_les_s:focs:fungicide = 0"
evA_mv_ctrl_beta = "bg_les_s + bg_les_s:foca = 0"
evA_mv_fung_beta = "bg_les_s + bg_les_s:fungicide + bg_les_s:foca + bg_les_s:foca:fungicide = 0"

# EvS background
evS_evS_ctrl_beta = "bg_les_s + bg_les_s:focs + bg_les_s:bgs + bg_les_s:focs:bgs = 0"
evS_evS_fung_beta = "bg_les_s + bg_les_s:focs + bg_les_s:bgs + bg_les_s:focs:bgs + bg_les_s:fungicide + bg_les_s:focs:fungicide + bg_les_s:bgs:fungicide + bg_les_s:focs:bgs:fungicide = 0"
mv_evS_ctrl_beta = "bg_les_s +  bg_les_s:bgs = 0"
mv_evS_fung_beta = "bg_les_s +  bg_les_s:bgs + bg_les_s:fungicide + bg_les_s:bgs:fungicide = 0"
evA_evS_ctrl_beta = "bg_les_s +  bg_les_s:bgs + bg_les_s:foca + bg_les_s:foca:bgs = 0"
evA_evS_fung_beta = "bg_les_s +  bg_les_s:bgs + bg_les_s:foca + bg_les_s:foca:bgs + bg_les_s:fungicide +  bg_les_s:bgs:fungicide + bg_les_s:foca:fungicide + bg_les_s:foca:bgs:fungicide = 0"

# EvA intra vs. inter
evA_evA_ctrl_beta = "bg_les_s + bg_les_s:foca + bg_les_s:bga + bg_les_s:foca:bga = 0"
evA_evA_fung_beta = "bg_les_s + bg_les_s:foca + bg_les_s:bga + bg_les_s:foca:bga + bg_les_s:fungicide + bg_les_s:foca:fungicide + bg_les_s:bga:fungicide + bg_les_s:foca:bga:fungicide = 0"
mv_evA_ctrl_beta = "bg_les_s +  bg_les_s:bga = 0"
mv_evA_fung_beta = "bg_les_s +  bg_les_s:bga + bg_les_s:fungicide + bg_les_s:bga:fungicide = 0"
evS_evA_ctrl_beta = "bg_les_s + bg_les_s:bga + bg_les_s:focs + bg_les_s:focs:bga = 0"
evS_evA_fung_beta = "bg_les_s + bg_les_s:bga + bg_les_s:focs + bg_les_s:focs:bga + bg_les_s:fungicide + bg_les_s:bga:fungicide + bg_les_s:focs:fungicide + bg_les_s:focs:bga:fungicide = 0"

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

# edit to save
betaDatSave <- betaDat %>%
  left_join(predDat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(treatment = fct_recode(treatment, "control" = "control (water)")) %>%
  select(year, focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper)

# save
write_csv(betaDatSave, "output/plot_transmission_coefficients_2018_2019_density_exp.csv")

# combine with preddat
predDat2 <- predDat %>%
  left_join(betaDat)

# raw data
figDat <- sevD1Dat %>%
  select(month, site, plot, treatment, focal, background, bg_les_s, foc_lesions_change) %>%
  mutate(year = "2018") %>%
  full_join(sevD2Dat %>%
              select(month, site, plot, treatment, focal, background, bg_les_s, foc_lesions_change) %>%
              mutate(year = "2019")) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# betas for figure
# sig beta values
betaDat2 <- betaDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, bg, background) %>%
              summarise(bg_les_s = max(bg_les_s))) %>%
  left_join(predDat2 %>%
              select(foc, focal) %>% # add focal names
              unique()) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(lower = min(CI.Lower)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(sev = min(foc_lesions_change))) %>%
              rowwise() %>%
              mutate(foc_lesions_change = min(c(sev, lower))) %>%
              ungroup() %>%
              select(-c(sev, lower))) %>%
  mutate(parm = round(Estimate, 2))

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.1, 0, 0.2, 2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

yearText <- tibble(year = c("2018", "2019"),
                   bg_les_s = c(-0.7, -0.2),
                   foc_lesions_change = c(4.2, 15),
                   background = "Ev adult",
                   focal = "Ev adult",
                   treatment = "fungicide")

textSize = 3

# water figure
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = bg_les_s, y = foc_lesions_change)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0, vjust = 1) +
  geom_text(data = filter(betaDat2, year == "2018"), 
            aes(label = paste("beta", " == ", parm, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 0, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab("Initial infected tissue") +
  ylab("Change in infected tissue") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# fungicide figure
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = bg_les_s, y = foc_lesions_change)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0, vjust = 1) +
  geom_text(data = filter(betaDat2, year == "2019"), 
            aes(label = paste("beta", " == ", parm, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 0, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab("Initial infected tissue") +
  ylab("Change in infected tissue") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg <- get_legend(pairD1Fig)

# combine plots
combFig <- plot_grid(pairD1Fig + theme(legend.position = "none"), pairD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.87),
                     label_x = c(0, -0.01))

# combine
pdf("output/plot_transmission_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()

