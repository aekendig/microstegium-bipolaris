##### info ####

# file: survival_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/7/21
# goal: analyses of growing season survival


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
sevD1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# severity data
sevD1Dat2 <- sevD1Dat %>%
  select(-lesions) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")

sevD2Dat2 <- sevD2Dat %>%
  select(-lesions) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")

# 2018 survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  left_join(plotsD) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  select(-c(month, field_notes, seeds_produced)) %>%
  filter(!is.na(survival)) %>%
  left_join(sevD1Dat2)

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
focD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3))) %>%
  mutate(sp = "Mv",
         age = "seedling") %>%
  full_join(plotsD %>%
              select(plot, treatment) %>%
              expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
              expand_grid(tibble(ID = c("1", "2", "3", "A"),
                                 age = c(rep("seedling", 3), "adult"))) %>%
              mutate(sp = "Ev"))

# 2019 focal survival
# 2018 survival starts in June because plants were replaced through May 24
survD2Dat2 <- survD2Dat %>%
  filter(focal == 1 & replace_date > 20190531) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(plantings = length(unique(replace_date)) + 1) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  left_join(plotsD) %>%
  mutate(plantings = replace_na(plantings, 1),
         survival = ifelse(plantings > 1, 0, 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = "")) %>%
  left_join(sevD2Dat2)

# winter survival 2018-2019
winSurvD1Dat <- survD1Dat %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  left_join(plotsD) %>%
  mutate(September = case_when(seeds_produced == 1 ~ 1, TRUE ~ September),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April) %>%
  left_join(sevD1Dat2)


#### fit models ####

survSevD1Mod <- brm(data = survD1Dat2, family = bernoulli,
                    survival ~ jul_severity * foc + (1|plotf),
                    prior <- c(prior(normal(0, 10), class = "Intercept"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(survSevD1Mod)

survSevD2Mod <- update(survSevD1Mod, newdata = survD2Dat2,
                       control = list(adapt_delta = 0.99)) 
mod_check_fun(survSevD2Mod)

winSurvSevD1Mod <- brm(data = winSurvD1Dat, family = bernoulli,
                       survival ~ sep_severity * foc + (1|plotf),
                       prior <- c(prior(normal(0, 10), class = "Intercept"),
                                  prior(normal(0, 10), class = "b")), # use default for sigma
                       iter = 6000, warmup = 1000, chains = 3, cores = 3) 
winSurvSevD1Mod <- update(winSurvSevD1Mod,
                          control = list(adapt_delta = 0.99))
mod_check_fun(winSurvSevD1Mod)

# save models
save(survSevD1Mod, file = "output/survival_severity_model_2018_density_exp.rda")
save(survSevD2Mod, file = "output/survival_severity_model_2019_density_exp.rda")
save(winSurvSevD1Mod, file = "output/winter_survival_severity_model_2018_density_exp.rda")


#### coefficients ####

mv_hyp = "jul_severity = 0"
evS_hyp = "jul_severity + jul_severity:focs = 0"
evA_hyp = "jul_severity + jul_severity:foca = 0"
evA_win_hyp = "sep_severity = 0"
evS_win_hyp = "sep_severity + sep_severity:focs = 0"


survD1Coef <- hypothesis(survSevD1Mod, c(mv_hyp, evS_hyp, evA_hyp))[[1]]
survD2Coef <- hypothesis(survSevD2Mod, c(mv_hyp, evS_hyp, evA_hyp))[[1]]
winSurvD1Coef <- hypothesis(winSurvSevD1Mod, c(evA_win_hyp, evS_win_hyp))[[1]]
# none are significantly affected by severity

survCoef <- survD1Coef %>%
  full_join(survD2Coef) %>%
  full_join(winSurvD1Coef) %>%
  mutate(year = c(rep(2018, 3), rep(2019, 3), rep(2018, 2)),
         focal = c("Mv", "Ev seedling", "Ev adult",
                   "Mv", "Ev seedling", "Ev adult",
                   "Ev seedling", "Ev adult"),
         season = c(rep("growing season", 6), rep("winter", 2))) %>%
  select(year, focal, season, Estimate:CI.Upper)

write_csv(survCoef, "output/survival_severity_coefficients_2018_2019_density_exp.csv")


#### 2018 Mv model ####

# initial visualization
ggplot(mvSurvD1Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# model
mvSurvD1Mod <- brm(data = mvSurvD1Dat, family = bernoulli,
                   survival ~ fungicide + (1|site/plotf),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                   control = list(adapt_delta = 0.99))
mod_check_fun(mvSurvD1Mod)

# posterior means
post_pred_fun2(mvSurvD1Mod)

# save
save(mvSurvD1Mod, file = "output/mv_growing_season_survival_model_2018_density_exp.rda")


#### 2018 Ev seedling model ####

# initial visualization
ggplot(evSSurvD1Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# model
evSSurvD1Mod <- update(mvSurvD1Mod, newdata = evSSurvD1Dat)
mod_check_fun(evSSurvD1Mod)

# posterior means
post_pred_fun(evSSurvD1Mod)

# save
save(evSSurvD1Mod, file = "output/evS_growing_season_survival_model_2018_density_exp.rda")


#### 2018 Ev adult model ####

# initial visualization
ggplot(evASurvD1Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# model
evASurvD1Mod <- update(mvSurvD1Mod, formula. = survival ~ fungicide + (1|site),
                       newdata = evASurvD1Dat)
mod_check_fun(evASurvD1Mod)

# posterior means
post_pred_fun(evASurvD1Mod)

# save
save(evASurvD1Mod, file = "output/evA_growing_season_survival_model_2018_density_exp.rda")


#### 2019 Mv model ####

# initial visualization
ggplot(mvSurvD2Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(mvSurvD2Dat, aes(x = plantings)) +
  geom_bar() +
  facet_wrap(~ treatment)

# model
mvSurvD2Mod <- update(mvSurvD1Mod, newdata = mvSurvD2Dat)
mod_check_fun(mvSurvD2Mod)

# posterior means
post_pred_fun3(mvSurvD2Mod)

# save
save(mvSurvD2Mod, file = "output/mv_growing_season_survival_model_2019_density_exp.rda")


#### 2019 Ev seedling model ####

# initial visualization
ggplot(evSSurvD2Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evSSurvD2Dat, aes(x = plantings)) +
  geom_bar() +
  facet_wrap(~ treatment)

# model
evSSurvD2Mod <- update(mvSurvD2Mod, newdata = evSSurvD2Dat,
                       control = list(adapt_delta = 0.999))
mod_check_fun(evSSurvD2Mod)

# posterior means
post_pred_fun(evSSurvD2Mod)

# save
save(evSSurvD2Mod, file = "output/evS_growing_season_survival_model_2019_density_exp.rda")


#### 2019 Ev adult model ####

# initial visualization
ggplot(evASurvD2Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evASurvD2Dat, aes(x = plantings)) +
  geom_bar() +
  facet_wrap(~ treatment)

# model
evASurvD2Mod <- update(evASurvD1Mod, newdata = evASurvD2Dat)
mod_check_fun(evASurvD2Mod)

# posterior means
post_pred_fun2(evASurvD2Mod)

# save
save(evASurvD2Mod, file = "output/evA_growing_season_survival_model_2019_density_exp.rda")


#### 2018 Ev seedling winter model ####

# initial visualization
ggplot(evSWinSurvD1Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# model
evSWinSurvD1Mod <- update(evSSurvD1Mod, newdata = evSWinSurvD1Dat,
                          control = list(adapt_delta = 0.999))
mod_check_fun(evSWinSurvD1Mod)

# posterior means
post_pred_fun(evSWinSurvD1Mod)

# save
save(evSWinSurvD1Mod, file = "output/evS_winter_survival_model_2018_density_exp.rda")


#### 2018 Ev adult winter model ####

# initial visualization
ggplot(evAWinSurvD1Dat, aes(treatment, survival, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# model
evAWinSurvD1Mod <- update(evASurvD1Mod, newdata = evAWinSurvD1Dat)
mod_check_fun(evAWinSurvD1Mod)

# posterior means
post_pred_fun(evAWinSurvD1Mod) 

# save
save(evAWinSurvD1Mod, file = "output/evA_winter_survival_model_2018_density_exp.rda")