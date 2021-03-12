##### info ####

# file: severity_change_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/11/21
# goal: analyses of severity change as a function of density


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import severity change data
sevChgDat <- read_csv("intermediate-data/severity_change_model_ran_slopes_2018_2019_density_exp.csv")
# severity_change_2018_2019_density_exp.R
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp.R

# biomass
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R

# function to transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

backtransform01 <- function(x) {
  (x * length(x) - 0.5) / (length(x) - 1)
}  

# function to transform logit to proportion
logit2prop <- function(x) {
  (exp(x) / (1 + exp(x)))
}


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(Mv_seedling_density = case_when(background == "Mv seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_seedling_density = case_when(background == "Ev seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_adult_density = case_when(background == "Ev adult" ~ background_density + 1, 
                                      TRUE ~ 1)) %>%
  select(plot, treatment, Mv_seedling_density, Ev_seedling_density, Ev_adult_density)

# plot biomass
# use average of other species in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         sp_age = paste(sp, age, sep = "_")) %>%
  pivot_wider(-c(sp, age),
              names_from = sp_age,
              names_glue = "{sp_age}_biomass",
              values_from = biomass.g) %>%
  mutate(Mv_seedling_biomass = replace_na(Mv_seedling_biomass, 0),
         Ev_seedling_biomass = replace_na(Ev_seedling_biomass, 0),
         Ev_adult_biomass = replace_na(Ev_adult_biomass, 0)) %>%
  select(-none_seedling_biomass) %>%
  left_join(mvBioD2Dat %>%
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g)) %>%
              group_by(site, plot, treatment, ) %>%
              summarise(Mv_seedling_biomass_foc = sum(biomass_weight.g)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID %in% c("1", "2", "3")) %>%
              group_by(site, plot, treatment) %>%
              mutate(weight_adj = mean(weight, na.rm = T)) %>%
              ungroup() %>%
              mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                        TRUE ~ weight)) %>%
              group_by(site, plot, treatment) %>%
              summarise(Ev_seedling_biomass_foc = sum(weight)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID == "A") %>%
              select(site, plot, treatment, weight) %>%
              rename(Ev_adult_biomass_foc = weight)) %>%
  mutate(Mv_seedling_biomass = Mv_seedling_biomass + Mv_seedling_biomass_foc,
         Ev_seedling_biomass = Ev_seedling_biomass + Ev_seedling_biomass_foc,
         Ev_adult_biomass = Ev_adult_biomass + Ev_adult_biomass_foc,
         total_biomass = Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass) %>%
  select(-c(Mv_seedling_biomass_foc, Ev_seedling_biomass_foc, Ev_adult_biomass_foc))

# raw severity
sevD1Dat2 <- sevD1Dat %>%
  filter(month == "sep") %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         severity_01 = transform01(severity)) %>%
  select(site, plot, treatment, sp, age, ID, severity, severity_01)

sevD2Dat2 <- sevD2Dat %>%
  filter(month == "late_aug") %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         severity_01 = transform01(severity)) %>%
  select(site, plot, treatment, sp, age, ID, severity, severity_01)

# combine and separate data
dat <- sevChgDat %>%
  rename(sev_chg = slope) %>%
  full_join(sevD1Dat2 %>%
              mutate(year = 2018)) %>%
  full_join(sevD2Dat2 %>%
              mutate(year = 2019)) %>%
  left_join(plotDens) %>%
  left_join(plotBioD2Dat) %>%
  mutate(plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plot_ID = paste(site, plot, substr(treatment, 1, 1), sep = "_"))

mvD2Dat <- dat %>%
  filter(plant_group == "Mv_seedling" & year == 2019 & !is.na(sev_chg)) %>%
  mutate(mv_bio_s = scale(Mv_seedling_biomass)[,1],
         evS_bio_s = scale(Ev_seedling_biomass)[,1],
         evA_bio_s = scale(Ev_adult_biomass)[,1],
         tot_bio_s = scale(total_biomass)[,1])

mvD1Dat <- dat %>%
  filter(plant_group == "Mv_seedling" & year == 2018 & !is.na(sev_chg)) %>%
  mutate(mv_bio_s = scale(Mv_seedling_biomass)[,1],
         evS_bio_s = scale(Ev_seedling_biomass)[,1],
         evA_bio_s = scale(Ev_adult_biomass)[,1],
         tot_bio_s = scale(total_biomass)[,1])

evSD2Dat <- dat %>%
  filter(plant_group == "Ev_seedling" & year == 2019 & !is.na(sev_chg)) %>%
  mutate(mv_bio_s = scale(Mv_seedling_biomass)[,1],
         evS_bio_s = scale(Ev_seedling_biomass)[,1],
         evA_bio_s = scale(Ev_adult_biomass)[,1],
         tot_bio_s = scale(total_biomass)[,1])

evAD2Dat <- dat %>%
  filter(plant_group == "Ev_adult" & year == 2019 & !is.na(sev_chg)) %>%
  mutate(mv_bio_s = scale(Mv_seedling_biomass)[,1],
         evS_bio_s = scale(Ev_seedling_biomass)[,1],
         evA_bio_s = scale(Ev_adult_biomass)[,1],
         tot_bio_s = scale(total_biomass)[,1])


#### initial visualizations ####

# histogram
ggplot(dat, aes(x = sev_chg)) +
  geom_histogram() +
  facet_grid(year ~ plant_group, scales = "free")

# plot-scale
ggplot(dat, aes(as.factor(plot), sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(year ~ plant_group, scales = "free")

# Mv density
dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(Mv_seedling_density, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(year ~ plant_group, scales = "free")

# Ev seedling density
dat %>%
  filter(plot %in% c(1, 5:7)) %>%
  ggplot(aes(Ev_seedling_density, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(year ~ plant_group, scales = "free")

# EvA density
dat %>%
  filter(plot %in% c(1, 8:10)) %>%
  ggplot(aes(Ev_adult_density, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(year ~ plant_group, scales = "free")

# Mv biomass
dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(Mv_seedling_biomass, sev_chg, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(year ~ plant_group, scales = "free")

# Ev seedling biomass
dat %>%
  filter(plot %in% c(1, 5:7)) %>%
  ggplot(aes(Ev_seedling_biomass, sev_chg, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(year ~ plant_group, scales = "free")

# EvA biomass
dat %>%
  filter(plot %in% c(1, 8:10)) %>%
  ggplot(aes(Ev_adult_biomass, sev_chg, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(year ~ plant_group, scales = "free")

# total biomass
dat %>%
  ggplot(aes(total_biomass, sev_chg, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(year ~ plant_group, scales = "free")

# 2019 site variation
ggplot(dat, aes(site, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(plotBioD2Dat %>% filter(plot %in% 1:4), aes(site, Mv_seedling_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### Mv 2019 model ####

# intercept prior
mvD2Dat %>%
  filter(treatment == "water") %>%
  summarise(sev_chg = mean(sev_chg))
# 0.57

# model
mvSevChgD2Mod <- brm(sev_chg ~ fungicide * (mv_bio_s + evS_bio_s + evA_bio_s) + (1|plot_ID), 
                data = mvD2Dat, family = "normal",
                prior <- c(prior(normal(0.6, 1), class = Intercept),
                           prior(normal(0, 1), class = b)), # use default priors for sigma and random effect
                iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(mvSevChgD2Mod)
pp_check(mvSevChgD2Mod, nsamples = 100)
plot(mvSevChgD2Mod)

# simulated data
mvMvD2Sim <- mvD2Dat %>%
  select(treatment, fungicide, mv_bio_s, Mv_seedling_biomass) %>%
  unique() %>%
  mutate(evS_bio_s = 0,
         evA_bio_s = 0) %>%
  mutate(pred = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

mvEvSD2Sim <- mvD2Dat %>%
  select(treatment, fungicide, evS_bio_s, Ev_seedling_biomass) %>%
  unique() %>%
  mutate(mv_bio_s = 0,
         evA_bio_s = 0) %>%
  mutate(pred = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

mvEvAD2Sim <- mvD2Dat %>%
  select(treatment, fungicide, evA_bio_s, Ev_adult_biomass) %>%
  unique() %>%
  mutate(mv_bio_s = 0,
         evS_bio_s = 0) %>%
  mutate(pred = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(mvSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

# figures
ggplot(mvD2Dat, aes(x = Mv_seedling_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = mvMvD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = mvMvD2Sim, aes(y = pred))

ggplot(mvD2Dat, aes(x = Ev_seedling_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = mvEvSD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = mvEvSD2Sim, aes(y = pred))

ggplot(mvD2Dat, aes(x = Ev_adult_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = mvEvAD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = mvEvAD2Sim, aes(y = pred))

# total biomass model
mvSevChgD2Mod2 <- brm(sev_chg ~ fungicide * tot_bio_s + (1|plot_ID), 
                     data = mvD2Dat, family = "normal",
                     prior <- c(prior(normal(0.6, 1), class = Intercept),
                                prior(normal(0, 1), class = b)), # use default priors for sigma and random effect
                     iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(mvSevChgD2Mod2)
pp_check(mvSevChgD2Mod2, nsamples = 100)
plot(mvSevChgD2Mod2)

# compare with separate biomass model
mvSevChgD2Loo <- loo(mvSevChgD2Mod, reloo = T)
mvSevChgD2Loo2 <- loo(mvSevChgD2Mod2, reloo = T)
loo_compare(mvSevChgD2Loo, mvSevChgD2Loo2)
# model 2 is a better fit

# simulated data
mvD2Sim <- mvD2Dat %>%
  select(treatment, fungicide, tot_bio_s, total_biomass) %>%
  unique() %>%
  mutate(pred = fitted(mvSevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(mvSevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(mvSevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])


#### Ev seedling 2019 model ####

# intercept prior
evSD2Dat %>%
  filter(treatment == "water") %>%
  summarise(sev_chg = mean(sev_chg))
# 0.62

# model
evSSevChgD2Mod <- update(mvSevChgD2Mod, newdata = evSD2Dat)
summary(evSSevChgD2Mod)
pp_check(evSSevChgD2Mod, nsamples = 100)
plot(evSSevChgD2Mod)

# simulated data
evSMvD2Sim <- evSD2Dat %>%
  select(treatment, fungicide, mv_bio_s, Mv_seedling_biomass) %>%
  unique() %>%
  mutate(evS_bio_s = 0,
         evA_bio_s = 0) %>%
  mutate(pred = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

evSEvSD2Sim <- evSD2Dat %>%
  select(treatment, fungicide, evS_bio_s, Ev_seedling_biomass) %>%
  unique() %>%
  mutate(mv_bio_s = 0,
         evA_bio_s = 0) %>%
  mutate(pred = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

evSEvAD2Sim <- evSD2Dat %>%
  select(treatment, fungicide, evA_bio_s, Ev_adult_biomass) %>%
  unique() %>%
  mutate(mv_bio_s = 0,
         evS_bio_s = 0) %>%
  mutate(pred = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evSSevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

# figures
ggplot(evSD2Dat, aes(x = Mv_seedling_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = evSMvD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = evSMvD2Sim, aes(y = pred))

ggplot(evSD2Dat, aes(x = Ev_seedling_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = evSEvSD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = evSEvSD2Sim, aes(y = pred))

ggplot(evSD2Dat, aes(x = Ev_adult_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = evSEvAD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = evSEvAD2Sim, aes(y = pred))

# total biomass model
evSSevChgD2Mod2 <- update(mvSevChgD2Mod2, newdata = evSD2Dat)
summary(evSSevChgD2Mod2)
pp_check(evSSevChgD2Mod2, nsamples = 100)
plot(evSSevChgD2Mod2)

# compare with separate biomass model
evSSevChgD2Loo <- loo(evSSevChgD2Mod, reloo = T)
evSSevChgD2Loo2 <- loo(evSSevChgD2Mod2, reloo = T)
loo_compare(evSSevChgD2Loo, evSSevChgD2Loo2)
# model 1 is a better fit

# simulated data
evSD2Sim <- evSD2Dat %>%
  select(treatment, fungicide, tot_bio_s, total_biomass) %>%
  unique() %>%
  mutate(pred = fitted(evSSevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evSSevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evSSevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])


#### Ev adult 2019 model ####

# intercept prior
evAD2Dat %>%
  filter(treatment == "water") %>%
  summarise(sev_chg = mean(sev_chg))
# 0.72

# model
evASevChgD2Mod <- brm(sev_chg ~ fungicide * (mv_bio_s + evS_bio_s + evA_bio_s), 
                      data = evAD2Dat, family = "normal",
                      prior <- c(prior(normal(0.6, 1), class = Intercept),
                                 prior(normal(0, 1), class = b)), # use default priors for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(evASevChgD2Mod)
pp_check(evASevChgD2Mod, nsamples = 100)
plot(evASevChgD2Mod)

# simulated data
evAMvD2Sim <- evAD2Dat %>%
  select(treatment, fungicide, mv_bio_s, Mv_seedling_biomass) %>%
  unique() %>%
  mutate(evS_bio_s = 0,
         evA_bio_s = 0) %>%
  mutate(pred = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

evAEvSD2Sim <- evAD2Dat %>%
  select(treatment, fungicide, evS_bio_s, Ev_seedling_biomass) %>%
  unique() %>%
  mutate(mv_bio_s = 0,
         evA_bio_s = 0) %>%
  mutate(pred = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

evAEvAD2Sim <- evAD2Dat %>%
  select(treatment, fungicide, evA_bio_s, Ev_adult_biomass) %>%
  unique() %>%
  mutate(mv_bio_s = 0,
         evS_bio_s = 0) %>%
  mutate(pred = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evASevChgD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])

# figures
ggplot(evAD2Dat, aes(x = Mv_seedling_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = evAMvD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = evAMvD2Sim, aes(y = pred))

ggplot(evAD2Dat, aes(x = Ev_seedling_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = evAEvSD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = evAEvSD2Sim, aes(y = pred))

ggplot(evAD2Dat, aes(x = Ev_adult_biomass, y = sev_chg, color = treatment)) +
  geom_point() +
  geom_ribbon(data = evAEvAD2Sim, aes(y = pred, ymin = pred_lower, ymax = pred_upper, fill = treatment), color = NA, alpha = 0.5) +
  geom_line(data = evAEvAD2Sim, aes(y = pred))

# total biomass model
evASevChgD2Mod2 <- brm(sev_chg ~ fungicide * tot_bio_s, 
                       data = evAD2Dat, family = "normal",
                       prior <- c(prior(normal(0.6, 1), class = Intercept),
                                  prior(normal(0, 1), class = b)), # use default priors for sigma
                       iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(evASevChgD2Mod2)
pp_check(evASevChgD2Mod2, nsamples = 100)
plot(evASevChgD2Mod2)

# compare with separate biomass model
evASevChgD2Loo <- loo(evASevChgD2Mod, reloo = T)
evASevChgD2Loo2 <- loo(evASevChgD2Mod2, reloo = T)
loo_compare(evASevChgD2Loo, evASevChgD2Loo2)
# model 2 is a better fit

# simulated data
evAD2Sim <- evAD2Dat %>%
  select(treatment, fungicide, tot_bio_s, total_biomass) %>%
  unique() %>%
  mutate(pred = fitted(evASevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Estimate"],
         pred_lower = fitted(evASevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Q2.5"],
         pred_upper = fitted(evASevChgD2Mod2, type = "response", newdata = ., re_formula = NA)[, "Q97.5"])


#### figure ####

# combine data
simDat <- mvMvD2Sim %>%
              mutate(plant_group = "Mv seedling",
                     bio_group = "Mv seedling") %>%
              rename(biomass = Mv_seedling_biomass) %>%
  full_join(evSMvD2Sim %>%
              mutate(plant_group = "Ev seedling",
                     bio_group = "Mv seedling") %>%
              rename(biomass = Mv_seedling_biomass)) %>%
  full_join(evAMvD2Sim %>%
              mutate(plant_group = "Ev adult",
                     bio_group = "Mv seedling") %>%
              rename(biomass = Mv_seedling_biomass)) %>%
  full_join(mvEvSD2Sim %>%
              mutate(plant_group = "Mv seedling",
                     bio_group = "Ev seedling") %>%
              rename(biomass = Ev_seedling_biomass)) %>%
  full_join(evSEvSD2Sim %>%
              mutate(plant_group = "Ev seedling",
                     bio_group = "Ev seedling") %>%
              rename(biomass = Ev_seedling_biomass)) %>%
  full_join(evAEvSD2Sim %>%
              mutate(plant_group = "Ev adult",
                     bio_group = "Ev seedling") %>%
              rename(biomass = Ev_seedling_biomass)) %>%
  full_join(mvEvAD2Sim %>%
              mutate(plant_group = "Mv seedling",
                     bio_group = "Ev adult") %>%
              rename(biomass = Ev_adult_biomass)) %>%
  full_join(evSEvAD2Sim %>%
              mutate(plant_group = "Ev seedling",
                     bio_group = "Ev adult") %>%
              rename(biomass = Ev_adult_biomass)) %>%
  full_join(evAEvAD2Sim %>%
              mutate(plant_group = "Ev adult",
                     bio_group = "Ev adult") %>%
              rename(biomass = Ev_adult_biomass)) %>%
  mutate(plant_group = fct_rev(plant_group),
         bio_group= fct_rev(bio_group))

# figure
pdf("output/severity_change_biomass_2019_density_exp.pdf", width = 10.5, height = 5)
simDat %>%
  filter(treatment == "water") %>%
ggplot(aes(biomass, pred, group = plant_group)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = plant_group), alpha = 0.5, color = NA) +
  geom_line(aes(color = plant_group), size = 1.5) +
  facet_wrap(~bio_group, scales = "free_x", strip.position = "bottom") +
  xlab(expression(paste("Biomass (g ", m^-2, ")", sep = ""))) +
  ylab("Change in disease severity") +
  scale_color_manual(values = c("#BEBEBE", "#A5BDBE", "#407879"), name = "Plant group") +
  scale_fill_manual(values = c("#BEBEBE", "#A5BDBE", "#407879"), name = "Plant group") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(0.075, 0.85),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        strip.placement = "outside")
dev.off()

simDat %>%
  filter(treatment == "fungicide") %>%
  ggplot(aes(biomass, pred, group = interaction(yearF, plant_group))) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = plant_group), alpha = 0.5, color = NA) +
  geom_line(aes(color = plant_group, linetype = yearF)) +
  facet_wrap(~bio_group, scales = "free_x") +
  xlab(expression(paste("Biomass (g ", m^-2, ")", sep = ""))) +
  ylab("Change in disease severity") +
  scale_color_manual(values = c("#BEBEBE", "#A5BDBE", "#407879"), name = "Plant group") +
  scale_fill_manual(values = c("#BEBEBE", "#A5BDBE", "#407879"), name = "Plant group") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Year") +
  theme_bw()
# type of biomass doesn't matter much

# total biomass figure
simTotDat <- mvD1Sim %>%
  mutate(year = 2018,
         plant_group = "Mv seedling") %>%
  full_join(mvD2Sim %>%
              mutate(year = 2019,
                     plant_group = "Mv seedling")) %>%
  full_join(evSD2Sim %>%
              mutate(year = 2019,
                     plant_group = "Ev seedling")) %>%
  full_join(evAD2Sim %>%
              mutate(year = 2019,
                     plant_group = "Ev adult")) %>%
  mutate(plant_group = fct_rev(plant_group),
         yearF = as.factor(year))

simTotDat %>%
  filter(treatment == "water") %>%
  ggplot(aes(total_biomass, pred, group = interaction(yearF, plant_group))) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper, fill = plant_group), alpha = 0.5, color = NA) +
  geom_line(aes(color = plant_group, linetype = yearF)) +
  xlab(expression(paste("Biomass (g ", m^-2, ")", sep = ""))) +
  ylab("Change in disease severity") +
  scale_color_manual(values = c("#BEBEBE", "#A5BDBE", "#407879"), name = "Plant group") +
  scale_fill_manual(values = c("#BEBEBE", "#A5BDBE", "#407879"), name = "Plant group") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Year") +
  theme_bw()


#### output ####

save(mvSevChgD2Mod, file = "output/mv_severity_change_group_biomass_model_2019_density_exp.rda")
save(mvSevChgD2Mod2, file = "output/mv_severity_change_total_biomass_model_2019_density_exp.rda")
save(evSSevChgD2Mod, file = "output/ev_seedling_severity_change_group_biomass_model_2019_density_exp.rda")
save(evSSevChgD2Mod2, file = "output/ev_seedling_severity_change_total_biomass_model_2019_density_exp.rda")
save(evASevChgD2Mod, file = "output/ev_adult_severity_change_group_biomass_model_2019_density_exp.rda")
save(evASevChgD2Mod2, file = "output/ev_adult_severity_change_total_biomass_model_2019_density_exp.rda")