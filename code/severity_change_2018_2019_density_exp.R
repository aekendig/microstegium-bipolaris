##### info ####

# file: severity_change_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/10/21
# goal: change in severity over growing season


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(car) # for logit
library(brms)
library(tidybayes) # for mean_hdi

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
bgSevD1Dat <- read_csv("intermediate-data/ev_background_leaf_scans_2018_density_exp.csv")
# using these because there are so few Ev from this year in the final figure
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") 
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R

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

sevD1Dat2 <- sevD1Dat %>%
  full_join(bgSevD1Dat) %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         logit_severity = logit(severity, adjust = 0.001),
         severity_01 = transform01(severity),
         prop_leaves = leaves_infec/leaves_tot, # for many of the bg plants, data are missing on leaves
         prop_area = lesion_area.pix/leaf_area.pix,
         prop_area = ifelse(prop_area > 1, 1, prop_area),
         month_num = case_when(month == "jul" ~ 0,
                               month == "late_aug" ~ 2,
                               month == "sep" ~ 3))  %>%
  mutate(plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         plant = paste(plant_group, site, plot, treatment, ID, sep = "_"),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  group_by(plant) %>%
  mutate(reps = n()) %>%
  ungroup() %>%
  filter(reps >= 3) # need at least 3 points to fit regressions

# remove plants missing leave counts
sevD1Dat3 <- sevD1Dat2 %>%
  filter(!is.na(leaves_infec)) %>%
  group_by(plant) %>%
  mutate(reps = n()) %>%
  ungroup() %>%
  filter(reps >= 3)

# format severity data
sevD2Dat2 <- sevD2Dat %>%
  filter(!(month %in% c("may", "sep"))) %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         logit_severity = logit(severity, adjust = 0.001),
         severity_01 = transform01(severity),
         prop_leaves = leaves_infec/leaves_tot,
         prop_area = lesion_area.pix/leaf_area.pix,
         prop_area = ifelse(prop_area > 1, 1, prop_area),
         month_num = case_when(month == "jun" ~ 0,
                               month == "jul" ~ 1,
                               month == "early_aug" ~ 2,
                               month == "late_aug" ~ 3))  %>%
  mutate(plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         plant = paste(plant_group, site, plot, treatment, ID, sep = "_"),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  group_by(plant) %>%
  mutate(reps = n()) %>%
  ungroup() %>%
  filter(reps >= 3)

# format edge severity
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(!(month %in% c("may", "sep"))) %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  select(month, site, plot, treatment, edge_severity)

# one site in one month is missing, use nearby plots
edgeSevD2Dat3 <- tibble(site = "D1", plot = 6, treatment = "fungicide", month = "early_aug") %>%
  mutate(edge_severity = edgeSevD2Dat2 %>% 
           filter(site == "D1" & plot %in% c(4, 5) & treatment == "water" & month == "early_aug") %>%
           summarise(sev = mean(edge_severity)) %>%
           pull(sev)) %>%
  full_join(edgeSevD2Dat2)

# plant IDs
plantD1 <- unique(sevD1Dat2$plant)
plantD2 <- unique(sevD2Dat2$plant)

# separate by plant group
mvSevD2Dat <- sevD2Dat2 %>%
  filter(plant_group == "Mv_seedling") %>%
  left_join(edgeSevD2Dat3) %>%
  mutate(plant = paste(site, plot, treatment, ID, sep = "_"),
         edge_sev_s = scale(edge_severity))

mvSevD1Dat <- sevD1Dat3 %>%
  filter(plant_group == "Mv_seedling") %>%
  mutate(plant = paste(site, plot, treatment, ID, sep = "_"))

evSSevD2Dat <- sevD2Dat2 %>%
  filter(plant_group == "Ev_seedling") %>%
  left_join(edgeSevD2Dat3) %>%
  mutate(plant = paste(site, plot, treatment, ID, sep = "_"),
         edge_sev_s = scale(edge_severity))

evSSevD1Dat <- sevD1Dat3 %>%
  filter(plant_group == "Ev_seedling") %>%
  mutate(plant = paste(site, plot, treatment, ID, sep = "_"))

evASevD2Dat <- sevD2Dat2 %>%
  filter(plant_group == "Ev_adult") %>%
  left_join(edgeSevD2Dat3) %>%
  mutate(plant = paste(site, plot, treatment, ID, sep = "_"),
         edge_sev_s = scale(edge_severity))

evASevD1Dat <- sevD1Dat3 %>%
  filter(plant_group == "Ev_adult") %>%
  mutate(plant = paste(site, plot, treatment, ID, sep = "_"))

# combine both datasets
sevDat <- sevD1Dat3 %>%
  mutate(year = 2018) %>%
  full_join(sevD2Dat2 %>%
              mutate(year = 2019) %>%
              left_join(edgeSevD2Dat3))


#### initial visualizations ####

# look at individual regressions
pdf("output/severity_change_by_plant_2018_density_exp.pdf")
for(i in plantD1){
  subDat <- filter(sevD1Dat2, plant == i)
  
  print(ggplot(subDat, aes(x = month_num, y = logit_severity)) +
          geom_point() +
          stat_smooth(method = "lm", formula = y ~ x) +
          ggtitle(i) +
          theme_bw())
}
dev.off()
# really high standard error, but linear model seems appropriate

pdf("output/severity_change_by_plant_2019_density_exp.pdf")
for(i in plantD2){
  subDat <- filter(sevD2Dat2, plant == i)
  
  print(ggplot(subDat, aes(x = month_num, y = logit_severity)) +
          geom_point() +
          stat_smooth(method = "lm", formula = y ~ x) +
          ggtitle(i) +
          theme_bw())
}
dev.off()
# trends vary in shape, so linear model seems most appropriate

pdf("output/leaves_change_by_plant_2018_density_exp.pdf")
for(i in plantD1){
  subDat <- filter(sevD1Dat2, plant == i)
  
  print(ggplot(subDat, aes(x = month_num, y = prop_leaves)) +
          geom_point() +
          geom_smooth(method = "glm", formula = y ~ x, 
                      method.args = list(family = "binomial"), 
                      aes(weight = leaves_tot)) +
          ggtitle(i) +
          theme_bw())
}
dev.off()
# really high standard errors, may be difficult to fit these with 3 points

pdf("output/leaves_change_by_plant_2019_density_exp.pdf")
for(i in plantD2){
  subDat <- filter(sevD2Dat2, plant == i)
  
  print(ggplot(subDat, aes(x = month_num, y = prop_leaves)) +
          geom_point() +
          geom_smooth(method = "glm", formula = y ~ x, 
                      method.args = list(family = "binomial"), 
                      aes(weight = leaves_tot)) +
          ggtitle(i) +
          theme_bw())
}
dev.off()
# seem like better fits for Ev than Mv

pdf("output/area_change_by_plant_2018_density_exp.pdf")
for(i in plantD1){
  subDat <- filter(sevD1Dat2, plant == i)
  
  print(ggplot(subDat, aes(x = month_num, y = prop_area)) +
          geom_point() +
          geom_smooth(method = "glm", formula = y ~ x, 
                      method.args = list(family = "binomial"), 
                      aes(weight = leaf_area.pix)) +
          ggtitle(i) +
          theme_bw())
}
dev.off()
# really low SE because leaf_area.pix values are so high
# okay fits, some of the middle measurements are really far

pdf("output/area_change_by_plant_2019_density_exp.pdf")
for(i in plantD2){
  subDat <- filter(sevD2Dat2, plant == i)
  
  print(ggplot(subDat, aes(x = month_num, y = prop_area)) +
          geom_point() +
          geom_smooth(method = "glm", formula = y ~ x, 
                      method.args = list(family = "binomial"), 
                      aes(weight = leaf_area.pix)) +
          ggtitle(i) +
          theme_bw())
}
dev.off()
# pretty consistent accumulation because so many final measurements are 1
# fit fails for some of the plants

# general trends in severity
ggplot(sevD1Dat3, aes(x = month_num, y = logit_severity, color = treatment)) +
  stat_smooth(method = "lm", se = F, size = 0.3, geom = "line", alpha = 0.3, formula = y ~ x,
              aes(group = plant)) +
  stat_smooth(method = "lm", se = F, formula = y ~ x) +
  facet_wrap(~ plant_group) +
  theme_bw()

ggplot(sevD2Dat2, aes(x = month_num, y = logit_severity, color = treatment)) +
  stat_smooth(method = "lm", se = F, size = 0.3, geom = "line", alpha = 0.2, formula = y ~ x,
              aes(group = plant)) +
  stat_smooth(method = "lm", se = F, formula = y ~ x) +
  facet_wrap(~ plant_group) +
  theme_bw()
# predicts negative severity values because logit_severity goes below minimum (which is 0 severity)
min(sevD2Dat2$logit_severity)

# general trends in leaves infected
ggplot(sevD1Dat3, aes(x = month_num, y = prop_leaves, color = treatment, weight = leaves_tot)) +
  geom_line(stat = "smooth", method = "glm", formula = y ~ x,
              size = 0.3, alpha = 0.3,
              method.args = list(family = "binomial"),
              aes(group = plant)) +
  geom_smooth(method = "glm", formula = y ~ x, se = F, 
              method.args = list(family = "binomial")) +
  facet_wrap(~ plant_group) +
  theme_bw()

ggplot(sevD2Dat2, aes(x = month_num, y = prop_leaves, color = treatment, weight = leaves_tot)) +
  geom_line(stat = "smooth", method = "glm", formula = y ~ x,
            size = 0.3, alpha = 0.3,
            method.args = list(family = "binomial"),
            aes(group = plant)) +
  geom_smooth(method = "glm", formula = y ~ x, se = F, 
              method.args = list(family = "binomial")) +
  facet_wrap(~ plant_group) +
  theme_bw()


# general trends in area infected
ggplot(sevD1Dat2, aes(x = month_num, y = prop_area, color = treatment, weight = leaf_area.pix)) +
  geom_line(stat = "smooth", method = "glm", formula = y ~ x,
            size = 0.3, alpha = 0.3,
            method.args = list(family = "binomial"),
            aes(group = plant)) +
  geom_smooth(method = "glm", formula = y ~ x, se = F, 
              method.args = list(family = "binomial")) +
  facet_wrap(~ plant_group) +
  theme_bw()

ggplot(sevD2Dat2, aes(x = month_num, y = prop_area, color = treatment, weight = leaf_area.pix)) +
  geom_line(stat = "smooth", method = "glm", formula = y ~ x,
            size = 0.3, alpha = 0.3,
            method.args = list(family = "binomial"),
            aes(group = plant)) +
  geom_smooth(method = "glm", formula = y ~ x, se = F, 
              method.args = list(family = "binomial")) +
  facet_wrap(~ plant_group) +
  theme_bw()
# can't see individual plants because fit fails for some

# edge severity
sevDat %>%
  filter(year == 2019) %>%
  ggplot(aes(edge_severity, severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm")


#### Mv models ####

# 2019 model
mvSevD2Mod <- brm(severity_01 ~ month_num * fungicide + edge_sev_s + (month_num||plant), 
                  data = mvSevD2Dat, family = "beta",
                  prior <- c(prior(normal(0, 1), class = Intercept),
                             prior(normal(0, 1), class = b)), # use default priors for phi and random effects
                  iter = 6000, warmup = 1000, chains = 3, cores = 2, 
                  control = list(adapt_delta = 0.99))
prior_summary(mvSevD2Mod)
summary(mvSevD2Mod)
pp_check(mvSevD2Mod, nsamples = 100)
plot(mvSevD2Mod)

# 2018 model
mvSevD1Mod <- brm(severity_01 ~ month_num * fungicide + (month_num||plant), 
                  data = mvSevD2Dat, family = "beta",
                  prior <- c(prior(normal(0, 1), class = Intercept),
                             prior(normal(0, 1), class = b)), # use default priors for phi and random effects
                  iter = 6000, warmup = 1000, chains = 3, cores = 2, 
                  control = list(adapt_delta = 0.99))
summary(mvSevD1Mod)
pp_check(mvSevD1Mod, nsamples = 100)
plot(mvSevD1Mod)

# predicted vs. measured
mvSevD1Dat %>%
  group_by(fungicide, month_num) %>%
  summarise(severity = mean(severity)) %>%
  ungroup() %>%
  mutate(edge_sev_s = 0) %>%
  mutate(pred_severity = fitted(mvSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# under estimates

mvSevD2Dat %>%
  mutate(pred_severity = fitted(mvSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# overestimates at low values and under estimates at high

# simulated data
mvSevD1Sim <- mvSevD1Dat %>%
  select(fungicide, treatment) %>%
  unique() %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(mvSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(mvSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(mvSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

mvSevD2Sim <- mvSevD2Dat %>%
  select(fungicide, treatment) %>%
  unique() %>%
  mutate(edge_sev_s = 0) %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(mvSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(mvSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(mvSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

# model fit
ggplot(mvSevD1Dat, aes(x = month_num, y = severity_01, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = mvSevD1Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = mvSevD1Sim, aes(y = pred_severity)) +
theme_bw()
# hits middle values well

ggplot(mvSevD2Dat, aes(x = month_num, y = severity_01, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = mvSevD2Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = mvSevD2Sim, aes(y = pred_severity)) +
  theme_bw()
# better fit than all data model

# Mv fungicide effects at 3 months
posterior_samples(mvSevD1Mod) %>%
  rename(b_month_num_fungicide = "b_month_num:fungicide") %>%
  transmute(water = logit2prop(b_Intercept + b_month_num * 3) %>%
              backtransform01,
            fungicide = logit2prop(b_Intercept + b_month_num * 3 + b_fungicide + b_month_num_fungicide * 3) %>%
              backtransform01) %>%
  mutate(fungicide_effect = (fungicide-water)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "severity") %>%
  group_by(treatment) %>%
  mean_hdi(severity)

posterior_samples(mvSevD2Mod) %>%
  rename(b_month_num_fungicide = "b_month_num:fungicide") %>%
  transmute(water = logit2prop(b_Intercept + b_month_num * 3) %>%
              backtransform01,
            fungicide = logit2prop(b_Intercept + b_month_num * 3 + b_fungicide + b_month_num_fungicide * 3) %>%
              backtransform01) %>%
  mutate(fungicide_effect = (fungicide-water)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "severity") %>%
  group_by(treatment) %>%
  mean_hdi(severity)


#### Ev seedling models ####

# 2019 model
evSSevD2Mod <- update(mvSevD2Mod, newdata = evSSevD2Dat)
summary(evSSevD2Mod)
pp_check(evSSevD2Mod, nsamples = 100)
plot(evSSevD2Mod)

# 2018 model
evSSevD1Mod <- update(mvSevD1Mod, newdata = evSSevD1Dat)
summary(evSSevD1Mod)
pp_check(evSSevD1Mod, nsamples = 100)
plot(evSSevD1Mod)

# predicted vs. measured
evSSevD1Dat %>%
  group_by(fungicide, month_num) %>%
  summarise(severity = mean(severity)) %>%
  ungroup() %>%
  mutate(pred_severity = fitted(evSSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# under estimates

evSSevD2Dat %>%
  mutate(pred_severity = fitted(evSSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# good estimates

# simulated data
evSSevD1Sim <- evSSevD1Dat %>%
  select(fungicide, treatment) %>%
  unique() %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(evSSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(evSSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(evSSevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

evSSevD2Sim <- evSSevD2Dat %>%
  select(fungicide, treatment) %>%
  unique() %>%
  mutate(edge_sev_s = 0) %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(evSSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(evSSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(evSSevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

# model fit
ggplot(evSSevD1Dat, aes(x = month_num, y = severity_01, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = evSSevD1Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = evSSevD1Sim, aes(y = pred_severity)) +
  theme_bw()
# no change in severity

ggplot(evSSevD2Dat, aes(x = month_num, y = severity_01, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = evSSevD2Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = evSSevD2Sim, aes(y = pred_severity)) +
  theme_bw()
# better fit than all data model


#### Ev adult models ####

# 2019 model
evASevD2Mod <- update(mvSevD2Mod, newdata = evASevD2Dat)
summary(evASevD2Mod)
pp_check(evASevD2Mod, nsamples = 100)
plot(evASevD2Mod)

# 2018 model
evASevD1Mod <- update(mvSevD1Mod, newdata = evASevD1Dat)
summary(evASevD1Mod)
pp_check(evASevD1Mod, nsamples = 100)
plot(evASevD1Mod)

# predicted vs. measured
evASevD1Dat %>%
  group_by(fungicide, month_num) %>%
  summarise(severity = mean(severity)) %>%
  ungroup() %>%
  mutate(pred_severity = fitted(evASevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# under estimates

evASevD2Dat %>%
  mutate(pred_severity = fitted(evASevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# generally good

# simulated data
evASevD1Sim <- evASevD1Dat %>%
  select(fungicide, treatment) %>%
  unique() %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(evASevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(evASevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(evASevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

evASevD2Sim <- evASevD2Dat %>%
  select(fungicide, treatment) %>%
  unique() %>%
  mutate(edge_sev_s = 0) %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(evASevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(evASevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(evASevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

# model fit
ggplot(evASevD1Dat, aes(x = month_num, y = severity_01, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = evASevD1Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = evASevD1Sim, aes(y = pred_severity)) +
  theme_bw()
# hits middle values well

ggplot(evASevD2Dat, aes(x = month_num, y = severity_01, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = evASevD2Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = evASevD2Sim, aes(y = pred_severity)) +
  theme_bw()


### extract plant-specific slopes ####
mvSevD1RanSlopes <- coef(mvSevD1Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(year = 2018,
         plant = rownames(coef(mvSevD1Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

mvSevD2RanSlopes <- coef(mvSevD2Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(year = 2019,
         plant = rownames(coef(mvSevD2Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

evSSevD1RanSlopes <- coef(evSSevD1Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(year = 2018,
         plant = rownames(coef(evSSevD1Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

evSSevD2RanSlopes <- coef(evSSevD2Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(year = 2019,
         plant = rownames(coef(evSSevD2Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

evASevD1RanSlopes <- coef(evASevD1Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(year = 2018,
         plant = rownames(coef(evASevD1Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

evASevD2RanSlopes <- coef(evASevD2Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(year = 2019,
         plant = rownames(coef(evASevD2Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

# combine
sevRanSlopes <- mvSevD1RanSlopes %>%
  full_join(mvSevD2RanSlopes) %>%
  mutate(plant_group = "Mv_seedling",
         sp = "Mv",
         age = "seedling") %>%
  full_join(evSSevD1RanSlopes %>%
              full_join(evSSevD2RanSlopes) %>%
              mutate(plant_group = "Ev_seedling",
                     sp = "Ev",
                     age = "seedling")) %>%
  full_join(evASevD1RanSlopes %>%
              full_join(evASevD2RanSlopes) %>%
              mutate(plant_group = "Ev_adult",
                     sp = "Ev",
                     age = "adult")) %>%
  rowwise() %>%
  mutate(site = strsplit(plant, "_")[[1]][1],
         plot = strsplit(plant, "_")[[1]][2] %>%
           as.numeric(),
         treatment = strsplit(plant, "_")[[1]][3],
         ID = strsplit(plant, "_")[[1]][4]) %>%
  ungroup() 


#### figure ####

# combine simulation data
sevSim <- mvSevD1Sim %>%
  mutate(year = 2018) %>%
  full_join(mvSevD2Sim %>%
              mutate(year = 2019)) %>%
  mutate(plant_group = "Mv seedling") %>%
  full_join(evSSevD1Sim %>%
              mutate(year = 2018) %>%
              full_join(evSSevD2Sim %>%
                          mutate(year = 2019)) %>%
              mutate(plant_group = "Ev seedling")) %>%
  full_join(evASevD1Sim %>%
              mutate(year = 2018) %>%
              full_join(evASevD2Sim %>%
                          mutate(year = 2019)) %>%
              mutate(plant_group = "Ev adult")) %>%
  mutate(plant_group = fct_rev(plant_group),
         yearF = as.factor(year),
         treatment = dplyr::recode(treatment, "water" = "water (control)") %>%
           fct_rev(),
         month_num = case_when(year == 2018 ~ month_num + 1,
                               TRUE ~ month_num))

# edit plant group in data
sevDat2 <- sevDat %>%
  mutate(plant_group = str_replace(plant_group, "_", " ") %>%
           fct_rev(),
         yearF = as.factor(year),
         treatment = dplyr::recode(treatment, "water" = "water (control)") %>%
           fct_rev(),
         month_num = case_when(year == 2018 ~ month_num + 1,
                               TRUE ~ month_num))

# figure
pdf("output/severity_change_2018_2019_density_exp.pdf", width = 10.5, height = 5)
ggplot(sevDat2, aes(month_num, severity_01, group = interaction(yearF, treatment))) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3), aes(color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2,  position = position_dodge(0.3), aes(color = treatment, shape = yearF)) +
  geom_ribbon(data = sevSim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = sevSim, aes(y = pred_severity, color = treatment, linetype = yearF)) +
  geom_text(x = 0, y = 0.85, check_overlap = T, hjust = 0, aes(label = plant_group), size = 5) +
  facet_wrap(~plant_group) +
  xlab("Month") +
  ylab("Disease severity") +
  scale_color_viridis_d(end = 0.6, name = "Treatment") +
  scale_fill_viridis_d(end = 0.6, name = "Treatment") +
  scale_shape_manual(values = c(17, 19), name = "Year") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Year") +
  scale_x_continuous(labels = c("Jun", "Jul", "early Aug", "late Aug", "Sep")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(0.1, 0.7),
        strip.background = element_blank(),
        strip.text = element_blank())
dev.off()


#### output ####
write_csv(sevRanSlopes, "intermediate-data/severity_change_model_ran_slopes_2018_2019_density_exp.csv")
save(mvSevD1Mod, file = "output/mv_severity_change_model_2018_density_exp.rda")
save(mvSevD2Mod, file = "output/mv_severity_change_model_2019_density_exp.rda")
save(evSSevD1Mod, file = "output/ev_seedling_severity_change_model_2018_density_exp.rda")
save(evSSevD2Mod, file = "output/ev_seedling_severity_change_model_2019_density_exp.rda")
save(evASevD1Mod, file = "output/ev_adult_severity_change_model_2018_density_exp.rda")
save(evASevD2Mod, file = "output/ev_adult_severity_change_model_2019_density_exp.rda")


#### old code ####

# models by plant group produced more accurate predictions

#### model with all data ####

# intercept prior
sevD2Dat2 %>%
  filter(month_num == 0 & treatment == "water" & plant_group == "Mv_seedling") %>%
  summarise(sev = mean(severity_01))
# 0.03

# 2019 model
sevD2Mod <- brm(severity_01 ~ month_num * fungicide * plant_group + (month_num||plant), 
                data = sevD2Dat2, family = "beta",
                prior <- c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b)), # use default priors for phi and random effects
                iter = 6000, warmup = 1000, chains = 3)
prior_summary(sevD2Mod)
summary(sevD2Mod)
pp_check(sevD2Mod, nsamples = 100)

# 2018 model
sevD1Mod <- update(sevD2Mod, newdata = sevD1Dat3)
prior_summary(sevD1Mod)
summary(sevD1Mod)
pp_check(sevD1Mod, nsamples = 100)
# overestimates low values


#### all data model output ####

# predicted vs. measured
sevD1Dat3 %>%
  group_by(plant_group, fungicide, month_num) %>%
  summarise(severity = mean(severity)) %>%
  ungroup() %>%
  mutate(pred_severity = fitted(sevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# poor fit at high severity

sevD2Dat2 %>%
  group_by(plant_group, fungicide, month_num) %>%
  summarise(severity = mean(severity)) %>%
  ungroup() %>%
  mutate(pred_severity = fitted(sevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01()) %>%
  ggplot(aes(severity, pred_severity, color = month_num)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
# poor fit at low month values

# Mv fungicide effects at 3 months
posterior_samples(sevD1Mod) %>%
  rename(b_month_num_fungicide = "b_month_num:fungicide") %>%
  transmute(water = logit2prop(b_Intercept + b_month_num * 3) %>%
              backtransform01,
            fungicide = logit2prop(b_Intercept + b_month_num * 3 + b_fungicide + b_month_num_fungicide * 3) %>%
              backtransform01) %>%
  mutate(fungicide_effect = (water-fungicide)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "severity") %>%
  group_by(treatment) %>%
  mean_hdi(severity)

posterior_samples(sevD2Mod) %>%
  rename(b_month_num_fungicide = "b_month_num:fungicide") %>%
  transmute(water = logit2prop(b_Intercept + b_month_num * 3) %>%
              backtransform01,
            fungicide = logit2prop(b_Intercept + b_month_num * 3 + b_fungicide + b_month_num_fungicide * 3) %>%
              backtransform01) %>%
  mutate(fungicide_effect = (water-fungicide)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "severity") %>%
  group_by(treatment) %>%
  mean_hdi(severity)

# compare to raw data
sevD1Dat3 %>%
  filter(plant_group == "Mv_seedling" & month_num == 3) %>%
  group_by(fungicide) %>%
  summarise(severity = mean(severity))

sevD2Dat2 %>%
  filter(plant_group == "Mv_seedling" & month_num == 3) %>%
  group_by(fungicide) %>%
  summarise(severity = mean(severity))

# simulated data
sevD1Sim <- sevD1Dat3 %>%
  select(plant_group, fungicide, treatment) %>%
  unique() %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(sevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(sevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(sevD1Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

sevD2Sim <- sevD2Dat2 %>%
  select(plant_group, fungicide, treatment) %>%
  unique() %>%
  expand_grid(tibble(month_num = seq(0, 3, length.out = 20))) %>%
  mutate(pred_severity = fitted(sevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Estimate"] %>%
           backtransform01(),
         pred_lower = fitted(sevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q2.5"] %>%
           backtransform01(),
         pred_upper = fitted(sevD2Mod, type = "response", newdata = ., re_formula = NA)[, "Q97.5"] %>%
           backtransform01())

# model fit
ggplot(sevD1Dat3, aes(x = month_num, y = severity, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = sevD1Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = sevD1Sim, aes(y = pred_severity)) +
  facet_wrap(~ plant_group) +
  theme_bw()
# less clear disease accumulation in Ev compared to 2019

ggplot(sevD2Dat2, aes(x = month_num, y = severity, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = sevD2Sim, aes(y = pred_severity, ymin = pred_lower, ymax = pred_upper), alpha = 0.5, color = NA) +
  geom_line(data = sevD2Sim, aes(y = pred_severity)) +
  facet_wrap(~ plant_group) +
  theme_bw()

# extract random slopes
sevD1RanSlopes <- ranef(sevD1Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(plant = rownames(ranef(sevD1Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)

sevD2RanSlopes <- ranef(sevD2Mod, pars = "month_num")[[1]] %>%
  as_tibble() %>%
  mutate(plant = rownames(ranef(sevD2Mod, pars = "month_num")[[1]])) %>%
  rename(slope = Estimate.month_num,
         error = Est.Error.month_num,
         lower = Q2.5.month_num,
         upper = Q97.5.month_num)


#### all data model output ####
save(sevD1Mod, file = "output/severity_change_model_2018_density_exp.rda")
write_csv(sevD1RanSlopes, "intermediate-data/severity_change_model_ran_slopes_2018_density_exp.csv")

save(sevD2Mod, file = "output/severity_change_model_2019_density_exp.rda")
write_csv(sevD2RanSlopes, "intermediate-data/severity_change_model_ran_slopes_2019_density_exp.csv")