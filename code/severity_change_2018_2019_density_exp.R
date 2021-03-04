##### info ####

# file: severity_change_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/3/21
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

# plant IDs
plantD1 <- unique(sevD1Dat2$plant)
plantD2 <- unique(sevD2Dat2$plant)


#### figure ####

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


#### statistics ####

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


#### model output ####

# fungicide effect
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

# Mvfungicide effects at 3 months
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
 

#### output ####
save(sevD1Mod, file = "output/severity_change_model_2018_density_exp.rda")
write_csv(sevD1RanSlopes, "intermediate-data/severity_change_model_ran_slopes_2018_density_exp.csv")

save(sevD2Mod, file = "output/severity_change_model_2019_density_exp.rda")
write_csv(sevD2RanSlopes, "intermediate-data/severity_change_model_ran_slopes_2019_density_exp.csv")
