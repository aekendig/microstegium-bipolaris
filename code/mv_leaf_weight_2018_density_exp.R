##### info ####

# file: mv_leaf_weight_2018_density_exp
# author: Amy Kendig
# date last edited: 3/28/21
# goal: analyses of leaf weight


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import data
mvLeafDat <- read_csv("data/mv_leaf_weight_2018_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")

# fungicide effect function
post_pred_fun <- function(mod, dat){
  
  posterior_samples(mod) %>%
    transmute(water = exp(b_Intercept),
              fung = exp(b_Intercept + b_fungicide),
              perc = 100 * (fung - water) / water) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% mean_hdi(effect)
  
}


#### edit data ####

unique(mvLeafDat$notes)

mvLeafD1Dat <- mvLeafDat %>%
  filter(!is.na(leaf_weight.g) & (notes != "Leaf was mutiliated" | is.na(notes))) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = ""),
         log_leaf_weight = log(leaf_weight.g))


#### 2018 Mv model ####

# initial visualization
ggplot(mvLeafD1Dat, aes(treatment, leaf_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(mvLeafD1Dat, aes(leaf_weight.g)) +
  geom_density() +
  facet_wrap(~ treatment)

ggplot(mvLeafD1Dat, aes(log_leaf_weight)) +
  geom_density() +
  facet_wrap(~ treatment)

# model
mvLeafD1Mod <- brm(data = mvLeafD1Dat, family = gaussian,
                   log_leaf_weight ~ fungicide + (1|site/plotf),
                   prior <- c(prior(normal(-5, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                   control = list(adapt_delta = 0.999))
mod_check_fun(mvLeafD1Mod)

# posterior means
post_pred_fun(mvLeafD1Mod)

# save
save(mvLeafD1Mod, file = "output/mv_leaf_weight_model_2018_density_exp.rda")

