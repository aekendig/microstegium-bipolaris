##### info ####

# file: fungicide_effect_biomass_simulations
# author: Amy Kendig
# date last edited: 3/12/21
# goal: how do changes in max biomass and competition coefficients affect total biomass?


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# load model for starting values
load("output/mv_biomass_mv_density_model_2019_density_exp.rda")


#### max biomass ####

# constant competition coefficient
alpha <- summary(mvMvD2Mod)$fixed["alpha_treatmentfungicide", "Estimate"] %>%
  round(2)

# max biomass with fungicide
mf <- summary(mvMvD2Mod)$fixed["b0_treatmentfungicide", "Estimate"] %>%
  round(2)

# strongest fungicide effect
maxFE <- 2

# make sure values don't go negative
mf / (1 + maxFE)

# max biomass with water (diseas effect)
mDat <- tibble(fungEff = seq(0, maxFE, by = 0.1)) %>%
  mutate(mf = mf,
         mw = mf / (1 + fungEff)) %>% # calculate max biomass with water based on different fungicide effects
  expand_grid(tibble(N = c(1, 10, 100, 1000))) %>% # range of density values
  mutate(bf = mf / (1 + alpha * N), #  calculate biomass
         bw = mw / (1 + alpha * N),
         bDiff = bf - bw, # fungicide effect on plant biomass
         bFE = bDiff / bw,
         tf = bf * N,
         tw = bw * N,
         totDiff = tf - tw, # fungicide effect on total biomass
         totFE = totDiff / tw) %>%
  pivot_longer(cols = c(mf, mw, bf, bw, tf, tw),
               names_to = c(".value", "treatment"),
               names_pattern = "(.)(.)") %>%
  mutate(treatment = dplyr::recode(treatment, "f" = "fungicide", "w" = "water")) %>%
  rename(maxBio = m, biomass = b, totBio = t)

# individual plant biomass
ggplot(mDat, aes(N, biomass, color = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ fungEff)

# total biomass
ggplot(mDat, aes(N, totBio, color = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ fungEff)


#### competition coefficient ####

# constant max biomass
m <- mf

# competition coefficient with fungicide
alphaf <- alpha

# max biomass with water (diseas effect)
alphaDat <- tibble(fungEff = seq(0, maxFE, by = 0.1)) %>%
  mutate(af = alphaf,
         aw = af * (1 + fungEff)) %>% # calculate max biomass with water based on different fungicide effects
  expand_grid(tibble(N = c(1, 10, 100, 1000))) %>% # range of density values
  mutate(bf = m / (1 + af * N), #  calculate biomass
         bw = m / (1 + aw * N),
         bDiff = bf - bw, # fungicide effect on plant biomass
         bFE = bDiff / bw,
         tf = bf * N,
         tw = bw * N,
         totDiff = tf - tw, # fungicide effect on total biomass
         totFE = totDiff / tw) %>%
  pivot_longer(cols = c(af, aw, bf, bw, tf, tw),
               names_to = c(".value", "treatment"),
               names_pattern = "(.)(.)") %>%
  mutate(treatment = dplyr::recode(treatment, "f" = "fungicide", "w" = "water")) %>%
  rename(alpha = a, biomass = b, totBio = t)

# individual plant biomass
ggplot(alphaDat, aes(N, biomass, color = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ fungEff)

# total biomass
ggplot(alphaDat, aes(N, totBio, color = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ fungEff)


#### both ####

bothDat <- tibble(mFE = seq(0, maxFE, by = 0.5)) %>%
  expand_grid(tibble(aFE = seq(0, maxFE, by = 0.5))) %>%
  mutate(mf = mf,
         af = alphaf,
         mw = mf / (1 + mFE),
         aw = af * (1 + aFE)) %>% # calculate max biomass with water based on different fungicide effects
  expand_grid(tibble(N = c(1, 10, 100, 1000))) %>% # range of density values
  mutate(bf = mf / (1 + af * N), #  calculate biomass
         bw = mw / (1 + aw * N),
         bDiff = bf - bw, # fungicide effect on plant biomass
         bFE = bDiff / bw,
         tf = bf * N,
         tw = bw * N,
         totDiff = tf - tw, # fungicide effect on total biomass
         totFE = totDiff / tw) %>%
  pivot_longer(cols = c(mf, mw, af, aw, bf, bw, tf, tw),
               names_to = c(".value", "treatment"),
               names_pattern = "(.)(.)") %>%
  mutate(treatment = dplyr::recode(treatment, "f" = "fungicide", "w" = "water")) %>%
  rename(maxBio = m, biomass = b, totBio = t)

# individual plant biomass
ggplot(bothDat, aes(N, biomass, color = treatment)) +
  geom_point() +
  geom_line() +
  facet_grid(mFE ~ aFE)

#### combine data ####

comDat <- alphaDat %>%
  mutate(param = "alpha") %>%
  full_join(mDat %>%
              mutate(param = "maxBio"))

# individual plant biomass
ggplot(comDat, aes(N, biomass, color = treatment)) +
  geom_point(aes(shape = param)) +
  geom_line(aes(linetype = param)) +
  facet_wrap(~ fungEff)

# total biomass
ggplot(comDat, aes(N, totBio, color = treatment)) +
  geom_point(aes(shape = param)) +
  geom_line(aes(linetype = param)) +
  facet_wrap(~ fungEff)

# summarise treatment effects by density
comSumDat <- comDat %>%
  select(param, N, fungEff, bDiff, bFE, totDiff, totFE) %>%
  unique()

# absolute difference in plant biomass
ggplot(comSumDat, aes(N, bDiff, color = as.factor(fungEff))) +
  geom_point(aes(shape = param)) +
  geom_line(aes(linetype = param)) +
  scale_x_log10()

# proportion change in plant biomass
ggplot(comSumDat, aes(N, bFE, color = as.factor(fungEff))) +
  geom_point(aes(shape = param)) +
  geom_line(aes(linetype = param)) +
  scale_x_log10()

# absolute difference in total biomass
ggplot(comSumDat, aes(N, totDiff, color = as.factor(fungEff))) +
  geom_point(aes(shape = param)) +
  geom_line(aes(linetype = param)) +
  scale_x_log10()

# proportion change in total biomass
ggplot(comSumDat, aes(N, totFE, color = as.factor(fungEff))) +
  geom_point(aes(shape = param)) +
  geom_line(aes(linetype = param)) +
  scale_x_log10()
