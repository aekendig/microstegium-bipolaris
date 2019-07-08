##### info ####

# file: mv-fungicide-effects-greenhouse-2019
# author: Amy Kendig
# date last edited: 7/2/19
# goal: see how fungicide treatment affects Mv biomass and seed head production


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

# import data
dat <- read_csv("./data/mv-biomass-seeds-height-jun-2019-fungicide-exp.csv")

# plot parameters
colpal = c("#3CBB75FF", "#39568CFF")
sm_txt = 12
lg_txt = 14
an_txt = 3


#### edit data ####

# remove the dead plant
dat <- dat %>%
  filter(is.na(notes)) %>%
  mutate(log_weight = log(weight.g),
         log_height = log(height.in))


#### visualize ####

# biomass
dat %>%
  ggplot(aes(x = treatment, y = weight.g)) +
  geom_point() +
  stat_summary(geom = "point", fun.y = "mean", color = "red", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "red", width = 0.1) 

# seeds
dat %>%
  ggplot(aes(x = treatment, y = seed_heads)) +
  geom_point() +
  stat_summary(geom = "point", fun.y = "mean", color = "red", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "red", width = 0.1) 

# height
dat %>%
  ggplot(aes(x = treatment, y = height.in)) +
  geom_point() +
  stat_summary(geom = "point", fun.y = "mean", color = "red", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "red", width = 0.1) 


#### statistical models ####

# biomass
mb <- brm(data = dat, family = gaussian,
          weight.g ~ treatment,
          prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b)),
          iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(mb)
plot(mb)
save(mb, file = "./output/mv-fungicide-effects-biomass.rda")
load("./output/mv-fungicide-effects-biomass.rda")

# log-transformed biomass
mbt <- update(mb, formula. = log_weight ~ treatment, newdata = dat)
summary(mbt)
plot(mbt)
save(mbt, file = "./output/mv-fungicide-effects-biomass-transformed.rda")
load("./output/mv-fungicide-effects-biomass-transformed.rda")

# height
mh <- update(mb, formula. = height.in ~ treatment, newdata = dat)
summary(mh)
plot(mh)
save(mh, file = "./output/mv-fungicide-effects-height.rda")
load("./output/mv-fungicide-effects-height.rda")

# log-transformed height
mht <- update(mb, formula. = log_height ~ treatment, newdata = dat)
summary(mht)
plot(mht)
save(mht, file = "./output/mv-fungicide-effects-height-transformed.rda")
load("./output/mv-fungicide-effects-height-transformed.rda")

# seed heads
mean(dat$seed_heads)
var(dat$seed_heads)

ms <- brm(data = dat, family = negbinomial(),
       seed_heads ~ treatment,
       prior <- c(prior(normal(0, 100), class = Intercept),
                  prior(normal(0, 10), class = b)),
       iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(ms)
plot(ms)
save(ms, file = "./output/mv-fungicide-effects-seed-heads.rda")
load("./output/mv-fungicide-effects-seed-heads.rda")


#### check models ####

# posterior predictive check
pp_check(mb, nsamples = 100)
pp_check(mbt, nsamples = 100)
pp_check(mh, nsamples = 100)
pp_check(mht, nsamples = 100)
pp_check(ms, nsamples = 100)

# compare transformed and untransformed
dat <- dat %>%
  mutate(pred_b = fitted(mb)[,1],
         pred_bt = fitted(mbt)[,1],
         pred_h = fitted(mh)[,1],
         pred_ht = fitted(mht)[,1],
         pred_s = fitted(ms)[,1])

cor.test(dat$weight.g, dat$pred_b)
cor.test(dat$log_weight, dat$pred_bt) # both are pretty bad

cor.test(dat$height.in, dat$pred_h)
cor.test(dat$log_height, dat$pred_ht) # both are pretty bad

cor.test(dat$seed_heads, dat$pred_s)


#### figure of coefficients ####

mcmc_intervals(as.array(mb), pars = c("b_Intercept", "b_treatmentwater", "sigma"))
mcmc_intervals(as.array(mbt), pars = c("b_Intercept", "b_treatmentwater", "sigma"))
mcmc_intervals(as.array(mh), pars = c("b_Intercept", "b_treatmentwater", "sigma"))
mcmc_intervals(as.array(mht), pars = c("b_Intercept", "b_treatmentwater", "sigma"))
mcmc_intervals(as.array(ms), pars = c("b_Intercept", "b_treatmentwater", "shape"))


#### final figures ####

#### start here ####

plot(marginal_effects(mbt), 
     points = T,
     errorbar_args = list(width = 0.1),
     point_args = list(alpha = 0.5)) +
  ggplot2::ylim(0, 10)
  
  
  ggplot2::ylab("ln(Microstegium biomass (g))") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.text.x = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.75, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
