##### outputs ####

# Tables S20-S22


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(broom.mixed)

# import data
mv_dat <- read_csv("data/mv_biomass_seeds_height_jun_2019_fungicide_exp.csv")
ev_dat <- read_csv("data/ev_biomass_dec_2019_fungicide_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}


#### edit data ####

# look at notes
unique(mv_dat$notes) # dead plant
unique(ev_dat$notes)

# remove the dead plant
# fungicide column
mv_dat2 <- mv_dat %>%
  filter(is.na(notes)) %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0"),
         log_bio.g = log(weight.g))

# combine live and dead weight
# fungicide column
ev_dat2 <- ev_dat %>%
  group_by(treatment, pot, sp) %>%
  summarise(weight.g = sum(weight.g)) %>%
  ungroup() %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0"),
         log_bio.g = log(weight.g))


#### visualize ####

# Mv biomass
mv_dat2 %>%
  ggplot(aes(x = treatment, y = weight.g)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " aboveground biomass (g)", sep = ""))) +
  xlab("Treatment")

# Ev biomass
ev_dat2 %>%
  ggplot(aes(x = treatment, y = weight.g)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Elymus), " aboveground biomass (g)", sep = ""))) +
  xlab("Treatment")

# seeds
mv_dat2 %>%
  ggplot(aes(x = treatment, y = seed_heads)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " seed heads", sep = ""))) +
  xlab("Treatment")

# height
mv_dat2 %>%
  ggplot(aes(x = treatment, y = height.in)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " height (in.)", sep = ""))) +
  xlab("Treatment")

# biomass/seeds
mv_dat2 %>%
  ggplot(aes(x = weight.g, y = seed_heads, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab(expression(paste(italic(Microstegium), " aboveground biomass (g)", sep = ""))) +
  ylab(expression(paste(italic(Microstegium), " seed heads", sep = "")))

# log-transformed values
mv_dat2 %>%
  ggplot(aes(x = treatment, y = log_bio.g)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " aboveground biomass (log-g)", sep = ""))) +
  xlab("Treatment")

ev_dat2 %>%
  ggplot(aes(x = treatment, y = log_bio.g)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Elymus), " aboveground biomass (log-g)", sep = ""))) +
  xlab("Treatment")


#### statistical models ####

# Mv biomass
mvBioGhMod <- brm(log_bio.g ~ fungicide,
                  data = mv_dat2, family = gaussian,
                  prior <- c(prior(normal(0, 10), class = Intercept),
                             prior(normal(0, 10), class = b)),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvBioGhMod)

# effect size
mv_dat2 %>%
  group_by(fungicide) %>%
  summarise(bio = mean(weight.g)) %>%
  ungroup()

hypothesis(mvBioGhMod, "exp(Intercept + fungicide1) - exp(Intercept) = 0")

# save
save(mvBioGhMod, file = "output/mv_biomass_fungicide_effects_greenhouse_2019.rda")
write_csv(tidy(mvBioGhMod), "output/mv_biomass_fungicide_effects_greenhouse_2019.csv")

# Ev biomass
evBioGhMod <- brm(log_bio.g ~ fungicide,
                  data = ev_dat2, family = gaussian,
                  prior <- c(prior(normal(0, 10), class = Intercept),
                             prior(normal(0, 10), class = b)),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(evBioGhMod)

# effect size
ev_dat2 %>%
  group_by(fungicide) %>%
  summarise(bio = mean(weight.g)) %>%
  ungroup()

hypothesis(evBioGhMod, "exp(Intercept + fungicide1) - exp(Intercept) = 0")

# save
save(evBioGhMod, file = "output/ev_biomass_fungicide_effects_greenhouse_2019.rda")
write_csv(tidy(evBioGhMod), "output/ev_biomass_fungicide_effects_greenhouse_2019.csv")
  
# Mv seed heads
mean(mv_dat2$seed_heads)
var(mv_dat2$seed_heads) # twice as large

mvSeedGhMod <- brm(seed_heads ~ fungicide,
                   data = mv_dat2, family = negbinomial,
                   prior <- c(prior(normal(0, 10), class = Intercept),
                              prior(normal(0, 10), class = b)),
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvSeedGhMod)

# effect size
mv_dat2 %>%
  group_by(fungicide) %>%
  summarise(seeds = mean(seed_heads)) %>%
  ungroup()

hypothesis(mvSeedGhMod, "exp(Intercept + fungicide1) - exp(Intercept) = 0")

# save
save(mvSeedGhMod, file = "output/mv_seed_fungicide_effects_greenhouse_2019.rda")
write_csv(tidy(mvSeedGhMod), "output/mv_seed_fungicide_effects_greenhouse_2019.csv")
