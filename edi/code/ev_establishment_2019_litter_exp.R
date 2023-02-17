##### outputs ####

# ev_litter_establishment_bh_model_2019_litter_exp.rda
# ev_litter_establishment_data_2018_litter_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
germ_jun <- read_csv("data/both_germination_disease_jun_2019_litter_exp.csv")
plots <- read_csv("data/litter_weight_apr_2019_litter_exp.csv")


#### edit data ####

# edit plot data
# put the removed litter into the addition litter
# convert units
# remove unnecessary variables
plots2 <- plots %>%
  spread(key = treatment, value = litter_weight.lb) %>%
  mutate(addition = addition + removal,
         removal = 0) %>%
  gather(key = "treatment", value = "litter_weight.lb", -c(date, site, block)) %>%
  mutate(litter_weight.g = litter_weight.lb * 453.592,
         litter.g.m2 = litter_weight.g / 4, # 4 because the plots are 2m^2
         litter.g.cm2 = litter.g.m2 / 10000,
         treatment = fct_relevel(treatment, "removal", "control"),
         plot = as.factor(paste(site, block, sep = "_"))) %>%
  select(-c(date))

# Ev planting data
plant <- tibble(treatment = c("removal", "control", "addition"),
                ev_tot = c(50, 26, 26)) 

# June germination data
# none of the germinants were infected
germ <- germ_jun %>%
  select(-c(date, flag_color, mv_germ, mv_infec)) %>%
  full_join(plots2) %>%
  full_join(plant) %>%
  mutate(ev_prop_germ = ev_germ / ev_tot,
         treatment = fct_relevel(treatment, "removal", "control"))


#### fit model ####

# priors
germ %>%
  filter(litter.g.m2 == 0) %>%
  summarize(mean = mean(ev_prop_germ))

x <- seq(0, 0.4, length.out = 20)
y <- 0.035/(1 + 10 * x)
plot(x, y, type = "l")

val <- seq(0, 20, length.out = 50)
dens <- dexp(val, 0.1)
plot(val, dens, type = "l")

# initial fit
evEstL2BhMod1 <- brm(bf(ev_prop_germ ~ maxEst / (1 + betaL * litter.g.cm2),
                        maxEst ~ 1,
                        betaL ~ 1,
                        nl = T),
                     data = germ, family = gaussian,
                     prior = c(prior(normal(0.04, 1), nlpar = "maxEst", class = "b", lb = 0),
                               prior(exponential(0.1), nlpar = "betaL", lb = 0)),
                     iter = 6000, warmup = 1000, chains = 1)
summary(evEstL2BhMod1)

# increase chains
evEstL2BhMod2 <- update(evEstL2BhMod1, chains = 3,
                        control = list(adapt_delta = 0.99))
summary(evEstL2BhMod2)
plot(evEstL2BhMod2)
pp_check(evEstL2BhMod1, ndraws = 50)


#### output ####

save(evEstL2BhMod2, file = "output/ev_litter_establishment_bh_model_2019_litter_exp.rda")
write_csv(germ, "intermediate-data/ev_litter_establishment_data_2018_litter_exp.csv")
