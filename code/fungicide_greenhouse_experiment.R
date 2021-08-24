#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(ggdist)
library(brms)
library(rladiesgnv)

# import data
dat <- read_csv("intermediate-data/fungicide_greenhouse_experiment.csv")


#### fit model ####

# model
mod <- brm(weight.g ~ fungicide,
           data = dat, family = gaussian,
           prior <- c(prior(normal(24, 10), class = Intercept),
                      prior(normal(0, 10), class = b)),
           iter = 6000, warmup = 1000, chains = 3, cores = 3)

# save model
save(mod, file = "output/fungicide_greenhouse_experiment.rda")

# load model
load("output/fungicide_greenhouse_experiment.rda")

# look at model
summary(mod)
plot(mod)


#### figure ####

# extract posterior samples
samps <- mod %>%
  posterior_samples() %>%
  transmute(control = b_Intercept,
            fungicide = b_Intercept + b_fungicide) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "biomass")

# figure
jpeg("output/fungicide_greenhouse_experiment.jpeg", width = 3.64, height = 2.73, units = "in", res = 600)
ggplot(samps, aes(x = biomass, y = treatment)) +
  stat_halfeye(aes(fill = treatment)) +
  geom_point(data = dat, aes(x = weight.g, y = treatment), position = position_nudge(y = 0.2)) +
  scale_fill_manual(values = gimme_color_codes()[c(6, 3)], guide = "none") +
  labs(x = "Biomass (g)", y = "Treatment") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"))
dev.off()