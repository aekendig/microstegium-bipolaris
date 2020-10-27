##### info ####

# file: elymus_litter_establishment_model_2018_greenhouse_exp
# author: Amy Kendig
# date last edited: 10/24/20
# goal: estimate the effect of litter on Elymus establishment using data from Benitez et al.


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(viridis)

# import data
evLitEstDat <- read_csv("../../microstegium-litter-reu/intermediate-data/ev_establishment_data.csv")


#### edit data ####

# add litter mass per m2
evLitEstDat2 <- evLitEstDat %>%
  mutate(litter.g.m2 = case_when(Litter.g == 0.91 ~ 50,
                                  Litter.g == 1.82 ~ 100,
                                  Litter.g == 3.64 ~ 200,
                                  TRUE ~ Litter.g),
         litter.g.cm2 = litter.g.m2/10000)


#### model ####

# initial fit
evLitEstMod1 <- brm(NewGermEv | trials(PlantedEv) ~ litter.g.cm2,
                    data = evLitEstDat2, family = binomial,
                    prior = c(prior(normal(0, 10), class = Intercept),
                              prior(normal(0, 10), class = b)),
                    iter = 6000, warmup = 1000, chains = 1)

summary(evLitEstMod1)

# increase chains
evLitEstMod2 <- update(evLitEstMod1, chains = 3)
summary(evLitEstMod2)
plot(evLitEstMod2)
pp_check(evLitEstMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(litter.g.cm2 = seq(0, 200/10000, length.out = 100)) %>%
  mutate(PlantedEv = 50) %>%
  mutate(Est = fitted(evLitEstMod2, newdata = ., type = "response")[, "Estimate"]/PlantedEv,
         Est_lower = fitted(evLitEstMod2, newdata = ., type = "response")[, "Q2.5"]/PlantedEv,
         Est_upper = fitted(evLitEstMod2, newdata = ., type = "response")[, "Q97.5"]/PlantedEv,
         litter.g.m2 = litter.g.cm2 * 10000)

# figure
evLitEstPlot <- ggplot(evLitEstDat2, aes(x = litter.g.m2, y = PropEstEv)) +
  geom_ribbon(data = fitDat, aes(y = Est, ymin = Est_lower, ymax = Est_upper), alpha = 0.6, color = NA, fill = viridis(0.6)) +
  geom_line(data = fitDat, aes(y = Est), color = viridis(0.6)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1), color = viridis(0.6)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1), color = viridis(0.6)) +
  theme_bw()


#### output ####

save(evLitEstMod2, file = "output/elymus_litter_establishment_model_2018_greenhouse_exp.rda")
save(evLitEstPlot, file = "output/elymus_litter_establishment_figure_2018_greenhouse_exp.rda")
