##### info ####

# file: elymus_litter_establishment_bh_model_2018_greenhouse_exp
# author: Amy Kendig
# date last edited: 10/27/20
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
evLitEstBhMod1 <- brm(bf(PropEstEv ~ maxEst / (1 + betaL * litter.g.cm2),
                       maxEst ~ 1,
                       betaL ~ 1,
                       nl = T),
                    data = evLitEstDat2, family = gaussian,
                    prior = c(prior(normal(1, 1), nlpar = "maxEst", ub = 1),
                              prior(exponential(0.5), nlpar = "betaL", lb = 0),
                              prior(cauchy(0, 1), class = "sigma")),
                    iter = 6000, warmup = 1000, chains = 1)
# 1 divergent transition
summary(evLitEstBhMod1)

# increase chains
evLitEstBhMod2 <- update(evLitEstBhMod1, chains = 3,
                         control = list(adapt_delta = 0.99))
summary(evLitEstBhMod2)
plot(evLitEstBhMod2)
pp_check(evLitEstBhMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(litter.g.cm2 = seq(0, 200/10000, length.out = 100)) %>%
  mutate(Est = fitted(evLitEstBhMod2, newdata = ., type = "response")[, "Estimate"],
         Est_lower = fitted(evLitEstBhMod2, newdata = ., type = "response")[, "Q2.5"],
         Est_upper = fitted(evLitEstBhMod2, newdata = ., type = "response")[, "Q97.5"],
         litter.g.m2 = litter.g.cm2 * 10000)

# figure
evLitEstBhPlot <- ggplot(evLitEstDat2, aes(x = litter.g.m2, y = PropEstEv)) +
  geom_ribbon(data = fitDat, aes(y = Est, ymin = Est_lower, ymax = Est_upper), alpha = 0.6, color = NA, fill = viridis(0.6)) +
  geom_line(data = fitDat, aes(y = Est), color = viridis(0.6)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1), color = viridis(0.6)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1), color = viridis(0.6)) +
  theme_bw()


#### output ####

save(evLitEstBhMod2, file = "output/elymus_litter_establishment_bh_model_2018_greenhouse_exp.rda")
save(evLitEstBhPlot, file = "output/elymus_litter_establishment_bh_figure_2018_greenhouse_exp.rda")
