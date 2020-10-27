##### info ####

# file: elymus_seedling_germination_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/17/20
# goal: analyze Elymus germination


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
germD2Dat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")


#### edit data ####

# check for unintuitive counts
filter(germD2Dat, week_2_emerg > seeds_planted | week_3_emerg > seeds_planted | week_4_emerg > seeds_planted)
filter(germD2Dat, week_4_emerg < week_3_emerg) # two cases, lost one, use week 3
filter(germD2Dat, week_3_cut_tops > week_2_emerg | week_4_cut_tops > week_2_emerg) # this cut-top count must be wrong (week 4)
filter(germD2Dat, week_4_cut_tops > week_3_cut_tops) # the one above, plus another that may have been mis-counted in week 3

# select data from year 2
# correct reduction in emergents between week 3 and 4
# correct the increase in cut-tops
# correct repair of cut tops between weeks 3 and 4
evGermD2Dat <- germD2Dat %>%
  filter(seeds_planted > 0 & year == 2019) %>%
  mutate(week_4_emerg = case_when(week_4_emerg < week_3_emerg ~ week_3_emerg,
                                  TRUE ~ week_4_emerg),
         week_3_cut_tops = case_when(week_4_cut_tops > week_3_cut_tops & week_4_cut_tops <= week_2_emerg ~ week_4_cut_tops,
                                     TRUE ~ week_3_cut_tops),
         week_4_cut_tops = case_when(week_4_cut_tops > week_2_emerg ~ week_3_cut_tops,
                                     week_4_cut_tops < week_3_cut_tops ~ week_3_cut_tops,
                                     TRUE ~ week_4_cut_tops),
         week_3_new_emerg = week_3_emerg - week_3_cut_tops,
         week_4_new_emerg = week_4_emerg - week_4_cut_tops,
         emerg = week_2_emerg + week_4_new_emerg + week_4_soil_germ,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot))

# check for issues in the data
filter(evGermD2Dat, week_4_new_emerg < week_3_new_emerg)
filter(evGermD2Dat, nonemerg < 0)

# sample sizes
evGermD2Dat %>%
  group_by(age, plot, treatment) %>%
  count()
# 3-4 reps each


#### initial visualizations ####

ggplot(evGermD2Dat, aes(x = as.factor(plot), y = emerg/seeds_planted)) +
  stat_summary(geom = "point", fun = "mean", size = 2, aes(color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(color = treatment)) +
  facet_wrap(~ age)
# fungicide increased germination among adults and decreased it among seedlings

ggplot(evGermD2Dat, aes(x = as.factor(plot), y = emerg/seeds_planted, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0)
# no effect of fungicide when seeds are combined by plant age

ggplot(evGermD2Dat, aes(age, emerg/seeds_planted)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0)
# germination of seeds from seedlings higher than from adults

ggplot(evGermD2Dat, aes(treatment, emerg/seeds_planted)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0)


#### model ####

# initial fit
evGermD2Mod1 <- brm(emerg | trials(seeds_planted) ~ fungicide + (1|site/plotr),
                     data = evGermD2Dat, family = binomial,
                     prior = c(prior(normal(0, 10), class = Intercept),
                               prior(normal(0, 10), class = b),
                               prior(cauchy(0, 1), class = sd)),
                     iter = 6000, warmup = 1000, chains = 1)
# 13 divergent transitions
summary(evGermD2Mod1)

# increase chains and adapt delta
evGermD2Mod2 <- update(evGermD2Mod1, chains = 3,
                        control = list(adapt_delta = 0.9999, max_treedepth = 15))
summary(evGermD2Mod2)
plot(evGermD2Mod2)
pp_check(evGermD2Mod2, nsamples = 50)

# simulate fit
fitDat <- tibble(fungicide = c(0, 1), seeds_planted = c(30, 30)) %>%
  mutate(germination = fitted(evGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"]/seeds_planted,
         germination_lower = fitted(evGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"]/seeds_planted,
         germination_upper = fitted(evGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]/seeds_planted,
         treatment = ifelse(fungicide == 0, "water", "fungicide"))

# summarize by site
vizDat <- evGermD2Dat %>%
  group_by(site, treatment) %>%
  summarise(germination = mean(emerg/seeds_planted))

# fit figure
evGermD2Plot <- ggplot(fitDat, aes(treatment, germination, color = treatment)) +
  geom_errorbar(aes(ymin = germination_lower, ymax = germination_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = vizDat, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  theme_bw()


#### output ####
save(evGermD2Mod2, file = "output/elymus_germination_model_2019_density_exp.rda")
save(evGermD2Plot, file = "output/elymus_germination_figure_2019_density_exp.rda")
