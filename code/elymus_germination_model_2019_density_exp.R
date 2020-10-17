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

# seedling data
evSGermD2Dat <- evGermD2Dat %>%
  filter(age == "seedling")

# adult data
evAGermD2Dat <- evGermD2Dat %>%
  filter(age == "adult")


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
evGermD2Mod1 <- brm(emerg | trials(seeds_planted) ~ fungicide + (1|site),
                     data = evGermD2Dat, family = binomial,
                     prior = c(prior(normal(0, 10), class = Intercept),
                               prior(normal(0, 10), class = b),
                               prior(cauchy(0, 1), class = sd)),
                     iter = 6000, warmup = 1000, chains = 1)
# 7 divergent transitions
summary(evGermD2Mod1)

# increase chains and adapt delta
evGermD2Mod2 <- update(evGermD2Mod1, chains = 3,
                        control = list(adapt_delta = 0.999))
summary(evGermD2Mod2)
plot(evGermD2Mod2)
pp_check(evGermD2Mod2, nsamples = 50)

# simulate fit
fitDat <- tibble(fungicide = c(0, 1), seeds_planted = c(30, 30)) %>%
  mutate(emerg = fitted(evGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         emerg_lower = fitted(evGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         emerg_upper = fitted(evGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
         treatment = ifelse(fungicide == 0, "water", "fungicide"))

# fit figure
ggplot(evGermD2Dat, aes(treatment, emerg/seeds_planted)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  geom_point(data = fitDat, size = 2, color = "red", alpha = 0.5) +
  geom_errorbar(data = fitDat, aes(ymin = emerg_lower/seeds_planted, ymax = emerg_upper/seeds_planted), color = "red", width = 0, alpha = 0.5) +
  theme_bw()

# did separate models by parent age (below), but seeds are grouped together in the model after they're produced, so the age of the parent is not tracked


#### seedling model ####

# initial fit
evSGermD2Mod1 <- brm(emerg | trials(seeds_planted) ~ fungicide + (1|site),
                    data = evSGermD2Dat, family = binomial,
                    prior = c(prior(normal(0, 10), class = Intercept),
                              prior(normal(0, 10), class = b),
                              prior(cauchy(0, 1), class = sd)),
                    iter = 6000, warmup = 1000, chains = 1)
# 7 divergent transitions
summary(evSGermD2Mod1)

# increase chains and adapt delta
evSGermD2Mod2 <- update(evSGermD2Mod1, chains = 3,
                        control = list(adapt_delta = 0.999))
summary(evSGermD2Mod2)
plot(evSGermD2Mod2)
pp_check(evSGermD2Mod2, nsamples = 50)

# simulate fit
fitDat1 <- tibble(fungicide = c(0, 1), seeds_planted = c(30, 30)) %>%
  mutate(emerg = fitted(evSGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         emerg_lower = fitted(evSGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         emerg_upper = fitted(evSGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
         treatment = ifelse(fungicide == 0, "water", "fungicide"))

# fit figure
ggplot(evSGermD2Dat, aes(treatment, emerg/seeds_planted)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  geom_point(data = fitDat1, size = 2, color = "red", alpha = 0.5) +
  geom_errorbar(data = fitDat1, aes(ymin = emerg_lower/seeds_planted, ymax = emerg_upper/seeds_planted), color = "red", width = 0, alpha = 0.5) +
  theme_bw()


#### adult model ####

# initial fit
evAGermD2Mod1 <- brm(emerg | trials(seeds_planted) ~ fungicide + (1|site),
                     data = evAGermD2Dat, family = binomial,
                     prior = c(prior(normal(0, 10), class = Intercept),
                               prior(normal(0, 10), class = b),
                               prior(cauchy(0, 1), class = sd)),
                     iter = 6000, warmup = 1000, chains = 1)
# 1 divergent transition
summary(evAGermD2Mod1)

# increase chains and adapt delta
evAGermD2Mod2 <- update(evAGermD2Mod1, chains = 3,
                        control = list(adapt_delta = 0.999))
summary(evAGermD2Mod2)
plot(evAGermD2Mod2)
pp_check(evAGermD2Mod2, nsamples = 50)

# simulate fit
fitDat2 <- tibble(fungicide = c(0, 1), seeds_planted = c(30, 30)) %>%
  mutate(emerg = fitted(evAGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         emerg_lower = fitted(evAGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         emerg_upper = fitted(evAGermD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
         treatment = ifelse(fungicide == 0, "water", "fungicide"))

# fit figure
ggplot(evAGermD2Dat, aes(treatment, emerg/seeds_planted)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  geom_point(data = fitDat2, size = 2, color = "red", alpha = 0.5) +
  geom_errorbar(data = fitDat2, aes(ymin = emerg_lower/seeds_planted, ymax = emerg_upper/seeds_planted), color = "red", width = 0, alpha = 0.5) +
  theme_bw()


#### output ####
save(evGermD2Mod2, file = "output/elymus_germination_model_2019_density_exp.rda")
save(evSGermD2Mod2, file = "output/elymus_seedling_germination_model_2019_density_exp.rda")
save(evAGermD2Mod2, file = "output/elymus_adult_germination_model_2019_density_exp.rda")
