##### info ####

# file: survival_analysis_2018_2019
# author: Amy Kendig
# date last edited: 10/4/20
# goal: estimate survival for each plant group


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
survL1Dat <- read_csv("intermediate-data/ev_processed_survival_2018_litter_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsV2D <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# check notes
unique(survD1Dat$field_notes)
unique(survL1Dat$field_notes)

# add plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# summer survival 2018
# remove background Microstegium (not individual plants)
# make survival 1 if the plant produced seeds in summer
# remove NA's 
sumSurvD1Dat <- survD1Dat %>%
  filter(month == "September" & !(sp == "Mv" & focal == 0)) %>%
  select(-month) %>%
  left_join(plotDens) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, TRUE ~ survival),
         plant_group = paste(sp, age, sep = " "),
         yearf = "year 1",
         experiment = "density",
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  filter(!is.na(survival))

# survival through winter given summer survival, 2018
# make summer survival 1 if the plant produced seeds in summer
# remove NA's 
# remove background Microstegium (not individual plants)
winSurvDat <- survD1Dat %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  mutate(September = case_when(seeds_produced == 1 ~ 1, TRUE ~ September)) %>%
  filter(September == 1 & !is.na(April) & !(sp == "Mv" & focal == 0)) %>%
  select(-September) %>%
  rename(survival = April)  %>%
  mutate(plant_group = paste(sp, age, sep = " "))

# 2019 survival data needs list of all plants (only replacements recorded)
# background plants
bgD2Dat <- plotsD %>% 
  select(plot, treatment, background, background_density) %>%
  filter(plot != 1) %>%
  mutate(start = 1,
         focal = 0) %>%
  group_by(plot, treatment, background, focal) %>%
  expand(start, background_density, ID = full_seq(start:background_density, 1)) %>%
  ungroup() %>%
  mutate(ID = as.character(ID),
         sp = substr(background, 1, 2),
         age = substr(background, 4, 11)) %>%
  select(-c(start, background, background_density)) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4"))

# all focal 2019 plants
# merge ID lists with plot
focD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(tibble(sp = c(rep("Mv", 3), rep("Ev", 4)),
                     age = c(rep("seedling", 6), "adult"),
                     ID = c(1, 2, 3, 1, 2, 3, "A"),
                     focal = 1))

# summer survival 2019
# add IDs to 2019 background plants
# use the latest replacement date
sumSurvD2Dat <- survD2Dat %>%
  group_by(site, plot, treatment, sp, age, focal) %>%
  mutate(Bg_ID = as.character(1:n()),
         ID = case_when(ID == "Bg" ~ Bg_ID,
                        TRUE ~ ID)) %>%
  ungroup() %>%
  select(-Bg_ID) %>%
  group_by(site, plot, treatment, sp, age, focal, ID) %>%
  summarise(replace_date = max(replace_date),
            survival = 0) %>%
  full_join(bgD2Dat) %>%
  full_join(focD2Dat) %>%
  left_join(plotDens) %>%
  mutate(survival = replace_na(survival, 1),
         plant_group = paste(sp, age, sep = " "),
         yearf = "year 2",
         experiment = "density",
         fungicide = ifelse(treatment == "fungicide", 1, 0))

# check that background ID's don't exceed density
sumSurvD2Dat %>%
  left_join(plotsD %>%
              select(plot, background_density)) %>%
  filter(focal == 0) %>%
  mutate(ID = as.numeric(ID)) %>%
  filter(ID > background_density)
# no

# check that number of plants is correct
nrow(sumSurvD2Dat)
2*4*(7*10 + 4 + 16 + 64 + 4 + 8 + 16 + 2 + 4 + 8)
# yes - matches 

# litter summer survival 2018
# make survival 1 if the plant produced seeds in summer
# remove NA's 
sumSurvL1Dat <- survL1Dat %>%
  filter(month == "September") %>%
  select(-month) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, TRUE ~ survival),
         plant_group = paste(sp, age, sep = " "),
         yearf = "year 1",
         experiment = "litter") %>%
  filter(!is.na(survival))


#### group and split data ####

# combine summer survival
sumSurvDat <- sumSurvD1Dat %>%
  select(experiment, yearf, site, plot, treatment, sp, age, plant_group, ID, focal, survival) %>%
  full_join(sumSurvD2Dat %>%
              select(experiment, yearf, site, plot, treatment, sp, age, plant_group, ID, focal, survival)) %>%
  full_join(sumSurvL1Dat %>%
              select(experiment, yearf, site, plot, sp, age, plant_group, ID, focal, survival)) %>%
  mutate(exp_yearf = paste(experiment, yearf, sep = "_"))

# divide data by plant group
evSSumSurvDat <- filter(sumSurvDat, plant_group == "Ev seedling")
evASumSurvDat <- filter(sumSurvDat, plant_group == "Ev adult")
mvSSumSurvDat <- filter(sumSurvDat, plant_group == "Mv seedling")

evSWinSurvDat <- filter(winSurvDat, plant_group == "Ev seedling")
evAWinSurvDat <- filter(winSurvDat, plant_group == "Ev adult")

# group data by site within experiment and year
grpSumSurvDat <-  sumSurvDat %>%
  group_by(plant_group, experiment, yearf, site) %>%
  summarise(surviving = sum(survival),
            plants = length(survival)) %>%
  ungroup()

grpWinSurvDat <-  winSurvDat %>%
  group_by(plant_group, site) %>%
  summarise(surviving = sum(survival),
            plants = length(survival)) %>%
  ungroup()

# divide data by plant group
evSGrpSumSurvDat <- filter(grpSumSurvDat, plant_group == "Ev seedling")
evAGrpSumSurvDat <- filter(grpSumSurvDat, plant_group == "Ev adult")
mvSGrpSumSurvDat <- filter(grpSumSurvDat, plant_group == "Mv seedling")

evSGrpWinSurvDat <- filter(grpWinSurvDat, plant_group == "Ev seedling")
evAGrpWinSurvDat <- filter(grpWinSurvDat, plant_group == "Ev adult")

# divide density data by plant group
evSSumSurvD1Dat <- filter(sumSurvD1Dat, plant_group == "Ev seedling" & focal == 1)
evASumSurvD1Dat <- filter(sumSurvD1Dat, plant_group == "Ev adult" & focal == 1)
mvSSumSurvD1Dat <- filter(sumSurvD1Dat, plant_group == "Mv seedling" & focal == 1)

evSSumSurvD2Dat <- filter(sumSurvD2Dat, plant_group == "Ev seedling" & focal == 1)
evASumSurvD2Dat <- filter(sumSurvD2Dat, plant_group == "Ev adult" & focal == 1)
mvSSumSurvD2Dat <- filter(sumSurvD2Dat, plant_group == "Mv seedling" & focal == 1)

# group density data by plot
grpSumSurvD1Dat <-sumSurvD1Dat %>%
  filter(focal == 1) %>%
  group_by(site, plot, mv_seedling_density, ev_seedling_density, ev_adult_density, treatment, fungicide, plant_group) %>%
  summarise(surviving = sum(survival),
            plants = length(survival)) %>%
  ungroup()

grpSumSurvD2Dat <-sumSurvD2Dat %>%
  filter(focal == 1) %>%
  group_by(site, plot, mv_seedling_density, ev_seedling_density, ev_adult_density, treatment, fungicide, plant_group) %>%
  summarise(surviving = sum(survival),
            plants = length(survival)) %>%
  ungroup()

# divide grouped density data by plant group
evSGrpSumSurvD1Dat <- filter(grpSumSurvD1Dat, plant_group == "Ev seedling")
evAGrpSumSurvD1Dat <- filter(grpSumSurvD1Dat, plant_group == "Ev adult")
mvSGrpSumSurvD1Dat <- filter(grpSumSurvD1Dat, plant_group == "Mv seedling")

evSGrpSumSurvD2Dat <- filter(grpSumSurvD2Dat, plant_group == "Ev seedling")
evAGrpSumSurvD2Dat <- filter(grpSumSurvD2Dat, plant_group == "Ev adult")
mvSGrpSumSurvD2Dat <- filter(grpSumSurvD2Dat, plant_group == "Mv seedling")


#### visualize ####

# mean summer survival by experiment and year
ggplot(sumSurvDat, aes(plant_group, survival, color = experiment, shape = yearf)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.2)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  theme_bw()

# mean winter survival
ggplot(winSurvDat, aes(plant_group, survival)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  theme_bw()

# density and fungicide effects from year 2
ggplot(sumSurvD2Dat, aes(background_density, survival)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.2)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  facet_grid(plant_group ~ background, scales = "free") +
  theme_bw()

ggplot(grpSumSurvD2Dat, aes(mv_seedling_density, surviving/plants, color = treatment)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~plant_group) +
  theme_bw()

ggplot(grpSumSurvD2Dat, aes(ev_seedling_density, surviving/plants, color = treatment)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~plant_group) +
  theme_bw()

ggplot(grpSumSurvD2Dat, aes(ev_adult_density, surviving/plants, color = treatment)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~plant_group) +
  theme_bw()


#### beta binomial regression settings ####

# source: https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

# define family
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]")

# provide Stan functions
stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"

# provide Stan variables
stanvars <- stanvar(scode = stan_funs, block = "functions")

# define log-likelihood
log_lik_beta_binomial2 <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

# define posterior prediction function
posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}


#### simulation data ####

simDat <- tibble(mv_seedling_density = c(seq(0, 64, length.out = 100), rep(0, 200)),
                 ev_seedling_density = c(rep(0, 100), seq(0, 16, length.out = 100), rep(0, 100)),
                 ev_adult_density = c(rep(0, 200), seq(0, 8, length.out = 100)),
                 background = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"),
         site = NA,
         background_density = case_when(background == "Mv seedling" ~ mv_seedling_density,
                                        background == "Ev seedling" ~ ev_seedling_density,
                                        TRUE ~ ev_adult_density))


#### Mv summer survival: all data ####

# initial model
mvSSumSurvMod1 <- brm(survival ~ 1 + (1|yearf) + (1|site/plot),
                      data = mvSSumSurvDat,
                      family = "bernoulli",
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(cauchy(0, 1), class = sd)),
                      iter = 6000, warmup = 1000, chains = 1)
# 118 divergent transitions
summary(mvSSumSurvMod1)
plot(mvSSumSurvMod1)

# look for issues that may be causing divergent transitions
group_by(mvSSumSurvDat, survival) %>% count()
# few zeros (will have same issue with Ev adults)

# model with grouped data
mvSGrpSumSurvMod1 <- brm(surviving | trials(plants) ~ 1 + (1|site),
                      data = mvSGrpSumSurvDat,
                      family = "binomial",
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(cauchy(0, 1), class = sd)),
                      iter = 6000, warmup = 1000, chains = 1)
# 6 divergent transitions
summary(mvSGrpSumSurvMod1)
plot(mvSGrpSumSurvMod1)

# increase adapt delta
mvSGrpSumSurvMod2 <- update(mvSGrpSumSurvMod1,
                            control = list(adapt_delta = 0.99))
# removed divergent transitions

# update number of chains
# need to increase adapt_delta to remove divergent transitions
mvSGrpSumSurvMod3 <- update(mvSGrpSumSurvMod2,
                            chains = 3,
                            control = list(adapt_delta = 0.999))
summary(mvSGrpSumSurvMod3)
plot(mvSGrpSumSurvMod3)

# try beta binomial function
mvSGrpSumSurvMod4 <- brm(surviving | vint(plants) ~ 1 + (1|site),
                         data = mvSGrpSumSurvDat,
                         family = beta_binomial2, 
                         stanvars = stanvars,
                         prior = c(prior(normal(0, 10), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                         iter = 6000, warmup = 1000, chains = 1)
# 1 divergent transition
summary(mvSGrpSumSurvMod4)

# increase chains and adapt delta
mvSGrpSumSurvMod5 <- update(mvSGrpSumSurvMod4,
                            chains = 3,
                            control = list(adapt_delta = 0.999))
summary(mvSGrpSumSurvMod5)
plot(mvSGrpSumSurvMod5)
expose_functions(mvSGrpSumSurvMod5, vectorize = TRUE)

# compare fits
mvSGrpSumSurvLoo <- loo(mvSGrpSumSurvMod3, mvSGrpSumSurvMod5, reloo = T)
mvSGrpSumSurvLoo
# beta-binomial provides a better fit

# examine fit
pp_check(mvSGrpSumSurvMod3, nsamples = 50)
pp_check(mvSGrpSumSurvMod5, nsamples = 50)


#### Mv summer 2 survival: density and fungicide effects ####

# beta-binomial
mvSGrpSumSurvD2Mod1 <- brm(surviving | vint(plants) ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site),
                         data = mvSGrpSumSurvD2Dat,
                         family = beta_binomial2, 
                         stanvars = stanvars,
                         prior = c(prior(normal(0, 10), class = Intercept),
                                   prior(normal(0, 10), class = b),
                                   prior(cauchy(0, 1), class = sd)),
                         iter = 6000, warmup = 1000, chains = 1)
# 3537 divergent transitions - the plants number may be too low
plot(mvSGrpSumSurvD2Mod1)
# big phi and sd values
summary(mvSGrpSumSurvD2Mod1)

# binomial
mvSSumSurvD2Mod1 <- brm(survival ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site/plot),
                        data = mvSSumSurvD2Dat,
                        family = bernoulli,
                        prior = c(prior(normal(0, 10), class = Intercept),
                                     prior(normal(0, 10), class = b),
                                     prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 1)
summary(mvSSumSurvD2Mod1)
plot(mvSSumSurvD2Mod1)

# increase chains
mvSSumSurvD2Mod2 <- update(mvSSumSurvD2Mod1, chains = 3)
summary(mvSSumSurvD2Mod2)
plot(mvSSumSurvD2Mod2)
# the interactions between fungicide and Ev seedling and Ev adult density is different from zero
pp_check(mvSSumSurvD2Mod2, nsamples = 50)

# simulate fit
mvSSumSurvD2Fit <- simDat %>%
  mutate(survival = fitted(mvSSumSurvD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(mvSSumSurvD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(mvSSumSurvD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# modify dataset so that zero plots are repeated
mvSGrpSumSurvD2Dat2 <- mvSGrpSumSurvD2Dat %>%
  full_join(plotsV2D)

# fit figure
ggplot(mvSSumSurvD2Fit, aes(background_density, survival, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  geom_point(data = mvSGrpSumSurvD2Dat2, aes(y = surviving/plants)) +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()
# survival is on average lower with fungicide
# this goes away with increasing density because survival maxes out at high density

# simplify model: remove interactions
mvSSumSurvD2Mod3 <- update(mvSSumSurvD2Mod2, formula = survival ~ fungicide + mv_seedling_density + ev_seedling_density + ev_adult_density + (1|site/plot),
                           control = list(adapt_delta = 0.999))

# compare model with and without interactions
loo(mvSSumSurvD2Mod2, mvSSumSurvD2Mod3)
# a lot of observations need to taken out to do LOO
mvSSumSurvD2Kfold <- kfold(mvSSumSurvD2Mod2, mvSSumSurvD2Mod3)
# run time: ~ 80 minutes (20 refits at 2-3 minutes each)
# full model is a better fit
waic(mvSSumSurvD2Mod2, mvSSumSurvD2Mod3)
# recommends loo instead

#### Mv summer 1 survival: density and fungicide effects ####

# binomial
mvSSumSurvD1Mod1 <- brm(survival ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site/plot),
                        data = mvSSumSurvD1Dat,
                        family = bernoulli,
                        prior = c(prior(normal(0, 10), class = Intercept),
                                  prior(normal(0, 10), class = b),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 1)
summary(mvSSumSurvD1Mod1)

# increase chains
mvSSumSurvD1Mod2 <- update(mvSSumSurvD1Mod1, chains = 3)
summary(mvSSumSurvD1Mod2)
plot(mvSSumSurvD1Mod2)
pp_check(mvSSumSurvD1Mod2, nsamples = 50)

# simulate fit
mvSSumSurvD1Fit <- simDat %>%
  mutate(survival = fitted(mvSSumSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(mvSSumSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(mvSSumSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# modify dataset so that zero plots are repeated
mvSGrpSumSurvD1Dat2 <- mvSGrpSumSurvD1Dat %>%
  full_join(plotsV2D)

# fit figure
ggplot(mvSSumSurvD1Fit, aes(background_density, survival, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  geom_point(data = mvSGrpSumSurvD1Dat2, aes(y = surviving/plants)) +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()
# survival is lower with fungicide and decreases slightly with density
# large error regions


#### Mv summer 2 survival: density and fungicide effects, with priors ####

# prior model
mvSSumSurvD1Mod3 <- brm(survival ~ 1,
                        data = mvSSumSurvD1Dat,
                        family = bernoulli,
                        prior = prior(normal(0, 10), class = Intercept),
                        iter = 6000, warmup = 1000, chains = 3)
summary(mvSSumSurvD1Mod3)

# binomial
mvSSumSurvD2Mod4 <- brm(survival ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site/plot),
                        data = mvSSumSurvD2Dat,
                        family = bernoulli,
                        prior = c(prior(normal(2.23, 0.22), class = Intercept),
                                  prior(normal(0, 10), class = b),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 1)
summary(mvSSumSurvD2Mod4)

# increase chains
mvSSumSurvD2Mod5 <- update(mvSSumSurvD2Mod4, chains = 3,
                           control = list(adapt_delta = 0.99))
summary(mvSSumSurvD2Mod5)
plot(mvSSumSurvD2Mod5)
pp_check(mvSSumSurvD2Mod2, nsamples = 50)

# simulate fit
mvSSumSurvD2Fit5 <- simDat %>%
  mutate(survival = fitted(mvSSumSurvD2Mod5, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(mvSSumSurvD2Mod5, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(mvSSumSurvD2Mod5, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# modify dataset so that zero plots are repeated
mvSGrpSumSurvD2Dat2 <- mvSGrpSumSurvD2Dat %>%
  full_join(plotsV2D)

# fit figure
ggplot(mvSSumSurvD2Fit5, aes(background_density, survival, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  geom_point(data = mvSGrpSumSurvD2Dat2, aes(y = surviving/plants)) +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()
# super low estimates

# allow sd to be larger
mvSSumSurvD2Mod6 <- update(mvSSumSurvD2Mod5,
                           prior = c(prior(normal(2.23, 10), class = Intercept),
                                     prior(normal(0, 10), class = b),
                                     prior(cauchy(0, 1), class = sd)))

summary(mvSSumSurvD2Mod6)

# simulate fit
mvSSumSurvD2Fit6 <- simDat %>%
  mutate(survival = fitted(mvSSumSurvD2Mod6, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(mvSSumSurvD2Mod6, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(mvSSumSurvD2Mod6, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# fit figure
ggplot(mvSSumSurvD2Fit6, aes(background_density, survival, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  geom_point(data = mvSGrpSumSurvD2Dat2, aes(y = surviving/plants)) +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()
# identical to one with general priors
  
  
#### output ####
save(mvSGrpSumSurvMod3, file = "output/mv_seedling_survival_all_data_binomial_model.rda")
save(mvSGrpSumSurvMod5, file = "output/mv_seedling_survival_all_data_beta_binomial_model.rda")
save(mvSGrpSumSurvLoo, file = "output/mv_seedling_survival_all_data_model_comparison.rda")

save(mvSSumSurvD1Mod2, file = "output/mv_seedling_survival_2018_density_exp_model.rda")
save(mvSSumSurvD2Mod2, file = "output/mv_seedling_survival_2019_density_exp_model.rda")
save(mvSSumSurvD2Kfold, file = "output/mv_seedling_survival_2019_density_exp_model_comparison.rda")
