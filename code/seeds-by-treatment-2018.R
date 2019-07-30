##### info ####

# file: seeds-by-treatment-2018
# author: Amy Kendig
# date last edited: 7/29/19
# goal: see how Mv and Ev seed production is affected by treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)

# run files
source("./code/bg-densities-data-processing-2018.R")
rm(list = setdiff(ls(), c("bgd", "esurv")))
source("./code/ev-seeds-data-processing-2018.R")
rm(list = setdiff(ls(), c("bgd", "esurv", "eseeds")))
source("./code/mv-seeds-data-processing-2018.R")
rm(list = setdiff(ls(), c("bgd", "esurv", "eseeds", "mseeds")))
source("./code/covariate-data-processing-2018.R")
rm(list = setdiff(ls(), c("bgd", "esurv", "eseeds", "mseeds", "covar")))

# import other data files
trt <- read_csv("./data/plot-treatments-for-figures-2018-density-exp.csv")
trt_s <- read_csv("./data/plot-treatments-2018-density-exp.csv")
til <- read_csv("./data/focal-size-disease-jul-2018-density-exp.csv")
mbseeds <- read_csv("./data/mv-biomass-oct-2018-density-exp.csv")

# non-linear species effects function
el_fun <- function(y0, x, a) {
  y = y0 * exp(a * log(x + 1))
  return(y)
}

# correlation between predicted and observed
cor_fun <- function(dat, mt, mn, sp){
  dat2 <- dat %>%
    mutate(transformed = predict(mt)[, 1],
           nonlinear = predict(mn)[, 1])
  if(sp == "ev"){
    print(cor.test(dat2$transformed, dat2$log_seeds))
    print(cor.test(dat2$nonlinear, dat2$seeds))
  }
  if(sp == "mv"){
    print(cor.test(dat2$transformed, dat2$log_seeds_pb))
    print(cor.test(dat2$nonlinear, dat2$seeds_per_bio))
  }
  
}

# function to create datasets for prediction
d_pred_fun = function(dat, mt, mn, max_val){
  
  d_pred <- tibble(
    background_density = rep(seq(0, max_val, length.out = 100), 2),
    counted_density = rep(seq(0, max_val, length.out = 100), 2),
    treatment = rep(c("water", "fungicide"), each = 100),
    sm_adj = mean(dat$sm_adj),
    cc_adj = mean(dat$cc_adj),
    pm_adj = mean(dat$pm_adj))
  
  d_pred2 <- d_pred %>%
    cbind(fitted(mt, newdata = d_pred, re_formula = NA, nsamples = 100)) %>%
    mutate(model = "linear") %>%
    full_join(d_pred %>%
                cbind(fitted(mn, newdata = d_pred, re_formula = NA, nsamples = 100)) %>%
                mutate(model = "nonlinear")) %>%
    rename(pred = Estimate,
           lower = Q2.5,
           upper = Q97.5) %>%
    as_tibble()
  
  return(d_pred2)
}

# plot parameters
colpal = c("#3CBB75FF", "#39568CFF")
sm_txt = 12
lg_txt = 14
an_txt = 3


#### edit data ####

# Mv tiller number
til2 <- til %>%
  filter(sp == "Mv" & !(site == "D4" & plot == 8 & treatment == "fungicide") & !(site == "D4" & plot == 10  & treatment == "fungicide")) %>%
  group_by(site, plot, treatment) %>%
  summarise(mean_till = mean(tillers, na.rm = T)) %>%
  ungroup()

# combine Ev across dates (remove ones with uncertain ID's)
eseeds2 <- eseeds %>%
  filter(ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, age, ID, focal) %>%
  summarise(
    seeds = sum(seeds, na.rm = T)
  ) %>%
  ungroup()

# edit soil/biomass labels and seed values 
mbseeds2 <- mbseeds %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment),
         seeds = seeds_bio + seeds_soil,
         seeds_per_bio = seeds / bio.g,
         log_seeds_pb = log(seeds_per_bio))

# September Ev survival data, leave in litter plots
esurv2 <- esurv %>%
  filter(month == "September" & !is.na(survival)) %>%
  select(-month)

# remove month
bgd2 <- bgd %>% select(-month)

# merge seed data with plot and other data
de <-  full_join(eseeds2, esurv2) %>% 
  filter(site %in% c("D1", "D2", "D3", "D4")) %>% 
  full_join(trt) %>% 
  full_join(bgd2) %>%
  left_join(covar) %>%
  mutate(density_level = factor(density_level, levels = c("none", "low", "medium", "high")))

le <- full_join(eseeds2, esurv2) %>% 
  filter(site %in% c("L1", "L2", "L3", "L4"))

dm <- full_join(mseeds, til2) %>% 
  full_join(trt)  %>% 
  full_join(bgd2) %>%
  left_join(covar) %>%
  mutate(density_level = factor(density_level, levels = c("none", "low", "medium", "high")))

dmb <- full_join(mbseeds2, trt)  %>% 
  full_join(bgd2) %>%
  left_join(covar) %>%
  mutate(density_level = factor(density_level, levels = c("none", "low", "medium", "high"))) 

# see if all have tiller counts
sum(is.na(dm$mean_till)) # yes

# estimate per capita seeds for Mv, assuming seeds were collected from 3 tillers
dm2 <- dm %>%
  mutate(
    seeds_percap = seeds / 3 * mean_till,
    log_seeds_percap = log(seeds_percap)
  )

# check for missing survival info
sum(is.na(de$survival)) # none missing (those with NA's had no seeds)
sum(is.na(le$survival))

# give Ev's 0's if no seeds were collected and they're still alive
filter(de, survival == 1 & is.na(seeds)) # 136 plants
filter(le, survival == 1 & is.na(seeds)) # 1 plant

de2 <- de %>%
  mutate(
    seeds = case_when(survival == 1 & is.na(seeds) ~ 0,
                      TRUE ~ seeds),
    log_seeds = log(seeds + 1)
  )

le2 <- le %>%
  mutate(
    seeds = case_when(survival == 1 & is.na(seeds) ~ 0,
                      TRUE ~ seeds),
    log_seeds = log(seeds + 1)
  )

# double check survival and seed relationships
filter(de2, survival_seeds == 1 & seeds == 0)
filter(le2, survival_seeds == 1 & seeds == 0) # looks good
filter(de2, survival == 0) %>% select(seeds) %>% unique()
filter(le2, survival == 0) %>% select(seeds) %>% unique() # looks good

# make sure all the data are there
de2 %>%
  filter(focal == 1) %>%
  group_by(site, plot, treatment, age) %>%
  summarise(n = n()) %>%
  filter((plot == 1 & age == "adult" & n != 3) |
           (plot == 1 & age == "seedling" & n != 9) |
           (plot != 1 & age == "adult" & n != 1) |
           (plot != 1 & age == "seedling" & n != 3))

de2 %>% select(site, plot, treatment) %>% unique() # 78

le2 %>% select(site, plot) %>% unique() # 27

dm2 %>% select(site, plot, treatment) %>% unique() # 78

dmb %>% select(site, plot, treatment) %>% unique() # 78


##### visualize #####

##  treatment effects

# Ev adult
de2 %>%
  filter(age == "adult" & !is.na(seeds) & survival == 1) %>%
  ggplot(aes(x = density_level, y = seeds, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) +
  ylab("Ev adult seeds") # not a big difference between fungicide and water, some benefits of Ev seedling background

# Ev adult, focal only
de2 %>%
  filter(age == "adult" & !is.na(seeds) & survival == 1 & focal == 1) %>%
  ggplot(aes(x = density_level, y = seeds, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) +
  ylab("Ev adult seeds") # similar to above

# Ev adult, focal only, log-transformed
de2 %>%
  filter(age == "adult" & !is.na(seeds) & survival == 1 & focal == 1) %>%
  ggplot(aes(x = density_level, y = log_seeds, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) +
  ylab("Ev adult seeds")

# Ev seedling
de2 %>%
  filter(age == "seedling" & !is.na(seeds) & survival == 1) %>%
  ggplot(aes(x = density_level, y = seeds, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) +
  ylab("Ev seedling seeds") # negative effect of fungicide with adults, similar to survival results

# Ev seedling, focal only
de2 %>%
  filter(age == "seedling" & !is.na(seeds) & survival == 1 & focal == 1) %>%
  ggplot(aes(x = density_level, y = seeds, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) +
  ylab("Ev seedling seeds") # similar to above

# Ev seedling, focal only, log-transformed
de2 %>%
  filter(age == "seedling" & !is.na(seeds) & survival == 1 & focal == 1) %>%
  ggplot(aes(x = density_level, y = log_seeds, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) +
  ylab("Ev seedling seeds") # similar to above

# Mv per capita
dm2 %>%
  filter(!is.na(seeds_percap)) %>%
  ggplot(aes(x = density_level, y = seeds_percap, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) # no clear fungicide pattern, lower with higher Mv under water treatment, different from biomass results

# Mv per biomass
dmb %>%
  filter(!is.na(seeds_per_bio)) %>%
  ggplot(aes(x = density_level, y = seeds_per_bio, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) # lower with water, especially with seedling backgrounds, lower with higher Mv density, different from biomass results

# Mv per biomass, log-transformed
dmb %>%
  filter(!is.na(seeds_per_bio)) %>%
  ggplot(aes(x = density_level, y = log_seeds_pb, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background) 


#### set up models ####

## trim datasets by seed and background
d_sa_ba <- filter(de2, age == "adult" & !is.na(seeds) & survival == 1 & focal == 1 & background == "Ev adult")
d_sa_bs <- filter(de2, age == "adult" & !is.na(seeds) & survival == 1 & focal == 1 & background == "Ev seedling")
d_sa_bm <- filter(de2, age == "adult" & !is.na(seeds) & survival == 1 & focal == 1 & background == "Mv seedling")

d_ss_ba <- filter(de2, age == "seedling" & !is.na(seeds) & survival == 1 & focal == 1 & background == "Ev adult")
d_ss_bs <- filter(de2, age == "seedling" & !is.na(seeds) & survival == 1 & focal == 1 & background == "Ev seedling")
d_ss_bm <- filter(de2, age == "seedling" & !is.na(seeds) & survival == 1 & focal == 1 & background == "Mv seedling")

d_sm_ba <-  filter(dmb, !is.na(seeds_per_bio) & background == "Ev adult")
d_sm_bs <- filter(dmb, !is.na(seeds_per_bio) & background == "Ev seedling")
d_sm_bm <- filter(dmb, !is.na(seeds_per_bio) & background == "Mv seedling")

## full linear log-transformed models

# Adult seeds, each background
# mt_sa_ba <- brm(data = d_sa_ba, family = gaussian,
#             log_seeds ~ background_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
#             prior <- c(prior(normal(0, 100), class = Intercept),
#                        prior(normal(0, 10), class = b),
#                        prior(cauchy(0, 1), class = sd)),
#             iter = 6000, warmup = 1000, chains = 3, cores = 2,
#             control = list(adapt_delta = 0.9999, max_treedepth = 15))
# summary(mt_sa_ba)
# save(mt_sa_ba, file = "./output/ev-adult-seeds-by-treatment-2018-log-transformed-ev-adult.rda")
load("./output/ev-adult-seeds-by-treatment-2018-log-transformed-ev-adult.rda")
 
# mt_sa_bs <- update(mt_sa_ba, 
#                    formula = log_seeds ~ counted_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
#                    newdata = d_sa_bs)
# summary(mt_sa_bs)
# save(mt_sa_bs, file = "./output/ev-adult-seeds-by-treatment-2018-log-transformed-ev-seedling.rda")
load("./output/ev-adult-seeds-by-treatment-2018-log-transformed-ev-seedling.rda")

# mt_sa_bm <- update(mt_sa_ba, newdata = d_sa_bm,
#                    control = list(adapt_delta = 0.99999))
# summary(mt_sa_bm)
# save(mt_sa_bm, file = "./output/ev-adult-seeds-by-treatment-2018-log-transformed-mv-seedling.rda")
load("./output/ev-adult-seeds-by-treatment-2018-log-transformed-mv-seedling.rda")

# Ev seedling seeds, each background
# mt_ss_ba <- update(mt_sa_ba, newdata = d_ss_ba)
# summary(mt_ss_ba)
# save(mt_ss_ba, file = "./output/ev-seedling-seeds-by-treatment-2018-log-transformed-ev-adult.rda")
load("./output/ev-seedling-seeds-by-treatment-2018-log-transformed-ev-adult.rda")

# mt_ss_bs <- update(mt_sa_bs, newdata = d_ss_bs)
# summary(mt_ss_bs)
# save(mt_ss_bs, file = "./output/ev-seedling-seeds-by-treatment-2018-log-transformed-ev-seedling.rda")
load("./output/ev-seedling-seeds-by-treatment-2018-log-transformed-ev-seedling.rda")

# mt_ss_bm <- update(mt_sa_ba, newdata = d_ss_bm,
#                    control = list(adapt_delta = 0.999999999))
# summary(mt_ss_bm)
# save(mt_ss_bm, file = "./output/ev-seedling-seeds-by-treatment-2018-log-transformed-mv-seedling.rda")
load("./output/ev-seedling-seeds-by-treatment-2018-log-transformed-mv-seedling.rda")

# Mv seedling seeds, each background
# mt_sm_ba <- update(mt_sa_ba, 
#                    formula = log_seeds_pb ~ background_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
#                    newdata = d_sm_ba)
# summary(mt_sm_ba)
# save(mt_sm_ba, file = "./output/mv-seedling-seeds-by-treatment-2018-log-transformed-ev-adult.rda")
load("./output/mv-seedling-seeds-by-treatment-2018-log-transformed-ev-adult.rda")

# mt_sm_bs <- update(mt_sa_bs, 
#                    formula = log_seeds_pb ~ counted_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
#                    newdata = d_sm_bs,
#                    control = list(adapt_delta = 0.99999))
# summary(mt_sm_bs)
# save(mt_sm_bs, file = "./output/mv-seedling-seeds-by-treatment-2018-log-transformed-ev-seedling.rda")
load("./output/mv-seedling-seeds-by-treatment-2018-log-transformed-ev-seedling.rda")

# mt_sm_bm <- update(mt_sa_bm, 
#                    formula = log_seeds_pb ~ background_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
#                    newdata = d_sm_bm)
# summary(mt_sm_bm)
# save(mt_sm_bm, file = "./output/mv-seedling-seeds-by-treatment-2018-log-transformed-mv-seedling.rda")
load("./output/mv-seedling-seeds-by-treatment-2018-log-transformed-mv-seedling.rda")

## full exponential untransformed model

# y0 is the value of seeds when density is zero
# just use one background dataset so that none plots aren't counted multiple times
d_sa_ba %>%
  filter(background_density == 0) %>%
  summarise(mean_seeds = mean(seeds), 
            se_seeds = sd(seeds)/sqrt(length(seeds)),
            min_seeds = min(seeds), 
            max_seeds = max(seeds))
# mean = 44.9, sd = 10.9
y0a = 44.9

d_ss_ba %>%
  filter(background_density == 0) %>%
  summarise(mean_seeds = mean(seeds), 
            se_seeds = sd(seeds)/sqrt(length(seeds)),
            min_seeds = min(seeds), 
            max_seeds = max(seeds))
# mean = 2.35, sd = 1.26
y0s = 2.35

d_sm_ba %>%
  filter(background_density == 0) %>%
  summarise(mean_seeds = mean(seeds_per_bio), 
            se_seeds = sd(seeds_per_bio)/sqrt(length(seeds_per_bio)),
            min_seeds = min(seeds_per_bio), 
            max_seeds = max(seeds_per_bio))
# mean = 60, sd = 16.6
y0m = 60

# a is the rate of increase or decrease - simulate data
sim_dat = expand.grid(0:64, c(seq(-2, 0, by = 0.5), seq(0.01, 0.5, by = 0.05)))
colnames(sim_dat) = c("density", "a")
sim_dat <- sim_dat %>%
  mutate(seeds = el_fun(y0m, density, a),
         pos = ifelse(a >= 0, 1, 0))

# figure
sim_dat %>%
  ggplot(aes(x = density, y = seeds, colour = as.factor(a))) + 
  geom_line() +
  facet_wrap(~pos, scales = "free")
# mean = 0, sd = 0.5 

# Adult seeds, each background
# mn_sa_ba <- brm(data = d_sa_ba, family = gaussian,
#                 bf(seeds ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T),
#                 prior <- c(prior(normal(45, 11), nlpar = "y0"),
#                            prior(normal(0, 0.5), nlpar = "a")),
#                 iter = 6000, warmup = 1000, chains = 3, cores = 2,
#                 control = list(adapt_delta = 0.999))
# summary(mn_sa_ba)
# save(mn_sa_ba, file = "./output/ev-adult-seeds-by-treatment-2018-nonlinear-ev-adult.rda")
load("./output/ev-adult-seeds-by-treatment-2018-nonlinear-ev-adult.rda")

# mn_sa_bs <- update(mn_sa_ba, newdata = d_sa_bs,
#                    bf(seeds ~ y0 * exp(a * log(counted_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T))
# summary(mn_sa_bs)
# save(mn_sa_bs, file = "./output/ev-adult-seeds-by-treatment-2018-nonlinear-ev-seedling.rda")
load("./output/ev-adult-seeds-by-treatment-2018-nonlinear-ev-seedling.rda")

# mn_sa_bm <- update(mn_sa_ba, newdata = d_sa_bm)
# summary(mn_sa_bm)
# save(mn_sa_bm, file = "./output/ev-adult-seeds-by-treatment-2018-nonlinear-mv-seedling.rda")
load("./output/ev-adult-seeds-by-treatment-2018-nonlinear-mv-seedling.rda")

# Ev seedling seeds, each background
# mn_ss_ba <- brm(data = d_ss_ba, family = gaussian,
#                 bf(seeds ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T),
#                 prior <- c(prior(normal(2, 1), nlpar = "y0"),
#                            prior(normal(0, 0.5), nlpar = "a")),
#                 iter = 6000, warmup = 1000, chains = 3, cores = 2,
#                 control = list(adapt_delta = 0.999))
# summary(mn_ss_ba)
# save(mn_ss_ba, file = "./output/ev-seedling-seeds-by-treatment-2018-nonlinear-ev-adult.rda")
load("./output/ev-seedling-seeds-by-treatment-2018-nonlinear-ev-adult.rda")

# mn_ss_bs <- update(mn_ss_ba, newdata = d_ss_bs,
#                    bf(seeds ~ y0 * exp(a * log(counted_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T))
# summary(mn_ss_bs)
# save(mn_ss_bs, file = "./output/ev-seedling-seeds-by-treatment-2018-nonlinear-ev-seedling.rda")
load("./output/ev-seedling-seeds-by-treatment-2018-nonlinear-ev-seedling.rda")

# mn_ss_bm <- update(mn_ss_ba, newdata = d_ss_bm)
# summary(mn_ss_bm)
# save(mn_ss_bm, file = "./output/ev-seedling-seeds-by-treatment-2018-nonlinear-mv-seedling.rda")
load("./output/ev-seedling-seeds-by-treatment-2018-nonlinear-mv-seedling.rda")

# Mv seedling seeds, each background
# mn_sm_ba <- brm(data = d_sm_ba, family = gaussian,
#                 bf(seeds_per_bio ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T),
#                 prior <- c(prior(normal(60, 17), nlpar = "y0"),
#                            prior(normal(0, 0.5), nlpar = "a")),
#                 iter = 6000, warmup = 1000, chains = 3, cores = 2,
#                 control = list(adapt_delta = 0.999))
# summary(mn_sm_ba)
# save(mn_sm_ba, file = "./output/mv-seedling-seeds-by-treatment-2018-nonlinear-ev-adult.rda")
load("./output/mv-seedling-seeds-by-treatment-2018-nonlinear-ev-adult.rda")

# mn_sm_bs <- update(mn_sm_ba, newdata = d_sm_bs,
#                    bf(seeds_per_bio ~ y0 * exp(a * log(counted_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T))
# summary(mn_sm_bs)
# save(mn_sm_bs, file = "./output/mv-seedling-seeds-by-treatment-2018-nonlinear-ev-seedling.rda")
load("./output/mv-seedling-seeds-by-treatment-2018-nonlinear-ev-seedling.rda")

# mn_sm_bm <- update(mn_sm_ba, newdata = d_sm_bm)
# summary(mn_sm_bm)
# save(mn_sm_bm, file = "./output/mv-seedling-seeds-by-treatment-2018-nonlinear-mv-seedling.rda")
load("./output/mv-seedling-seeds-by-treatment-2018-nonlinear-mv-seedling.rda")


#### check models ####

# convergence of chains
plot(mt_sa_ba)
plot(mt_sa_bs)
plot(mt_sa_bm)

plot(mt_ss_ba)
plot(mt_ss_bs)
plot(mt_ss_bm)

plot(mt_sm_ba)
plot(mt_sm_bs)
plot(mt_sm_bm)

plot(mn_sa_ba)
plot(mn_sa_bs)
plot(mn_sa_bm)

plot(mn_ss_ba)
plot(mn_ss_bs)
plot(mn_ss_bm)

plot(mn_sm_ba)
plot(mn_sm_bs)
plot(mn_sm_bm)

# compare pp_check between model types
pp_check(mt_sa_ba, nsamples = 50)
pp_check(mn_sa_ba, nsamples = 50) # t looks slightly better
pp_check(mt_sa_bs, nsamples = 50)
pp_check(mn_sa_bs, nsamples = 50) # n looks slightly better
pp_check(mt_sa_bm, nsamples = 50)
pp_check(mn_sa_bm, nsamples = 50) # similar

pp_check(mt_ss_ba, nsamples = 50)
pp_check(mn_ss_ba, nsamples = 50) # n slightly better
pp_check(mt_ss_bs, nsamples = 50)
pp_check(mn_ss_bs, nsamples = 50) # t clearly better
pp_check(mt_ss_bm, nsamples = 50)
pp_check(mn_ss_bm, nsamples = 50) # t clearly better

pp_check(mt_sm_ba, nsamples = 50)
pp_check(mn_sm_ba, nsamples = 50) # similar
pp_check(mt_sm_bs, nsamples = 50)
pp_check(mn_sm_bs, nsamples = 50) # similar
pp_check(mt_sm_bm, nsamples = 50)
pp_check(mn_sm_bm, nsamples = 50) # similar

# correlation between predicted and observed
cor_fun(d_sa_ba, mt_sa_ba, mn_sa_ba, "ev") # transformed
cor_fun(d_sa_bs, mt_sa_bs, mn_sa_bs, "ev") # non-linear, by 0.03
cor_fun(d_sa_bm, mt_sa_bm, mn_sa_bm, "ev") # definitely transformed

cor_fun(d_ss_ba, mt_ss_ba, mn_ss_ba, "ev") # tranformed
cor_fun(d_ss_bs, mt_ss_bs, mn_ss_bs, "ev") # tranformed
cor_fun(d_ss_bm, mt_ss_bm, mn_ss_bm, "ev") # tranformed by 0.02

cor_fun(d_sm_ba, mt_sm_ba, mn_sm_ba, "mv") # definitely tranformed
cor_fun(d_sm_bs, mt_sm_bs, mn_sm_bs, "mv") # tranformed
cor_fun(d_sm_bm, mt_sm_bm, mn_sm_bm, "mv") # tranformed


#### simplify models ####

# come back to this later: in mv-biomass-by-treatment-2018 and ev-survival-by-treatment-2018 all of the models are very similar and it takes a long time to run the analysis


#### visualize ####

# Prediction datasets
pd_sa_bm <- d_pred_fun(d_sa_bm, mt_sa_bm, mn_sa_bm, 64)
pd_sa_bs <- d_pred_fun(d_sa_bs, mt_sa_bs, mn_sa_bs, 10)
pd_sa_ba <- d_pred_fun(d_sa_ba, mt_sa_ba, mn_sa_ba, 8)

pd_ss_bm <- d_pred_fun(d_ss_bm, mt_ss_bm, mn_ss_bm, 64)
pd_ss_bs <- d_pred_fun(d_ss_bs, mt_ss_bs, mn_ss_bs, 10)
pd_ss_ba <- d_pred_fun(d_ss_ba, mt_ss_ba, mn_ss_ba, 8)

pd_sm_bm <- d_pred_fun(d_sm_bm, mt_sm_bm, mn_sm_bm, 64)
pd_sm_bs <- d_pred_fun(d_sm_bs, mt_sm_bs, mn_sm_bs, 10)
pd_sm_ba <- d_pred_fun(d_sm_ba, mt_sm_ba, mn_sm_ba, 8)

# non-linear models (just to examine, manually updated datasets)
pd_ss_bm %>%
  filter(model == "nonlinear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_ss_bm, aes(y = seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(data = d_ss_bm, aes(y = seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(1))

pd_sm_bs %>%
  filter(model == "nonlinear") %>%
  ggplot(aes(x = counted_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sm_bs, aes(y = seeds_per_bio), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(data = d_sm_bs, aes(y = seeds_per_bio), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(1))

# Mv seeds, each background
p_sm_bm <- pd_sm_bm %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sm_bm, aes(y = log_seeds_pb), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(data = d_sm_bm, aes(y = log_seeds_pb), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(1)) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylab("ln(Microstegium seeds)") +
  ylim(1.5, 5.5)

p_sm_bs <- pd_sm_bs %>%
  filter(model == "linear") %>%
  ggplot(aes(x = counted_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sm_bs, aes(y = log_seeds_pb), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.2)) +
  stat_summary(data = d_sm_bs, aes(y = log_seeds_pb), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylim(1.5, 5.5)

p_sm_ba <- pd_sm_ba %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sm_ba, aes(y = log_seeds_pb), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.2)) +
  stat_summary(data = d_sm_ba, aes(y = log_seeds_pb), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylim(1.5, 5.5)

# Ev seedling seeds, each background
p_ss_bm <- pd_ss_bm %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_ss_bm, aes(y = log_seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(data = d_ss_bm, aes(y = log_seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(1)) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylab("ln(Elymus seedling seeds)") +
  ylim(-4, 3.7)

p_ss_bs <- pd_ss_bs %>%
  filter(model == "linear") %>%
  ggplot(aes(x = counted_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_ss_bs, aes(y = log_seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.2)) +
  stat_summary(data = d_ss_bs, aes(y = log_seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylim(-4, 3.7)

p_ss_ba <- pd_ss_ba %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_ss_ba, aes(y = log_seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.2)) +
  stat_summary(data = d_ss_ba, aes(y = log_seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylim(-4, 3.7)

# Ev adult seeds, each background
p_sa_bm <- pd_sa_bm %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sa_bm, aes(y = log_seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(data = d_sa_bm, aes(y = log_seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(1)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylab("ln(Elymus adult seeds)") +
  xlab("Microstegium density") +
  ylim(0, 6)

p_sa_bs <- pd_sa_bs %>%
  filter(model == "linear") %>%
  ggplot(aes(x = counted_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sa_bs, aes(y = log_seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.2)) +
  stat_summary(data = d_sa_bs, aes(y = log_seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling density") +
  ylim(0, 6)

p_sa_ba <- pd_sa_ba %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_sa_ba, aes(y = log_seeds), width = 0.1, geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.2)) +
  stat_summary(data = d_sa_ba, aes(y = log_seeds), size = 2, shape = 21, color = "black", geom = "point", fun.y = "mean", position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus adult density") +
  ylim(0, 6)

# combine
pdf("./output/seeds-by-treatment-2018-full-models.pdf", height = 9, width = 9)
plot_grid(p_sm_bm, p_sm_bs, p_sm_ba, p_ss_bm, p_ss_bs, p_ss_ba, p_sa_bm, p_sa_bs, p_sa_ba, nrow = 3)
dev.off()
