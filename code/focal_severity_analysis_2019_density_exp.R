##### info ####

# file: focal_severity_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 6/17/20
# goal: analyze focal disease severity


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(glmmTMB)
library(brms)
library(cowplot)

# import all raw data files
dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
edge_dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv") # didn't use this yet...


#### edit data ####

# edge data
edge_dat2 <- edge_dat %>%
  filter(month != "sep") %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         month = recode(month, "may" = "jun", "jun" = "jul", "jul" = "early_aug", "early_aug" = "late_aug", "late_aug" = "sep")) %>%
  select(month, site, plot, treatment, edge_severity)

# calculate severity
# add month columns
# add background disease data
# transform severity
dat2 <- dat %>%
  mutate(leaf_severity = lesion_area.pix / leaf_area.pix,
         plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         plant_severity_adjusted = case_when(plant_severity == 0 ~ 0.001,
                                             plant_severity >= 1 ~ 0.999,
                                             TRUE ~ plant_severity),
         plant_severity_transformed = asin(sqrt(plant_severity_adjusted)),
         Month = recode(month, early_aug = "Early August", jul = "July", 
                        jun = "June", late_aug = "Late August", may = "May", sep = "September") %>%
           fct_relevel("May", "June", "July", "Early August", "Late August", "September"),
         Month_num = as.numeric(Month),
         plant_type = paste(sp, age, sep = " "),
         fungicide = case_when(treatment == "water" ~ 0,
                               TRUE ~ 1),
         exp_plot = paste0(plot, toupper(substr(treatment, 1, 1))),
         sp_trt = paste0(sp, toupper(substr(treatment, 1, 1)))) %>%
  left_join(edge_dat2)

# separate out May data
may_dat <- dat2 %>%
  filter(Month == "May")

# separate out focal data
foc_dat <- dat2 %>%
  filter(focal == 1)

# add treatment plots (repeats no neighbor plots)
foc_dat_plots <- foc_dat %>%
  left_join(plots) %>%
  mutate(background_pres = case_when(density_level == "none" ~ 0,
                                     TRUE ~ 1))

# add treatment plots (repeats no neighbor plots)
may_dat_plots <- may_dat %>%
  left_join(plots)

# only allow for sites with both treatment in September (some lost in the mail)
# information obtained from examine data section
foc_dat2 <- foc_dat %>%
  filter(!(sp == "Ev" & site != "D2" & Month == "September")) %>%
  filter(!(sp == "Mv" & site == "D2" & Month == "September"))

foc_dat_plots2 <- foc_dat_plots %>%
  filter(!(sp == "Ev" & site != "D2" & Month == "September")) %>%
  filter(!(sp == "Mv" & site == "D2" & Month == "September"))

# separate focal data by month and species
may_ev_dat <- filter(foc_dat2, Month == "May" & sp == "Ev")
jun_ev_dat <- filter(foc_dat2, Month == "June" & sp == "Ev")
jul_ev_dat <- filter(foc_dat2, Month == "July" & sp == "Ev")
eau_ev_dat <- filter(foc_dat2, Month == "Early August" & sp == "Ev")
lau_ev_dat <- filter(foc_dat2, Month == "Late August" & sp == "Ev")
sep_ev_dat <- filter(foc_dat2, Month == "September" & sp == "Ev")

jun_mv_dat <- filter(foc_dat2, Month == "June" & sp == "Mv")
jul_mv_dat <- filter(foc_dat2, Month == "July" & sp == "Mv")
eau_mv_dat <- filter(foc_dat2, Month == "Early August" & sp == "Mv")
lau_mv_dat <- filter(foc_dat2, Month == "Late August" & sp == "Mv")
sep_mv_dat <- filter(foc_dat2, Month == "September" & sp == "Mv")


#### visualize ####

# all data
(all_dat_viz <- foc_dat2 %>%
  filter(!is.na(leaf_severity)) %>%
  ggplot(aes(Month_num, leaf_severity, color = as.factor(plot))) +
  stat_summary(fun = "mean", geom = "point", aes(shape = treatment)) +
  stat_summary(fun = "mean", geom = "line", aes(linetype = treatment)) +
  facet_grid(site ~ sp))
# late August and September data probably aren't usable (late August data mostly 1's, September data missing from D1 and D2 because of shipping issue)
# some clear and unexpected separations between fungicide and water treatments that change over time

# remove late data
all_dat_viz %+%
  filter(foc_dat2, !is.na(leaf_severity) & Month_num < 5)
# Why does Mv drop so low in July? generally healthy looking leaves
# generally healthy-looking leaves in June too. The high values seem to be from tip senescence
# lesions aren't super common in early August even. seems like senescence is a big contributor

# look at relationship with edge severity
(edge_viz <- foc_dat2 %>%
  filter(!is.na(leaf_severity) & !is.na(edge_severity)) %>%
  ggplot(aes(edge_severity, leaf_severity, color = treatment)) +
  stat_summary(fun = "mean", geom = "point") +
  facet_grid(sp ~ Month, scales = "free_x"))
# maybe most positive for Mv in early August (makes sense)

edge_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(edge_severity)) %+%
  aes(y = plant_severity)
# similar to leaf metric

# treatment effect
(treat_viz <- foc_dat2 %>%
  filter(!is.na(leaf_severity)) %>%
  ggplot(aes(Month_num, leaf_severity, color = treatment, shape = sp, linetype = sp)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1))

treat_viz %+%
  filter(foc_dat2, !is.na(plant_severity))  %+%
  aes(y = plant_severity)
# much more intuitive pattern

# treatment effect in May for Ev from previous year
(may_treat_viz <- may_dat %>%
  ggplot(aes(treatment, leaf_severity)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1))
# legacy fungicide effect

may_treat_viz %+%
  aes(y = plant_severity)
# similar to leaf metric, but smaller values

# May density effect
(month_viz <- may_dat_plots %>%
  ggplot(aes(background_density, leaf_severity, color = treatment)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
    facet_grid(plant_type ~ background, scales = "free"))
# no clear density effect and a lot of variation within density value

month_viz %+%
  aes(y = plant_severity)
# no clear density effect and a lot of variation within density value

# June density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "June")
# no strong density effects

month_viz %+%
  filter(foc_dat_plots2, Month == "June") %+%
  aes(y = plant_severity)

# July density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "July")
# no strong density effects

month_viz %+%
  filter(foc_dat_plots2, Month == "July") %+%
  aes(y = plant_severity)

# early August density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "Early August")
# mv background Ev seedling is what I would expect

month_viz %+%
  filter(foc_dat_plots2, Month == "Early August") %+%
  aes(y = plant_severity)
# fungicide effects seem stronger when the background is Mv

month_viz %+%
  filter(foc_dat_plots2, Month == "Early August") %+%
  aes(y = plant_severity) +
  facet_grid(sp ~ background, scales = "free")
# mv density almost suppressing Ev disease

# late August density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "Late August")
# no strong density effects

month_viz %+%
  filter(foc_dat_plots2, Month == "Late August") %+%
  aes(y = plant_severity)
# slight Mv on Mv density effect, maybe best example

month_viz %+%
  filter(foc_dat_plots2, Month == "Late August") %+%
  aes(y = plant_severity) +
  facet_grid(sp ~ background, scales = "free")
# no effect on Ev

# September density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "September")
# strong density effect of Mv on Mv 

month_viz %+%
  filter(foc_dat_plots2, Month == "September") %+%
  aes(y = plant_severity)
# strong density effect of Mv on Mv 


#### examine data ####

# September data for Mv
foc_dat %>%
  filter(sp == "Mv" & Month == "September") %>%
  group_by(plot, treatment) %>%
  summarise(sites = length(unique(site)),
            plants = n(),
            site_names = paste(unique(site), collapse = ",")) %>%
  data.frame()
# one more site for all water plots (D2)

# September data for Ev
foc_dat %>%
  filter(sp == "Ev" & Month == "September") %>%
  group_by(plot, treatment) %>%
  summarise(sites = length(unique(site)),
            plants = n(),
            site_names = paste(unique(site), collapse = ",")) %>%
  data.frame()
# D2 is the only site used throughout

# high lesion values
dat2 %>%
  filter(plant_severity > 1)

# make edits in edit data section


#### edge analysis ####

# divide data
eau_ev_edge_dat <- filter(eau_ev_dat, treatment == "water")
lau_ev_edge_dat <- filter(lau_ev_dat, treatment == "water")
lau_mv_edge_dat <- filter(lau_mv_dat, treatment == "water")

# early August ev
eau_ev_edge_mod <- glmmTMB(plant_severity_adjusted ~ edge_severity, data = eau_ev_edge_dat, family = "beta_family")
summary(eau_ev_edge_mod)
# positive, not sig

# late August mv
lau_mv_edge_mod <- glmmTMB(plant_severity_adjusted ~ edge_severity, data = lau_mv_edge_dat, family = "beta_family")
summary(lau_mv_edge_mod)
# positive, significant

# late August ev
lau_ev_edge_mod <- glmmTMB(plant_severity_adjusted ~ edge_severity, data = lau_ev_edge_dat, family = "beta_family")
summary(lau_ev_edge_mod)
# positive, significant


#### treatment analysis ####

# edge severity is not significantly different between fungicide and water plots (edge_severity_analysis_2019_density_exp.R)

# all Ev May
may_ev_all_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site/exp_plot), data = may_dat, family = "beta_family")
summary(may_ev_all_mod)

# significance codes dataset
foc_trt_sig = tibble(Month = rep(levels(dat2$Month), 2),
                     sp = rep(c("Ev", "Mv"), each = 6)) %>%
  mutate(sig = "",
         treatment = NA_character_) %>%
  filter(!(Month == "May" & sp == "Mv"))

# focal Ev May
may_ev_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
summary(may_ev_mod)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "May" & sp == "Ev" ~ "*",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "May" & sp == "Ev" ~ "water",
                                           TRUE ~ treatment))

# focal Ev June
jun_ev_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = jun_ev_dat, family = "beta_family")
summary(jun_ev_mod)
jun_ev_mod2 <- update(jun_ev_mod, .~. - edge_severity)
summary(jun_ev_mod2)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "June" & sp == "Ev" ~ "*",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "June" & sp == "Ev" ~ "water",
                                           TRUE ~ treatment))

# focal Ev July
jul_ev_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = jul_ev_dat, family = "beta_family")
summary(jul_ev_mod)
jul_ev_mod2 <- update(jul_ev_mod, .~. - edge_severity)
summary(jul_ev_mod2)

# focal Ev Early August
eau_ev_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = eau_ev_dat, family = "beta_family")
summary(eau_ev_mod)
eau_ev_mod2 <- update(eau_ev_mod, .~. - edge_severity)
summary(eau_ev_mod2)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "Early August" & sp == "Ev" ~ "**",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "Early August" & sp == "Ev" ~ "water",
                                           TRUE ~ treatment))

# focal Ev Late August
lau_ev_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = lau_ev_dat, family = "beta_family")
summary(lau_ev_mod)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "Late August" & sp == "Ev" ~ ".",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "Late August" & sp == "Ev" ~ "fungicide",
                                           TRUE ~ treatment))

# focal Ev September
sep_ev_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = sep_ev_dat, family = "beta_family")
summary(sep_ev_mod)
sep_ev_mod2 <- update(sep_ev_mod, .~. - edge_severity)
summary(sep_ev_mod2)

# focal Mv June
jun_mv_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = jun_mv_dat, family = "beta_family")
summary(jun_mv_mod)
jun_mv_mod2 <- update(jun_mv_mod, .~. - edge_severity)
summary(jun_mv_mod2)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "June" & sp == "Mv" ~ "*",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "June" & sp == "Mv" ~ "water",
                                           TRUE ~ treatment))

# focal Mv July
jul_mv_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = jul_mv_dat, family = "beta_family")
summary(jul_mv_mod)
jul_mv_mod2 <- update(jul_mv_mod, .~. - edge_severity)
summary(jul_mv_mod2)

# focal Mv Early August
eau_mv_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = eau_mv_dat, family = "beta_family")
summary(eau_mv_mod)
eau_mv_mod2 <- update(eau_mv_mod, .~. - edge_severity)
summary(eau_mv_mod2)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "Early August" & sp == "Mv" ~ "***",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "Early August" & sp == "Mv" ~ "water",
                                           TRUE ~ treatment))

# focal Mv Late August
lau_mv_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = lau_mv_dat, family = "beta_family")
summary(lau_mv_mod)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "Late August" & sp == "Mv" ~ "***",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "Late August" & sp == "Mv" ~ "water",
                                           TRUE ~ treatment))

# focal Mv September
sep_mv_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + edge_severity + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")
summary(sep_mv_mod)
sep_mv_mod2 <- update(sep_mv_mod, .~. - edge_severity)
summary(sep_mv_mod2)
foc_trt_sig = mutate(foc_trt_sig,
                     sig = case_when(Month == "September" & sp == "Mv" ~ "***",
                                     TRUE ~ sig),
                     treatment = case_when(Month == "September" & sp == "Mv" ~ "water",
                                           TRUE ~ treatment))

#### density analysis ####

# select Mv background data from late August
lau_mv_mv_dat <- foc_dat_plots2 %>%
  filter(background == "Mv seedling" &
           Month == "Late August" & 
           sp == "Mv")

# select Ev with Mv background for early August
eau_ev_mv_dat <- foc_dat_plots2 %>%
  filter(background == "Mv seedling" &
           Month == "Early August" & 
           sp == "Ev") 

# fit model
lau_mv_mv_den_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * background_density + edge_severity + (1|site), data = lau_mv_mv_dat, family = "beta_family")
summary(lau_mv_mv_den_mod)
# only fungicide and edge sig

eau_ev_mv_den_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * background_density + edge_severity + (1|site), data = eau_ev_mv_dat, family = "beta_family")
summary(eau_ev_mv_den_mod)
# fungicide marginally sig

# check model fit
lau_mv_mv_dat %>%
  mutate(pred = predict(lau_mv_mv_den_mod, newdata = lau_mv_mv_dat, type = "response", re.form = NA)) %>%
  ggplot(aes(background_density, plant_severity_adjusted, color = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(aes(y = pred), fun = "mean", geom = "line")

eau_ev_mv_dat %>%
  mutate(pred = predict(eau_ev_mv_den_mod, newdata = eau_ev_mv_dat, type = "response", re.form = NA)) %>%
  ggplot(aes(background_density, plant_severity_adjusted, color = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(aes(y = pred), fun = "mean", geom = "line")
# general negative effect of density

# # tried a logistic regression with glmer (cbind(tiller_area_infec.pix, tiller_area_healthy.pix))
# # could not address convergence and large eigenvalue issues
# # tried a negative binomial regression with tiller_area_infec.pix ~ offset(tiller_area.pix)
# # also received convergence warnings
# 
# # non-linear model with transformed severity
# 
# # priors
# filter(lau_mv_mv_dat, background_density == 0) %>%
#   group_by(treatment) %>%
#   summarise(R0 = mean(plant_severity_transformed))
# 
# filter(lau_mv_mv_dat, background_density > 10) %>%
#   group_by(treatment) %>%
#   summarise(asym = mean(plant_severity_transformed))
# 
# # distributions
# x <- seq(0, 10, length.out = 100)
# y <- dgamma(x, shape = 1, scale = 1) # note that this scale is 1/(stan scale)
# plot(x, y, type = "l")
# 
# # model
# lau_mv_mv_mod2 <- brm(data = lau_mv_mv_dat, family = gaussian,
#                  bf(plant_severity_transformed ~ asym + (R0 - asym) * exp(-exp(lrc) * background_density), 
#                     asym ~ 0 + treatment, 
#                     R0 ~ 0 + treatment,
#                     lrc ~ 0 + treatment,
#                     nl = T),
#                  prior <- c(prior(gamma(1, 1), nlpar = "asym", lb = 0),
#                             prior(gamma(1, 1), nlpar = "R0", lb = 0),
#                             prior(gamma(1, 1), nlpar = "lrc", lb = 0),
#                             prior(cauchy(0, 1), class = sigma)),
#                  iter = 6000, warmup = 1000, chains = 1, cores = 1)
# summary(lau_mv_mv_mod2)
# 
# # simulate data
# lau_mv_mv_sim_dat <- tibble(background_density = seq(0, 64, length.out = 300)) %>%
#   merge(tibble(treatment = c("water", "fungicide")),
#         all = T) %>%
#   as_tibble() %>%
#   mutate(plant_severity_transformed = fitted(lau_mv_mv_mod2, newdata = .)[, "Estimate"],
#          lower = fitted(lau_mv_mv_mod2, newdata = .)[, "Q2.5"],
#          upper = fitted(lau_mv_mv_mod2, newdata = .)[, "Q97.5"]) 
# 
# # plot model
# ggplot(lau_mv_mv_dat, aes(x = background_density, y = plant_severity_transformed, color = treatment, fill = treatment)) +
#   stat_summary(geom = "point", fun = "mean", size = 2) +
#   stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
#   geom_line(data = lau_mv_mv_sim_dat)+
#   geom_ribbon(data = lau_mv_mv_sim_dat, aes(ymin = lower, ymax = upper), alpha = 0.3)
# # not a great fit, mostly showing differences at high density

# significance codes dataset
foc_den_sig = tibble(background_density = rep(c(0, 4, 16, 64), 2),
                     sp = rep(c("Ev", "Mv"), each = 4)) %>%
  mutate(sig = "",
         treatment = NA_character_)

# # divide data by density treatment
# lau_mv_mv_0_dat <- filter(lau_mv_mv_dat, background_density == 0)
# lau_mv_mv_4_dat <- filter(lau_mv_mv_dat, background_density == 4)
# lau_mv_mv_16_dat <- filter(lau_mv_mv_dat, background_density == 16)
# lau_mv_mv_64_dat <- filter(lau_mv_mv_dat, background_density == 64)
# 
# eau_ev_mv_0_dat <- filter(eau_ev_mv_dat, background_density == 0)
# eau_ev_mv_4_dat <- filter(eau_ev_mv_dat, background_density == 4)
# eau_ev_mv_16_dat <- filter(eau_ev_mv_dat, background_density == 16)
# eau_ev_mv_64_dat <- filter(eau_ev_mv_dat, background_density == 64)
# 
# # 0 density Mv
# lau_mv_mv_0_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = lau_mv_mv_0_dat, family = "beta_family")
# summary(lau_mv_mv_0_mod)
# foc_den_sig = mutate(foc_den_sig,
#                      sig = case_when(background_density == 0 & sp == "Mv" ~ ".",
#                                      TRUE ~ sig),
#                      treatment = case_when(background_density == 0 & sp == "Mv" ~ "water",
#                                            TRUE ~ treatment))
# 
# # 4 density Mv
# lau_mv_mv_4_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = lau_mv_mv_4_dat, family = "beta_family")
# summary(lau_mv_mv_4_mod)
# foc_den_sig = mutate(foc_den_sig,
#                      sig = case_when(background_density == 4 & sp == "Mv" ~ "**",
#                                      TRUE ~ sig),
#                      treatment = case_when(background_density == 4 & sp == "Mv" ~ "water",
#                                            TRUE ~ treatment))
# 
# # 16 density Mv
# lau_mv_mv_16_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = lau_mv_mv_16_dat, family = "beta_family")
# summary(lau_mv_mv_16_mod)
# foc_den_sig = mutate(foc_den_sig,
#                      sig = case_when(background_density == 16 & sp == "Mv" ~ "**",
#                                      TRUE ~ sig),
#                      treatment = case_when(background_density == 16 & sp == "Mv" ~ "water",
#                                            TRUE ~ treatment))
# 
# # 64 density Mv
# lau_mv_mv_64_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = lau_mv_mv_64_dat, family = "beta_family")
# summary(lau_mv_mv_64_mod)
# 
# # 0 density Ev
# eau_ev_mv_0_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = eau_ev_mv_0_dat, family = "beta_family")
# summary(eau_ev_mv_0_mod)
# 
# # 4 density Ev
# eau_ev_mv_4_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = eau_ev_mv_4_dat, family = "beta_family")
# summary(eau_ev_mv_4_mod)
# foc_den_sig = mutate(foc_den_sig,
#                      sig = case_when(background_density == 4 & sp == "Ev" ~ "*",
#                                      TRUE ~ sig),
#                      treatment = case_when(background_density == 4 & sp == "Ev" ~ "water",
#                                            TRUE ~ treatment))
# 
# # 16 density
# eau_ev_mv_16_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = eau_ev_mv_16_dat, family = "beta_family")
# summary(eau_ev_mv_16_mod)
# 
# # 64 density
# eau_ev_mv_64_mod <- glmmTMB(plant_severity_adjusted ~ fungicide + (1|site), data = eau_ev_mv_64_dat, family = "beta_family")
# summary(eau_ev_mv_64_mod)
# foc_den_sig = mutate(foc_den_sig,
#                      sig = case_when(background_density == 64 & sp == "Ev" ~ "*",
#                                      TRUE ~ sig),
#                      treatment = case_when(background_density == 64 & sp == "Ev" ~ "water",
#                                            TRUE ~ treatment))


#### Microstegium analysis ####

# mv
lau_mv_mv_pres_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * background_pres + edge_severity + (1|site), data = lau_mv_mv_dat, family = "beta_family")
summary(lau_mv_mv_pres_mod)
# presence not sig

# ev
eau_ev_mv_pres_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * background_pres + edge_severity + (1|site), data = eau_ev_mv_dat, family = "beta_family")
summary(eau_ev_mv_pres_mod)
# none sig


#### figures ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "italic"),
        strip.placement = "outside")

# colors
col_pal = c("#018571", "#a6611a")

# shapes
shape_pal = c(21, 24)

# label size
label_size = 4

# join treatment sig to rest of treatment data
# capitalize columns
# calculate mean and 95% CI
# order month again
foc_trt_dat <- foc_dat2 %>%
  left_join(foc_trt_sig) %>%
  rename(Treatment = treatment, Species = sp) %>%
  filter(!is.na(plant_severity_adjusted))%>%
  group_by(Month, Species, Treatment, sp_trt, sig) %>%
  summarise(Severity = mean_cl_boot(plant_severity_adjusted)$y,
            lower = mean_cl_boot(plant_severity_adjusted)$ymin,
            upper = mean_cl_boot(plant_severity_adjusted)$ymax) %>%
  ungroup() %>%
  mutate(Month = recode(Month, "Early August" = "Early\nAugust", "Late August" = "Late\nAugust") %>% 
           fct_relevel("May", "June", "July", "Early\nAugust", "Late\nAugust", "September"),
         Treatment = recode(Treatment, water = "control (water)"),
         September = case_when(Month == "September" ~ "yes",
                               TRUE ~ "no"),
         Species = recode(Species, Ev = "Elymus", Mv = "Microstegium"))

# treatment effect over time
foc_trt_plot <- foc_trt_dat %>%
  ggplot(aes(Month, Severity, group = sp_trt)) +
  geom_line(aes(linetype = Species, color = Treatment)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, alpha = September), width = 0.1) +
  geom_point(aes(shape = Species, fill = Treatment, alpha = September), size = 3) +
  geom_text(aes(y = upper + 0.02, label = sig, alpha = September), check_overlap = T, size = 4) +
  annotate("text", label = "A", x = 0.7, y = 0.8, size = label_size) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  scale_alpha_manual(values = c(1, 0.5), guide = F) +
  guides(fill = guide_legend(override.aes = list(shape = shape_pal[1])), 
         shape = guide_legend(label.theme = element_text(size = 8, face = "italic")), 
         linetype = guide_legend(label.theme = element_text(size = 8, face = "italic"))) +
  xlab("Four-week interval") +
  ylab("Proportion leaf area brown") +
  temp_theme +
  theme(legend.position = c(0.25, 0.7))
# foc_trt_plot

# join edge data with predictions and labels
# capitalize columns
# calculate mean and 95% CI
foc_edge_dat <- lau_ev_edge_dat %>%
  mutate(pred = predict(lau_ev_edge_mod, newdata = lau_ev_edge_dat, re.form = NA, type = "response"),
         pred.se = predict(lau_ev_edge_mod, newdata = lau_ev_edge_dat, re.form = NA, type = "response", se.fit = T)$se.fit) %>%
  full_join(lau_mv_edge_dat %>%
              mutate(pred = predict(lau_mv_edge_mod, newdata = lau_mv_edge_dat, re.form = NA, type = "response"),
                     pred.se = predict(lau_mv_edge_mod, newdata = lau_mv_edge_dat, re.form = NA, type = "response", se.fit = T)$se.fit)) %>%
  mutate(label = case_when(sp == "Ev" ~ "B",
                           sp == "Mv" ~ "C"),
         label_x = case_when(sp == "Ev" ~ 0.001,
                             sp == "Mv" ~ 0.01)) %>%
  rename(Treatment = treatment, Species = sp, Edge = edge_severity, Severity = plant_severity_adjusted) %>%
  mutate(Treatment = recode(Treatment, water = "control (water)"),
         Species = recode(Species, Ev = "Elymus", Mv = "Microstegium"))

# treatment effect over density
foc_edge_plot <- foc_edge_dat %>%
  ggplot(aes(Edge, Severity)) +
  geom_point(aes(shape = Species), fill = col_pal[1]) +
  geom_line(aes(y = pred)) +
  geom_ribbon(aes(ymin = pred - pred.se, ymax = pred + pred.se), alpha = 0.5, fill = col_pal[1]) +
  geom_text(aes(label = label, x = label_x), y = 1.05, check_overlap = T, size = label_size) +
  facet_wrap(~ Species, nrow = 1, scales = "free_x") +
  scale_shape_manual(values = shape_pal) +
  xlab("Environmental disease pressure") +
  ylab("Proportion leaf area brown") +
  ylim(0, 1.05) +
  temp_theme +
  theme(legend.position = "none")
# foc_edge_plot

# join density data with labels
# capitalize columns
# calculate mean and 95% CI
foc_den_dat <- foc_dat_plots2 %>%
  filter(background == "Mv seedling" & ((Month == "Early August" & sp == "Ev") | (Month == "Late August" & sp == "Mv"))) %>%
  mutate(label = case_when(sp == "Ev" ~ "D",
                           sp == "Mv" ~ "E")) %>%
  rename(Treatment = treatment, Species = sp, Density = background_density) %>%
  filter(!is.na(plant_severity_adjusted))%>%
  group_by(Density, Species, Treatment, sp_trt, label) %>%
  summarise(Severity = mean_cl_boot(plant_severity_adjusted)$y,
            lower = mean_cl_boot(plant_severity_adjusted)$ymin,
            upper = mean_cl_boot(plant_severity_adjusted)$ymax) %>%
  ungroup() %>%
  mutate(Treatment = recode(Treatment, water = "control (water)"))

# treatment effect over density
foc_den_plot <- foc_den_dat %>%
  ggplot(aes(Density, Severity, group = sp_trt)) +
  geom_line(aes(linetype = Species, color = Treatment)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 3) +
  geom_point(aes(shape = Species, fill = Treatment), size = 3) +
  geom_text(aes(label = label), x = -0.3, y = 0.65, check_overlap = T, size = label_size) +
  facet_wrap(~ Species, nrow = 1) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  guides(fill = guide_legend(override.aes = list(shape = shape_pal[1]))) +
  xlab(expression(paste(italic("Microstegium"), " density", sep = ""))) +
  ylab("Proportion leaf area brown") +
  temp_theme +
  theme(legend.position = "none",
        strip.text = element_blank())
# foc_den_plot

# combine figures
foc_right_plots <- plot_grid(foc_edge_plot, foc_den_plot,
                             nrow = 2)
foc_plot <- plot_grid(foc_trt_plot, foc_right_plots,
                      nrow = 1,
                      rel_widths = c(0.85, 1))

#### output ####
pdf("./output/focal_severity_analysis_2019_density_exp.pdf", width = 7, height = 5)
foc_plot
dev.off()
