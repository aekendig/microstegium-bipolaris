##### info ####

# file: mv-biomass-by-treatment-2018
# author: Amy Kendig
# date last edited: 6/4/19
# goal: see how treatments and disease severity affected Mv biomass


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)
library(loo)
library(cowplot)

# run survival files
source("./code/ev-survival-data-processing-2018.R")
rm(list = setdiff(ls(), "esurv"))
source("./code/mv-survival-data-processing-2018.R")
rm(list = setdiff(ls(), c("esurv", "msurv")))
source("./code/bg-densities-data-processing-2018.R")
rm(list = setdiff(ls(), c("esurv", "msurv", "bgd")))

# run leaf scan files
source("./code/mv-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), c("mleaf", "esurv", "msurv", "bgd")))

# run covariate file
source("./code/covariate-data-processing-2018.R")
rm(list = setdiff(ls(), c("mleaf", "covar", "esurv", "msurv", "bgd")))

# import data
mb <- read_csv("./data/mv-biomass-oct-2018-density-exp.csv")
trt <- read_csv("./data/plot-treatments-for-figures-2018-density-exp.csv")
trt_s <- read_csv("./data/plot-treatments-2018-density-exp.csv")
til <- read_csv("./data/mv-disease-sep-2018-density-exp.csv")

# non-linear species effects function
el_fun <- function(y0, x, a) {
  y = y0 * exp(a * log(x + 1))
  return(y)
}

# function to create datasets for prediction
d_pred_fun = function(max_val, length_out, bg_group){
  
  d_pred <- tibble(background_density = rep(seq(0, max_val, length.out = length_out), 2), 
                   treatment = rep(c("water", "fungicide"), each = length_out),
                   sm_adj = mean(filter(dat, background == bg_group)$sm_adj),
                   cc_adj = mean(filter(dat, background == bg_group)$cc_adj),
                   pm_adj = mean(filter(dat, background == bg_group)$pm_adj))
  
  return(d_pred)
}

# plot parameters
colpal = c("#3CBB75FF", "#39568CFF")
sm_txt = 12
lg_txt = 14
an_txt = 3


#### edit data ####

# edit leaf scan data for focals
ml_f <- mleaf %>%
  ungroup() %>%
  filter(month == "September" & focal == 1) %>%
  select(site, plot, treatment, ID, lesion_area.pix, green_area.pix, leaf_area.pix) %>%
    mutate(plot = as.numeric(plot)) %>%
  full_join(filter(til, ID %in% c("1", "2", "3"))) %>%
  mutate(leaves_infec2 = ifelse(leaves_infec == 0 & !is.na(leaf_area.pix), 1, leaves_infec),
         ind_severity = (lesion_area.pix - green_area.pix) * leaves_infec2 / (leaf_area.pix * leaves_tot)) %>%
  group_by(site, plot, treatment) %>%
  summarise(severity = mean(ind_severity, na.rm = T),
            n = sum(!is.na(ind_severity)),
            severity_se = sd(ind_severity, na.rm = T)/n)

# edit leaf scan data for background (skip tillers - only have background numbers for Mv background plots)
ml_b <- mleaf %>%
  ungroup() %>%
  filter(month == "September" & focal == 0) %>%
  select(site, plot, treatment, ID, lesion_area.pix, green_area.pix, leaf_area.pix) %>%
  mutate(ind_severity = (lesion_area.pix - green_area.pix) / leaf_area.pix,
         plot = as.numeric(plot)) %>%
  group_by(site, plot, treatment) %>%
  summarise(bg_severity = mean(ind_severity, na.rm = T),
            bg_n = sum(!is.na(ind_severity)),
            bg_severity_se = sd(ind_severity, na.rm = T) / bg_n)

# check notes
unique(mb$processing_notes)

# edit values and merge with other data
dat <- mb %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment),
         log_bio.g = log(bio.g)) %>%
  left_join(trt) %>%
  left_join(covar) %>%
  mutate(density_level = factor(density_level, levels = c("none", "low", "medium", "high")),
         disease = ifelse(treatment == "fungicide", 0, 1))

# edit values and merge with simplified treatment data
dat_s <- mb %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment),
         log_bio.g = log(bio.g)) %>%
  left_join(trt_s) %>%
  left_join(ml_f) %>%
  left_join(ml_b) %>%
  left_join(covar) %>%
  mutate(density_level = factor(density_level, levels = c("none", "low", "medium", "high")),
         disease = ifelse(treatment == "fungicide", 0, 1))

# subset by background group
d_ea = filter(dat, background == "Ev adult")
d_es = filter(dat, background == "Ev seedling")
d_ms = filter(dat, background == "Mv seedling")


##### visualize #####

# check for extreme values
dat_s %>%
  ggplot(aes(x = bio.g)) +
  geom_histogram()

filter(dat_s, bio.g > 30) # these values match the hand-written weights

# treatment effects
dat %>%
  ggplot(aes(x = density_level, y = bio.g, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background)

dat %>%
  ggplot(aes(x = density_level, y = log_bio.g, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background)

dat %>%
  ggplot(aes(x = background_density, y = log_bio.g, colour = treatment)) + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background, scales = "free")

dat_s %>%
  ggplot(aes(x = treatment, y = bio.g, colour = density_level))  + 
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~background)

# severity
dat_s %>%
  ggplot(aes(x = severity, y = bio.g, colour = treatment)) + 
  geom_point()

dat_s %>%
  ggplot(aes(x = bg_severity, y = bio.g, colour = treatment)) + 
  geom_point()


#### set up models ####

# full linear log-transformed models
# mt_ea <- brm(data = d_ea, family = gaussian,
#            log_bio.g ~ background_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
#            prior <- c(prior(normal(0, 100), class = Intercept),
#                       prior(normal(0, 10), class = b),
#                       prior(cauchy(0, 1), class = sd)),
#            iter = 6000, warmup = 1000, chains = 3, cores = 2,
#            control = list(adapt_delta = 0.999))
# summary(mt_ea)
# save(mt_ea, file = "./output/mv-biomass-by-treatment-2018-log-transformed-ev-adult.rda")
load("./output/mv-biomass-by-treatment-2018-log-transformed-ev-adult.rda")

# mt_es <- update(mt_ea, newdata = d_es,
#                 control = list(adapt_delta = 0.9999))
# summary(mt_es)
# save(mt_es, file = "./output/mv-biomass-by-treatment-2018-log-transformed-ev-seedling.rda")
load("./output/mv-biomass-by-treatment-2018-log-transformed-ev-seedling.rda")

# mt_ms <- update(mt_ea, newdata = d_ms,
#                 control = list(adapt_delta = 0.9999))
# summary(mt_ms)
# save(mt_ms, file = "./output/mv-biomass-by-treatment-2018-log-transformed-mv.rda")
load("./output/mv-biomass-by-treatment-2018-log-transformed-mv.rda")

# exponential untransformed model

# y0 is the value of biomass when density is zero
# just use one dataset so that none plots aren't counted multiple times
d_ea %>%
  filter(background_density == 0) %>%
  summarise(mean_bio = mean(bio.g), 
            min_bio = min(bio.g), 
            max_bio = max(bio.g))
# mean = 13.8, sd = 10
y0 = 13.8

# a is the rate of increase or decrease - simulate data
sim_dat = expand.grid(0:64, c(seq(-2, 0, by = 0.5), seq(0.01, 0.1, by = 0.01)))
colnames(sim_dat) = c("density", "a")
sim_dat <- sim_dat %>%
  mutate(bio = el_fun(y0, density, a),
         pos = ifelse(a >= 0, 1, 0))

# figure
sim_dat %>%
  ggplot(aes(x = density, y = bio, colour = as.factor(a))) + 
  geom_line() +
  facet_wrap(~pos, scales = "free")
# mean = 0, sd = 0.5 
  
# models
# mn_ea <- brm(data = d_ea, family = gaussian,
#              bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T),
#              prior <- c(prior(normal(13.8, 10), nlpar = "y0"),
#                         prior(normal(0, 0.5), nlpar = "a")),
#              iter = 6000, warmup = 1000, chains = 3, cores = 2,
#              control = list(adapt_delta = 0.99))
# summary(mn_ea)
# save(mn_ea, file = "./output/mv-biomass-by-treatment-2018-nonlinear-ev-adult.rda")
load("./output/mv-biomass-by-treatment-2018-nonlinear-ev-adult.rda")

# mn_es <- update(mn_ea, newdata = d_es)
# summary(mn_es)
# save(mn_es, file = "./output/mv-biomass-by-treatment-2018-nonlinear-ev-seedling.rda")
load("./output/mv-biomass-by-treatment-2018-nonlinear-ev-seedling.rda")

# mn_ms <- update(mn_ea, newdata = d_ms)
# summary(mn_ms)
# save(mn_ms, file = "./output/mv-biomass-by-treatment-2018-nonlinear-mv.rda")
load("./output/mv-biomass-by-treatment-2018-nonlinear-mv.rda")


#### check models ####

# convergence of chains
plot(mt_ea)
plot(mt_es)
plot(mt_ms)

plot(mn_ea)
plot(mn_es)
plot(mn_ms)

# compare pp_check
pp_check(mt_ea, nsamples = 50)
pp_check(mn_ea, nsamples = 50) # similar

pp_check(mt_es, nsamples = 50)
pp_check(mn_es, nsamples = 50) # similar

pp_check(mt_ms, nsamples = 50)
pp_check(mn_ms, nsamples = 50) # similar

# correlation between predicted and observed
d_ea_comp <- d_ea %>%
  mutate(transformed = predict(mt_ea)[, 1],
         nonlinear = predict(mn_ea)[, 1])
cor.test(d_ea_comp$transformed, d_ea_comp$log_bio.g)
cor.test(d_ea_comp$nonlinear, d_ea_comp$bio.g) # both close to 0.7, transformed a little higher

d_es_comp <- d_es %>%
  mutate(transformed = predict(mt_es)[, 1],
         nonlinear = predict(mn_es)[, 1])
cor.test(d_es_comp$transformed, d_es_comp$log_bio.g)
cor.test(d_es_comp$nonlinear, d_es_comp$bio.g) # transformed higher (0.45 vs. 0.33) & sig

d_ms_comp <- d_ms %>%
  mutate(transformed = predict(mt_ms)[, 1],
         nonlinear = predict(mn_ms)[, 1])
cor.test(d_ms_comp$transformed, d_ms_comp$log_bio.g)
cor.test(d_ms_comp$nonlinear, d_ms_comp$bio.g) # transformed higher (0.78 vs. 0.69)

# visualize

# Ev adult
d_plot_ea <- d_ea %>%
  select(background_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g) %>%
  mutate(model = "linear") %>%
  rbind(d_ea %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "nonlinear"))
  
d_pred_ea <- d_pred_fun(8, 100, "Ev adult") %>%
  mutate(pred = fitted(mt_ea, newdata = d_pred_fun(8, 100, "Ev adult"), re_formula = NA, nsamples = 100)[,1],
         lower = fitted(mt_ea, newdata = d_pred_fun(8, 100, "Ev adult"), re_formula = NA, nsamples = 100)[,3],
         upper = fitted(mt_ea, newdata = d_pred_fun(8, 100, "Ev adult"), re_formula = NA, nsamples = 100)[,4],
         model = "linear") %>%
  full_join(d_pred_fun(8, 100, "Ev adult") %>%
              mutate(pred = fitted(mn_ea, newdata = d_pred_fun(8, 100, "Ev adult"), re_formula = NA, nsamples = 100)[,1],
                     lower = fitted(mn_ea, newdata = d_pred_fun(8, 100, "Ev adult"), re_formula = NA, nsamples = 100)[,3],
                     upper = fitted(mn_ea, newdata = d_pred_fun(8, 100, "Ev adult"), re_formula = NA, nsamples = 100)[,4],
                     model = "nonlinear"))

pdf("./output/mv-biomass-by-treatment-2018-ev-adult-full-models.pdf", height = 5, width = 7)
d_pred_ea %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_plot_ea, aes(y = bio.g), fun.data = "mean_se", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = d_plot_ea, aes(y = bio.g), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_wrap(~model, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
           axis.text = element_text(color = "black", size = sm_txt),
           strip.text = element_text(color = "black", size = lg_txt),
           legend.title = element_text(color = "black", size = sm_txt),
           legend.text = element_text(color = "black", size = sm_txt),
           legend.position = c(0.9, 0.9),
           legend.background = element_blank(),
           legend.key = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Adult Elymus density") +
  ylab("Microstegium biomass (g)")
dev.off()

# Ev seedling
d_plot_es <- d_es %>%
  select(background_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g) %>%
  mutate(model = "linear") %>%
  rbind(d_es %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "nonlinear"))

d_pred_es <- d_pred_fun(16, 100, "Ev seedling") %>%
  mutate(pred = fitted(mt_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,1],
         lower = fitted(mt_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,3],
         upper = fitted(mt_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,4],
         model = "linear") %>%
  full_join(d_pred_fun(16, 100, "Ev seedling") %>%
              mutate(pred = fitted(mn_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,1],
                     lower = fitted(mn_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,3],
                     upper = fitted(mn_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,4],
                     model = "nonlinear"))

pdf("./output/mv-biomass-by-treatment-2018-ev-seedling-full-models.pdf", height = 5, width = 7)
d_pred_es %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_plot_es, aes(y = bio.g), fun.data = "mean_se", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = d_plot_es, aes(y = bio.g), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_wrap(~model, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling density") +
  ylab("Microstegium biomass (g)")
dev.off()

# Mv seedling
d_plot_ms <- d_ms %>%
  select(background_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g) %>%
  mutate(model = "linear") %>%
  rbind(d_ms %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "nonlinear"))

d_pred_ms <- d_pred_fun(64, 100, "Mv seedling") %>%
  mutate(pred = fitted(mt_ms, newdata = d_pred_fun(64, 100, "Mv seedling"), re_formula = NA, nsamples = 100)[,1],
         lower = fitted(mt_ms, newdata = d_pred_fun(64, 100, "Mv seedling"), re_formula = NA, nsamples = 100)[,3],
         upper = fitted(mt_ms, newdata = d_pred_fun(64, 100, "Mv seedling"), re_formula = NA, nsamples = 100)[,4],
         model = "linear") %>%
  full_join(d_pred_fun(64, 100, "Mv seedling") %>%
              mutate(pred = fitted(mn_ms, newdata = d_pred_fun(64, 100, "Mv seedling"), re_formula = NA, nsamples = 100)[,1],
                     lower = fitted(mn_ms, newdata = d_pred_fun(64, 100, "Mv seedling"), re_formula = NA, nsamples = 100)[,3],
                     upper = fitted(mn_ms, newdata = d_pred_fun(64, 100, "Mv seedling"), re_formula = NA, nsamples = 100)[,4],
                     model = "nonlinear"))

pdf("./output/mv-biomass-by-treatment-2018-mv-full-models.pdf", height = 5, width = 7)
d_pred_ms %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = d_plot_ms, aes(y = bio.g), fun.data = "mean_se", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = d_plot_ms, aes(y = bio.g), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_wrap(~model, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Microstegium density") +
  ylab("Microstegium biomass (g)")
dev.off()


#### models with actual density values ####

# survival rates for tracked individuals
msurv %>%
  filter(month == "September" & !is.na(survival)) %>%
  summarise(sum(survival) / length(survival))
# high survival and poor estimates for background

esurv %>%
  filter(month == "September" & !is.na(survival)) %>%
  group_by(age, focal) %>%
  summarise(sum(survival) / length(survival))
# adult background survival is 100%
# only really need to do the Ev seedling plots

bgd %>%
  group_by(plot) %>%
  summarise(mean(counted_density, na.rm = T))

bgd %>%
  filter(plot == 9) # one adult died

# merge Ev seedling data with counted densities
d_es_2 <- d_es %>%
  left_join(bgd)

# update models with counted density
# mt_es_2 <- update(mt_es, 
#                   formula = log_bio.g ~ counted_density * treatment + sm_adj + cc_adj + pm_adj + (1 | site),
#                   newdata = d_es_2,
#                   control = list(adapt_delta = 0.999))
# save(mt_es_2, file = "./output/mv-biomass-by-treatment-2018-log-transformed-ev-seedling-counted.rda")
load("./output/mv-biomass-by-treatment-2018-log-transformed-ev-seedling-counted.rda")

# mn_es_2 <- update(mn_es, 
#                   formula = bf(bio.g ~ y0 * exp(a * log(counted_density + 1)), 
#                                y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), 
#                                a ~ treatment, nl = T), 
#                   newdata = d_es_2)
# save(mn_es_2, file = "./output/mv-biomass-by-treatment-2018-nonlinear-ev-seedling-counted.rda")
load("./output/mv-biomass-by-treatment-2018-nonlinear-ev-seedling-counted.rda")

# PP check with planted density models
pp_check(mt_es, nsamples = 50)
pp_check(mt_es_2, nsamples = 50) # similar, maybe more underestimates

pp_check(mn_es, nsamples = 50)
pp_check(mn_es_2, nsamples = 50) # similar

# loo comparison with planted density models
loo_es_t = list(loo(mt_es, reloo = T), loo(mt_es_2, reloo = T)) 
loo_es_t
loo_compare(loo_es_t) # counted density model

loo_es_n = list(loo(mn_es, reloo = T), loo(mn_es_2, reloo = T)) 
loo_es_n
loo_compare(loo_es_n) # equal

# correlations
d_es_comp_2 <- d_es_2 %>%
  mutate(transformed = predict(mt_es_2)[, 1],
         nonlinear = predict(mn_es_2)[, 1])
cor.test(d_es_comp_2$transformed, d_es_comp_2$log_bio.g)
cor.test(d_es_comp_2$nonlinear, d_es_comp_2$bio.g) # transformed better

# visual comparison with planted density models
d_plot_es_2 <- d_es_2 %>%
  select(counted_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g,) %>%
  mutate(model = "linear") %>%
  rbind(d_es_2 %>%
          select(counted_density, bio.g, treatment) %>%
          mutate(model = "nonlinear")) %>%
  mutate(density = "counted") %>%
  rename(background_density = counted_density) %>%
  rbind(d_plot_es %>%
          mutate(density = "planted"))

d_pred_es_2 <- d_pred_fun(16, 100, "Ev seedling") %>%
  mutate(pred = fitted(mt_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,1],
         lower = fitted(mt_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,3],
         upper = fitted(mt_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,4],
         model = "linear") %>%
  full_join(d_pred_fun(16, 100, "Ev seedling") %>%
              mutate(pred = fitted(mn_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,1],
                     lower = fitted(mn_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,3],
                     upper = fitted(mn_es, newdata = d_pred_fun(16, 100, "Ev seedling"), re_formula = NA, nsamples = 100)[,4],
                     model = "nonlinear")) %>%
  mutate(density = "counted") %>%
  full_join(d_pred_es %>%
              mutate(density = "planted"))

pdf("./output/mv-biomass-by-treatment-2018-ev-seedling-full-models-counted.pdf", height = 7, width = 7)
d_pred_es_2 %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  geom_point(data = d_plot_es_2, aes(y = bio.g), size = 2, shape = 21, color = "black") +
  facet_grid(model~density, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.4, 0.4),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling density") +
  ylab("Microstegium biomass (g)")
dev.off()


#### simplify models ####

# Ev adult, nonlinear sub-models
mt_ea_c <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + cc_adj + (1|site), control = list(adapt_delta = 0.9999))
mt_ea_s <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + sm_adj + (1|site), control = list(adapt_delta = 0.99999))
mt_ea_p <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + pm_adj + (1|site), control = list(adapt_delta = 0.99999))
mt_ea_cs <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + sm_adj + cc_adj + (1|site), control = list(adapt_delta = 0.999999))
mt_ea_sp <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + sm_adj + pm_adj + (1|site), control = list(adapt_delta = 0.9999999))
mt_ea_pc <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + cc_adj + pm_adj + (1|site), control = list(adapt_delta = 0.9999))
mt_ea_n <- update(mt_ea, formula. = log_bio.g ~ background_density * treatment + (1|site), control = list(adapt_delta = 0.9999))

# model comparison
lc_ea = list(loo(mt_ea, reloo = T), loo(mt_ea_c, reloo = T), loo(mt_ea_s, reloo = T), loo(mt_ea_p, reloo = T), loo(mt_ea_cs, reloo = T), loo(mt_ea_sp, reloo = T), loo(mt_ea_pc, reloo = T), loo(mt_ea_n, reloo = T))
lc_ea # check that k values are low and that p_loo values are near number of parameters (7 in full)
loo_compare(lc_ea) # they're all pretty similar, with plot-edge mv included in the lowest models

# save model comparison
save(lc_ea, file = "./output/mv-biomass-by-treatment-2018-submodel-loo-list-ev-adult.rda")

# Ev seedling, nonlinear sub-models
mt_es_c <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + cc_adj + (1|site), newdata = d_es_2, control = list(adapt_delta = 0.9999))
mt_es_s <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + sm_adj + (1|site), newdata = d_es_2, control = list(adapt_delta = 0.99999))
mt_es_p <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + pm_adj + (1|site), newdata = d_es_2, control = list(adapt_delta = 0.9999))
mt_es_cs <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + sm_adj + cc_adj + (1|site), newdata = d_es_2, control = list(adapt_delta = 0.9999))
mt_es_sp <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + sm_adj + pm_adj + (1|site), newdata = d_es_2)
mt_es_pc <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + cc_adj + pm_adj + (1|site), newdata = d_es_2, control = list(adapt_delta = 0.9999))
mt_es_n <- update(mt_es_2, formula. = log_bio.g ~ background_density * treatment + (1|site), newdata = d_es_2, control = list(adapt_delta = 0.9999))

# model comparison
lc_es = list(loo(mt_es_2, reloo = T), loo(mt_es_c, reloo = T), loo(mt_es_s, reloo = T), loo(mt_es_p, reloo = T), loo(mt_es_cs, reloo = T), loo(mt_es_sp, reloo = T), loo(mt_es_pc, reloo = T), loo(mt_es_n, reloo = T))
lc_es # check that k values are low and that p_loo values are near number of parameters (7 in full)
loo_compare(lc_es) # all very similar - fewer covariates tended to be better

# save model comparison
save(lc_es, file = "./output/mv-biomass-by-treatment-2018-submodel-loo-list-ev-seedling.rda")

# Mv, nonlinear sub-models
mt_ms_c <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + cc_adj + (1|site), control = list(adapt_delta = 0.999999))
mt_ms_s <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + sm_adj + (1|site), control = list(adapt_delta = 0.999999))
mt_ms_p <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + pm_adj + (1|site))
mt_ms_cs <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + sm_adj + cc_adj + (1|site))
mt_ms_sp <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + sm_adj + pm_adj + (1|site))
mt_ms_pc <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + cc_adj + pm_adj + (1|site))
mt_ms_n <- update(mt_ms, formula. = log_bio.g ~ background_density * treatment + (1|site))

# model comparison
lc_ms = list(loo(mt_ms, reloo = T), loo(mt_ms_c, reloo = T), loo(mt_ms_s, reloo = T), loo(mt_ms_p, reloo = T), loo(mt_ms_cs, reloo = T), loo(mt_ms_sp, reloo = T), loo(mt_ms_pc, reloo = T), loo(mt_ms_n, reloo = T))
lc_ms # check that k values are low and that p_loo values are near number of parameters (7 in full)
loo_compare(lc_ms) # all very similar, full model and models with multiple seem useful

# save model comparison
save(lc_ms, file = "./output/mv-biomass-by-treatment-2018-submodel-loo-list-mv.rda")


#### final model figure ####

# Mv
p_m <- d_pred_ms %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  geom_point(data = filter(d_plot_ms, model == "linear"), aes(y = bio.g), size = 2, shape = 21, color = "black") +
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
  xlab("Microstegium density") +
  ylab("ln(Microstegium biomass)")

# Ev adult
p_ea <- d_pred_ea %>%
  filter(model == "linear") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  geom_point(data = filter(d_plot_ea, model == "linear"), aes(y = bio.g), size = 2, shape = 21, color = "black") +
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
  xlab("Adult Elymus density")

# Ev seedling
p_es <- d_pred_es_2 %>%
  filter(model == "linear" & density == "counted") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  geom_point(data = filter(d_plot_es_2, model == "linear" & density == "counted"), aes(y = bio.g), size = 2, shape = 21, color = "black") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.8, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling density") +
  xlim(0, 10)

# combine
pdf("./output/mv-biomass-by-treatment-2018-full-models.pdf", height = 5, width = 9)
plot_grid(p_m, p_es, p_ea, nrow = 1)
dev.off()


#### interaction coefficients ####

# posterior samples
post_mt_ea <- posterior_samples(mt_ea)
post_mt_es <- posterior_samples(mt_es_2)
post_mt_ms <- posterior_samples(mt_ms)

# rename columns
colnames(post_mt_ea) <- colnames(post_mt_es) <- colnames(post_mt_ms) <-  c("int", "density", "trt", "sm", "cc", "pm", "density_trt", "site", "sigma", "site_D1", "site_D2", "site_D3", "site_D4", "lp")

# posterior distributions
d_den_ea <- post_mt_ea %>%
  transmute(fungicide = density,
            water = density + density_trt) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(background = "Elymus adult") %>%
  as_tibble()

d_den_es <- post_mt_es %>%
  transmute(fungicide = density,
            water = density + density_trt) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(background = "Elymus seedling") %>%
  as_tibble()

d_den_ms <- post_mt_ms %>%
  transmute(fungicide = density,
            water = density + density_trt) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(background = "Microstegium") %>%
  as_tibble()

d_den <- d_den_ea %>%
  full_join(d_den_es) %>%
  full_join(d_den_ms) 

# visualize
pdf("./output/mv-biomass-by-treatment-2018-interaction-effects.pdf", width = 5, height = 5)
d_den %>%
  ggplot(aes(x = effect, y = background, color = treatment, fill = treatment)) +
  stat_density_ridges(rel_min_height = 0.001, scale = 1, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  scale_fill_manual(values = alpha(colpal, 0.7)) +
  scale_color_manual(values = colpal) +
  xlab("Group effect") +
  ylab("Interacting group") +
  ggtitle("Microstegium biomass") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.8, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.2)))
dev.off()


#### severity correlation ####

# correlation
cor.test(dat_s$severity, dat_s$bio.g) # weak

cor.test(filter(dat_s, treatment == "water")$severity, filter(dat_s, treatment == "water")$bio.g) # weak

# figure
pdf("./output/mv-biomass-by-treatment-2018-severity.pdf", width = 5, height = 5)
dat_s %>%
  ggplot(aes(x = severity, y = bio.g, colour = treatment)) + 
  geom_point() +
  scale_color_manual(values = colpal) +
  xlab("Proportion of leaf area with lesions") +
  ylab("Microstegium biomass (g)") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.8, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()