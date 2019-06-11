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

# run leaf scan files
source("./code/mv-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), "mleaf"))
source("./code/covariate-data-processing-2018.R")
rm(list = setdiff(ls(), c("mleaf", "covar")))

# import data
mb <- read_csv("./data/mv-biomass-oct-2018-density-exp.csv")
trt <- read_csv("./data/plot-treatments-for-figures-2018-density-exp.csv")
til <- read_csv("./data/mv-disease-sep-2018-density-exp.csv")

# non-linear species effects function
el_fun <- function(y0, x, a) {
  y = y0 * exp(a * log(x + 1))
  return(y)
}


#### edit data ####

# edit leaf scan data
ml <- mleaf %>%
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

# check notes
unique(mb$processing_notes)

# edit values and merge with other data
dat <- mb %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment),
         log_bio.g = log(bio.g)) %>%
  full_join(trt) %>%
  left_join(ml) %>%
  left_join(covar) %>%
  mutate(density_level = factor(density_level, levels = c("none", "low", "medium", "high")),
         disease = ifelse(treatment == "fungicide", 0, 1))


##### visualize #####

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

# severity
dat %>%
  ggplot(aes(x = severity, y = bio.g, colour = treatment)) + 
  geom_point()

dat %>%
  ggplot(aes(x = severity, y = log_bio.g, colour = treatment)) + 
  geom_point()


#### set up models ####

# check data
filter(dat, background == "Ev adult") %>%
  select(log_bio.g, background_density, treatment, sm_adj, cc_adj, pm_adj, site)

# subset
d_ea = filter(dat, background == "Ev adult")
d_es = filter(dat, background == "Ev seedling")
d_ms = filter(dat, background == "Mv seedling")

# full linear log-transformed models
mt_ea <- brm(data = d_ea, family = gaussian,
           log_bio.g ~ background_density * treatment + sm_adj + cc_adj + pm_adj + (1|site),
           prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(cauchy(0, 1), class = sd)),
           iter = 6000, warmup = 1000, chains = 3, cores = 2,
           control = list(adapt_delta = 0.999))
summary(mt_ea)
save(mt_ea, file = "./output/mv-biomass-by-treatment-2018-log-transformed-ev-adult.rda")

mt_es <- update(mt_ea, newdata = d_es,
                control = list(adapt_delta = 0.9999))
summary(mt_es)
save(mt_es, file = "./output/mv-biomass-by-treatment-2018-log-transformed-ev-seedling.rda")

mt_ms <- update(mt_ea, newdata = d_ms,
                control = list(adapt_delta = 0.9999))
summary(mt_ms)
save(mt_ms, file = "./output/mv-biomass-by-treatment-2018-log-transformed-mv.rda")

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
mn_ea <- brm(data = filter(dat, background == "Ev adult"), family = gaussian,
             bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + cc_adj + pm_adj + (1|site), a ~ treatment, nl = T),
             prior <- c(prior(normal(13.8, 10), nlpar = "y0"),
                        prior(normal(0, 0.5), nlpar = "a")),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
summary(mn_ea)
save(mn_ea, file = "./output/mv-biomass-by-treatment-2018-nonlinear-ev-adult.rda")

mn_es <- update(mn_ea, newdata = d_es)
summary(mn_es)
save(mn_es, file = "./output/mv-biomass-by-treatment-2018-nonlinear-ev-seedling.rda")

mn_ms <- update(mn_ea, newdata = d_ms)
summary(mn_ms)
save(mn_ms, file = "./output/mv-biomass-by-treatment-2018-nonlinear-mv.rda")


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
pp_check(mn_ea, nsamples = 50)

pp_check(mt_es, nsamples = 50)
pp_check(mn_es, nsamples = 50)

pp_check(mt_ms, nsamples = 50)
pp_check(mn_ms, nsamples = 50)

# data and model fits
d_pred <- tibble(background_density = rep(0:64, 4), 
                 treatment = rep(rep(c("water", "fungicide"), each = 65), 2),
                 sm_adj = mean(filter(dat, background == "Ev adult")$sm_adj),
                 cc_adj = mean(filter(dat, background == "Ev adult")$cc_adj),
                 pm_adj = mean(filter(dat, background == "Ev adult")$pm_adj),
                 model = rep(c("transformed", "nonlinear"), each = 65*2))

# Ev adult
d_plot_ea <- d_ea %>%
  select(background_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g) %>%
  mutate(model = "transformed") %>%
  rbind(d_ea %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "nonlinear"))
  
d_pred_ea <- d_pred %>%
  mutate(pred = c(fitted(mt_ea, newdata = d_pred[1:130,], re_formula = NA)[,1], fitted(mn_ea, newdata = d_pred[131:260,], re_formula = NA)[,1]),
         lower = c(fitted(mt_ea, newdata = d_pred[1:130,], re_formula = NA)[,3], fitted(mn_ea, newdata = d_pred[131:260,], re_formula = NA)[,3]),
         upper = c(fitted(mt_ea, newdata = d_pred[1:130,], re_formula = NA)[,4], fitted(mn_ea, newdata = d_pred[131:260,], re_formula = NA)[,4])) %>%
  full_join(d_plot_ea)

d_pred_ea %>%
  filter(background_density <= 8) %>%
  ggplot(aes(x = background_density, y = bio.g, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, colour = NA)+
  geom_smooth(aes(y = pred), size = 1.5) +
  facet_wrap(~model, scales = "free") +
  theme_bw()

# Ev seedling
d_plot_es <- d_es %>%
  select(background_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g) %>%
  mutate(model = "transformed") %>%
  rbind(d_es %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "nonlinear"))

d_pred_es <- d_pred %>%
  mutate(pred = c(fitted(mt_es, newdata = d_pred[1:130,], re_formula = NA)[,1], fitted(mn_es, newdata = d_pred[131:260,], re_formula = NA)[,1]),
         lower = c(fitted(mt_es, newdata = d_pred[1:130,], re_formula = NA)[,3], fitted(mn_es, newdata = d_pred[131:260,], re_formula = NA)[,3]),
         upper = c(fitted(mt_es, newdata = d_pred[1:130,], re_formula = NA)[,4], fitted(mn_es, newdata = d_pred[131:260,], re_formula = NA)[,4])) %>%
  full_join(d_plot_es)

d_pred_es %>%
  filter(background_density <= 16) %>%
  ggplot(aes(x = background_density, y = bio.g, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, colour = NA)+
  geom_smooth(aes(y = pred), size = 1.5) +
  facet_wrap(~model, scales = "free") +
  theme_bw()

# Mv seedling
d_plot_ms <- d_ms %>%
  select(background_density, log_bio.g, treatment) %>%
  rename(bio.g = log_bio.g) %>%
  mutate(model = "transformed") %>%
  rbind(d_ms %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "nonlinear"))

d_pred_ms <- d_pred %>%
  mutate(pred = c(fitted(mt_ms, newdata = d_pred[1:130,], re_formula = NA)[,1], fitted(mn_ms, newdata = d_pred[131:260,], re_formula = NA)[,1]),
         lower = c(fitted(mt_ms, newdata = d_pred[1:130,], re_formula = NA)[,3], fitted(mn_ms, newdata = d_pred[131:260,], re_formula = NA)[,3]),
         upper = c(fitted(mt_ms, newdata = d_pred[1:130,], re_formula = NA)[,4], fitted(mn_ms, newdata = d_pred[131:260,], re_formula = NA)[,4])) %>%
  full_join(d_plot_ms)

d_pred_ms %>%
  ggplot(aes(x = background_density, y = bio.g, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, colour = NA)+
  geom_smooth(aes(y = pred), size = 1.5) +
  facet_wrap(~model, scales = "free") +
  theme_bw()


#### simplify models ####

# Ev adult, nonlinear sub-models
mn_ea_c <- update(mn_ea, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + cc_adj + (1|site), a ~ treatment, nl = T))
mn_ea_s <- update(mn_ea, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + (1|site), a ~ treatment, nl = T))
mn_ea_p <- update(mn_ea, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + pm_adj + (1|site), a ~ treatment, nl = T))
mn_ea_cs <- update(mn_ea, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + cc_adj + sm_adj + (1|site), a ~ treatment, nl = T))
mn_ea_sp <- update(mn_ea, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + pm_adj + (1|site), a ~ treatment, nl = T))
mn_ea_pc <- update(mn_ea, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + pm_adj + cc_adj + (1|site), a ~ treatment, nl = T))

# model comparison
lc_ea = list(loo(mn_ea, reloo = T), loo(mn_ea_c, reloo = T), loo(mn_ea_s, reloo = T), loo(mn_ea_p, reloo = T), loo(mn_ea_cs, reloo = T), loo(mn_ea_sp, reloo = T), loo(mn_ea_pc, reloo = T))
lc_ea # check that k values are low and that p_loo values are near number of parameters (7 in full)
loo_model_weights(lc_ea) # split between full model and model with pm_adj. Adding cc and sm improve the model fit a little when added together, but not when added individually

# Ev seedling, nonlinear sub-models
mn_es_c <- update(mn_es, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + cc_adj + (1|site), a ~ treatment, nl = T))
mn_es_s <- update(mn_es, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + (1|site), a ~ treatment, nl = T))
mn_es_p <- update(mn_es, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + pm_adj + (1|site), a ~ treatment, nl = T))
mn_es_cs <- update(mn_es, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + cc_adj + sm_adj + (1|site), a ~ treatment, nl = T))
mn_es_sp <- update(mn_es, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + pm_adj + (1|site), a ~ treatment, nl = T))
mn_es_pc <- update(mn_es, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + pm_adj + cc_adj + (1|site), a ~ treatment, nl = T))

# model comparison
lc_es = list(loo(mn_es, reloo = T), loo(mn_es_c, reloo = T), loo(mn_es_s, reloo = T), loo(mn_es_p, reloo = T), loo(mn_es_cs, reloo = T), loo(mn_es_sp, reloo = T), loo(mn_es_pc, reloo = T))
lc_es # check that k values are low and that p_loo values are near number of parameters (7 in full)
loo_model_weights(lc_es) # the canopy cover model is the highest, followed by the soil one

# save best model
save(mn_es_c, file = "./output/mv-biomass-by-treatment-2018-nonlinear-canopy-ev-seedling.rda")

# Mv, nonlinear sub-models
mn_ms_c <- update(mn_ms, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + cc_adj + (1|site), a ~ treatment, nl = T))
mn_ms_s <- update(mn_ms, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + (1|site), a ~ treatment, nl = T))
mn_ms_p <- update(mn_ms, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + pm_adj + (1|site), a ~ treatment, nl = T))
mn_ms_cs <- update(mn_ms, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + cc_adj + sm_adj + (1|site), a ~ treatment, nl = T))
mn_ms_sp <- update(mn_ms, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + sm_adj + pm_adj + (1|site), a ~ treatment, nl = T))
mn_ms_pc <- update(mn_ms, formula. = bf(bio.g ~ y0 * exp(a * log(background_density + 1)), y0 ~ treatment + pm_adj + cc_adj + (1|site), a ~ treatment, nl = T))

# model comparison
lc_ms = list(loo(mn_ms, reloo = T), loo(mn_ms_c, reloo = T), loo(mn_ms_s, reloo = T), loo(mn_ms_p, reloo = T), loo(mn_ms_cs, reloo = T), loo(mn_ms_sp, reloo = T), loo(mn_ms_pc, reloo = T))
lc_ms # check that k values are low and that p_loo values are near number of parameters (7 in full)
loo_model_weights(lc_ms) # the soil and plot-edge mv model is the best, followed by soil only

# save best model
save(mn_ms_sp, file = "./output/mv-biomass-by-treatment-2018-nonlinear-soil-plot-edge-mv.rda")

# ppc's of sub-models
pp_check(mn_ea, nsamples = 50)
pp_check(mn_ea_s, nsamples = 50)

pp_check(mn_es_c, nsamples = 50)
pp_check(mn_es_s, nsamples = 50)

pp_check(mn_ms_sp, nsamples = 50)
pp_check(mn_ms_s, nsamples = 50)

# data and model fits

# Ev adult
d_plot_ea_2 <- d_ea %>%
  select(background_density, bio.g, treatment) %>%
  mutate(model = "full") %>%
  rbind(d_ea %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "soil"))

d_pred_ea_2 <- d_pred %>%
  mutate(model = case_when(model == "transformed" ~ "full",
                           model == "nonlinear" ~ "soil"),
    pred = c(fitted(mn_ea, newdata = d_pred[1:130,], re_formula = NA)[,1], fitted(mn_ea_s, newdata = d_pred[131:260,], re_formula = NA)[,1]),
         lower = c(fitted(mn_ea, newdata = d_pred[1:130,], re_formula = NA)[,3], fitted(mn_ea_s, newdata = d_pred[131:260,], re_formula = NA)[,3]),
         upper = c(fitted(mn_ea, newdata = d_pred[1:130,], re_formula = NA)[,4], fitted(mn_ea_s, newdata = d_pred[131:260,], re_formula = NA)[,4])) %>%
  full_join(d_plot_ea_2)

d_pred_ea_2 %>%
  filter(background_density <= 8) %>%
  ggplot(aes(x = background_density, y = bio.g, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, colour = NA)+
  geom_smooth(aes(y = pred), size = 1.5) +
  facet_wrap(~model, scales = "free") +
  theme_bw()

# Ev seedling
d_plot_es_2 <- d_es %>%
  select(background_density, bio.g, treatment) %>%
  mutate(model = "canopy") %>%
  rbind(d_es %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "soil"))

d_pred_es_2 <- d_pred %>%
  mutate(model = case_when(model == "transformed" ~ "canopy",
                           model == "nonlinear" ~ "soil"),
         pred = c(fitted(mn_es_c, newdata = d_pred[1:130,], re_formula = NA)[,1], fitted(mn_es_s, newdata = d_pred[131:260,], re_formula = NA)[,1]),
         lower = c(fitted(mn_es_c, newdata = d_pred[1:130,], re_formula = NA)[,3], fitted(mn_es_s, newdata = d_pred[131:260,], re_formula = NA)[,3]),
         upper = c(fitted(mn_es_c, newdata = d_pred[1:130,], re_formula = NA)[,4], fitted(mn_es_s, newdata = d_pred[131:260,], re_formula = NA)[,4])) %>%
  full_join(d_plot_es_2)

d_pred_es_2 %>%
  filter(background_density <= 16) %>%
  ggplot(aes(x = background_density, y = bio.g, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, colour = NA)+
  geom_smooth(aes(y = pred), size = 1.5) +
  facet_wrap(~model, scales = "free") +
  theme_bw()

# Mv seedling
d_plot_ms_2 <- d_ms %>%
  select(background_density, bio.g, treatment) %>%
  mutate(model = "soil+mv") %>%
  rbind(d_ms %>%
          select(background_density, bio.g, treatment) %>%
          mutate(model = "soil"))

d_pred_ms_2 <- d_pred %>%
  mutate(model = case_when(model == "transformed" ~ "soil+mv",
                           model == "nonlinear" ~ "soil"),
         pred = c(fitted(mn_ms_sp, newdata = d_pred[1:130,], re_formula = NA)[,1], fitted(mn_ms_s, newdata = d_pred[131:260,], re_formula = NA)[,1]),
         lower = c(fitted(mn_ms_sp, newdata = d_pred[1:130,], re_formula = NA)[,3], fitted(mn_ms_s, newdata = d_pred[131:260,], re_formula = NA)[,3]),
         upper = c(fitted(mn_ms_sp, newdata = d_pred[1:130,], re_formula = NA)[,4], fitted(mn_ms_s, newdata = d_pred[131:260,], re_formula = NA)[,4])) %>%
  full_join(d_plot_ms_2)

d_pred_ms_2 %>%
  ggplot(aes(x = background_density, y = bio.g, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, colour = NA)+
  geom_smooth(aes(y = pred), size = 1.5) +
  facet_wrap(~model, scales = "free") +
  theme_bw()


#### interaction coefficients ####

# posterior samples
post_mn_ea <- posterior_samples(mn_ea)
post_mn_es <- posterior_samples(mn_es_c)
post_mn_ms <- posterior_samples(mn_ms_sp)

# rename columns
colnames(post_mn_ea) <-  c("y0_int", "y0_trt", "sm", "cc", "pm", "a_int", "a_trt", "site", "sigma", "site_D1", "site_D2", "site_D3", "site_D4", "lp")
colnames(post_mn_es) <- c("y0_int", "y0_trt", "cc", "a_int", "a_trt", "site", "sigma", "site_D1", "site_D2", "site_D3", "site_D4", "lp")
colnames(post_mn_ms) <- c("y0_int", "y0_trt", "sm", "pm", "a_int", "a_trt", "site", "sigma", "site_D1", "site_D2", "site_D3", "site_D4", "lp")

# posterior distributions
d_trt_ea <- post_mn_ea %>%
  transmute(fungicide = a_int,
            water = a_int + a_trt) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(background = "Elymus adult") %>%
  as_tibble()

d_trt_es <- post_mn_es %>%
  transmute(fungicide = a_int,
            water = a_int + a_trt) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(background = "Elymus seedling") %>%
  as_tibble()

d_trt_ms <- post_mn_ms %>%
  transmute(fungicide = a_int,
            water = a_int + a_trt) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(background = "Microstegium") %>%
  as_tibble()

d_trt <- d_trt_ea %>%
  full_join(d_trt_es) %>%
  full_join(d_trt_ms) 

# visualize
pdf("./output/mv-biomass-by-treatment-2018-interaction-effects.pdf")
d_trt %>%
  ggplot(aes(x = effect, y = background, color = treatment, fill = treatment)) +
  geom_halfeyeh() +
  scale_fill_manual(values = alpha(c("thistle1", "lightsteelblue1"), 0.7)) +
  scale_color_manual(values = c("thistle3", "lightsteelblue3")) +
  xlab("Group effect") +
  ylab("Interacting group") +
  ggtitle("Microstegium biomass") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()