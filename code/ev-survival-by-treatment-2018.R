##### info ####

# file: ev-survival-by-treatment-2018
# author: Amy Kendig
# date last edited: 6/24/19
# goal: evaluate treatment effects on Ev survival


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)
library(cowplot)

# run data files and clear processing data
source("./code/bg-densities-data-processing-2018.R")
rm(list = setdiff(ls(), c("esurv", "bgd")))
source("./code/covariate-data-processing-2018.R")
rm(list = setdiff(ls(), c("esurv", "bgd", "covar")))

# import data
trt <- read_csv("./data/plot-treatments-for-figures-2018-density-exp.csv")
trt_s <- read_csv("./data/plot-treatments-2018-density-exp.csv")

# function to create datasets for prediction
d_pred_fun = function(max_val, length_out, dat){
  
  d_pred <- tibble(background_density = rep(seq(0, max_val, length.out = length_out), 4), 
                   counted_density = rep(seq(0, max_val, length.out = length_out), 4),
                   treatment = rep(rep(c("water", "fungicide"), each = length_out), 2),
                   age = rep(c("seedling", "adult"), each = length_out * 2),
                   sm_adj = mean(dat$sm_adj),
                   cc_adj = mean(dat$cc_adj),
                   pm_adj = mean(dat$pm_adj))
  
  return(d_pred)
}

# plot parameters
colpal = c("#3CBB75FF", "#39568CFF")
sm_txt = 12
lg_txt = 14
an_txt = 3


#### edit data ####

# merge with treatment, covariates, and counted densities
d <- left_join(esurv, trt) %>%
  left_join(covar) %>%
  left_join(select(bgd, -month))

# select September data (summer survival, remove litter experiment with the treatment column)
ds <- d %>%
  filter(month == "September" & !is.na(survival) & !is.na(treatment) & focal == 1)

# check that it makes sense
ds %>%
  group_by(site, plot, treatment, age) %>%
  summarise(n = length(survival)) %>%
  data.frame()

# select April data (winter survival given summer survival)
da <- d %>%
  filter(month %in% c("April", "September") & !is.na(treatment) & focal == 1) %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April)

# check that it makes sense
da %>%
  group_by(site, plot, treatment, age) %>%
  summarise(n = length(survival)) %>%
  data.frame()

# data subsets
dsa <- filter(ds, background == "Ev adult")
dss <- filter(ds, background == "Ev seedling")
dsm <- filter(ds, background == "Mv seedling")

daa <- filter(da, background == "Ev adult")
das <- filter(da, background == "Ev seedling")
dam <- filter(da, background == "Mv seedling")


#### visualize data ####

# proportion September
ds %>%
  group_by(age, treatment, background, background_density) %>%
  summarise(p = sum(survival) / length(survival)) %>%
  ggplot(aes(x = background_density, y = p, colour = treatment)) +
  geom_point(position = position_dodge(0.1)) +
  facet_grid(age~background, scales = "free") +
  theme_bw()

# proportion September, counted densities
ds %>%
  group_by(age, treatment, background, counted_density) %>%
  summarise(p = sum(survival) / length(survival)) %>%
  ggplot(aes(x = counted_density, y = p, colour = treatment)) +
  geom_point(position = position_dodge(0.1)) +
  facet_grid(age~background, scales = "free") +
  theme_bw()

# proportion April
da %>%
  group_by(age, treatment, background, background_density) %>%
  summarise(p = sum(survival) / length(survival)) %>%
  ggplot(aes(x = background_density, y = p, colour = treatment)) +
  geom_point(position = position_dodge(0.1)) +
  facet_grid(age~background, scales = "free") +
  theme_bw()

# proportion April, counted density
da %>%
  group_by(age, treatment, background, counted_density) %>%
  summarise(p = sum(survival) / length(survival)) %>%
  ggplot(aes(x = counted_density, y = p, colour = treatment)) +
  geom_point(position = position_dodge(0.1)) +
  facet_grid(age~background, scales = "free") +
  theme_bw()


#### statistical models ####

# adult plots, summer
msa <- brm(data = dsa, family = bernoulli,
             survival ~ background_density * treatment * age + sm_adj + cc_adj + pm_adj + (1|site),
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b),
                       prior(cauchy(0, 1), class = sd)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(msa)
save(msa, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer.rda")

# seedling plots, summer
mss <- update(msa, newdata = dss, control = list(max_treedepth = 15))
summary(mss)
save(mss, file = "./output/ev-survival-by-treatment-2018-ev-seedling-summer.rda")

# microstegium plots, summer
msm <- update(msa, newdata = dsm, control = list(max_treedepth = 15))
summary(msm)
save(msm, file = "./output/ev-survival-by-treatment-2018-mv-summer.rda")

# seedling plots, summer, counted
mss_d <- update(mss, formula = survival ~ counted_density * treatment * age + sm_adj + cc_adj + pm_adj + (1|site), newdata = dss)
summary(mss_d)
save(mss_d, file = "./output/ev-survival-by-treatment-2018-ev-seedling-summer-counted.rda")

# adult plots, winter
maa <- update(msa, newdata = daa)
summary(maa)
save(maa, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter.rda")

# seedling plots, winter
mas <- update(maa, newdata = das)
summary(mas)
save(mas, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter.rda")

# microstegium plots, winter
mam <- update(maa, newdata = dam, control = list(adapt_delta = 0.9999999, max_treedepth = 15))
# can't converge (too many 1's?)
sum(dam$survival)
nrow(dam) # only 2 0's

# seedling plots, winter, counted
mas_d <- update(mas, formula = survival ~ counted_density * treatment * age + sm_adj + cc_adj + pm_adj + (1|site), newdata = das)
summary(mas_d)
save(mas_d, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-counted.rda")


#### check models ####

# convergence of chains
plot(msa)
plot(mss)
plot(msm)
plot(mss_d)

plot(maa)
plot(mas)
plot(mas_d)

# posterior predictive check
pp_check(msa, nsamples = 50)
pp_check(mss, nsamples = 50)
pp_check(mss_d, nsamples = 50) # less spread
pp_check(msm, nsamples = 50)

pp_check(maa, nsamples = 50)
pp_check(mas, nsamples = 50)
pp_check(mas_d, nsamples = 50) # more spread

# loo comparison
loo_mss = list(loo(mss, reloo = T), loo(mss_d, reloo = T)) # very slow, but too many high k values when reloo = F
save(loo_mss, file = "./output/ev-survival-by-treatment-2018-density-loo-list-ev-seedling-summer.rda")
loo_mss 
loo::loo_compare(loo_mss) # pretty equal

loo_mas = list(loo(mas, reloo = T), loo(mas_d, reloo = T))
save(loo_mas, file = "./output/ev-survival-by-treatment-2018-density-loo-list-ev-seedling-winter.rda")
loo_mas
loo::loo_compare(loo_mas) # pretty equal

# correlations
dss_comp <- dss %>%
  mutate(planted = predict(mss)[, 1],
         counted = predict(mss_d)[, 1])
cor.test(dss_comp$planted, dss_comp$survival)
cor.test(dss_comp$counted, dss_comp$survival) # pretty much equal
cor.test(dss_comp$planted, dss_comp$counted)

das_comp <- das %>%
  mutate(planted = predict(mas)[, 1],
         counted = predict(mas_d)[, 1])
cor.test(das_comp$planted, das_comp$survival)
cor.test(das_comp$counted, das_comp$survival) # counted a little lower, but not by much
cor.test(das_comp$planted, das_comp$counted)


#### visualize models ####

# adult plots, summer
dsa_pred <- d_pred_fun(8, 100, dsa) %>%
  cbind(fitted(msa, newdata = d_pred_fun(8, 100, dsa), re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5) %>%
  as_tibble()

pdf("./output/ev-survival-by-treatment-2018-ev-adult-summer-full-models.pdf", height = 5, width = 7)
dsa_pred %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = dsa, aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = dsa, aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_wrap(~age, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.3, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Adult Elymus background density") +
  ylab("Focal Elymus summer survival")
dev.off()

# seedling plots, summer
dss_pred <- d_pred_fun(16, 100, dss) %>%
  cbind(fitted(mss, newdata = d_pred_fun(16, 100, dss), re_formula = NA, nsamples = 100)) %>%
  mutate(density = "planted") %>%
  select(-counted_density) %>%
  full_join(
    d_pred_fun(10, 100, dss) %>%
      cbind(fitted(mss_d, newdata = d_pred_fun(10, 100, dss), re_formula = NA, nsamples = 100)) %>%
      mutate(density = "counted") %>%
      select(-background_density) %>%
      rename(background_density = counted_density)
  ) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5) %>%
  as_tibble()

dss_plot <- dss %>%
  select(background_density, treatment, age, survival) %>%
  mutate(density = "planted") %>%
  full_join(dss %>%
              select(counted_density, treatment, age, survival) %>%
              rename(background_density = counted_density) %>%
              mutate(density = "counted"))

pdf("./output/ev-survival-by-treatment-2018-ev-seedling-summer-full-models.pdf", height = 7, width = 7)
dss_pred %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = dss_plot, aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = dss_plot, aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_grid(age ~ density, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.4, 0.6),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling background density") +
  ylab("Focal Elymus summer survival")
dev.off()

# microstegium plots, summer
dsm_pred <- d_pred_fun(64, 100, dsm) %>%
  cbind(fitted(msm, newdata = d_pred_fun(64, 100, dsa), re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5) %>%
  as_tibble()

pdf("./output/ev-survival-by-treatment-2018-mv-summer-full-models.pdf", height = 5, width = 7)
dsm_pred %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = dsm, aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = dsm, aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_wrap(~age, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.2, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Microstegium background density") +
  ylab("Focal Elymus summer survival")
dev.off()

# adult plots, winter
daa_pred <- d_pred_fun(8, 100, daa) %>%
  cbind(fitted(maa, newdata = d_pred_fun(8, 100, daa), re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5) %>%
  as_tibble() 

pdf("./output/ev-survival-by-treatment-2018-ev-adult-winter-full-models.pdf", height = 5, width = 7)
daa_pred %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = daa, aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = daa, aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_wrap(~age, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.3, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Adult Elymus background density") +
  ylab("Focal Elymus winter survival")
dev.off()

# seedling plots, winter
das_pred <- d_pred_fun(16, 100, dss) %>%
  cbind(fitted(mas, newdata = d_pred_fun(16, 100, das), re_formula = NA, nsamples = 100)) %>%
  mutate(density = "planted") %>%
  select(-counted_density) %>%
  full_join(
    d_pred_fun(10, 100, das) %>%
      cbind(fitted(mas_d, newdata = d_pred_fun(10, 100, das), re_formula = NA, nsamples = 100)) %>%
      mutate(density = "counted") %>%
      select(-background_density) %>%
      rename(background_density = counted_density)
  ) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5) %>%
  as_tibble() 

das_plot <- das %>%
  select(background_density, treatment, age, survival) %>%
  mutate(density = "planted") %>%
  full_join(das %>%
              select(counted_density, treatment, age, survival) %>%
              rename(background_density = counted_density) %>%
              mutate(density = "counted"))

pdf("./output/ev-survival-by-treatment-2018-ev-seedling-winter-full-models.pdf", height = 7, width = 7)
das_pred %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = das_plot, aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = das_plot, aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  facet_grid(age ~ density, scales = "free") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.3, 0.6),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling background density") +
  ylab("Focal Elymus winter survival")
dev.off()


#### simplify models ####

# adult plots, summer
msa_n <- update(msa, formula. = survival ~ background_density * treatment * age + (1|site))
msa_s <- update(msa, formula. = survival ~ background_density * treatment * age + sm_adj + (1|site))
msa_c <- update(msa, formula. = survival ~ background_density * treatment * age + cc_adj + (1|site))
msa_p <- update(msa, formula. = survival ~ background_density * treatment * age + pm_adj + (1|site))
msa_sc <- update(msa, formula. = survival ~ background_density * treatment * age + sm_adj + cc_adj + (1|site))
msa_cp <- update(msa, formula. = survival ~ background_density * treatment * age + cc_adj + pm_adj + (1|site))
msa_sp <- update(msa, formula. = survival ~ background_density * treatment * age + sm_adj + pm_adj + (1|site))

save(msa_n, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-n.rda")
save(msa_s, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-s.rda")
save(msa_c, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-c.rda")
save(msa_p, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-p.rda")
save(msa_sc, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-sc.rda")
save(msa_cp, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-cp.rda")
save(msa_sp, file = "./output/ev-survival-by-treatment-2018-ev-adult-summer-sp.rda")

load("./output/ev-survival-by-treatment-2018-ev-adult-summer-n.rda")
load("./output/ev-survival-by-treatment-2018-ev-adult-summer-s.rda")
load("./output/ev-survival-by-treatment-2018-ev-adult-summer-c.rda")
load("./output/ev-survival-by-treatment-2018-ev-adult-summer-p.rda")
load("./output/ev-survival-by-treatment-2018-ev-adult-summer-sc.rda")
load("./output/ev-survival-by-treatment-2018-ev-adult-summer-cp.rda")
load("./output/ev-survival-by-treatment-2018-ev-adult-summer-sp.rda")

loo_msa <- list(loo(msa, reloo = T), loo(msa_n, reloo = T), loo(msa_s, reloo = T), loo(msa_c, reloo = T), loo(msa_p, reloo = T), loo(msa_sc, reloo = T), loo(msa_cp, reloo = T), loo(msa_sp, reloo = T))
save(loo_msa, file = "./output/ev-survival-by-treatment-2018-submodel-loo-list-ev-adult-summer.rda")
loo::loo_compare(loo_msa)
# no covariates is the best, all are pretty similar

# seedling plots, summer (counted)
mss_n <- update(mss_d, formula. = survival ~ counted_density * treatment * age + (1|site))
mss_s <- update(mss_d, formula. = survival ~ counted_density * treatment * age + sm_adj + (1|site))
mss_c <- update(mss_d, formula. = survival ~ counted_density * treatment * age + cc_adj + (1|site))
mss_p <- update(mss_d, formula. = survival ~ counted_density * treatment * age + pm_adj + (1|site))
mss_sc <- update(mss_d, formula. = survival ~ counted_density * treatment * age + sm_adj + cc_adj + (1|site))
mss_cp <- update(mss_d, formula. = survival ~ counted_density * treatment * age + cc_adj + pm_adj + (1|site))
mss_sp <- update(mss_d, formula. = survival ~ counted_density * treatment * age + sm_adj + pm_adj + (1|site))

save(mss_n, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-n.rda")
save(mss_s, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-s.rda")
save(mss_c, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-c.rda")
save(mss_p, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-p.rda")
save(mss_sc, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-sc.rda")
save(mss_cp, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-cp.rda")
save(mss_sp, "./output/ev-survival-by-treatment-2018-ev-seedling-summer-sp.rda")

load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-n.rda")
load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-s.rda")
load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-c.rda")
load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-p.rda")
load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-sc.rda")
load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-cp.rda")
load("./output/ev-survival-by-treatment-2018-ev-seedling-summer-sp.rda")

loo_mss_2 <- list(loo(mss_n, reloo = T), loo(mss_s, reloo = T), loo(mss_c, reloo = T), loo(mss_p, reloo = T), loo(mss_sc, reloo = T), loo(mss_cp, reloo = T), loo(mss_sp, reloo = T))
save(loo_mss_2, file = "./output/ev-survival-by-treatment-2018-submodel-loo-list-ev-seedling-summer.rda")
loo_mss_2[[8]] <- loo_mss[[2]]
loo::loo_compare(loo_mss_2)
# no covariates is the best, all are pretty similar

# microstegium plots, summer
msm_n <- update(msm, formula. = survival ~ background_density * treatment * age + (1|site))
msm_s <- update(msm, formula. = survival ~ background_density * treatment * age + sm_adj + (1|site))
msm_c <- update(msm, formula. = survival ~ background_density * treatment * age + cc_adj + (1|site))
msm_p <- update(msm, formula. = survival ~ background_density * treatment * age + pm_adj + (1|site))
msm_sc <- update(msm, formula. = survival ~ background_density * treatment * age + sm_adj + cc_adj + (1|site))
msm_cp <- update(msm, formula. = survival ~ background_density * treatment * age + cc_adj + pm_adj + (1|site))
msm_sp <- update(msm, formula. = survival ~ background_density * treatment * age + sm_adj + pm_adj + (1|site))

save(msm_n, "./output/ev-survival-by-treatment-2018-mv-summer-n.rda")
save(msm_s, "./output/ev-survival-by-treatment-2018-mv-summer-s.rda")
save(msm_c, "./output/ev-survival-by-treatment-2018-mv-summer-c.rda")
save(msm_p, "./output/ev-survival-by-treatment-2018-mv-summer-p.rda")
save(msm_sc, "./output/ev-survival-by-treatment-2018-mv-summer-sc.rda")
save(msm_cp, "./output/ev-survival-by-treatment-2018-mv-summer-cp.rda")
save(msm_sp, "./output/ev-survival-by-treatment-2018-mv-summer-sp.rda")

load("./output/ev-survival-by-treatment-2018-mv-summer-n.rda")
load("./output/ev-survival-by-treatment-2018-mv-summer-s.rda")
load("./output/ev-survival-by-treatment-2018-mv-summer-c.rda")
load("./output/ev-survival-by-treatment-2018-mv-summer-p.rda")
load("./output/ev-survival-by-treatment-2018-mv-summer-sc.rda")
load("./output/ev-survival-by-treatment-2018-mv-summer-cp.rda")
load("./output/ev-survival-by-treatment-2018-mv-summer-sp.rda")

loo_msm <- list(loo(msm, reloo = T), loo(msm_n, reloo = T), loo(msm_s, reloo = T), loo(msm_c, reloo = T), loo(msm_p, reloo = T), loo(msm_sc, reloo = T), loo(msm_cp, reloo = T), loo(msm_sp, reloo = T))
save(loo_msm, file = "./output/ev-survival-by-treatment-2018-submodel-loo-list-mv-summer.rda")
loo::loo_compare(loo_msm)
# no covariates is the best, all are pretty similar

# adult plots, winter
maa_n <- update(maa, formula. = survival ~ background_density * treatment * age + (1|site))
maa_s <- update(maa, formula. = survival ~ background_density * treatment * age + sm_adj + (1|site))
maa_c <- update(maa, formula. = survival ~ background_density * treatment * age + cc_adj + (1|site))
maa_p <- update(maa, formula. = survival ~ background_density * treatment * age + pm_adj + (1|site))
maa_sc <- update(maa, formula. = survival ~ background_density * treatment * age + sm_adj + cc_adj + (1|site))
maa_cp <- update(maa, formula. = survival ~ background_density * treatment * age + cc_adj + pm_adj + (1|site))
maa_sp <- update(maa, formula. = survival ~ background_density * treatment * age + sm_adj + pm_adj + (1|site))

save(maa_n, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-n.rda")
save(maa_s, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-s.rda")
save(maa_c, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-c.rda")
save(maa_p, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-p.rda")
save(maa_sc, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-sc.rda")
save(maa_cp, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-cp.rda")
save(maa_sp, file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-sp.rda")

load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-n.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-s.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-c.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-p.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-sc.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-cp.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-adult-winter-sp.rda")

loo_maa <- list(loo(maa, reloo = T), loo(maa_n, reloo = T), loo(maa_s, reloo = T), loo(maa_c, reloo = T), loo(maa_p, reloo = T), loo(maa_sc, reloo = T), loo(maa_cp, reloo = T), loo(maa_sp, reloo = T))
save(loo_maa, file = "./output/ev-survival-by-treatment-2018-submodel-loo-list-ev-adult-winter.rda")
loo::loo_compare(loo_maa)
# no covariates is the best, all are pretty similar

# seedling plots, winter (counted)
mas_n <- update(mas_d, formula. = survival ~ counted_density * treatment * age + (1|site))
mas_s <- update(mas_d, formula. = survival ~ counted_density * treatment * age + sm_adj + (1|site))
mas_c <- update(mas_d, formula. = survival ~ counted_density * treatment * age + cc_adj + (1|site))
mas_p <- update(mas_d, formula. = survival ~ counted_density * treatment * age + pm_adj + (1|site))
mas_sc <- update(mas_d, formula. = survival ~ counted_density * treatment * age + sm_adj + cc_adj + (1|site))
mas_cp <- update(mas_d, formula. = survival ~ counted_density * treatment * age + cc_adj + pm_adj + (1|site))
mas_sp <- update(mas_d, formula. = survival ~ counted_density * treatment * age + sm_adj + pm_adj + (1|site))

save(mas_n, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-n.rda")
save(mas_s, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-s.rda")
save(mas_c, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-c.rda")
save(mas_p, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-p.rda")
save(mas_sc, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-sc.rda")
save(mas_cp, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-cp.rda")
save(mas_sp, file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-sp.rda")

load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-n.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-s.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-c.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-p.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-sc.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-cp.rda")
load(file = "./output/ev-survival-by-treatment-2018-ev-seedling-winter-sp.rda")

loo_mas_2 <- list(loo(mas_n, reloo = T), loo(mas_s, reloo = T), loo(mas_c, reloo = T), loo(mas_p, reloo = T), loo(mas_sc, reloo = T), loo(mas_cp, reloo = T), loo(mas_sp, reloo = T))
save(loo_mas_2, file = "./output/ev-survival-by-treatment-2018-submodel-loo-list-ev-seedling-winter.rda")
loo_mas_2[[8]] <- loo_mas[[2]]
loo::loo_compare(loo_mas_2)
# no covariates is the best, all are pretty similar


#### final model figure ####

# Mv, summer, seedling
psms <- dsm_pred %>%
  filter(age == "seedling") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(dsm, age == "seedling"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(1)) +
  stat_summary(data = filter(dsm, age == "seedling"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(1), shape = 21, color = "black") +
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
  ylab("Elymus seedling summer survival") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Mv, summer, adult
psma <- dsm_pred %>%
  filter(age == "adult") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(dsm, age == "adult"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(1)) +
  stat_summary(data = filter(dsm, age == "adult"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(1), shape = 21, color = "black") +
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
  ylab("Adult Elymus summer survival") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Ev seedling, summer, seedling
psss <- dss_pred %>%
  filter(age == "seedling" & density == "counted") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(dss_plot, age == "seedling" & density == "counted"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(dss_plot, age == "seedling" & density == "counted"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
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
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Ev seedling, summer, adult
pssa <- dss_pred %>%
  filter(age == "adult" & density == "counted") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(dss_plot, age == "adult" & density == "counted"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(dss_plot, age == "adult" & density == "counted"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.7, 0.3),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling density") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Ev adult, summer, seedling
psas <- dsa_pred %>%
  filter(age == "seedling") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(dsa, age == "seedling"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(dsa, age == "seedling"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
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
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Ev adult, summer, adult
psaa <- dsa_pred %>%
  filter(age == "adult") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(dsa, age == "adult"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(dsa, age == "adult"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
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
  xlab("Adult Elymus density") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Ev seedling, winter, seedling
pass <- das_pred %>%
  filter(age == "seedling" & density == "counted") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(das_plot, age == "seedling" & density == "counted"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(das_plot, age == "seedling" & density == "counted"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
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
  ylab("Elymus seedling winter survival") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Ev seedling, winter, adult
pasa <- das_pred %>%
  filter(age == "adult" & density == "counted") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(das_plot, age == "adult" & density == "counted"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(das_plot, age == "adult" & density == "counted"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.3, 0.3),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  xlab("Elymus seedling density") +
  ylab("Elymus adult winter survival") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

paas <- daa_pred %>%
  filter(age == "seedling") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(daa, age == "seedling"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(daa, age == "seedling"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
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
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

paaa <- daa_pred %>%
  filter(age == "adult") %>%
  ggplot(aes(x = background_density, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_smooth(aes(y = pred), size = 1.5) +
  stat_summary(data = filter(daa, age == "adult"), aes(y = survival), fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  stat_summary(data = filter(daa, age == "adult"), aes(y = survival), fun.y = "mean", geom = "point", size = 2, position = position_dodge(0.2), shape = 21, color = "black") +
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
  xlab("Adult Elymus density") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# combine
pdf("./output/ev-survival-by-treatment-2018-summer-full-models.pdf", height = 6.5, width = 9)
plot_grid(psms, psss, psas, psma, pssa, psaa, nrow = 2)
dev.off()

pdf("./output/ev-survival-by-treatment-2018-winter-full-models.pdf", height = 6, width = 6)
plot_grid(pass, paas, pasa, paaa, nrow = 2)
dev.off()

