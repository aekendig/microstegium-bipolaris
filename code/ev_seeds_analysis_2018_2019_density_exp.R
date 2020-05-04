##### info ####

# file: ev_seeds_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/28/20
# goal: evaluate the effects of density treatments and environmental covariates on the seed production of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(bayesplot)
library(invgamma)

# import data
spike18 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
spike19 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
surv18 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# focal plants 2019 (all were replaced and tracked)
foc19 <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 4),
                ID = rep(c("1", "2", "3", "A"), 4),
                age = rep(c(rep("seedling", 3), "adult"), 4)) %>%
  mutate(sp = "Ev",
         focal = 1) %>%
  merge(plots, all = T) %>%
  as_tibble()

# 2018 plants
plants18 <- surv18 %>%
  filter(month == "September" & !is.na(survival) & !(sp == "Mv" & focal == 0))

# check that all focal are there
plants18 %>%
  filter(focal == 1) %>%
  group_by(sp, age, ID) %>%
  summarise(plots = n())

# focal plants 2018
ev18 <- plants18 %>%
  filter(sp == "Ev" & survival == 1) %>%
  select(site, plot, treatment, sp, ID, age, focal)

# check 2018 data
unique(spike18$ID_unclear)
unique(spike18$spikelet_notes)

# 2018 seed
seed18 <- spike18 %>%
  filter(ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, ID, focal, age) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g),
            seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(ev18) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         seeds_round = round(seeds))

# separate focal
fseed18 <- seed18 %>%
  filter(focal == 1)

# check 2019 data
unique(spike19$spikelet_notes)
filter(spike19, spikelet_notes == "D2?") # does this need to be removed?
unique(spike19$processing_notes)

# 2019 seed
seed19 <- spike19 %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g),
            seeds = sum(seeds)) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1) %>%
  ungroup() %>%
  full_join(foc19) %>%
  left_join(covar) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         seeds_round = round(seeds))


#### check datasets ####

# sample sizes
samps_19 <- seed19 %>%
  group_by(site, plot, treatment) %>%
  count() %>%
  mutate(seeds = 405)

samps_18 <- fseed18 %>%
  group_by(site, plot, treatment) %>%
  count() %>%
  mutate(seeds = 350)

# visualize
samp_plot_19 <- ggplot(seed19, aes(x = plot, y = seeds)) +
  geom_point(alpha = 0.5, aes(color = age)) +
  facet_grid(treatment ~ site)

samp_plot_19 +
  geom_text(data = samps_19, aes(label = n), size = 2)

samp_plot_19 %+%
  fseed18 +
  geom_text(data = samps_18, aes(label = n), size = 2)


#### visualize treatment effects ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside")

# template figure
temp_fig <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

# save treatment figures
pdf("./output/ev_seeds_visualize_treatment_2018_2019_density_exp.pdf")

# focal plants 2018
plot_grid(temp_fig %+%
            fseed18 %+%
            aes(y = spikelet_weight.g) +
            ylab("Year 1 focal spikelet weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig %+%
            fseed18 %+%
            aes(y = seeds) +
            ylab("Year 1 focal seeds"),
          nrow = 2)
  
# focal plants 2019
plot_grid(temp_fig %+%
            seed19 %+%
            aes(y = spikelet_weight.g) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig %+%
            seed19 %+%
            aes(y = seeds) +
            ylab("Year 2 focal seeds"),
          nrow = 2)

dev.off()

# all background treatments
ggplot(filter(seed19, fungicide == 0), 
       aes(x = background_density, y = seeds, color = background)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(1)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(1)) +
  facet_wrap(~ plant_type, scales = "free", nrow = 2) +
  temp_theme +
  xlab("Background density")


#### visualize covariate effects ####

# covariate template figure
temp_fig_cov <- ggplot(fseed18, aes(x = soil_moisture_jun.prop, y = spikelet_weight.g)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "loess") +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  temp_theme

# save covariate figures
pdf("./output/ev_seeds_visualize_covariates_2018_2019_density_exp.pdf")

# June soil moisture focal year 1
plot_grid(temp_fig_cov +
            ylab("Year 1 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(y = seeds) +
            xlab("June soil moisture") +
            ylab("Year 1 focal seeds"),
          nrow = 2)

# June soil moisture focal year 2
plot_grid(temp_fig_cov %+%
            seed19 +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            seed19 %+%
            aes(y = seeds) +
            xlab("June soil moisture") +
            ylab("Year 2 focal seeds"),
          nrow = 2)

# October soil moisture focal year 1
plot_grid(temp_fig_cov %+%
            aes(x = soil_moisture_oct.prop) +
            ylab("Year 1 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(x = soil_moisture_oct.prop, y = seeds) +
            xlab("October soil moisture") +
            ylab("Year 1 focal seeds"),
          nrow = 2)

# October soil moisture focal year 2
plot_grid(temp_fig_cov %+%
            seed19 %+%
            aes(x = soil_moisture_oct.prop) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            seed19 %+%
            aes(x = soil_moisture_oct.prop, y = seeds) +
            xlab("October soil moisture") +
            ylab("Year 2 focal seeds"),
          nrow = 2)

# canopy cover focal year 1
plot_grid(temp_fig_cov %+%
            aes(x = canopy_cover.prop) +
            ylab("Year 1 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(x = canopy_cover.prop, y = seeds) +
            xlab("Canopy cover") +
            ylab("Year 1 focal seeds"),
          nrow = 2)

# canopy cover focal year 2
plot_grid(temp_fig_cov %+%
            seed19 %+%
            aes(x = canopy_cover.prop) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            seed19 %+%
            aes(x = canopy_cover.prop, y = seeds) +
            xlab("Canopy cover") +
            ylab("Year 2 focal seeds"),
          nrow = 2)

# September biomass focal year 1
plot_grid(temp_fig_cov %+%
            aes(x = mv_sep.g) +
            ylab("Year 1 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(x = mv_sep.g, y = seeds) +
            xlab("September biomass") +
            ylab("Year 1 focal seeds"),
          nrow = 2)

# September biomass focal year 2
plot_grid(temp_fig_cov %+%
            seed19 %+%
            aes(x = mv_sep.g) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            seed19 %+%
            aes(x = mv_sep.g, y = seeds) +
            xlab("September biomass") +
            ylab("Year 2 focal seeds"),
          nrow = 2)

dev.off()


#### separate datasets ####

# data
# scale seeds
# scale density
ea_ea_dat <- filter(seed19, age == "adult" & background == "Ev adult") %>%
  mutate(seeds_scale = scale(seeds),
         density_scale = background_density / max(background_density))
ea_es_dat <- filter(seed19, age == "adult" & background == "Ev seedling")
ea_ms_dat <- filter(seed19, age == "adult" & background == "Mv seedling")

# histograms
ggplot(ea_ea_dat, aes(x = seeds_round)) +
  geom_histogram(binwidth = 3)
ggplot(ea_es_dat, aes(x = seeds_round)) +
  geom_histogram(binwidth = 3)
ggplot(ea_ms_dat, aes(x = seeds_round)) +
  geom_histogram(binwidth = 3)

# summary statistics
mean(ea_ea_dat$seeds_round)
var(ea_ea_dat$seeds_round)


#### manual model fitting ####

# sigmoid Beverton-Holt
sbh_fun_simple <- function(dat_in, r, a, d){
  
  # extract values
  xmin = min(dat_in$density_scale)
  xmax = max(dat_in$density_scale)
  y0 = filter(dat_in,
              density_scale == 0 & fungicide == 0) %>%
          summarise(mean_seed = mean(seeds_scale)) %>%
          as.numeric() %>%
  round(2)
  print(y0)
  
  # create data
  dat <- tibble(x = rep(seq(xmin, xmax, length.out = 100), 2),
                treatment = rep(c("water", "fungicide"), each = 100)) %>%
    mutate(y = y0 + r * x^(d-1) / (1 + r * a * x^d))
  
  # plot
  print(ggplot(dat_in, aes(x = density_scale, y = seeds_scale)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment))) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
    geom_line(data = dat, aes(x = x, y = y))
}

sbh_fun_treatment <- function(dat_in, rw, rf, aw, af, d, y0w, y0f){
  
  # extract values
  xmin = min(dat_in$density_scale)
  xmax = max(dat_in$density_scale)
  
  # create data
  dat <- tibble(x = rep(seq(xmin, xmax, length.out = 100), 2),
                treatment = rep(c("water", "fungicide"), each = 100)) %>%
    mutate(y = case_when(treatment == "water" ~ y0w + rw * x^(d-1) / (1 + rw * aw * x^d),
                         treatment == "fungicide" ~ y0f + rf * x^(d-1) / (1 + rf * af * x^d)))
  
  # plot
  print(ggplot(dat_in, aes(x = density_scale, y = seeds_scale, color = treatment)) +
          stat_summary(geom = "point", fun = "mean")) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
          geom_line(data = dat, aes(x = x, y = y))
}

# try values
sbh_fun_simple(ea_ea_dat, 10, 2.5, 2) # alpha needs to be positive
sbh_fun_simple(ea_ea_dat, 2, 5, 2)
sbh_fun_simple(ea_ea_dat, 0.3, 0.0005, 3)
sbh_fun_simple(ea_es_dat, 10, 0.1, 2)
sbh_fun_simple(ea_ms_dat, 0.3, 0.0005, 3)

# distributions
x <- seq(0, 10, length.out = 200)
y <- dgamma(x, shape = 2.5, scale = 1)
# y <- dlnorm(x, 0.25, 7)
# y <- dinvgamma(x, 3, 0.5)
# y <- dnorm(x, 0, 100)
plot(x, y, type = "l")


#### model building, testing priors ####

# priors
# ea_seeds_ea_bg_2019_prior <- brm(data = ea_ea_dat, family = negbinomial,
#                                  bf(seeds_round ~ 107 + background_density / ((1/k) + alpha * background_density^2),
#                                     k + alpha ~ 1,
#                                     nl = T),
#                                  prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
#                                             prior(lognormal(0.05, 5), nlpar = "alpha")),
#                                  sample_prior = "only", seed = 1234)

# alpha priors
# gamma(0.005, 1): values way too small
# gamma(0.05, 1): larger numbers, but diverget transitions and funny looking distribution
# gamma(0.5, 1): normal looking distributions, but numbers likely way too large for fit. lots of divergent transitions
# lognormal(0.05, 5): divergent transitions and funny looking distribution

# re-paramaterized priors
ea_seeds_ea_bg_2019_prior <- brm(data = ea_ea_dat, family = gaussian(),
                                 bf(seeds_scale ~ -0.32 + k * density_scale / (1 + k * alpha * density_scale^2),
                                    k + alpha ~ 1,
                                    nl = T),
                                 prior <- c(prior(normal(10, 10), nlpar = "k"),
                                            prior(gamma(2.5, 1), nlpar = "alpha", lb = 0)),
                                 sample_prior = "only", seed = 1234)
# 2829 divergent transitions with family = negbinomial

# figures of priors
mcmc_areas(as.array(ea_seeds_ea_bg_2019_prior),
           pars = c("b_alpha_Intercept"),
           prob = 0.8,
           prob_outer = 0.99,
           point_est = "mean") +
  ggplot2::labs(title = "Prior parameter distribution",
                subtitle = "with means and 80% intervals")
pairs(ea_seeds_ea_bg_2019_prior)
plot(ea_seeds_ea_bg_2019_prior)

# k intercept model
ea_seeds_ea_bg_2019_mod_1 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ 107 + k * background_density / (1 + 0.25 * background_density^2),
                                    k ~ 1,
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0)),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1) 
summary(ea_seeds_ea_bg_2019_mod_1)
plot(ea_seeds_ea_bg_2019_mod_1) # super small estimates when family = poisson and y = seeds_round

# k and alpha intercept model
ea_seeds_ea_bg_2019_mod_2 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ 107 + k * background_density / (1 + alpha * background_density^2),
                                    k + alpha ~ 1,
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
                                            prior(inv_gamma(2.5, 0.7), nlpar = "alpha", lb = 0)),
                               iter = 6000, warmup = 1000, chains = 1, cores = 1)

summary(ea_seeds_ea_bg_2019_mod_2)
plot(ea_seeds_ea_bg_2019_mod_2)

# k, alpha, and y0 intercept model
ea_seeds_ea_bg_2019_mod_3 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k + alpha + y0 ~ 1,
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
                                            prior(inv_gamma(2.5, 0.7), nlpar = "alpha", lb = 0),
                                            prior(gamma(107, 1), nlpar = "y0", lb = 0)),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1)

summary(ea_seeds_ea_bg_2019_mod_3)
plot(ea_seeds_ea_bg_2019_mod_3)

# k, alpha intercept model, y0 random effect
ea_seeds_ea_bg_2019_mod_4 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k + alpha ~ 1,
                                    y0 ~ 1 + (1|site),
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
                                            prior(inv_gamma(2.5, 0.7), nlpar = "alpha", lb = 0),
                                            prior(gamma(107, 1), nlpar = "y0", lb = 0)),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1)
# 8 divergent transitions
ea_seeds_ea_bg_2019_mod_4 <- update(ea_seeds_ea_bg_2019_mod_4, 
                                    control = list(adapt_delta = 0.99))
summary(ea_seeds_ea_bg_2019_mod_4)
plot(ea_seeds_ea_bg_2019_mod_4)

# k, alpha intercept model, y0 random effect and treatment
ea_seeds_ea_bg_2019_mod_5 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k + alpha ~ 1,
                                    y0 ~ fungicide + (1|site),
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
                                            prior(inv_gamma(2.5, 0.7), nlpar = "alpha", lb = 0),
                                            prior(gamma(107, 1), nlpar = "y0", coef = "Intercept"),
                                            prior(normal(0, 1), nlpar = "y0", coef = "fungicide")),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1, 
                                 control = list(adapt_delta = 0.99))

summary(ea_seeds_ea_bg_2019_mod_5)
plot(ea_seeds_ea_bg_2019_mod_5)

# k intercept model, y0 random effect and treatment, alpha treatment
ea_seeds_ea_bg_2019_mod_6 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k ~ 1,
                                    alpha ~ fungicide,
                                    y0 ~ fungicide + (1|site),
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
                                            prior(inv_gamma(2.5, 0.7), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 1), nlpar = "alpha", coef = "fungicide"),
                                            prior(gamma(107, 1), nlpar = "y0", lb = 0),
                                            prior(normal(0, 1), nlpar = "y0", coef = "fungicide")),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1, 
                                 control = list(adapt_delta = 0.99))

prior_summary(ea_seeds_ea_bg_2019_mod_6)
summary(ea_seeds_ea_bg_2019_mod_6)

# k and alpha treatment, y0 random effect and treatment
ea_seeds_ea_bg_2019_mod_7 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k ~ fungicide,
                                    alpha ~ fungicide,
                                    y0 ~ fungicide + (1|site),
                                    nl = T),
                                 prior <- c(prior(gamma(50, 1), nlpar = "k", lb = 0),
                                            prior(normal(0, 1), nlpar = "k", coef = "fungicide"),
                                            prior(inv_gamma(2.5, 0.7), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 1), nlpar = "alpha", coef = "fungicide"),
                                            prior(gamma(107, 1), nlpar = "y0", lb = 0),
                                            prior(normal(0, 1), nlpar = "y0", coef = "fungicide")),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1, 
                                 control = list(adapt_delta = 0.99))

prior_summary(ea_seeds_ea_bg_2019_mod_7)
summary(ea_seeds_ea_bg_2019_mod_7)
sbh_fun_treatment(ea_ea_dat, 49.69, (49.69 + 0.8), 0.57, (0.57 + 0.73), 2, 106.17, (106.17 + 0.82))


#### Ev adult seeds, Ev adult background ####

# increase number of chains, refine priors based on treatments, use scaled values
ea_seeds_ea_bg_2019_mod_8 <- brm(data = ea_ea_dat, family = gaussian,
                                 bf(seeds_scale ~ y0 + k * density_scale / (1 + k * alpha * density_scale^2),
                                    k ~ fungicide,
                                    alpha ~ fungicide,
                                    y0 ~ fungicide + (1|site),
                                    nl = T),
                                 prior <- c(prior(normal(10, 10), nlpar = "k"),
                                            prior(normal(0, 10), nlpar = "k", coef = "fungicide"),
                                            prior(gamma(2.5, 1), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 10), nlpar = "alpha", coef = "fungicide"),
                                            prior(normal(0, 1), nlpar = "y0"),
                                            prior(normal(0, 1), nlpar = "y0", coef = "fungicide")),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1, 
                                 control = list(adapt_delta = 0.99))
summary(ea_seeds_ea_bg_2019_mod_8)
sbh_fun_treatment(ea_ea_dat, 99.81, (99.81 + 0.8), 0.53, (0.53 + 0.83), 2, 79.5, (79.5 + 50.42))


#### Ev adult seeds, Ev seedling background ####

ea_seeds_es_bg_2019_mod_8 <- brm(data = ea_es_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k ~ fungicide,
                                    alpha ~ fungicide,
                                    y0 ~ fungicide + (1|site),
                                    nl = T),
                                 prior <- c(prior(gamma(10, 1), nlpar = "k", lb = 0),
                                            prior(normal(0, 1), nlpar = "k", coef = "fungicide"),
                                            prior(inv_gamma(3, 0.7), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 1), nlpar = "alpha", coef = "fungicide"),
                                            prior(gamma(81, 1), nlpar = "y0", lb = 0),
                                            prior(normal(0, 100), nlpar = "y0", coef = "fungicide")),
                                 iter = 6000, warmup = 1000, chains = 3, cores = 2, 
                                 control = list(adapt_delta = 0.99))
summary(ea_seeds_es_bg_2019_mod_8)
sbh_fun_treatment(ea_es_dat, 10.17, (10.17 + 0.8), 0.33, (0.33 + 0.82), 2, 80.61, (80.61 + 20.75))


#### Ev adult seeds, Mv seedling background ####

ea_seeds_ms_bg_2019_mod_8 <- brm(data = ea_ms_dat, family = gaussian,
                                 bf(seeds ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k ~ fungicide,
                                    alpha ~ fungicide,
                                    y0 ~ fungicide + (1|site),
                                    nl = T),
                                 prior <- c(prior(gamma(10, 1), nlpar = "k", lb = 0),
                                            prior(normal(0, 10), nlpar = "k", coef = "fungicide"),
                                            prior(inv_gamma(3, 0.5), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 10), nlpar = "alpha", coef = "fungicide"),
                                            prior(gamma(81, 1), nlpar = "y0", lb = 0),
                                            prior(normal(0, 100), nlpar = "y0", coef = "fungicide")),
                                 iter = 6000, warmup = 1000, chains = 3, cores = 2, 
                                 control = list(adapt_delta = 0.99))
summary(ea_seeds_ms_bg_2019_mod_8)
sbh_fun_treatment(ea_ms_dat, 9.98, (9.98 + 7.97), 0.26, (0.26 + 7.99), 2, 80.4, (80.4 + 37.47))


#### integrating spikelet weight ####
# bf(spikelet_weight ~ beta*seeds,
#    seeds ~ s0 + background...,
   # priors: log-normal, exponential, gamma, inverse gamma - strictly positive, doesn't make a big difference which one you use - if results vary a lot - look into those
   # look up latent variable automation - Simon will look up
   # can go with regression approach or specify in Stan
   # start with minimum model and add on treatment-specific models
   # maybe pick simpler

