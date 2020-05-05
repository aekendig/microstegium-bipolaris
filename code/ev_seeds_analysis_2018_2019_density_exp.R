##### info ####

# file: ev_seeds_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 5/5/20
# goal: evaluate the effects of density treatments and environmental covariates on the seed production of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
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

foc19b <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 4),
                ID = rep(c("1", "2", "3", "A"), 4),
                age = rep(c(rep("seedling", 3), "adult"), 4)) %>%
  mutate(sp = "Ev",
         focal = 1) %>%
  merge(treat, all = T) %>%
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

# 2019 seed with "none" treatment repeated across backgrounds
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
         fungicide = recode(treatment, water = 0, fungicide = 1))

# 2019 seed with single "none" treatment
seed19b <- spike19 %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g),
            seeds = sum(seeds)) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1) %>%
  ungroup() %>%
  full_join(foc19b) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1))


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
        legend.title = element_text(size = 12),
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

# all background treatments
plot_grid(ggplot(filter(seed19b, fungicide == 0), 
                 aes(x = background_density, y = seeds, color = background)) +
            stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(1)) +
            stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(1)) +
            facet_wrap(~ plant_type, scales = "free", nrow = 2, strip.position = "right") +
            temp_theme +
            xlab("Background density") +
            ggtitle("Water treatment") +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.position = "none"),
          ggplot(filter(seed19b, fungicide == 1), 
                 aes(x = background_density, y = seeds, color = background)) +
            stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(1)) +
            stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(1)) +
            facet_wrap(~ plant_type, scales = "free", nrow = 2, strip.position = "right") +
            temp_theme +
            xlab("Background density") +
            ggtitle("Fungicide treatment") +
            theme(plot.title = element_text(hjust = 0.5)),
          nrow = 2)

dev.off()


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


#### separate data ####

# non-repeat data
ea_dat_19b <- filter(seed19b, age == "adult")
es_dat_19b <- filter(seed19b, age == "seedling")

# data
# use mean and sd from non-repeat data to scale
ea_dat_19 <- filter(seed19, age == "adult") %>%
  mutate(seeds_scale = (seeds - mean(ea_dat_19b$seeds)) / sd(ea_dat_19b$seeds))
es_dat_19 <- filter(seed19, age == "seedling") %>%
  mutate(seeds_scale = (seeds - mean(es_dat_19b$seeds)) / sd(es_dat_19b$seeds),
         background = fct_relevel(background, "Ev seedling"))

# histograms
ggplot(ea_dat_19, aes(x = seeds_scale)) +
  geom_histogram(binwidth = 0.1)
ggplot(es_dat_19, aes(x = seeds_scale)) +
  geom_histogram(binwidth = 0.1)

# treatment
ea_dat_fungicide_19 <- filter(ea_dat_19, fungicide == 1)
ea_dat_water_19 <- filter(ea_dat_19, fungicide == 0)
es_dat_fungicide_19 <- filter(es_dat_19, fungicide == 1)
es_dat_water_19 <- filter(es_dat_19, fungicide == 0)


#### manual model fitting ####

# sigmoid Beverton-Holt
sbh_fun_simple <- function(dat_in, r, a, d){
  
  # extract values
  xmin = min(dat_in$background_density)
  xmax = max(dat_in$background_density)
  y0 = filter(dat_in,
              background_density == 0 ) %>%
          summarise(mean_seed = mean(seeds_scale)) %>%
          as.numeric() %>%
  round(2)
  print(y0)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 300)) %>%
    mutate(y = y0 + r * x^(d-1) / (1 + a * x^d))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds_scale)) +
          stat_summary(geom = "point", fun = "mean", shape = 21, aes(color = background)) +
          stat_summary(geom = "line", fun = "mean", aes(color = background)) +
    geom_line(data = dat, aes(x = x, y = y)))
}

# try values
sbh_fun_simple(ea_dat_water_19, 1, 1, 2) # alpha needs to be positive
sbh_fun_simple(ea_dat_fungicide_19, 0.5, 1, 2)
sbh_fun_simple(es_dat_water_19, 0.4, 0.1, 2)

# distributions
x <- seq(0, 10, length.out = 200)
y <- dgamma(x, shape = 1.5, scale = 1)
# y <- dlnorm(x, 0.25, 7)
# y <- dinvgamma(x, 3, 0.5)
# y <- dnorm(x, 0, 100)
plot(x, y, type = "l")

# function for trying results
sbh_fun_treatment <- function(dat_in, mod){
  
  # extract values
  xmin = min(dat_in$background_density)
  xmax = max(dat_in$background_density)
  
  # create data
  dat <- tibble(background_density = rep(seq(xmin, xmax, length.out = 300), 3),
                background = rep(c("Ev adult", "Ev seedling", "Mv seedling"), each = 300))
  dat$seeds_scale = fitted(mod, newdata = dat, re_formula = NA)[,"Estimate"]
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds_scale, color = background)) +
          stat_summary(geom = "point", fun = "mean", shape = 21) +
          geom_line(data = dat))
}


#### model building, testing priors ####

# chose gaussian because I couldn't fit the simplest negative binomial model without errors and I don't understand how the non-normal priors work with Poisson model (super small estimates when there shouldn't have been)

mean(filter(ea_dat_19, background_density == 0)$seeds_scale)

# priors
ea_seeds_2019_prior <- brm(data = ea_dat_19, family = gaussian,
                                 bf(seeds_scale ~ 0.07 + k * background_density / (1 + alpha * background_density^2),
                                    k + alpha ~ 1,
                                    nl = T),
                                 prior <- c(prior(normal(1, 10), nlpar = "k"),
                                            prior(gamma(1.5, 1), nlpar = "alpha", lb = 0)),
                                 sample_prior = "only", seed = 1234)
plot(ea_seeds_2019_prior)

# k intercept model
ea_seeds_2019_mod_1 <- brm(data = ea_dat_19, family = gaussian,
                                 bf(seeds_scale ~ 0.07 + k * background_density / (1 + background_density^2),
                                    k ~ 1,
                                    nl = T),
                                 prior <- c(prior(normal(1, 10), nlpar = "k")),
                                 iter = 6000, warmup = 1000, chains = 1, cores = 1) 
summary(ea_seeds_2019_mod_1)
plot(ea_seeds_2019_mod_1)

# k and alpha intercept model
ea_seeds_2019_mod_2 <- brm(data = ea_dat_19, family = gaussian,
                           bf(seeds_scale ~ 0.07 + k * background_density / (1 + alpha * background_density^2),
                              k + alpha ~ 1,
                              nl = T),
                           prior <- c(prior(normal(1, 10), nlpar = "k"),
                                      prior(gamma(1.5, 1), nlpar = "alpha", lb = 0)),
                           iter = 6000, warmup = 1000, chains = 1, cores = 1) 
summary(ea_seeds_2019_mod_2)
plot(ea_seeds_2019_mod_2)
# 53 divergent transitions
pairs(ea_seeds_2019_mod_2)
ea_seeds_2019_mod_2 <- update(ea_seeds_2019_mod_2,
                              control = list(adapt_delta = 0.99))

# k, alpha, and y0 intercept model
ea_seeds_2019_mod_3 <- brm(data = ea_dat_19, family = gaussian,
                                 bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    k + alpha + y0 ~ 1,
                                    nl = T),
                           prior <- c(prior(normal(1, 10), nlpar = "k"),
                                            prior(gamma(1.5, 1), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 1), nlpar = "y0")),
                           iter = 6000, warmup = 1000, chains = 1, cores = 1,
                           control = list(adapt_delta = 0.99))
# 46 divergent transisions
ea_seeds_2019_mod_3 <- update(ea_seeds_2019_mod_3, control = list(adapt_delta = 0.9999))
summary(ea_seeds_2019_mod_3)
plot(ea_seeds_2019_mod_3)

# k, alpha intercept model, y0 random effect
ea_seeds_2019_mod_4 <- brm(data = ea_dat_19, family = gaussian,
                           bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                              k + alpha ~ 1,
                              y0 ~ 1 + (1|site),
                              nl = T),
                           prior <- c(prior(normal(1, 10), nlpar = "k"),
                                      prior(gamma(1.5, 1), nlpar = "alpha", lb = 0),
                                      prior(normal(0, 1), nlpar = "y0")),
                           iter = 6000, warmup = 1000, chains = 1, cores = 1,
                           control = list(adapt_delta = 0.9999))
summary(ea_seeds_2019_mod_4)
plot(ea_seeds_2019_mod_4)

# k, alpha intercept model, y0 random effect and treatment
ea_seeds_2019_mod_5 <- brm(data = ea_dat_19, family = gaussian,
                           bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                              k + alpha ~ 1,
                              y0 ~ fungicide + (1|site),
                              nl = T),
                           prior <- c(prior(normal(1, 10), nlpar = "k"),
                                      prior(gamma(1.5, 1), nlpar = "alpha", lb = 0),
                                      prior(normal(0, 1), nlpar = "y0")),
                           iter = 6000, warmup = 1000, chains = 1, cores = 1,
                           control = list(adapt_delta = 0.9999))
ea_seeds_2019_mod_5 <- update(ea_seeds_2019_mod_5,
                              control = list(adapt_delta = 0.9999, max_treedepth = 15))
summary(ea_seeds_2019_mod_5)

# k background, y0 random effect and treatment, alpha intercept
ea_seeds_2019_mod_6 <- brm(data = ea_dat_19, family = gaussian,
                           bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                              alpha ~ 1,
                              k ~ background,
                              y0 ~ fungicide + (1|site),
                              nl = T),
                           prior <- c(prior(normal(1, 10), nlpar = "k"),
                                      prior(normal(0, 1), nlpar = "k", coef = "backgroundEvseedling"),
                                      prior(normal(0, 1), nlpar = "k", coef = "backgroundMvseedling"),
                                      prior(gamma(1.5, 1), nlpar = "alpha", lb = 0),
                                      prior(normal(0, 1), nlpar = "y0")),
                           iter = 6000, warmup = 1000, chains = 1, cores = 1,
                           control = list(adapt_delta = 0.9999, max_treedepth = 15))
summary(ea_seeds_2019_mod_6)

# k and alpha background, y0 random effect, water-only
ea_seeds_water_2019_mod_7 <- brm(data = ea_dat_water_19, family = gaussian,
                           bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                              alpha + k ~ background,
                              y0 ~ 1 + (1|site),
                              nl = T),
                           prior <- c(prior(normal(1, 10), nlpar = "k"),
                                      prior(normal(0, 1), nlpar = "k", coef = "backgroundEvseedling"),
                                      prior(normal(0, 1), nlpar = "k", coef = "backgroundMvseedling"),
                                      prior(gamma(1.5, 1), nlpar = "alpha", lb = 0),
                                      prior(normal(0, 1), nlpar = "alpha", coef = "backgroundEvseedling"),
                                      prior(normal(0, 1), nlpar = "alpha", coef = "backgroundMvseedling"),
                                      prior(normal(0, 1), nlpar = "y0")),
                           iter = 6000, warmup = 1000, chains = 1, cores = 1,
                           control = list(adapt_delta = 0.99999, max_treedepth = 15))
summary(ea_seeds_water_2019_mod_7)
sbh_fun_treatment(ea_dat_water_19, ea_seeds_water_2019_mod_7)


#### Ev adult seed models 2019 ####

# water treatment
ea_seeds_water_2019_mod <- update(ea_seeds_water_2019_mod_7,
                           chains = 3, cores = 2)
summary(ea_seeds_water_2019_mod)
sbh_fun_treatment(ea_dat_water_19, ea_seeds_water_2019_mod)
pp_check(ea_seeds_water_2019_mod, nsamples = 50)

# fungicide treatment
ea_seeds_fungicide_2019_mod <- brm(data = ea_dat_fungicide_19, family = gaussian,
                                 bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                    alpha + k ~ background,
                                    y0 ~ 1 + (1|site),
                                    nl = T),
                                 prior <- c(prior(normal(0.5, 10), nlpar = "k"),
                                            prior(normal(0, 10), nlpar = "k", coef = "backgroundEvseedling"),
                                            prior(normal(0, 10), nlpar = "k", coef = "backgroundMvseedling"),
                                            prior(gamma(1.5, 1), nlpar = "alpha", lb = 0),
                                            prior(normal(0, 1), nlpar = "alpha", coef = "backgroundEvseedling"),
                                            prior(normal(0, 1), nlpar = "alpha", coef = "backgroundMvseedling"),
                                            prior(normal(0, 1), nlpar = "y0")),
                                 iter = 6000, warmup = 1000, chains = 3, cores = 2,
                                 control = list(adapt_delta = 0.99999, max_treedepth = 15))
summary(ea_seeds_fungicide_2019_mod)
sbh_fun_treatment(ea_dat_fungicide_19, ea_seeds_fungicide_2019_mod)
pp_check(ea_seeds_fungicide_2019_mod, nsamples = 50)


#### Ev seedling seed models 2019 ####

# water treatment
es_seeds_water_2019_mod <- brm(data = es_dat_water_19, family = gaussian,
                                   bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                      alpha + k ~ background,
                                      y0 ~ 1 + (1|site),
                                      nl = T),
                                   prior <- c(prior(normal(0.4, 10), nlpar = "k"),
                                              prior(normal(0, 10), nlpar = "k", coef = "backgroundEvadult"),
                                              prior(normal(0, 10), nlpar = "k", coef = "backgroundMvseedling"),
                                              prior(gamma(1.1, 1), nlpar = "alpha", lb = 0),
                                              prior(normal(0, 1), nlpar = "alpha", coef = "backgroundEvadult"),
                                              prior(normal(0, 1), nlpar = "alpha", coef = "backgroundMvseedling"),
                                              prior(normal(0, 1), nlpar = "y0")),
                                   iter = 6000, warmup = 1000, chains = 3, cores = 2,
                                   control = list(adapt_delta = 0.99999, max_treedepth = 15))
summary(es_seeds_water_2019_mod)
sbh_fun_treatment(es_dat_water_19, es_seeds_water_2019_mod)
pp_check(es_seeds_water_2019_mod, nsamples = 50)

# fungicide treatment
es_seeds_fungicide_2019_mod <- brm(data = es_dat_fungicide_19, family = gaussian,
                               bf(seeds_scale ~ y0 + k * background_density / (1 + alpha * background_density^2),
                                  alpha + k ~ background,
                                  y0 ~ 1 + (1|site),
                                  nl = T),
                               prior <- c(prior(normal(0.4, 10), nlpar = "k"),
                                          prior(normal(0, 10), nlpar = "k", coef = "backgroundEvadult"),
                                          prior(normal(0, 10), nlpar = "k", coef = "backgroundMvseedling"),
                                          prior(gamma(1.1, 1), nlpar = "alpha", lb = 0),
                                          prior(normal(0, 1), nlpar = "alpha", coef = "backgroundEvadult"),
                                          prior(normal(0, 1), nlpar = "alpha", coef = "backgroundMvseedling"),
                                          prior(normal(0, 1), nlpar = "y0")),
                               iter = 6000, warmup = 1000, chains = 3, cores = 2,
                               control = list(adapt_delta = 0.99999, max_treedepth = 15))
summary(es_seeds_fungicide_2019_mod)
sbh_fun_treatment(es_dat_fungicide_19, es_seeds_fungicide_2019_mod)
pp_check(es_seeds_fungicide_2019_mod, nsamples = 50)


#### estimate figure ####

# posterior draws
ea_seeds_water_2019_draws <- posterior_samples(ea_seeds_water_2019_mod) %>%
  as_tibble() %>%
  mutate(treatment = "water")
ea_seeds_fungicide_2019_draws <- posterior_samples(ea_seeds_fungicide_2019_mod) %>%
  as_tibble() %>%
  mutate(treatment = "fungicide")
es_seeds_water_2019_draws <- posterior_samples(es_seeds_water_2019_mod) %>%
  as_tibble() %>%
  mutate(treatment = "water")
es_seeds_fungicide_2019_draws <- posterior_samples(es_seeds_fungicide_2019_mod) %>%
  as_tibble() %>%
  mutate(treatment = "fungicide")

# combine adult
ea_seeds_2019_draws <- ea_seeds_water_2019_draws %>%
  full_join(ea_seeds_fungicide_2019_draws) %>%
  transmute(age = "adult",
            treatment = treatment,
            y0 = b_y0_Intercept,
            alpha_Ev_adult = b_alpha_Intercept,
            alpha_Ev_seedling = b_alpha_Intercept + b_alpha_backgroundEvseedling,
            alpha_Mv_seedling = b_alpha_Intercept + b_alpha_backgroundMvseedling,
            r_Ev_adult = b_k_Intercept,
            r_Ev_seedling = b_k_Intercept + b_k_backgroundEvseedling,
            r_Mv_seedling = b_k_Intercept + b_k_backgroundMvseedling)

# combine seedling
es_seeds_2019_draws <- es_seeds_water_2019_draws %>%
  full_join(es_seeds_fungicide_2019_draws) %>%
  transmute(age = "seedling",
            treatment = treatment,
            y0 = b_y0_Intercept,
            alpha_Ev_adult = b_alpha_Intercept + b_alpha_backgroundEvadult,
            alpha_Ev_seedling = b_alpha_Intercept,
            alpha_Mv_seedling = b_alpha_Intercept + b_alpha_backgroundMvseedling,
            r_Ev_adult = b_k_Intercept + b_k_backgroundEvadult,
            r_Ev_seedling = b_k_Intercept,
            r_Mv_seedling = b_k_Intercept + b_k_backgroundMvseedling)

# combine ages
seeds_2019_draws <- ea_seeds_2019_draws %>%
  full_join(es_seeds_2019_draws)  %>%
  mutate(comp_Ev_adult = alpha_Ev_adult / r_Ev_adult,
         comp_Ev_seedling = alpha_Ev_seedling / r_Ev_seedling,
         comp_Mv_seedling = alpha_Mv_seedling / r_Mv_seedling) %>%
  gather(key = "coefficient", value = "value", -c(age, treatment)) %>%
  mutate(Response = recode(age, adult = "Adult Ev seeds", seedling = "First-year Ev seeds"))

# all coefficient figure
ggplot(seeds_2019_draws, aes(y = coefficient, x = value, color = treatment)) +
  stat_pointintervalh(point_interval = median_hdi) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ age)

# alpha values
alpha_seeds_2019_draws <- seeds_2019_draws %>%
  filter(grepl("alpha", coefficient) == T) %>%
  mutate(background = gsub("alpha_", "", coefficient) %>%
           gsub("_", " ", .))

# alpha figure
ggplot(alpha_seeds_2019_draws, aes(x = background, y = value, color = treatment)) +
  stat_pointinterval(point_interval = median_hdi, position = position_dodge(0.3)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  facet_wrap(~ Response) +
  temp_theme +
  theme(strip.text = element_text(size = 12))  +
  ylab("alpha value") +
  xlab("Neighbors") +
  coord_flip()

# r values
r_seeds_2019_draws <- seeds_2019_draws %>%
  filter(grepl("r", coefficient) == T) %>%
  mutate(background = gsub("r_", "", coefficient) %>%
           gsub("_", " ", .))

# r figure
ggplot(r_seeds_2019_draws, aes(x = background, y = value, color = treatment)) +
  stat_pointinterval(point_interval = median_hdi, position = position_dodge(0.3)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  facet_wrap(~ Response) +
  temp_theme +
  theme(strip.text = element_text(size = 12))  +
  ylab("r value") +
  xlab("Neighbors") +
  coord_flip()

# comp values
comp_seeds_2019_draws <- seeds_2019_draws %>%
  filter(grepl("comp", coefficient) == T) %>%
  mutate(background = gsub("comp_", "", coefficient) %>%
           gsub("_", " ", .))

# comp figure
ggplot(comp_seeds_2019_draws, aes(x = treatment, y = value, color = background)) +
  stat_pointinterval(point_interval = median_hdi, position = position_dodge(0.3)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  facet_wrap(~ Response) +
  temp_theme +
  theme(strip.text = element_text(size = 12))  +
  ylab("comp value") +
  xlab("Neighbors")


#### combined figure ####

# simulation dataset
sim_dat <- tibble(background = rep(c("Ev adult", "Ev seedling", "Mv seedling"), each = 300),
                  background_density = c(seq(0, 8, length.out = 300), seq(0, 16, length.out = 300), seq(0, 64, length.out = 300)))

# predictions
ea_seeds_water_2019_pred <- sim_dat %>%
  mutate(treatment = "water",
         age = "adult",
         response = "Adult Ev seeds",
         seeds_scale = fitted(ea_seeds_water_2019_mod, newdata = sim_dat, re_formula = NA)[,"Estimate"])

ea_seeds_fungicide_2019_pred <- sim_dat %>%
  mutate(treatment = "fungicide",
         age = "adult",
         response = "Adult Ev seeds",
         seeds_scale = fitted(ea_seeds_fungicide_2019_mod, newdata = sim_dat, re_formula = NA)[,"Estimate"])

es_seeds_water_2019_pred <- sim_dat %>%
  mutate(treatment = "water",
         age = "seedling",
         response = "First-year Ev seeds",
         seeds_scale = fitted(es_seeds_water_2019_mod, newdata = sim_dat, re_formula = NA)[,"Estimate"])

es_seeds_fungicide_2019_pred <- sim_dat %>%
  mutate(treatment = "fungicide",
         age = "seedling",
         response = "First-year Ev seeds",
         seeds_scale = fitted(es_seeds_fungicide_2019_mod, newdata = sim_dat, re_formula = NA)[,"Estimate"])

# combine by age
ea_seeds_2019_pred <- ea_seeds_water_2019_pred %>%
  full_join(ea_seeds_fungicide_2019_pred) %>%
  mutate(seeds = seeds_scale * sd(ea_dat_19b$seeds) + mean(ea_dat_19b$seeds))

es_seeds_2019_pred <- es_seeds_water_2019_pred %>%
  full_join(es_seeds_fungicide_2019_pred) %>%
  mutate(seeds = seeds_scale * sd(es_dat_19b$seeds) + mean(es_dat_19b$seeds))

# coefficient datasets
ea_seeds_2019_coef <- r_seeds_2019_draws %>%
  mutate(coef = "r") %>%
  full_join(alpha_seeds_2019_draws %>%
              mutate(coef = "alpha")) %>%
  filter(age == "adult") %>%
  mutate(Treatment = recode(treatment, fungicide = "Fungicide", water = "Control (water)"))

es_seeds_2019_coef <- r_seeds_2019_draws %>%
  mutate(coef = "r") %>%
  full_join(alpha_seeds_2019_draws %>%
              mutate(coef = "alpha")) %>%
  filter(age == "seedling") %>%
  mutate(Treatment = recode(treatment, fungicide = "Fungicide", water = "Control (water)"))

# color palette
col_pal = c("#0072B2", "#56B4E9", "#D55E00")
dodge_value = 1

# simulation figures
ea_seeds_2019_sim_fig <- ea_dat_19 %>%
  mutate(Treatment = recode(treatment, fungicide = "Fungicide", water = "Control (water)")) %>%
  ggplot(aes(x = background_density, y = seeds, color = background)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(dodge_value)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.3, position = position_dodge(dodge_value)) +
  geom_line(data = ea_seeds_2019_pred %>%
              mutate(Treatment = recode(treatment, fungicide = "Fungicide", water = "Control (water)"))) +
  facet_wrap(~ Treatment) +
  scale_color_manual(values = col_pal) +
  temp_theme +
  theme(strip.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust= 0.5),
        legend.position = "none")  +
  ylab("Seeds") +
  xlab(expression(paste("Neighbor density (plants ", m^-2, ")", sep = ""))) +
  ggtitle("Adult Ev Seeds")

es_seeds_2019_sim_fig <- es_dat_19 %>%
  mutate(Treatment = recode(treatment, fungicide = "Fungicide", water = "Control (water)"),
         background = fct_relevel(background, "Ev adult")) %>%
  ggplot(aes(x = background_density, y = seeds, color = background)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(dodge_value)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.3, position = position_dodge(dodge_value)) +
  geom_line(data = es_seeds_2019_pred %>%
              mutate(Treatment = recode(treatment, fungicide = "Fungicide", water = "Control (water)"))) +
  facet_wrap(~ Treatment) +
  scale_color_manual(values = col_pal) +
  temp_theme +
  theme(strip.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust= 0.5),
        legend.position = "none")  +
  ylab("Seeds") +
  xlab(expression(paste("Neighbor density (plants ", m^-2, ")", sep = ""))) +
  ggtitle("First-year Ev Seeds")

# coefficient figures
ea_seeds_2019_coef_fig <- ggplot(ea_seeds_2019_coef, aes(x = coef, y = value, color = background)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  stat_pointinterval(point_interval = median_hdi, position = position_dodge(0.5)) +
  facet_wrap(~ Treatment) +
  scale_color_manual(values = col_pal, name = "Neighbors") +
  scale_x_discrete(labels = c("alpha" = expression(alpha),
                              "r" = "r")) +
  temp_theme +
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", face = "italic"))  +
  ylab("Estimate")

es_seeds_2019_coef_fig <- ggplot(es_seeds_2019_coef, aes(x = coef, y = value, color = background)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  stat_pointinterval(point_interval = median_hdi, position = position_dodge(0.5)) +
  facet_wrap(~ Treatment) +
  scale_color_manual(values = col_pal, name = "Neighbors") +
  scale_x_discrete(labels = c("alpha" = expression(alpha),
                              "r" = "r")) +
  temp_theme +
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", face = "italic"),
        legend.position = "none")  +
  ylab("Estimate")

# legend
leg <- get_legend(ea_seeds_2019_coef_fig +
                    theme(legend.margin = margin(0, 0, 0, 0, unit="cm")))

# combine figures
figs <- plot_grid(ea_seeds_2019_sim_fig,
                  es_seeds_2019_sim_fig,
                  ea_seeds_2019_coef_fig +
                    theme(legend.position = "none"),
                  es_seeds_2019_coef_fig,
                  nrow = 2,
                  rel_heights = c(1, 0.8),
                  labels = c("a", "c", "b", "d"),
                  label_y = c(0.9, 0.9, 1, 1))

figs_leg <- plot_grid(figs, leg,
                      nrow = 2,
                      rel_heights = c(1, 0.08))


#### output ####
pdf("./output/ev_seeds_analysis_models_2019_density_exp.pdf")
figs_leg
dev.off()

save(ea_seeds_water_2019_mod, file = "./output/ev_seeds_adult_water_2019_density_exp.rda")
save(ea_seeds_fungicide_2019_mod, file = "./output/ev_seeds_adult_fungicide_2019_density_exp.rda")
save(es_seeds_water_2019_mod, file = "./output/ev_seeds_seedling_water_2019_density_exp.rda")
save(es_seeds_fungicide_2019_mod, file = "./output/ev_seeds_seedling_fungicide_2019_density_exp.rda")

  #### integrating spikelet weight ####
# bf(spikelet_weight ~ beta*seeds,
#    seeds ~ s0 + background...,
   # priors: log-normal, exponential, gamma, inverse gamma - strictly positive, doesn't make a big difference which one you use - if results vary a lot - look into those
   # look up latent variable automation - Simon will look up
   # can go with regression approach or specify in Stan
   # start with minimum model and add on treatment-specific models
   # maybe pick simpler

