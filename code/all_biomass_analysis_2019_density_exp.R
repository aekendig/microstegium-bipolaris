##### info ####

# file: all_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/3/20
# goal: analyze total biomass


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
# library(DHARMa)
# library(MASS)
# library(MuMIn)
library(tidyverse)
library(brms)
library(cowplot)
library(lubridate)

# import all raw data files
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


#### edit data ####

# modified from focal_severity_analysis_2019_density_exp.R to include separate values for Ev seedling and adult

# background biomass data
filter(bg_bio, is.na(biomass.g))
# add species and age
# remove no background plots (biomass = 0)
bg_bio2 <- bg_bio %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         focal = 0) %>%
  filter(sp != "none")

# mv biomass data
unique(mv_bio$process_notes) # all issues addressed
group_by(mv_bio, site, plot, treatment) %>%
  count() %>%
  filter(n != 3) # all plots have 3 plants
filter(mv_bio, is.na(biomass_weight.g))
filter(mv_bio, (site == "D3" & plot == 7 & treatment == "fungicide") | (site == "D4" & plot == 8 & treatment == "water"))
# because two individuals are missing, multiply the average individual weight by 3
# vegetative biomass from individuals within plots
mv_bio2 <- mv_bio %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(biomass.g = mean(biomass_weight.g, na.rm = T) * 3) %>%
  ungroup() %>%
  mutate(age = "seedling",
         focal = 1)

# ev biomass data
unique(ev_bio$processing_notes) # all issues addressed
group_by(ev_bio, site, plot, treatment) %>%
  count() %>%
  filter(n != 4) # all plots have 4 plants
filter(ev_bio, is.na(weight)) # one seedling missing
# because one individual is missing, multiply the average individual weight by 3
# vegetative biomass from individuals within plots
ev_bio2 <- ev_bio %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         TRUE ~ "seedling")) %>%
  group_by(site, plot, treatment, sp, age) %>%
  summarise(biomass.g = mean(weight, na.rm = T)) %>%
  ungroup() %>%
  mutate(biomass.g = case_when(age == "seedling" ~ biomass.g * 3,
                               TRUE ~ biomass.g),
         focal = 1)

# combine all biomass
bio <- full_join(bg_bio2, mv_bio2) %>%
  full_join(ev_bio2) %>%
  mutate(sp_age = paste(sp, age, sep = "_")) %>%
  group_by(site, plot, treatment, sp_age) %>%
  summarise(biomass.g = sum(biomass.g)) %>%
  ungroup() %>%
  spread(key = sp_age, value = biomass.g) %>%
  mutate(Ev_biomass.g = Ev_adult + Ev_seedling,
         total_biomass.g = Ev_biomass.g + Mv_seedling) %>%
  rename(Ev_adult_biomass.g = Ev_adult, Ev_seedling_biomass.g = Ev_seedling, Mv_biomass.g = Mv_seedling)

# density
dens_dat <- plots_simple %>%
  select(plot, background, background_density) %>%
  unique() %>%
  filter(background_density > 0) %>%
  spread(key = background, value = background_density) %>%
  rename(Ev_adult_density = "Ev adult", Ev_seedling_density = "Ev seedling", Mv_density = "Mv seedling") %>%
  full_join(tibble(plot = 1, Ev_adult_density = 0, Ev_seedling_density = 0, Mv_density = 0)) %>%
  mutate(Ev_adult_density = replace_na(Ev_adult_density, 0) + 1,
         Ev_seedling_density  = replace_na(Ev_seedling_density , 0) + 3,
         Ev_density = Ev_adult_density + Ev_seedling_density,
         Mv_density  = replace_na(Mv_density , 0) + 3,
         total_density = Ev_adult_density + Ev_seedling_density + Mv_density,
         Ev_present = case_when(plot > 4 ~ 1,
                                TRUE ~ 0),
         Mv_present = case_when(plot > 1 & plot < 5 ~ 1,
                                TRUE ~ 0),
         bg_present = case_when(plot > 1 ~ 1,
                                TRUE ~ 0))

# combine biomass and density data
dat <- left_join(bio, dens_dat) %>%
  mutate(Treatment = recode(treatment, water = "control (water)"),
         fungicide = recode(treatment, water = 0, fungicide = 1))

# make long
dat_bio_long <- dat %>%
  select(site:total_biomass.g) %>%
  pivot_longer(cols = ends_with("biomass.g"),
               names_to = "biomass_plants",
               values_to = "biomass.g") %>%
  mutate(biomass_plants = str_replace(biomass_plants, "_biomass.g", ""),
         biomass_plants = str_replace(biomass_plants, "_", " "))

dat_dens_long <- dat %>%
  select(c(site:treatment, Ev_adult_density:total_density)) %>%
  pivot_longer(cols = ends_with("density"),
               names_to = "density_plants",
               values_to = "density") %>%
  mutate(density_plants = str_replace(density_plants, "_density", ""),
         density_plants = str_replace(density_plants, "_", " "))

dat_long <- full_join(dat_bio_long, dat_dens_long)


#### visualize ####

# all combinations 
dat_long %>%
  ggplot(aes(x = density, y = biomass.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(biomass_plants ~ density_plants, scales = "free")
# look at interspecific interactions using the focal biomass only
# note that lowest densities are averaged across all plots where the density species is low
# the intraspecific/group effects are most interesting

# Mv density and biomass
mvd_mvb_viz <- dat %>%
  filter(plot < 5) %>%
  ggplot(aes(Mv_density, Mv_biomass.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 3)
mvd_mvb_viz
# fungicide has a stronger effect at higher density/biomass
# individual Mv biomass doesn't totally compensate for unused resources, but the individuals in lower density plots are probably bigger because the relationship isn't linear
# higher Mv density will increase biomass, up to a certain point?

# Ev density and Ev biomass
evd_evb_viz <- mvd_mvb_viz %+%
  filter(dat, !(plot %in% c(2, 3, 4))) %+%
  aes(x = Ev_density, y = Ev_biomass.g)
evd_evb_viz
# Ev biomass asymptotes (larger individuals in lower density plots)

# Ev seedling density and Ev seedling biomass
evsd_evsb_viz <- mvd_mvb_viz %+%
  filter(dat, plot %in% c(1, 5, 6, 7)) %+%
  aes(x = Ev_seedling_density, y = Ev_seedling_biomass.g)
evsd_evsb_viz
# Ev seedling biomass increases with density

# Ev adult density and Ev adult biomass
evad_evab_viz <- mvd_mvb_viz %+%
  filter(dat, plot == 1 | plot > 7) %+%
  aes(x = Ev_adult_density, y = Ev_adult_biomass.g)
evad_evab_viz
# Ev adult biomass increases across adult gradient

# combine intra-group figures
plot_grid(mvd_mvb_viz, evd_evb_viz, evad_evab_viz, evsd_evsb_viz,
          nrow = 2)


#### Mv biomass model ####

# select plots
mv_dat <- dat %>%
  filter(plot < 5)

# initial asymptote
mv_dat %>%
  filter(Mv_density > 60) %>%
  summarise(mean = mean(Mv_biomass.g),
            sd = sd(Mv_biomass.g))
# 391, 115

# initial lrc
391/2
# looks like it would be where density = 10
log(log(2)/10)
log(2)/exp(-2.7)
# -2.7
# this parameter seems to be able to take on any value
# very high -> time goes to 0
# very low -> time goes to infinity

# model
mv_bio_mod <- brm(data = mv_dat, family = gaussian,
                  bf(Mv_biomass.g ~ asym * (1 - exp(-exp(lrc) * Mv_density)),
                     asym ~ fungicide + (1|site),
                     lrc ~ fungicide,
                     nl = T),
                  prior <- c(prior(normal(391, 100), nlpar = "asym", coef = "Intercept"),
                             prior(normal(0, 100), nlpar = "asym", coef = "fungicide"),
                             prior(cauchy(0, 1), nlpar = "asym", class = sd),
                             prior(normal(-2.7, 1), nlpar = "lrc", coef = "Intercept"),
                             prior(normal(0, 1), nlpar = "lrc", coef = "fungicide"),
                             prior(cauchy(0, 1), class = sigma)),
                  iter = 6000, warmup = 1000, chains = 1, cores = 1)

# check model and increase chains
summary(mv_bio_mod)
prior_summary(mv_bio_mod)
mv_bio_mod <- update(mv_bio_mod, chains = 3, control = list(adapt_delta = 0.99))
plot(mv_bio_mod)
pp_check(mv_bio_mod, nsamples = 100)

# include priors informed by greenhouse fungicide experiment

# values in fungicide_effects_greenhouse_2019.R
21.05 / 23.80 # fungicide effect = 0.88 * control
1 - 0.88

mv_bio_fung_mod <- brm(data = mv_dat, family = gaussian,
                  bf(Mv_biomass.g ~ (asym * fung) * (1 - exp(-exp(lrc) * Mv_density)),
                     asym ~ fungicide + (1|site),
                     fung ~ fungicide,
                     lrc ~ fungicide,
                     nl = T),
                  prior <- c(prior(normal(391, 100), nlpar = "asym", coef = "Intercept"),
                             prior(normal(0, 100), nlpar = "asym", coef = "fungicide"),
                             prior(cauchy(0, 1), nlpar = "asym", class = sd),
                             prior(normal(1, 0.001), nlpar = "fung", coef = "Intercept"),
                             prior(normal(-0.12, 0.001), nlpar = "fung", coef = "fungicide"),
                             prior(normal(-2.7, 1), nlpar = "lrc", coef = "Intercept"),
                             prior(normal(0, 1), nlpar = "lrc", coef = "fungicide"),
                             prior(cauchy(0, 1), class = sigma)),
                  iter = 6000, warmup = 1000, chains = 1, cores = 1, 
                  control = list(adapt_delta = 0.99))

# check model and increase chains
summary(mv_bio_fung_mod)
prior_summary(mv_bio_fung_mod)
mv_bio_fung_mod <- update(mv_bio_fung_mod, chains = 3)
plot(mv_bio_fung_mod)
pp_check(mv_bio_fung_mod, nsamples = 100)


#### Ev biomass model ####

# select plots
ev_dat <- dat %>%
  filter(plot == 1 | plot > 4)

# initial asymptote
ev_dat %>%
  filter(Ev_density > 15) %>%
  summarise(mean = mean(Ev_biomass.g),
            sd = sd(Ev_biomass.g))
# 56, 23

# initial lrc
56/2
# looks like it would be where density = 6
log(log(2)/6)
log(2)/exp(-2.2)
# -2.2

# model
ev_bio_mod <- brm(data = ev_dat, family = gaussian,
                  bf(Ev_biomass.g ~ asym * (1 - exp(-exp(lrc) * Ev_density)),
                     asym ~ fungicide + (1|site),
                     lrc ~ fungicide,
                     nl = T),
                  prior <- c(prior(normal(56, 10), nlpar = "asym", coef = "Intercept"),
                             prior(normal(0, 10), nlpar = "asym", coef = "fungicide"),
                             prior(cauchy(0, 1), nlpar = "asym", class = sd),
                             prior(normal(-2.2, 1), nlpar = "lrc", coef = "Intercept"),
                             prior(normal(0, 1), nlpar = "lrc", coef = "fungicide"),
                             prior(cauchy(0, 1), class = sigma)),
                  iter = 6000, warmup = 1000, chains = 1, cores = 1,
                  control = list(adapt_delta = 0.99))

# check model and increase chains
summary(ev_bio_mod)
ev_bio_mod <- update(ev_bio_mod, chains = 3)
plot(ev_bio_mod)
pp_check(ev_bio_mod, nsamples = 100)


# include priors informed by greenhouse fungicide experiment

# values in fungicide_effects_greenhouse_2019.R
13.84 / 13.39 # fungicide effect = 1.03 * control

ev_bio_fung_mod <- brm(data = ev_dat, family = gaussian,
                  bf(Ev_biomass.g ~ (asym * fung) * (1 - exp(-exp(lrc) * Ev_density)),
                     asym ~ fungicide + (1|site),
                     fung ~ fungicide,
                     lrc ~ fungicide,
                     nl = T),
                  prior <- c(prior(normal(56, 10), nlpar = "asym", coef = "Intercept"),
                             prior(normal(0, 10), nlpar = "asym", coef = "fungicide"),
                             prior(cauchy(0, 1), nlpar = "asym", class = sd),
                             prior(normal(1, 0.001), nlpar = "fung", coef = "Intercept"),
                             prior(normal(0.03, 0.001), nlpar = "fung", coef = "fungicide"),
                             prior(normal(-2.2, 1), nlpar = "lrc", coef = "Intercept"),
                             prior(normal(0, 1), nlpar = "lrc", coef = "fungicide"),
                             prior(cauchy(0, 1), class = sigma)),
                  iter = 6000, warmup = 1000, chains = 1, cores = 1,
                  control = list(adapt_delta = 0.99))

# check model and increase chains
summary(ev_bio_fung_mod)
prior_summary(ev_bio_fung_mod)
ev_bio_fung_mod <- update(ev_bio_fung_mod, chains = 3)
plot(ev_bio_fung_mod)
pp_check(ev_bio_fung_mod, nsamples = 100)


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
        legend.position = "bottom", 
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        plot.title = element_text(size = 10, hjust = 0.5))

# colors
col_pal = c("#a6611a", "#018571")

# simulate data
mv_sim_dat <- tibble(Mv_density = rep(seq(0, 67, length.out = 300), 2),
                     fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Mv_biomass.g = fitted(mv_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
         Mv_biomass_lower = fitted(mv_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         Mv_biomass_upper = fitted(mv_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) 

mv_fung_sim_dat <- tibble(Mv_density = rep(seq(0, 67, length.out = 300), 2),
                     fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Mv_biomass.g = fitted(mv_bio_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         Mv_biomass_lower = fitted(mv_bio_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         Mv_biomass_upper = fitted(mv_bio_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) 

ev_sim_dat <- tibble(Ev_density = rep(seq(0, 20, length.out = 300), 2),
                     fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Ev_biomass.g = fitted(ev_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
         Ev_biomass_lower = fitted(ev_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         Ev_biomass_upper = fitted(ev_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) 

ev_fung_sim_dat <- tibble(Ev_density = rep(seq(0, 20, length.out = 300), 2),
                     fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Ev_biomass.g = fitted(ev_bio_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         Ev_biomass_lower = fitted(ev_bio_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         Ev_biomass_upper = fitted(ev_bio_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) 

# coefficients
mv_coef <- fixef(mv_bio_mod)[c(2, 4), ] %>%
  as_tibble() %>%
  mutate(Parameter = c("asymptote", "rate"))
mv_fung_coef <- fixef(mv_bio_fung_mod)[c(2, 6), ] %>%
  as_tibble() %>%
  mutate(Parameter = c("asymptote", "rate"))
ev_coef <- fixef(ev_bio_mod)[c(2, 4), ] %>%
  as_tibble() %>%
  mutate(Parameter = c("asymptote", "rate"))
ev_fung_coef <- fixef(ev_bio_fung_mod)[c(2, 6), ] %>%
  as_tibble() %>%
  mutate(Parameter = c("asymptote", "rate"))

# coefficient plots
mv_coef_fig <- ggplot(mv_coef, aes(x = Parameter, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.01) +
  geom_point(size = 2) +
  ylab("Fungicide effect") +
  temp_theme

mv_fung_coef_fig <- ggplot(mv_fung_coef, aes(x = Parameter, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.01) +
  geom_point(size = 2) +
  ylab("Fungicide effect") +
  temp_theme

ev_coef_fig <- ggplot(ev_coef, aes(x = Parameter, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.01) +
  geom_point(size = 2) +
  ylab("Fungicide effect") +
  temp_theme

ev_fung_coef_fig <- ggplot(ev_fung_coef, aes(x = Parameter, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.01) +
  geom_point(size = 2) +
  ylab("Fungicide effect") +
  temp_theme

# Mv plot
mv_fig <- ggplot(mv_dat, aes(x = Mv_density, y = Mv_biomass.g)) + 
  geom_ribbon(data = mv_sim_dat, aes(ymin = Mv_biomass_lower, ymax = Mv_biomass_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 1, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegiun"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " biomass (g ", m^-1, ")", sep = ""))) +
  temp_theme

mv_fung_fig <- ggplot(mv_dat, aes(x = Mv_density, y = Mv_biomass.g)) + 
  geom_ribbon(data = mv_fung_sim_dat, aes(ymin = Mv_biomass_lower, ymax = Mv_biomass_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_fung_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 1, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegiun"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " biomass (g ", m^-1, ")", sep = ""))) +
  temp_theme

# Ev plot
ev_fig <- ggplot(ev_dat, aes(x = Ev_density, y = Ev_biomass.g)) + 
  geom_ribbon(data = ev_sim_dat, aes(ymin = Ev_biomass_lower, ymax = Ev_biomass_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = ev_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.3, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Elymus"), " biomass (g ", m^-1, ")", sep = ""))) +
  temp_theme

ev_fung_fig <- ggplot(ev_dat, aes(x = Ev_density, y = Ev_biomass.g)) + 
  geom_ribbon(data = ev_fung_sim_dat, aes(ymin = Ev_biomass_lower, ymax = Ev_biomass_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = ev_fung_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.3, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Elymus"), " biomass (g ", m^-1, ")", sep = ""))) +
  temp_theme


#### output ####
pdf("output/all_biomass_analysis_2019_density_exp.pdf")
plot_grid(mv_fig, mv_coef_fig,
          nrow = 1,
          rel_widths = c(1, 0.4))
plot_grid(ev_fig, ev_coef_fig,
          nrow = 1,
          rel_widths = c(1, 0.4))
dev.off()

pdf("output/all_biomass_analysis_greenhouse_fungicide_2019_density_exp.pdf")
plot_grid(mv_fung_fig, mv_fung_coef_fig,
          nrow = 1,
          rel_widths = c(1, 0.4))
plot_grid(ev_fung_fig, ev_fung_coef_fig,
          nrow = 1,
          rel_widths = c(1, 0.4))
dev.off()

mv_all_bio_mod <- mv_bio_mod
save(mv_all_bio_mod, file ="output/Mv_all_biomass_model_2019_density_exp.rda")
ev_all_bio_mod <- ev_bio_mod
save(ev_all_bio_mod, file ="output/Ev_all_biomass_model_2019_density_exp.rda")