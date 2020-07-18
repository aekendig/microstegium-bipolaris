##### info ####

# file: ev_seeds_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/15/20
# goal: evaluate the effects of density treatments and environmental covariates on the seed production of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(glmmTMB)

# import data
spike18 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
spike19 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
surv18 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
bgbio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")
sevy1_0 <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
sevy1_bgev_0 <- read_csv("intermediate-data/ev_background_leaf_scans_2018_density_exp.csv")
sevy2_0 <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")

# import seed model
load("./output/ev_seeds_data_processing_2018_2019_seed_spikelet_weight_mod.rda")


#### edit data ####

# disease severity
sevy1_1 <- sevy1_0 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Ev")

sevy1_bgev_1 <- sevy1_bgev_0 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}")

sevy2_1 <- sevy2_0 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Ev")

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
  left_join(sevy1_1) %>%
  left_join(sevy1_bgev_1) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"),
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
  left_join(bgbio %>%
              rename(background_biomass = biomass.g)) %>%
  left_join(sevy2_1) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"),
         background_pc_bio = background_biomass / background_density)

# combine years
seed <- fseed18 %>%
  mutate(yearf = "Year 1") %>%
  full_join(seed19 %>%
              mutate(yearf = "Year 2"))

# remove repeating 0 plots
seed2 <- seed %>%
  filter(plot > 1 | (plot == 1 & background == "Mv seedling"))



#### check datasets ####

# sample sizes
samps_19 <- seed19 %>%
  group_by(site, plot, treatment) %>%
  summarise(n = sum(!is.na(seeds))) %>%
  mutate(seeds = 405)

samps_18 <- fseed18 %>%
  group_by(site, plot, treatment) %>%
  summarise(n = sum(!is.na(seeds))) %>%
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

# template figures
temp_fig <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

temp_fig2 <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

temp_fig3 <- ggplot(plots, aes(x = treatment, y = plot)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  temp_theme +
  xlab("Treatment")

# save treatment figures
pdf("./output/ev_seeds_visualize_treatment_2018_2019_density_exp.pdf")

# density 2018
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
  
# density 2019
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

# biomass 2019
plot_grid(temp_fig2 %+%
            seed19 %+%
            aes(x = background_biomass, y = spikelet_weight.g) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig2 %+%
            seed19 %+%
            aes(x = background_biomass, y = seeds) +
            xlab("Background biomass (g)") +
            ylab("Year 2 focal seeds"),
          nrow = 2)

# Per capita biomass 2019
plot_grid(temp_fig2 %+%
            seed19 %+%
            aes(x = background_pc_bio, y = spikelet_weight.g) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig2 %+%
            seed19 %+%
            aes(x = background_pc_bio, y = seeds) +
            xlab("Per capita background biomass (g)") +
            ylab("Year 2 focal seeds"),
          nrow = 2)

# treatment within background 2018
plot_grid(temp_fig3 %+%
            fseed180 %+%
            aes(y = spikelet_weight.g) +
            ylab("Year 1 focal spikelet weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig3 %+%
            fseed180 %+%
            aes(y = seeds) +
            ylab("Year 1 focal seeds"),
          nrow = 2)

plot_grid(temp_fig3 %+%
            seed190 %+%
            aes(y = spikelet_weight.g) +
            ylab("Year 2 focal spikelet weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig3 %+%
            seed190 %+%
            aes(y = seeds) +
            ylab("Year 2 focal seeds"),
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


#### direct disease models ####

# no background plants
no_back_plant_dat = filter(seed2, plot == 1)

# divide by age
evs_no_back_plant_dat = filter(no_back_plant_dat, age == "seedling")
eva_no_back_plant_dat = filter(no_back_plant_dat, age == "adult")

# control seeds
filter(evs_no_back_plant_dat, fungicide == 0) %>%
  summarise(bio = mean(seeds))
filter(eva_no_back_plant_dat, fungicide == 0) %>%
  summarise(bio = mean(seeds))

# Ev seedling model
evs_no_back_seed_mod <- brm(data = evs_no_back_plant_dat, family = gaussian,
                           seeds ~ fungicide + (1|site) + (1|yearf),
                           prior <- c(prior(normal(6.6, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd),
                                      prior(cauchy(0, 1), class = sigma)),
                           iter = 6000, warmup = 1000, chains = 1,
                           control = list(adapt_delta = 0.9999))

# check model and add chains
summary(evs_no_back_seed_mod)
prior_summary(evs_no_back_seed_mod)
evs_no_back_seed_mod <- update(evs_no_back_seed_mod, chains = 3)
plot(evs_no_back_seed_mod)
pp_check(evs_no_back_seed_mod, nsamples = 100)

# Ev adult model
eva_no_back_seed_mod <- brm(data = eva_no_back_plant_dat, family = gaussian,
                           seeds ~ fungicide + (1|site) + (1|yearf),
                           prior <- c(prior(normal(57.7, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd),
                                      prior(cauchy(0, 1), class = sigma)),
                           iter = 6000, warmup = 1000, chains = 1,
                           control = list(adapt_delta = 0.9999))

# check model and add chains
summary(eva_no_back_seed_mod)
prior_summary(eva_no_back_seed_mod)
eva_no_back_seed_mod <- update(eva_no_back_seed_mod, chains = 3)
plot(eva_no_back_seed_mod)
pp_check(eva_no_back_seed_mod, nsamples = 100)


#### Ev seedling disease and density models ####

# separate data
evs_mv_seed_dat <- seed19 %>%
  filter(background == "Mv seedling" & plot != 1 & age == "seedling")

evs_evs_seed_dat <- seed19 %>%
  filter(background == "Ev seedling" & plot != 1  & age == "seedling")

evs_eva_seed_dat <- seed19 %>%
  filter(background == "Ev adult" & plot != 1  & age == "seedling")

# check for complete values
sum(is.na(evs_mv_seed_dat$seeds))
sum(is.na(evs_evs_seed_dat$seeds))
sum(is.na(evs_eva_seed_dat$seeds))

# Beverton-Holt function
bh_fun <- function(dat_in, a, y0){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$background_density)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# try values
mean(filter(evs_mv_seed_dat, density_level == "low")$seeds)
bh_fun(evs_mv_seed_dat, 0.04, 25)
mean(evs_evs_seed_dat$seeds)
bh_fun(evs_evs_seed_dat, 0, 13)
mean(filter(evs_eva_seed_dat, density_level == "low")$seeds)
bh_fun(evs_eva_seed_dat, 0.2, 25)

# mv background model
evs_mv_seed_mod <- brm(data = evs_mv_seed_dat, family = gaussian,
                      bf(seeds ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(25, 10), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(evs_mv_seed_mod)
prior_summary(evs_mv_seed_mod)
evs_mv_seed_mod <- update(evs_mv_seed_mod, chains = 3)
plot(evs_mv_seed_mod)
pp_check(evs_mv_seed_mod, nsamples = 100)

# ev seedling background model
evs_evs_seed_mod <- brm(data = evs_evs_seed_dat, family = gaussian,
                       bf(seeds ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(13, 10), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 6000, warmup = 1000, chains = 1,
                       control = list(adapt_delta = 0.999))

# check model and increase chains
summary(evs_evs_seed_mod)
evs_evs_seed_mod <- update(evs_evs_seed_mod, chains = 3)
plot(evs_evs_seed_mod)
pp_check(evs_evs_seed_mod, nsamples = 100)

# ev adult background model
evs_eva_seed_mod <- brm(data = evs_eva_seed_dat, family = gaussian,
                       bf(seeds ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(25, 10), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 6000, warmup = 1000, chains = 1,
                       control = list(adapt_delta = 0.999))

# check model and increase chains
summary(evs_eva_seed_mod)
evs_eva_seed_mod <- update(evs_eva_seed_mod, chains = 3)
plot(evs_eva_seed_mod)
pp_check(evs_eva_seed_mod, nsamples = 100)


#### Ev adult disease and density models ####

# separate data
eva_mv_seed_dat <- seed19 %>%
  filter(background == "Mv seedling" & plot != 1 & age == "adult")

eva_evs_seed_dat <- seed19 %>%
  filter(background == "Ev seedling" & plot != 1  & age == "adult")

eva_eva_seed_dat <- seed19 %>%
  filter(background == "Ev adult" & plot != 1  & age == "adult")

# check for complete values
sum(is.na(eva_mv_seed_dat$seeds))
sum(is.na(eva_evs_seed_dat$seeds))
sum(is.na(eva_eva_seed_dat$seeds))

# try values
mean(filter(eva_mv_seed_dat, density_level == "low")$seeds)
bh_fun(eva_mv_seed_dat, 0.03, 200)
mean(filter(eva_evs_seed_dat, density_level == "low")$seeds)
bh_fun(eva_evs_seed_dat, 0.3, 200)
mean(eva_eva_seed_dat$seeds)
bh_fun(eva_eva_seed_dat, 0, 122)

# mv background model
eva_mv_seed_mod <- brm(data = eva_mv_seed_dat, family = gaussian,
                      bf(seeds ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(200, 100), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(eva_mv_seed_mod)
prior_summary(eva_mv_seed_mod)
eva_mv_seed_mod <- update(eva_mv_seed_mod, chains = 3,
                          control = list(adapt_delta = 0.99999))
plot(eva_mv_seed_mod)
pp_check(eva_mv_seed_mod, nsamples = 100)

# ev seedling background model
eva_evs_seed_mod <- brm(data = eva_evs_seed_dat, family = gaussian,
                       bf(seeds ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(200, 100), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 6000, warmup = 1000, chains = 1,
                       control = list(adapt_delta = 0.999))

# check model and increase chains
summary(eva_evs_seed_mod)
eva_evs_seed_mod <- update(eva_evs_seed_mod, chains = 3)
plot(eva_evs_seed_mod)
pp_check(eva_evs_seed_mod, nsamples = 100)

# ev adult background model
eva_eva_seed_mod <- brm(data = eva_eva_seed_dat, family = gaussian,
                       bf(seeds ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(122, 100), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 6000, warmup = 1000, chains = 1,
                       control = list(adapt_delta = 0.999))

# check model and increase chains
summary(eva_eva_seed_mod)
eva_eva_seed_mod <- update(eva_eva_seed_mod, chains = 3)
plot(eva_eva_seed_mod)
pp_check(eva_eva_seed_mod, nsamples = 100)


#### model fit figure ####

# model fits over simulated data
seed_mv_sim_dat <- tibble(background_density = rep(seq(0, 64, length.out = 300), 2),
                         fungicide = rep(c(0, 1), each = 300)) %>%
  merge(tibble(age = c("seedling", "adult")), all = T) %>%
  as_tibble() %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = case_when(age == "seedling" ~ fitted(evs_mv_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
                                  age == "adult" ~ fitted(eva_mv_seed_mod, newdata = ., re_formula = NA)[, "Estimate"]),
         seed_lower = case_when(age == "seedling" ~ fitted(evs_mv_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
                               age == "adult" ~ fitted(eva_mv_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"]),
         seed_upper = case_when(age == "seedling" ~ fitted(evs_mv_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"],
                               age == "adult" ~ fitted(eva_mv_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"]))

seed_evs_sim_dat <- tibble(background_density = rep(seq(0, 16, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  merge(tibble(age = c("seedling", "adult")), all = T) %>%
  as_tibble() %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = case_when(age == "seedling" ~ fitted(evs_evs_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
                                  age == "adult" ~ fitted(eva_evs_seed_mod, newdata = ., re_formula = NA)[, "Estimate"]),
         seed_lower = case_when(age == "seedling" ~ fitted(evs_evs_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
                               age == "adult" ~ fitted(eva_evs_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"]),
         seed_upper = case_when(age == "seedling" ~ fitted(evs_evs_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"],
                               age == "adult" ~ fitted(eva_evs_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"]))

seed_eva_sim_dat <- tibble(background_density = rep(seq(0, 8, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  merge(tibble(age = c("seedling", "adult")), all = T) %>%
  as_tibble() %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = case_when(age == "seedling" ~ fitted(evs_eva_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
                                  age == "adult" ~ fitted(eva_eva_seed_mod, newdata = ., re_formula = NA)[, "Estimate"]),
         seed_lower = case_when(age == "seedling" ~ fitted(evs_eva_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
                               age == "adult" ~ fitted(eva_eva_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"]),
         seed_upper = case_when(age == "seedling" ~ fitted(evs_eva_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"],
                               age == "adult" ~ fitted(eva_eva_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"]))

# combine simulated data
seed_sim_dat <- seed_mv_sim_dat %>%
  mutate(background = "Mv seedling") %>%
  full_join(seed_evs_sim_dat %>%
              mutate(background = "Ev seedling")) %>%
  full_join(seed_eva_sim_dat %>%
              mutate(background = "Ev adult"))

# colors
col_pal = c("#a6611a", "#018571")

# figure
pdf("output/ev_seed_analysis_model_fits_2019_density_exp.pdf")
ggplot(seed19, aes(x = background_density, y = seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = Treatment)) +
  geom_ribbon(data = seed_sim_dat, alpha = 0.5, aes(ymin = seed_lower, ymax = seed_upper, fill = Treatment)) +
  geom_line(data = seed_sim_dat, aes(color = Treatment)) +
  facet_grid(age~background, scales = "free", switch = "both") +
  xlab("Background density") +
  ylab(expression(paste(italic(Elymus), " seeds (", plant^-1, ")", sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  temp_theme
dev.off()


#### severity and seeds ####

# divide data
# remove extra 1 plots
evs_sev_2018_dat <- filter(seed18, treatment == "water" & plant_type == "Ev seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_")) %>%
  as.data.frame()
eva_sev_2018_dat <- filter(seed18, treatment == "water" & plant_type == "Ev adult" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_")) %>%
  as.data.frame()
evs_sev_2019_dat <- filter(seed19, treatment == "water" & plant_type == "Ev seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_")) %>%
  as.data.frame()
eva_sev_2019_dat <- filter(seed19, treatment == "water" & plant_type == "Ev adult" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_")) %>%
  as.data.frame()

# check for density-severity correlations
for(i in 20:22){
  evs_sev_2018_dat$severity = evs_sev_2018_dat[,i]
  print(cor.test(evs_sev_2018_dat$background_density, evs_sev_2018_dat$severity))
  
  eva_sev_2018_dat$severity = eva_sev_2018_dat[,i]
  print(cor.test(eva_sev_2018_dat$background_density, eva_sev_2018_dat$severity))
}
# none significantly correlated

for(i in 22:25){
  evs_sev_2019_dat$severity = evs_sev_2019_dat[,i]
  print(cor.test(evs_sev_2019_dat$background_density, evs_sev_2019_dat$severity))
  
  eva_sev_2019_dat$severity = eva_sev_2019_dat[,i]
  print(cor.test(eva_sev_2019_dat$background_density, eva_sev_2019_dat$severity))
}

# biomass-severity relationships
evs_seed_sev_jul_18_mod <- glmmTMB(seeds ~ severity_jul + (1|site/plot), family = gaussian, data = evs_sev_2018_dat)
summary(evs_seed_sev_jul_18_mod) # not sig
evs_seed_sev_lau_18_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site/plot), family = gaussian, data = evs_sev_2018_dat)
summary(evs_seed_sev_lau_18_mod) # not sig

eva_seed_sev_jul_18_mod <- glmmTMB(seeds ~ severity_jul + (1|site/plot), family = gaussian, data = eva_sev_2018_dat) # convergence issue
eva_seed_sev_jul_18_mod <- glmmTMB(seeds ~ severity_jul + (1|site_plot), family = gaussian, data = eva_sev_2018_dat) 
summary(eva_seed_sev_jul_18_mod) # not sig
eva_seed_sev_lau_18_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site/plot), family = gaussian, data = eva_sev_2018_dat) # convergence issue
eva_seed_sev_lau_18_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site_plot), family = gaussian, data = eva_sev_2018_dat)
summary(eva_seed_sev_lau_18_mod) # not sig

evs_seed_sev_jun_19_mod <- glmmTMB(seeds ~ severity_jun + (1|site/plot), family = gaussian, data = evs_sev_2019_dat)
summary(evs_seed_sev_jun_19_mod) # not sig
evs_seed_sev_jul_19_mod <- glmmTMB(seeds ~ severity_jul + (1|site/plot), family = gaussian, data = evs_sev_2019_dat)
summary(evs_seed_sev_jul_19_mod) # not sig
evs_seed_sev_eau_19_mod <- glmmTMB(seeds ~ severity_early_aug + (1|site/plot), family = gaussian, data = evs_sev_2019_dat)
summary(evs_seed_sev_eau_19_mod) # not sig
evs_seed_sev_lau_19_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site/plot), family = gaussian, data = evs_sev_2019_dat)
summary(evs_seed_sev_lau_19_mod) # not sig

eva_seed_sev_jun_19_mod <- glmmTMB(seeds ~ severity_jun + (1|site/plot), family = gaussian, data = eva_sev_2019_dat) # convergence issue
eva_seed_sev_jun_19_mod <- glmmTMB(seeds ~ severity_jun + (1|site_plot), family = gaussian, data = eva_sev_2019_dat) 
summary(eva_seed_sev_jun_19_mod) # not sig
eva_seed_sev_jul_19_mod <- glmmTMB(seeds ~ severity_jul + (1|site/plot), family = gaussian, data = eva_sev_2019_dat) # convergence issue
eva_seed_sev_jul_19_mod <- glmmTMB(seeds ~ severity_jul + (1|site_plot), family = gaussian, data = eva_sev_2019_dat) # convergence issue
eva_seed_sev_eau_19_mod <- glmmTMB(seeds ~ severity_early_aug + (1|site/plot), family = gaussian, data = eva_sev_2019_dat) # convergence issue
eva_seed_sev_eau_19_mod <- glmmTMB(seeds ~ severity_early_aug + (1|site_plot), family = gaussian, data = eva_sev_2019_dat) # convergence issue
eva_seed_sev_lau_19_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site/plot), family = gaussian, data = eva_sev_2019_dat) # convergence issue
eva_seed_sev_lau_19_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site_plot), family = gaussian, data = eva_sev_2019_dat)
summary(eva_seed_sev_lau_19_mod) # not sig


#### output ####

save(evs_no_back_seed_mod, file = "output/ev_seedling_seeds_no_background_model_2018_2019_density_exp.rda")
save(eva_no_back_seed_mod, file = "output/ev_adult_seeds_no_background_model_2018_2019_density_exp.rda")
save(evs_mv_seed_mod, file = "output/ev_seedling_seeds_mv_background_model_2019_density_exp.rda")
save(evs_evs_seed_mod, file = "output/ev_seedling_seeds_ev_seedling_background_model_2019_density_exp.rda")
save(evs_eva_seed_mod, file = "output/ev_seedling_seeds_ev_adult_background_model_2019_density_exp.rda")
save(eva_mv_seed_mod, file = "output/ev_adult_seeds_mv_background_model_2019_density_exp.rda")
save(eva_evs_seed_mod, file = "output/ev_adult_seeds_ev_seedling_background_model_2019_density_exp.rda")
save(eva_eva_seed_mod, file = "output/ev_adult_seeds_ev_adult_background_model_2019_density_exp.rda")
