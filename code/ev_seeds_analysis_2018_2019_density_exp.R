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

# import data
spike18 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
spike19 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
surv18 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")


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
         fungicide = recode(treatment, water = 0, fungicide = 1))

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


#### 2019 seed models ####


### manual model fitting ###

# sigmoid Beverton-Holt
sbh_fun <- function(dat_in, r, a, d){
  
  # extract values
  xmin = min(dat_in$background_density)
  xmax = max(dat_in$background_density)
  y0 = filter(dat_in, fungicide == 0) %>%
          summarise(mean_seed = mean(seeds)) %>%
          as.numeric() %>%
          round()
  print(y0)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 + x^(d-1) / (1/r + a * x^d))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# separate data
ea_ea_dat <- filter(seed19, age == "adult" & background == "Ev adult")

# try values
sbh_fun(ea_ea_dat, 30, 0.005, 2)


### Ev adult seeds, Ev adult background ###

# ea_seeds_ea_bg_mod <- brm(data = ea_ea_dat, family = gaussian,
#                           bf(seeds ~ s0 + background_density / (k + alpha * background_density^2),
#                              s0 ~ fungicide + (1|site),
#                              k ~ fungicide,
#                              alpha ~ fungicide,
#                              nl = T),
#                           prior <- c(prior(normal(97, 10), nlpar = "s0"),
#                                      prior(normal(0.03, 1), nlpar = "k"),
#                                      prior(normal(0.005, 1), nlpar = "alpha")),
#                           iter = 6000, warmup = 1000, chains = 3, cores = 2,
#                           control = list(adapt_delta = 0.9999, max_treedepth = 15))
# 1178 divergent transitions

# remove site random effect

ea_seeds_ea_bg_mod <- brm(data = ea_ea_dat, family = gaussian,
                          bf(seeds ~ s0 + background_density / (k + alpha * background_density^2),
                             s0 ~ fungicide,
                             k ~ fungicide,
                             alpha ~ fungicide,
                             nl = T),
                          prior <- c(prior(normal(97, 10), nlpar = "s0"),
                                     prior(normal(0.03, 1), nlpar = "k"),
                                     prior(normal(0.005, 1), nlpar = "alpha")),
                          iter = 6000, warmup = 1000, chains = 3, cores = 2,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15))
# 1006 divergent transitions




