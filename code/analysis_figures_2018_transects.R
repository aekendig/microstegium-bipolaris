
##### info ####

# file: analysis_figures_2018_transects.R
# author: Amy Kendig
# date last edited: 10/27/21
# goal: effect of distance on infection


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)
library(glmmTMB)
library(DHARMa)

# import data
sev <- read_csv("intermediate-data/mv_transect_leaf_scans_2018_density_exp.csv")
prev_eau <- read_csv("data/prevalence_early_aug_2018_transects.csv")
prev_lau <- read_csv("data/prevalence_late_aug_2018_transects.csv")
prev_sep <- read_csv("data/prevalence_sep_2018_transects.csv")
soil <- read_csv("data/soil_moisture_oct_2018_transects.csv")
canopy <- read_csv("data/canopy_cover_nov_2018_transects.csv")
mv_edge <- read_csv("data/plot_edge_mv_weight_jul_2018_transects.csv")


#### edit data ####

# canopy cover
canopy2 <- canopy %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                       type == "c" ~ (count * 1.04) / 100)) %>%
  select(site, transect, distance.m, canopy_cover.prop)

# soil moisture
soil2 <- soil %>%
  mutate(soil_moisture.prop = soil_moisture.vwc/100) %>%
  select(site, transect, distance.m, soil_moisture.prop)

# mv_edge
unique(mv_edge$process_notes)

mv_edge2 <- mv_edge %>%
  mutate(mv_tot.g = mv.g + mv_inf.g,
         mv_inf.prop = mv_inf.g / mv_tot.g) %>%
  select(site, transect, distance.m, mv_tot.g, mv_inf.prop)

# severity
sev2 <- sev %>%
  mutate(severity = lesion_area.pix / leaf_area.pix,
         month = "September") %>%
  select(month, site, transect, distance.m, severity)

# early aug prev
filter(prev_eau, infected > live)

prev_eau2 <- prev_eau %>%
  select(-date) %>%
  mutate(prev = infected / live,
         uninfected = live - infected,
         month = "early August")

# late aug prev
filter(prev_lau, infected > live)

prev_lau2 <- prev_lau %>%
  select(-date) %>%
  mutate(prev = infected / live,
         uninfected = live - infected,
         month = "late August")

# sep prev
prev_sep2 <- prev_sep %>%
  filter(!is.na(infected)) %>% # remove dead plants
  group_by(site, transect, distance.m) %>%
  summarize(live = n_distinct(plant),
            infected = sum(infected),
            leaves_infected = sum(leaves_infected),
            leaves = sum(leaves_total)) %>%
  ungroup() %>%
  mutate(prev = infected / live,
         uninfected = live - infected,
         month = "September")

# combine data
transects <- prev_eau2 %>%
  full_join(prev_lau2) %>%
  full_join(prev_sep2) %>%
  full_join(sev2) %>%
  full_join(canopy2) %>%
  full_join(soil2) %>%
  full_join(mv_edge2) %>%
  mutate(site_transect = paste(site, transect, sep = "_"))


#### initial visualizations ####

# correlations among covariates
transects %>%
  select(canopy_cover.prop, soil_moisture.prop, mv_tot.g, mv_inf.prop) %>%
  unique() %>%
  ggpairs()
# higher canopy cover --> lower mv total biomass

# prevalence
ggplot(transects, aes(x = distance.m, y = prev, color = month)) +
  geom_point()

# severity
ggplot(transects, aes(x = distance.m, y = severity)) +
  geom_point()


#### prevalence model ####

# fit model without covariates
prev_mod <- glmmTMB(cbind(infected, uninfected) ~ month * distance.m + (1|site_transect),
                    data = transects,
                    family = "binomial")
summary(prev_mod)
plot(simulateResiduals(prev_mod))


#### prevalence change figure ####

# simulate data
simdat <- transects %>%
  select(month) %>%
  unique() %>%
  expand_grid(tibble(distance.m = seq(0, 10, length.out = 50))) %>%
  mutate(site_transect = NA) %>%
  mutate(prev = predict(prev_mod, newdata = ., re.form = NA, type = "response"),
         prev_se = predict(prev_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

# figure
pdf("output/mv_infection_prevalence_2018_transects.pdf", width = 3, height = 3)
ggplot(simdat, aes(x = distance.m, y = prev, color = month)) +
  geom_ribbon(aes(ymin = prev - prev_se, ymax = prev + prev_se, fill = month), color = NA, alpha = 0.5) +
  geom_line() +
  stat_summary(data = transects, geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(data = transects, geom = "point", size = 2, fun = "mean") +
  labs(x = expression(paste("Distance from infected ", italic("M. vimineum"), " (m)", sep = "")),
       y = expression(paste("Proportion of sentinel ", italic("M. vimineum"), " infected", sep = ""))) +
  scale_fill_viridis_d(direction = -1, labels = c("early\nAugust", "late\nAugust", "September")) +
  scale_color_viridis_d(direction = -1, labels = c("early\nAugust", "late\nAugust", "September")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10, hjust = 1),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.18, 0.16),
        legend.key.size = unit(0.5, 'cm'))
dev.off()
