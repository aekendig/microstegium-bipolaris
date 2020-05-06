##### info ####

# file: ev_seeds_data_processing_2019
# author: Chris Wojan, Amy Kendig
# date last edited: 4/24/20
# goal: edit Ev seed data and check for errors
# background: spikelet counts for all samples, seed counts for a subset

#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(MASS)
library(tidyverse)
library(lubridate)

# import data
spikelets <- read_csv("./data/ev_spikelets_2019_density_exp.csv")
biomass <- read_csv("./data/ev_biomass_seeds_oct_2019_density_exp.csv")
seeds <- read_csv("./data/ev_seed_subset_2019_density_exp.csv")
plants_jun <- read_csv("./data/focal_disease_jun_2019_density_exp.csv")
plants_jul <- read_csv("./data/focal_disease_jul_2019_density_exp.csv")
plants_eau <- read_csv("./data/focal_disease_early_aug_2019_density_exp.csv")
plants_lau <- read_csv("./data/focal_disease_late_aug_2019_density_exp.csv")
plants_sep <- read_csv("./data/focal_disease_sep_2019_density_exp.csv")
dat18 <- read_csv("./intermediate-data/ev_spikelets_seeds_2018.csv")


#### combine spikelet and plant data ####

# notes
unique(spikelets$spikelet_notes)
unique(plants_jun$field_notes) # no seed collection data
unique(plants_jul$field_notes)
unique(plants_eau$field_notes)
unique(plants_lau$field_notes)
unique(plants_sep$field_notes)

# seed collection dates
sort(unique(spikelets$collect_date))
# biomass data (Oct dates) have no dates

# spikelet columns
# combine by collection date (sometimes multiple bags)
spikelets2 <- spikelets %>%
  mutate(treatment = str_extract(plot, "[aA-zZ]+") %>% recode("W" = "water", "F" = "fungicide"),
         plot = str_extract(plot, "[0-9]+") %>% as.numeric(),
         sp = "Ev",
         ID = str_remove(plant, "EV"),
         seeds_green = case_when(grepl("green", spikelet_notes) == T ~ 1,
                                 grepl("GREEN", spikelet_notes) == T ~ 1,
                                 TRUE ~ 0),
         collect_date = ifelse(collect_date > 20190928, 20191022, collect_date)) %>%
  group_by(site, plot, treatment, sp, ID, collect_date) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g),
            spikelet_notes = paste(spikelet_notes, collapse = "; "),
            spikelet_bags_green = sum(seeds_green))

# collection dates
collect <- select(plants_jul, date:ID, seeds_collected) %>%
  full_join(select(plants_eau, date:ID, seeds_collected)) %>%
  full_join(select(plants_lau, date:ID, seeds_collected)) %>%
  full_join(select(plants_sep, date:ID, seeds_collected)) %>%
  full_join(select(biomass, site:ID, seeds_collected) %>%
              mutate(date = 20191022)) %>%
  filter(seeds_collected == 1) %>%
  rename(collect_date = date)

# combine data
spike_collect <- full_join(spikelets2, collect)  


#### check data for mislabel ####

# check for plants where seeds were recorded as collected, but no bag was weighed
filter(spike_collect, is.na(spikelet_weight.g))
# D3 1W Ev A on 7/29/19
# D3 7F Ev 2 on 9/26/19
# D4 4W Ev 2 on 8/28/19

# check for bags that were weighed, but seeds weren't recorded as collected
filter(spike_collect, is.na(seeds_collected))
# D2 10W Ev 2 on 10/22/19
# D3 2F Ev A on 7/3/19
# D3 7F Ev 2 on 9/25/19 - this is the same as above, the date is off one day
# D4 6F Ev 2 on 9/24/19
# none of the other seem like mis-labels

# update spikelet data
# add in dummy variables for plotting sampling density
spikelets3 <- spikelets2 %>%
  mutate(collect_date = case_when(site == "D3" & plot == 7 & treatment == "fungicide" & ID == "2" & collect_date == 20190925 ~ 20190926,
                                  TRUE ~ collect_date),
         lower = -7,
         upper = -3,
         Year = "2019")


#### combine seed and spikelet data ####

# notes
unique(seeds$processing_notes)
sort(unique(seeds$collect_date))

# select rows with dummy seeds
# change collect dates for October
# update age column
seeds2 <- seeds %>%
  filter(is.na(processing_notes) | !(processing_notes %in% c("lost two seeds; this row has weight of reduced seed count", "lost one seed; this row has weight of reduced seed count"))) %>%
  mutate(collect_date = ifelse(collect_date > 20190928, 20191022, collect_date),
         age = recode(age, A = "adult", S = "seedling")) %>%
  group_by(site, plot, treatment, ID, age, collect_date) %>%
  summarise(seeds = sum(seeds),
            seed_weight = sum(seed_weight)) %>%
  ungroup()

# check seed data for this plant
filter(seeds2, site == "D3" & plot == 7 & treatment == "fungicide" & ID == 2)
# not included

# combine data
# extract month
spike_seed = seeds2 %>%
  left_join(spikelets3) %>%
  mutate(collect_month = as.character(collect_date) %>% 
           as.Date(., format = "%Y%m%d") %>% 
           month(label = T))

# check for missing data
filter(spike_seed, is.na(spikelet_weight.g))


#### visualize seed relationships ####

# seeds and spikelet_weight
spike_seed %>%
  ggplot(aes(x = spikelet_weight.g, y = seeds)) +
  geom_point(size = 2, aes(colour = age, shape = collect_month)) +
  stat_smooth(method = "lm")

# seed_weight and spikelet_weight
spike_seed %>%
  ggplot(aes(x = spikelet_weight.g, y = seed_weight)) +
  geom_point(size = 2, aes(colour = age, shape = collect_month)) +
  stat_smooth(method = "lm")

# spikelet weight distribution
spikelets3 %>%
  ggplot(aes(x = spikelet_weight.g)) +
  geom_histogram()
# a handful are that large


#### regressions ####

# examine data
mean(spike_seed$seeds)
var(spike_seed$seeds)
ggplot(spike_seed, aes(x = seeds)) +
  geom_histogram(binwidth = 1)

# poisson
mod_19a <- glm(seeds ~ spikelet_weight.g, data = spike_seed, family = poisson)
summary(mod_19a)
plot(mod_19a)

# neg binomial
mod_19b <- glm.nb(seeds ~ spikelet_weight.g, data = spike_seed)
summary(mod_19b)
# dispersion parameter is close to 1
plot(mod_19b)
# much smaller std deviance residuals

# compare
AIC(mod_19a, mod_19b)
# neg binomial better

# simulation
sim_dat = tibble(spikelet_weight.g = seq(0, max(spike_seed$spikelet_weight.g), length.out = 100))
sim_dat$yp = predict(mod_19a, newdata = sim_dat, type = "response")
sim_dat$yn = predict(mod_19b, newdata = sim_dat, type = "response")

ggplot(spike_seed, aes(x = spikelet_weight.g)) +
  geom_point(aes(y = seeds)) +
  geom_line(data = sim_dat, aes(y = yp), color = "blue") +
  geom_line(data = sim_dat, aes(y = yn), color = "green")
# Poisson looks way better

# Poisson, origin
mod_19c <- glm(seeds ~ 0 + spikelet_weight.g, data = spike_seed, family = poisson)
summary(mod_19c)
plot(mod_19c)

AIC(mod_19a, mod_19c)
# a is a better fit

# simulation
sim_dat$yp0 = predict(mod_19c, newdata = sim_dat, type = "response")

ggplot(spike_seed, aes(x = spikelet_weight.g)) +
  geom_point(aes(y = seeds)) +
  geom_line(data = sim_dat, aes(y = yp), color = "blue") +
  geom_line(data = sim_dat, aes(y = yp0), color = "purple")
# the 0 intercept forces it too low

# linear models
mod_19d <- lm(seeds ~ spikelet_weight.g, data = spike_seed)
mod_19e <- lm(seeds ~ 0 + spikelet_weight.g, data = spike_seed)

# compare
AIC(mod_19a, mod_19b, mod_19d, mod_19e)
# linear are better, no difference with intercept

# simulation
sim_dat$yl = predict(mod_19d, newdata = sim_dat)
sim_dat$yl0 = predict(mod_19e, newdata = sim_dat)

ggplot(spike_seed, aes(x = spikelet_weight.g)) +
  geom_point(aes(y = seeds)) +
  geom_line(data = sim_dat, aes(y = yp), color = "blue") +
  geom_line(data = sim_dat, aes(y = yp0), color = "purple") +
  geom_line(data = sim_dat, aes(y = yl), color = "red") +
  geom_line(data = sim_dat, aes(y = yl0), color = "orange")

# predicted values
pred_dat <- sim_dat %>%
  select(spikelet_weight.g, yl0) %>%
  rename(seeds = yl0) %>%
  mutate(seeds.se = predict(mod_19e, newdata = ., se.fit = T)$se.fit)

# figure
pred_plot <- spike_seed %>%
  ggplot(aes(x = spikelet_weight.g, y = seeds)) +
  geom_ribbon(data = pred_dat, alpha = 0.5, aes(ymin = seeds - seeds.se, ymax = seeds + seeds.se)) +
  geom_line(data = pred_dat, linetype = "dashed") +
  geom_point(size = 2, aes(colour = age, shape = treatment)) +
  theme_bw() +
  xlab("spikelet weight (g)")



#### combine data from both years ####

# edit 2018 spikelet data
spikelets18 <- dat18 %>%
  filter(experiment == "density") %>%
  mutate(spikelet_weight.g = spikelet_weight.mg / 1000,
         lower = -7,
         upper = -3,
         Year = "2018")

# edit 2018 seed data  
dat18b <- spikelets18 %>%
  filter(!is.na(seeds)) %>%
  select(site, plot, treatment, sp, ID, age, collect_date, spikelet_weight.g, seeds)

# linear model for 2018
mod_18 <- lm(seeds ~ 0 + spikelet_weight.g, data = dat18b)
summary(mod_18)

# combined data
dat <- spike_seed %>%
  select(site, plot, treatment, sp, ID, age, collect_date, spikelet_weight.g, seeds, spikelet_bags_green) %>%
  full_join(dat18b) %>%
  mutate(year = as.character(collect_date) %>% 
           as.Date(., format = "%Y%m%d") %>% 
           year(),
         Year = as.factor(year))

# prediction data
pred_dat_comb <- pred_dat %>%
  mutate(Year = "2019") %>%
  full_join(tibble(spikelet_weight.g = seq(min(dat18b$spikelet_weight.g), max(dat18b$spikelet_weight.g), length.out = 100),
                        Year = rep("2018", each = 100)) %>%
  mutate(seeds = predict(mod_18, newdata = .),
         seeds.se = predict(mod_18, newdata = ., se.fit = T)$se.fit))

# visualize
pred_plot_comb <- dat %>%
  ggplot(aes(x = spikelet_weight.g, y = seeds, color = Year)) +
  geom_ribbon(data = pred_dat_comb, alpha = 0.5, aes(ymin = seeds - seeds.se, ymax = seeds + seeds.se)) +
  geom_line(data = pred_dat_comb, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  theme_bw() +
  xlab("spikelet weight (g)")

# coefficients from linear models
summary(mod_18)
summary(mod_19b)

# fit model to both years
mod_both <- lm(seeds ~ 0 + spikelet_weight.g, data = dat)
summary(mod_both)

# predict values
pred_dat_both <- tibble(spikelet_weight.g = seq(0, max(spikelets18$spikelet_weight.g), length.out = 100)) %>%
  mutate(seeds = predict(mod_both, newdata = .),
         seeds.se = predict(mod_both, newdata = ., se.fit = T)$se.fit)

# visualize
pred_plot_both <- pred_dat_both %>%
  ggplot(aes(x = spikelet_weight.g)) +
  geom_ribbon(aes(y = seeds, ymin = seeds - seeds.se, ymax = seeds + seeds.se), alpha = 0.5) +
  geom_line(aes(y = seeds)) +
  geom_point(data = dat, alpha = 0.75, aes(y = seeds, color = Year)) +
  geom_linerange(data = spikelets3, aes(ymin = lower, ymax = upper, color = Year), alpha = 0.5) +
  geom_linerange(data = spikelets18, aes(ymin = lower, ymax = upper, color = Year), alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8)) +
  xlab("spikelet weight (g)")


#### predict seeds ####

# 2019 data
spikelets3_out <- spikelets3 %>%
  select(-c(lower:Year)) %>%
  left_join(seeds2) %>%
  mutate(seeds = case_when(is.na(seeds) ~ spikelet_weight.g * coef(mod_both),
                           TRUE ~ seeds))

# check
pred_plot_both +
  geom_point(data = spikelets3_out, aes(y = seeds))

# 2018 data
spikelets18_out <- spikelets18 %>%
  select(-c(lower:Year)) %>%
  mutate(seeds = case_when(is.na(seeds) ~ spikelet_weight.g * coef(mod_both),
                           TRUE ~ seeds))

# check
pred_plot_both +
  geom_point(data = spikelets18_out, aes(y = seeds))


#### save output ####

# seeds prediction relationship
pdf("./output/ev_seeds_data_processing_2019_seeds_spikelet_weight.pdf", width = 5, height = 4)
pred_plot
dev.off()

pdf("./output/ev_seeds_data_processing_2018_2019_seeds_spikelet_weight.pdf", width = 5, height = 4)
pred_plot_comb
dev.off()

pdf("./output/ev_seeds_data_processing_2019_seeds_spikelet_weight_combined.pdf", width = 4, height = 4)
pred_plot_both
dev.off()

# model
save(mod_both, file = "./output/ev_seeds_data_processing_2018_2019_seed_spikelet_weight_mod.rda")

# data
write_csv(spikelets3_out, "./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
write_csv(spikelets18_out, "./intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
