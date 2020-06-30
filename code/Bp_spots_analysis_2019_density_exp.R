##### info ####

# file: Bp_spots_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 6/26/20
# goal: analyze occurrence of Bp spots


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(DHARMa)
library(MASS)
library(tidyverse)
library(glmmTMB)

# import all raw data files
dt_may <- read_csv("data/ev_disease_may_2019_density_exp.csv")
dt_jun <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
dt_jul <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
dt_early_aug <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
dt_late_aug <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")
dt_sep <- read_csv("data/focal_disease_sep_2019_density_exp.csv")
edge_dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
temp_hum <- read_csv("intermediate-data/temp_humidity_daily_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


#### edit data ####

# follow the same protocols used in focal_severity_analysis_2019_density_exp for temp_hum, bio, and edge dat

# temp/hum data
# take average of daily values over the month leading up to sampling
temp_hum2 <- temp_hum %>%
  mutate(month = case_when(day < as.Date("2019-07-29") ~ "early_aug",
                           day >= as.Date("2019-07-29") & day < as.Date("2019-08-28") ~ "late_aug",
                           day >= as.Date("2019-08-28") & day < as.Date("2019-09-24") ~ "sep",
                           TRUE ~ "oct")) %>%
  group_by(site, plot, treatment, month) %>%
  summarise(temp_avg = mean(temp_avg),
            temp_min = mean(temp_min),
            temp_max = mean(temp_max),
            hum_avg = mean(hum_avg),
            hum_min = mean(hum_min),
            hum_max = mean(hum_max),
            hum_dur = mean(hum_dur)) %>%
  ungroup() %>%
  filter(month != "oct")

# edge data
# use data collected one month before leaf severity sampling
edge_dat2 <- edge_dat %>%
  filter(month != "sep") %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         month = recode(month, "may" = "jun", "jun" = "jul", "jul" = "early_aug", "early_aug" = "late_aug", "late_aug" = "sep")) %>%
  select(month, site, plot, treatment, edge_severity)

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
  group_by(site, treatment, plot, sp) %>%
  summarise(biomass.g = sum(biomass.g)) %>%
  ungroup() %>%
  spread(key = sp, value = biomass.g) %>%
  mutate(total_biomass.g = Ev + Mv) %>%
  rename("Ev_biomass.g" = "Ev", "Mv_biomass.g" = "Mv")

# density
dens_dat <- plots_simple %>%
  select(plot, background_sp, background_density) %>%
  unique() %>%
  filter(background_density > 0) %>%
  spread(key = background_sp, value = background_density) %>%
  full_join(tibble(plot = 1, Ev = 0, Mv = 0)) %>%
  mutate(Ev = replace_na(Ev, 0) + 4,
         Mv = replace_na(Mv, 0) + 3,
         total_density = Ev + Mv,
         Ev_present = case_when(plot > 4 ~ 1,
                                TRUE ~ 0),
         Mv_present = case_when(plot > 1 & plot < 5 ~ 1,
                                TRUE ~ 0),
         bg_present = case_when(plot > 1 ~ 1,
                                TRUE ~ 0)) %>%
  rename(Ev_density = Ev, Mv_density = Mv)

# dataset function
dat_fun <- function(dat_in, mth){
  
  dat_out <- dat_in %>%
    filter(sp == "Ev" & !is.na(Bp_spots)) %>%
    mutate(month = mth,
           fungicide = case_when(treatment == "water" ~ 0,
                                 TRUE ~ 1),
           exp_plot = paste0(plot, toupper(substr(treatment, 1, 1))),
           site_exp_plot = paste(site, exp_plot, sep = " ")) %>%
    left_join(temp_hum2) %>%
    left_join(edge_dat2) %>%
    left_join(bio) %>%
    left_join(dens_dat)
  
  return(dat_out)
}

# datasets
may_dat <- dat_fun(dt_may, "may")
may_ev_dat <- dat_fun(filter(dt_may, ID %in% c("1", "2", "3", "A")), "may")
jun_ev_dat <- dat_fun(dt_jun, "jun")
jul_ev_dat <- dat_fun(dt_jul, "jul")
eau_ev_dat <- dat_fun(dt_early_aug, "early_aug")
lau_ev_dat <- dat_fun(dt_late_aug, "late_aug")
sep_ev_dat <- dat_fun(dt_sep, "sep")


#### model selection for water treatment ####

# include humidity data (only available for control plots)

# divide data
eau_ev_water_dat <- filter(eau_ev_dat, treatment == "water")
lau_ev_water_dat <- filter(lau_ev_dat, treatment == "water")
sep_ev_water_dat <- filter(sep_ev_dat, treatment == "water")

## early August ##
eau_ev_totb_water_mod <- glmmTMB(Bp_spots ~ total_biomass.g + edge_severity + hum_avg + (1|site/plot), data = eau_ev_water_dat, family = "binomial")
eau_ev_sepb_water_mod <- glmmTMB(Bp_spots ~ Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg + (1|site/plot), data = eau_ev_water_dat, family = "binomial")
eau_ev_totd_water_mod <- glmmTMB(Bp_spots ~ total_density + edge_severity + hum_avg + (1|site/plot), data = eau_ev_water_dat, family = "binomial")
eau_ev_sepd_water_mod <- glmmTMB(Bp_spots ~ Ev_density + Mv_density + edge_severity + hum_avg + (1|site/plot), data = eau_ev_water_dat, family = "binomial")
eau_ev_totp_water_mod <- glmmTMB(Bp_spots ~ bg_present + edge_severity + hum_avg + (1|site/plot), data = eau_ev_water_dat, family = "binomial")
eau_ev_sepp_water_mod <- glmmTMB(Bp_spots ~ Ev_present + Mv_present + edge_severity + hum_avg + (1|site/plot), data = eau_ev_water_dat, family = "binomial")

# compare different plant models
AIC(eau_ev_totb_water_mod, eau_ev_sepb_water_mod, eau_ev_totd_water_mod, eau_ev_sepd_water_mod, eau_ev_totp_water_mod, eau_ev_sepp_water_mod)
# sep biomass is the lowest
summary(eau_ev_sepb_water_mod)
# all significant
stepAIC(eau_ev_sepb_water_mod)
# keep full model
plot(simulateResiduals(eau_ev_sepb_water_mod))


## late August ##
lau_ev_totb_water_mod <- glmmTMB(Bp_spots ~ total_biomass.g + edge_severity + hum_avg + (1|site/plot), data = lau_ev_water_dat, family = "binomial")
lau_ev_sepb_water_mod <- glmmTMB(Bp_spots ~ Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg + (1|site/plot), data = lau_ev_water_dat, family = "binomial")
lau_ev_totd_water_mod <- glmmTMB(Bp_spots ~ total_density + edge_severity + hum_avg + (1|site/plot), data = lau_ev_water_dat, family = "binomial")
lau_ev_sepd_water_mod <- glmmTMB(Bp_spots ~ Ev_density + Mv_density + edge_severity + hum_avg + (1|site/plot), data = lau_ev_water_dat, family = "binomial")
lau_ev_totp_water_mod <- glmmTMB(Bp_spots ~ bg_present + edge_severity + hum_avg + (1|site/plot), data = lau_ev_water_dat, family = "binomial")
lau_ev_sepp_water_mod <- glmmTMB(Bp_spots ~ Ev_present + Mv_present + edge_severity + hum_avg + (1|site/plot), data = lau_ev_water_dat, family = "binomial")

# compare different plant models
AIC(lau_ev_totb_water_mod, lau_ev_sepb_water_mod, lau_ev_totd_water_mod, lau_ev_sepd_water_mod, lau_ev_totp_water_mod, lau_ev_sepp_water_mod)
# tot biomass is the lowest
summary(lau_ev_totb_water_mod)
# none significant
stepAIC(lau_ev_totb_water_mod)
# keep humidity only
lau_ev_totb_water_mod2 <- update(lau_ev_totb_water_mod, .~. - total_biomass.g - edge_severity)
summary(lau_ev_totb_water_mod2)
# not sig
plot(simulateResiduals(lau_ev_totb_water_mod2))


## September ##
sep_ev_totb_water_mod <- glmmTMB(Bp_spots ~ total_biomass.g + edge_severity + hum_avg + (1|site/plot), data = sep_ev_water_dat, family = "binomial")
# can't converge
sep_ev_totb_water_mod <- glmmTMB(Bp_spots ~ total_biomass.g + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_sepb_water_mod <- glmmTMB(Bp_spots ~ Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_totd_water_mod <- glmmTMB(Bp_spots ~ total_density + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_sepd_water_mod <- glmmTMB(Bp_spots ~ Ev_density + Mv_density + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_totp_water_mod <- glmmTMB(Bp_spots ~ bg_present + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_sepp_water_mod <- glmmTMB(Bp_spots ~ Ev_present + Mv_present + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "binomial")

# compare different plant models
AIC(sep_ev_totb_water_mod, sep_ev_sepb_water_mod, sep_ev_totd_water_mod, sep_ev_sepd_water_mod, sep_ev_totp_water_mod, sep_ev_sepp_water_mod)
# tot biomass is the lowest
summary(sep_ev_totb_water_mod)
# none significant
stepAIC(sep_ev_totb_water_mod)
# keep humidity only
sep_ev_totb_water_mod2 <- update(sep_ev_totb_water_mod, .~. - total_biomass.g - edge_severity)
summary(sep_ev_totb_water_mod2)
# not sig
plot(simulateResiduals(sep_ev_totb_water_mod2))


#### model selection for both treatments ####

## May all Ev ##
may_ev_all_totb_mod <- glmmTMB(Bp_spots ~ fungicide * total_biomass.g + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g) + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_totd_mod <- glmmTMB(Bp_spots ~ fungicide * total_density + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density) + (1|site/exp_plot), data = may_dat, family = "binomial")
# convergence issue
summary(may_ev_all_sepd_mod) # remove site
may_ev_all_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density) + (1|site_exp_plot), data = may_dat, family = "binomial")
# convergence issue
filter(may_dat, Mv_present == 1 & treatment == "fungicide")
# only 0's
may_ev_all_totp_mod <- glmmTMB(Bp_spots ~ fungicide * bg_present + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present) + (1|site/exp_plot), data = may_dat, family = "binomial")
# convergence issue

# compare different plant models
AIC(may_ev_all_totb_mod, may_ev_all_sepb_mod, may_ev_all_totd_mod, may_ev_all_totp_mod)
# separate biomass is the lowest
summary(may_ev_all_sepb_mod)
stepAIC(may_ev_all_sepb_mod)
# remove Mv_biomass and interactions
may_ev_all_sepb_mod2 <- update(may_ev_all_sepb_mod, .~. - Mv_biomass.g - fungicide:Mv_biomass.g)
summary(may_ev_all_sepb_mod2)
# none sig
plot(simulateResiduals(may_ev_all_sepb_mod2))


## May Ev ##
may_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * total_biomass.g + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
summary(may_ev_totb_mod)
may_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * total_biomass.g + (1|site_exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
may_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g) + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
may_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * total_density + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
may_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density) + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
may_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * bg_present + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
may_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present) + (1|site/exp_plot), data = may_ev_dat, family = "binomial")

# compare different plant models
AIC(may_ev_sepb_mod, may_ev_totd_mod)
# total density
summary(may_ev_totd_mod)
# none
stepAIC(may_ev_totd_mod)
# keep full model
plot(simulateResiduals(may_ev_totd_mod))
# sig dev


## June Ev ##
jun_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) + (1|site/exp_plot), data = jun_ev_dat, family = "binomial")
jun_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|site/exp_plot), data = jun_ev_dat, family = "binomial")
jun_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) + (1|site/exp_plot), data = jun_ev_dat, family = "binomial")
jun_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|site/exp_plot), data = jun_ev_dat, family = "binomial")
jun_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) + (1|site/exp_plot), data = jun_ev_dat, family = "binomial")
jun_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|site/exp_plot), data = jun_ev_dat, family = "binomial")

# compare different plant models
AIC(jun_ev_totb_mod, jun_ev_sepb_mod, jun_ev_totd_mod, jun_ev_sepd_mod, jun_ev_totp_mod, jun_ev_sepp_mod)
# total biomass
summary(jun_ev_totb_mod)
# none
stepAIC(jun_ev_totb_mod)
# remove total biomass
jun_ev_totb_mod2 <- update(jun_ev_totb_mod, .~. - total_biomass.g - fungicide:total_biomass.g)
summary(jun_ev_totb_mod2)
# none
plot(simulateResiduals(jun_ev_totb_mod2))


## July Ev ##
jul_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) + (1|site/exp_plot), data = jul_ev_dat, family = "binomial")
jul_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|site/exp_plot), data = jul_ev_dat, family = "binomial")
jul_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) + (1|site/exp_plot), data = jul_ev_dat, family = "binomial")
jul_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|site/exp_plot), data = jul_ev_dat, family = "binomial")
jul_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) + (1|site/exp_plot), data = jul_ev_dat, family = "binomial")
jul_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|site/exp_plot), data = jul_ev_dat, family = "binomial")

# compare different plant models
AIC(jul_ev_totb_mod, jul_ev_sepb_mod, jul_ev_totd_mod, jul_ev_sepd_mod, jul_ev_totp_mod, jul_ev_sepp_mod)
# total density
summary(jul_ev_totd_mod)
# none
stepAIC(jul_ev_totd_mod)
# keep full model
plot(simulateResiduals(jul_ev_totd_mod))


## early Aug Ev ##
eau_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) + (1|site/exp_plot), data = eau_ev_dat, family = "binomial")
eau_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|site/exp_plot), data = eau_ev_dat, family = "binomial")
eau_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) + (1|site/exp_plot), data = eau_ev_dat, family = "binomial")
eau_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|site/exp_plot), data = eau_ev_dat, family = "binomial")
eau_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) + (1|site/exp_plot), data = eau_ev_dat, family = "binomial")
eau_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|site/exp_plot), data = eau_ev_dat, family = "binomial")

# compare different plant models
AIC(eau_ev_totb_mod, eau_ev_sepb_mod, eau_ev_totd_mod, eau_ev_sepd_mod, eau_ev_totp_mod, eau_ev_sepp_mod)
# total presence
summary(eau_ev_totp_mod)
# bg present (*), edge severity (**)
stepAIC(eau_ev_totp_mod)
# keep full model
plot(simulateResiduals(eau_ev_totp_mod))


## late Aug Ev ##
lau_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) + (1|site/exp_plot), data = lau_ev_dat, family = "binomial")
lau_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|site/exp_plot), data = lau_ev_dat, family = "binomial")
lau_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) + (1|site/exp_plot), data = lau_ev_dat, family = "binomial")
lau_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|site/exp_plot), data = lau_ev_dat, family = "binomial")
lau_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) + (1|site/exp_plot), data = lau_ev_dat, family = "binomial")
lau_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|site/exp_plot), data = lau_ev_dat, family = "binomial")

# compare different plant models
AIC(lau_ev_totb_mod, lau_ev_sepb_mod, lau_ev_totd_mod, lau_ev_sepd_mod, lau_ev_totp_mod, lau_ev_sepp_mod)
# total density
summary(lau_ev_totd_mod)
# fungicide:total density (*) fungicide:edge (*)
stepAIC(lau_ev_totd_mod)
# keep full model
plot(simulateResiduals(lau_ev_totd_mod))


## September Ev ##
sep_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")

# compare different plant models
AIC(sep_ev_totb_mod, sep_ev_sepb_mod, sep_ev_totd_mod, sep_ev_sepd_mod, sep_ev_totp_mod, sep_ev_sepp_mod)
# separate presence
summary(sep_ev_sepp_mod)
# none
stepAIC(sep_ev_sepp_mod)
# remove fung:Ev and edge severity
sep_ev_sepp_mod2 <- update(sep_ev_sepp_mod, .~. - fungicide:Ev_present - edge_severity - fungicide:edge_severity)
summary(sep_ev_sepp_mod2)
# fung:Mv present (*)
plot(simulateResiduals(sep_ev_sepp_mod2))


#### treatment by month ####

trt_sig_dat = tibble(month = c("may", "jun", "jul", "early_aug", "late_aug", "sep"),
                     sig = c("", "", "", "", "*", "*"),
                     treatment = c(NA_character_, NA_character_, NA_character_, NA_character_, "fungicide", "water")) %>%
  mutate(sp = "Ev",
         metric = "spots")


#### output ####

write_csv(trt_sig_dat, "output/Bp_spots_treatment_sig_2019_density_exp.csv")
bp_spots_water_mod <- eau_ev_sepb_water_mod
save(bp_spots_water_mod, file ="output/Bp_spots_water_model_2019_density_exp.rda")
bp_spots_mod <- eau_ev_totp_mod
save(bp_spots_mod, file ="output/Bp_spots_model_2019_density_exp.rda")
