##### info ####

# file: focal_severity_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/2/20
# goal: analyze focal disease severity


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(DHARMa)
library(MASS)
library(MuMIn)
library(tidyverse)
library(glmmTMB)
library(cowplot)
library(lubridate)
library(GGally)

# import all raw data files
dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
edge_dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
temp_hum <- read_csv("intermediate-data/temp_humidity_daily_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


#### edit data ####

# severity dates
dat %>%
  group_by(month) %>%
  summarise(first_day = min(date),
            last_day = max(date))

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

# combine temp/hum, edge, and biomass
cov_dat <- full_join(edge_dat2, temp_hum2) %>%
  full_join(bio) %>%
  full_join(dens_dat)

# calculate severity
# add month columns
# add background disease data
# transform severity
dat2 <- dat %>%
  mutate(leaf_severity = lesion_area.pix / leaf_area.pix,
         plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         plant_severity_adjusted = case_when(plant_severity == 0 ~ 0.001,
                                             plant_severity >= 1 ~ 0.999,
                                             TRUE ~ plant_severity),
         plant_severity_transformed = asin(sqrt(plant_severity_adjusted)),
         Month = recode(month, early_aug = "Early August", jul = "July", 
                        jun = "June", late_aug = "Late August", may = "May", sep = "September") %>%
           fct_relevel("May", "June", "July", "Early August", "Late August", "September"),
         Month_num = as.numeric(Month),
         plant_type = paste(sp, age, sep = " "),
         fungicide = case_when(treatment == "water" ~ 0,
                               TRUE ~ 1),
         exp_plot = paste0(plot, toupper(substr(treatment, 1, 1))),
         sp_trt = paste0(sp, toupper(substr(treatment, 1, 1))),
         Treatment = recode(treatment, water = "control (water)")) %>%
  left_join(edge_dat2) %>%
  left_join(temp_hum2) %>%
  left_join(bio) %>%
  left_join(dens_dat)

# separate out May data
may_dat <- dat2 %>%
  filter(Month == "May")

# separate out focal data
foc_dat <- dat2 %>%
  filter(focal == 1)

# add treatment plots (repeats no neighbor plots)
may_dat_plots <- may_dat %>%
  left_join(plots)

# add treatment plots (repeats no neighbor plots)
foc_dat_plots <- foc_dat %>%
  left_join(plots) %>%
  mutate(background_pres = case_when(density_level == "none" ~ 0,
                                     TRUE ~ 1))

# only allow for sites with both treatment in September (some lost in the mail)
# information obtained from examine data section
foc_dat2 <- foc_dat %>%
  filter(!(sp == "Ev" & site != "D2" & Month == "September")) %>%
  filter(!(sp == "Mv" & site == "D2" & Month == "September"))

foc_dat_plots2 <- foc_dat_plots %>%
  filter(!(sp == "Ev" & site != "D2" & Month == "September")) %>%
  filter(!(sp == "Mv" & site == "D2" & Month == "September"))


#### visualize ####

# all data
(all_dat_viz <- foc_dat2 %>%
  filter(!is.na(leaf_severity)) %>%
  ggplot(aes(Month_num, leaf_severity, color = as.factor(plot))) +
  stat_summary(fun = "mean", geom = "point", aes(shape = treatment)) +
  stat_summary(fun = "mean", geom = "line", aes(linetype = treatment)) +
  facet_grid(site ~ sp))
# late August and September data probably aren't usable (late August data mostly 1's, September data missing from D1 and D2 because of shipping issue)
# some clear and unexpected separations between fungicide and water treatments that change over time

# remove late data
all_dat_viz %+%
  filter(foc_dat2, !is.na(leaf_severity) & Month_num < 5)
# Why does Mv drop so low in July? generally healthy looking leaves
# generally healthy-looking leaves in June too. The high values seem to be from tip senescence
# lesions aren't super common in early August even. seems like senescence is a big contributor

# treatment effect
(treat_viz <- foc_dat2 %>%
  filter(!is.na(leaf_severity)) %>%
  ggplot(aes(Month_num, leaf_severity, color = treatment, shape = sp, linetype = sp)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1))

treat_viz %+%
  filter(foc_dat2, !is.na(plant_severity))  %+%
  aes(y = plant_severity)
# much more intuitive pattern

# for comparison to Vida's results
treat_viz %+%
  aes(y = (leaves_infec/leaves_tot))

# treatment effect in May for Ev from previous year
(may_treat_viz <- may_dat %>%
  ggplot(aes(treatment, leaf_severity)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1))
# legacy fungicide effect

may_treat_viz %+%
  aes(y = plant_severity)
# similar to leaf metric, but smaller values

# May density effect
(month_viz <- may_dat_plots %>%
  ggplot(aes(background_density, leaf_severity, color = treatment)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
    facet_grid(plant_type ~ background, scales = "free"))
# no clear density effect and a lot of variation within density value

month_viz %+%
  aes(y = plant_severity)
# no clear density effect and a lot of variation within density value

# June density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "June")
# no strong density effects

month_viz %+%
  filter(foc_dat_plots2, Month == "June") %+%
  aes(y = plant_severity)

# July density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "July")
# no strong density effects

month_viz %+%
  filter(foc_dat_plots2, Month == "July") %+%
  aes(y = plant_severity)

# early August density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "Early August")
# mv background Ev seedling is what I would expect

month_viz %+%
  filter(foc_dat_plots2, Month == "Early August") %+%
  aes(y = plant_severity)
# fungicide effects seem stronger when the background is Mv

month_viz %+%
  filter(foc_dat_plots2, Month == "Early August") %+%
  aes(y = plant_severity) +
  facet_grid(sp ~ background, scales = "free")
# mv density almost suppressing Ev disease

# late August density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "Late August")
# no strong density effects

month_viz %+%
  filter(foc_dat_plots2, Month == "Late August") %+%
  aes(y = plant_severity)
# slight Mv on Mv density effect, maybe best example

month_viz %+%
  filter(foc_dat_plots2, Month == "Late August") %+%
  aes(y = plant_severity) +
  facet_grid(sp ~ background, scales = "free")
# no effect on Ev

# September density effect
month_viz %+%
  filter(foc_dat_plots2, Month == "September")
# strong density effect of Mv on Mv 

month_viz %+%
  filter(foc_dat_plots2, Month == "September") %+%
  aes(y = plant_severity)
# strong density effect of Mv on Mv 

# correlations between covariates
filter(cov_dat, treatment == "water" & month == "early_aug") %>%
  select(temp_avg:total_density) %>%
  ggpairs()
# can only pick one temp or humidity metric, all are correlated except max temp and humidity
# total biomass and Mv biomass are correlated
# density and biomass are correlated
# Mv and Ev biomass are negatively correlated at 0.36

filter(cov_dat, treatment == "water" & month == "late_aug") %>%
  select(temp_avg:total_density) %>%
  ggpairs()
# same as above
# humidity average seems like a representative variable because it's correlated with all three temperature metrics and all three humidity metrics in both months

# focal severity and edge severity
(edge_viz <- foc_dat2 %>%
    filter(!is.na(leaf_severity) & !is.na(edge_severity)) %>%
    ggplot(aes(edge_severity, leaf_severity, color = treatment)) +
    stat_summary(fun = "mean", geom = "point") +
    facet_grid(sp ~ Month, scales = "free_x"))
# maybe most positive for Mv in early August (makes sense)

edge_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(edge_severity)) %+%
  aes(y = plant_severity)
# similar to leaf metric

# focal severity and biomass
(bio_viz <- foc_dat2 %>%
    filter(!is.na(leaf_severity) & !is.na(total_biomass.g)) %>%
    ggplot(aes(total_biomass.g, leaf_severity, color = treatment)) +
    stat_summary(fun = "mean", geom = "point") +
    facet_grid(sp ~ Month, scales = "free_x"))

bio_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(total_biomass.g)) %+%
  aes(y = plant_severity)

# focal severity and Mv biomass
(mbio_viz <- foc_dat2 %>%
    filter(!is.na(leaf_severity) & !is.na(Mv_biomass.g)) %>%
    ggplot(aes(Mv_biomass.g, leaf_severity, color = treatment)) +
    stat_summary(fun = "mean", geom = "point") +
    facet_grid(sp ~ Month, scales = "free_x"))

mbio_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(Mv_biomass.g)) %+%
  aes(y = plant_severity)

# focal severity and Ev biomass
(ebio_viz <- foc_dat2 %>%
    filter(!is.na(leaf_severity) & !is.na(Ev_biomass.g)) %>%
    ggplot(aes(Ev_biomass.g, leaf_severity, color = treatment)) +
    stat_summary(fun = "mean", geom = "point") +
    facet_grid(sp ~ Month, scales = "free_x"))

ebio_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(Ev_biomass.g)) %+%
  aes(y = plant_severity)

# focal severity and temperature
(temp_viz <- foc_dat2 %>%
    filter(!is.na(leaf_severity) & !is.na(temp_avg)) %>%
    ggplot(aes(temp_avg, leaf_severity, color = treatment)) +
    geom_point() +
    facet_grid(sp ~ Month, scales = "free_x"))

temp_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(temp_avg)) %+%
  aes(y = plant_severity)

filter(foc_dat2, !is.na(plant_severity) & !is.na(temp_avg)) %>%
  ggplot(aes(temp_avg, plant_severity)) +
  geom_point() +
  facet_wrap(~sp, scales = "free")

# focal severity and humidity
(hum_viz <- foc_dat2 %>%
    filter(!is.na(leaf_severity) & !is.na(hum_avg)) %>%
    ggplot(aes(hum_avg, leaf_severity, color = treatment)) +
    geom_point() +
    facet_grid(sp ~ Month, scales = "free_x"))

hum_viz %+%
  filter(foc_dat2, !is.na(plant_severity) & !is.na(hum_avg)) %+%
  aes(y = plant_severity)

filter(foc_dat2, !is.na(plant_severity) & !is.na(hum_avg)) %>%
  ggplot(aes(hum_avg, plant_severity)) +
  geom_point() +
  facet_wrap(~sp, scales = "free")


#### examine data ####

# September data for Mv
foc_dat %>%
  filter(sp == "Mv" & Month == "September") %>%
  group_by(plot, treatment) %>%
  summarise(sites = length(unique(site)),
            plants = n(),
            site_names = paste(unique(site), collapse = ",")) %>%
  data.frame()
# one more site for all water plots (D2)

# September data for Ev
foc_dat %>%
  filter(sp == "Ev" & Month == "September") %>%
  group_by(plot, treatment) %>%
  summarise(sites = length(unique(site)),
            plants = n(),
            site_names = paste(unique(site), collapse = ",")) %>%
  data.frame()
# D2 is the only site used throughout

# high lesion values
dat2 %>%
  filter(plant_severity > 1)

# make edits in edit data section


#### model selection for water treatment ####

# include humidity data (only available for control plots)

# water treatment by species
ev_water_dat <- filter(foc_dat2, sp == "Ev" &  treatment == "water" & month %in% c("early_aug", "late_aug"))
sep_ev_water_dat <- filter(foc_dat2, sp == "Ev" &  treatment == "water" & month == "sep")
mv_water_dat <- filter(foc_dat2, sp == "Mv" &  treatment == "water" & month %in% c("early_aug", "late_aug"))
sep_mv_water_dat <- filter(foc_dat2, sp == "Mv" &  treatment == "water" & month == "sep")


## early and late August Ev ##
ev_totb_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (total_biomass.g + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "beta_family")
# can ignore warning "NA/NaN function evaluation" if model converges
ev_sepb_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "beta_family")
ev_totd_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (total_density + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "beta_family")
ev_sepd_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (Ev_density + Mv_density + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "beta_family")
ev_totp_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (bg_present + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "beta_family")
ev_sepp_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (Ev_present + Mv_present + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "beta_family")

# compare different plant models
AIC(ev_totb_water_mod, ev_sepb_water_mod, ev_totd_water_mod, ev_sepd_water_mod, ev_totp_water_mod, ev_sepp_water_mod)
# total density is the lowest
summary(ev_totd_water_mod)
# none significant
stepAIC(ev_totd_water_mod)
# keep full model
plot(simulateResiduals(ev_totd_water_mod))
# sig deviation


## September Ev ##
sep_ev_totb_water_mod <- glmmTMB(plant_severity_adjusted ~ total_biomass.g + edge_severity + hum_avg + (1|plot), data = sep_ev_water_dat, family = "beta_family")
# could not converge
summary(sep_ev_totb_water_mod)
sep_ev_totb_water_mod <- glmmTMB(plant_severity_adjusted ~ total_biomass.g + edge_severity + hum_avg, data = sep_ev_water_dat, family = "beta_family")
sep_ev_sepb_water_mod <- glmmTMB(plant_severity_adjusted ~ Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg, data = sep_ev_water_dat, family = "beta_family")
sep_ev_totd_water_mod <- glmmTMB(plant_severity_adjusted ~ total_density + edge_severity + hum_avg, data = sep_ev_water_dat, family = "beta_family")
sep_ev_sepd_water_mod <- glmmTMB(plant_severity_adjusted ~ Ev_density + Mv_density + edge_severity + hum_avg, data = sep_ev_water_dat, family = "beta_family")
sep_ev_totp_water_mod <- glmmTMB(plant_severity_adjusted ~ bg_present + edge_severity + hum_avg, data = sep_ev_water_dat, family = "beta_family")
sep_ev_sepp_water_mod <- glmmTMB(plant_severity_adjusted ~ Ev_present + Mv_present + edge_severity + hum_avg, data = sep_ev_water_dat, family = "beta_family")

# compare different plant models
AIC(sep_ev_totb_water_mod, sep_ev_sepb_water_mod, sep_ev_totd_water_mod, sep_ev_sepd_water_mod, sep_ev_totp_water_mod, sep_ev_sepp_water_mod)
# total density is the lowest
summary(sep_ev_totd_water_mod)
# edge severity is significant
stepAIC(sep_ev_totd_water_mod)
# remove humidity
sep_ev_totd_water_mod2 <- update(sep_ev_totd_water_mod, .~. -hum_avg)
summary(sep_ev_totd_water_mod2)
# none sig
plot(simulateResiduals(sep_ev_totd_water_mod2))


## early and late August Mv ##
mv_totb_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (total_biomass.g + edge_severity + hum_avg) + (1|site/plot), data = mv_water_dat, family = "beta_family")
mv_sepb_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg) + (1|site/plot), data = mv_water_dat, family = "beta_family")
mv_totd_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (total_density + edge_severity + hum_avg) + (1|site/plot), data = mv_water_dat, family = "beta_family")
mv_sepd_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (Ev_density + Mv_density + edge_severity + hum_avg) + (1|site/plot), data = mv_water_dat, family = "beta_family")
mv_totp_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (bg_present + edge_severity + hum_avg) + (1|site/plot), data = mv_water_dat, family = "beta_family")
mv_sepp_water_mod <- glmmTMB(plant_severity_adjusted ~ Month * (Ev_present + Mv_present + edge_severity + hum_avg) + (1|site/plot), data = mv_water_dat, family = "beta_family")

# compare different plant models
AIC(mv_totb_water_mod, mv_sepb_water_mod, mv_totd_water_mod, mv_sepd_water_mod, mv_totp_water_mod, mv_sepp_water_mod)
# total biomass is the lowest
summary(mv_totb_water_mod)
# none significant
stepAIC(mv_totb_water_mod)
# keep full model
plot(simulateResiduals(mv_totb_water_mod))


## September Mv ##
sep_mv_totb_water_mod <- glmmTMB(plant_severity_adjusted ~ total_biomass.g + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")
sep_mv_totb_water_mod <- glmmTMB(plant_severity_adjusted ~ total_biomass.g + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")
sep_mv_sepb_water_mod <- glmmTMB(plant_severity_adjusted ~ Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")
sep_mv_totd_water_mod <- glmmTMB(plant_severity_adjusted ~ total_density + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")
sep_mv_sepd_water_mod <- glmmTMB(plant_severity_adjusted ~ Ev_density + Mv_density + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")
sep_mv_totp_water_mod <- glmmTMB(plant_severity_adjusted ~ bg_present + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")
sep_mv_sepp_water_mod <- glmmTMB(plant_severity_adjusted ~ Ev_present + Mv_present + edge_severity + hum_avg + (1|site/plot), data = sep_mv_water_dat, family = "beta_family")

# compare different plant models
AIC(sep_mv_totb_water_mod, sep_mv_sepb_water_mod, sep_mv_totd_water_mod, sep_mv_sepd_water_mod, sep_mv_totp_water_mod, sep_mv_sepp_water_mod)
# total presence is the lowest
summary(sep_mv_totp_water_mod)
# humidity is significant
stepAIC(sep_mv_totp_water_mod)
# remove everything except humidity
sep_mv_totp_water_mod2 <- update(sep_mv_totp_water_mod, .~. -bg_present - edge_severity)
summary(sep_mv_totp_water_mod2)
# humidity (***)
plot(simulateResiduals(sep_mv_totp_water_mod2))


#### model selection for both treatments ####

# edge severity is not significantly different between fungicide and water plots (edge_severity_analysis_2019_density_exp.R)

# split data by month
jun_lau_ev_dat <- filter(foc_dat2, sp == "Ev" & month %in% c("jun", "jul", "early_aug", "late_aug"))
may_ev_dat <- filter(foc_dat2, sp == "Ev" & month == "may")
sep_ev_dat <- filter(foc_dat2, sp == "Ev" & month == "sep")

jun_lau_mv_dat <- filter(foc_dat2, sp == "Mv" & month %in% c("jun", "jul", "early_aug", "late_aug"))
may_mv_dat <- filter(foc_dat2, sp == "Mv" & month == "may")
sep_mv_dat <- filter(foc_dat2, sp == "Mv" & month == "sep")


## May all Ev ##
may_ev_all_totb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * total_biomass.g + (1|site/exp_plot), data = may_dat, family = "beta_family")
may_ev_all_sepb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_biomass.g + Mv_biomass.g) + (1|site/exp_plot), data = may_dat, family = "beta_family")
may_ev_all_totd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * total_density + (1|site/exp_plot), data = may_dat, family = "beta_family")
may_ev_all_sepd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_density + Mv_density) + (1|site/exp_plot), data = may_dat, family = "beta_family")
may_ev_all_totp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * bg_present + (1|site/exp_plot), data = may_dat, family = "beta_family")
may_ev_all_sepp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_present + Mv_present) + (1|site/exp_plot), data = may_dat, family = "beta_family")

# compare different plant models
AIC(may_ev_all_totb_mod, may_ev_all_sepb_mod, may_ev_all_totd_mod, may_ev_all_sepd_mod, may_ev_all_totp_mod, may_ev_all_sepp_mod)
# separate biomass is the lowest
summary(may_ev_all_sepb_mod)
# Ev biomass sig
stepAIC(may_ev_all_sepb_mod)
# remove interactions
may_ev_all_sepb_mod2 <- update(may_ev_all_sepb_mod, .~. - fungicide:Ev_biomass.g - fungicide:Mv_biomass.g)
summary(may_ev_all_sepb_mod2)
# fungicide (*) and Ev biomass sig
plot(simulateResiduals(may_ev_all_sepb_mod2))


## May Ev ##
may_ev_totb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * total_biomass.g + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
may_ev_sepb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_biomass.g + Mv_biomass.g) + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
may_ev_totd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * total_density + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
may_ev_sepd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_density + Mv_density) + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
may_ev_totp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * bg_present + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
may_ev_sepp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_present + Mv_present) + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")

# compare different plant models
AIC(may_ev_totb_mod, may_ev_sepb_mod, may_ev_totd_mod, may_ev_sepd_mod, may_ev_totp_mod, may_ev_sepp_mod)
# sep presence
summary(may_ev_sepp_mod)
# none
stepAIC(may_ev_sepp_mod)
# remove everything except fungicide and Mv presence
may_ev_sepp_mod2 <- glmmTMB(plant_severity_adjusted ~ fungicide + Mv_present + (1|site/exp_plot), data = may_ev_dat, family = "beta_family")
summary(may_ev_sepp_mod2)
# fungicide (*) and Mv presence (**)
plot(simulateResiduals(may_ev_sepb_mod2))


## June-Late August Ev ##
ev_totb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_biomass.g + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "beta_family")
ev_sepb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "beta_family")
ev_totd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_density + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "beta_family")
ev_sepd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_density + Mv_density + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "beta_family")
ev_totp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (bg_present + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "beta_family")
ev_sepp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_present + Mv_present + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "beta_family")

# compare different plant models
AIC(ev_totb_mod, ev_sepb_mod, ev_totd_mod, ev_sepd_mod, ev_totp_mod, ev_sepp_mod)
# sep biomass
summary(ev_sepb_mod)
# fungicide:edge severity:month (***)
# Mv biomass:month (*)
# Ev biomass:month (*)
# fungicide:Ev biomass (**)
stepAIC(ev_sepb_mod)
# keep full model
plot(simulateResiduals(ev_sepb_mod))
# sig deviance


## June-late August Mv ##
mv_totb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_biomass.g + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_mv_dat, family = "beta_family")
mv_sepb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_mv_dat, family = "beta_family")
mv_totd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_density + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_mv_dat, family = "beta_family")
mv_sepd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_density + Mv_density + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_mv_dat, family = "beta_family")
mv_totp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (bg_present + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_mv_dat, family = "beta_family")
mv_sepp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_present + Mv_present + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_mv_dat, family = "beta_family")

# compare different plant models
AIC(mv_totb_mod, mv_sepb_mod, mv_totd_mod, mv_sepd_mod, mv_totp_mod, mv_sepp_mod)
# total biomass
summary(mv_totb_mod)
# fungicide:edge severity:month
# total biomass:month
stepAIC(mv_totb_mod)
# keep full model
plot(simulateResiduals(mv_totb_mod))
# sig deviance


## September Ev ##
sep_ev_totb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_biomass.g + edge_severity), data = sep_ev_dat, family = "beta_family")
sep_ev_sepb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity), data = sep_ev_dat, family = "beta_family")
sep_ev_totd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_density + edge_severity), data = sep_ev_dat, family = "beta_family")
sep_ev_sepd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_density + Mv_density + edge_severity), data = sep_ev_dat, family = "beta_family")
sep_ev_totp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (bg_present + edge_severity), data = sep_ev_dat, family = "beta_family")
sep_ev_sepp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_present + Mv_present + edge_severity), data = sep_ev_dat, family = "beta_family")

# compare different plant models
AIC(sep_ev_totb_mod, sep_ev_sepb_mod, sep_ev_totd_mod, sep_ev_sepd_mod, sep_ev_totp_mod, sep_ev_sepp_mod)
# sep biomass is the lowest
summary(sep_ev_sepb_mod)
# fungicide:Mv bio (*)
stepAIC(sep_ev_sepb_mod)
# keep full model
plot(simulateResiduals(sep_ev_sepb_mod))


## September Mv ##
sep_mv_totb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_biomass.g + edge_severity) + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")
sep_mv_sepb_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")
sep_mv_totd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (total_density + edge_severity) + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")
sep_mv_sepd_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")
sep_mv_totp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (bg_present + edge_severity) + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")
sep_mv_sepp_mod <- glmmTMB(plant_severity_adjusted ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|site/exp_plot), data = sep_mv_dat, family = "beta_family")

# compare different plant models
AIC(sep_mv_totb_mod, sep_mv_sepb_mod, sep_mv_totd_mod, sep_mv_sepd_mod, sep_mv_totp_mod, sep_mv_sepp_mod)
# total biomass is the lowest
summary(sep_mv_totb_mod)
# total biomass
stepAIC(sep_mv_totb_mod)
# remove fungicide:edge
sep_mv_totb_mod2 <- update(sep_mv_totb_mod, .~. - fugicide:edge_severity)
summary(sep_mv_totb_mod2)
# total biomass (**)
plot(simulateResiduals(sep_mv_totb_mod2))


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


## Humidity Mv September ##
sep_mv_water_dat_sim <- sep_mv_water_dat %>%
  mutate(pred = predict(sep_mv_totp_water_mod2, newdata = sep_mv_water_dat, re.form = NA, type = "response"),
         pred_se = predict(sep_mv_totp_water_mod2, newdata = sep_mv_water_dat, re.form = NA, type = "response", se.fit = T)$se.fit)

hum_mv_fig <- ggplot(sep_mv_water_dat_sim, aes(x = hum_avg, y = plant_severity_adjusted)) +
  geom_point(alpha = 0.3, aes(color = Treatment)) +
  geom_ribbon(alpha = 0.5, aes(fill = Treatment, ymin = pred - pred_se, ymax = pred + pred_se)) +
  geom_line(aes(y = pred, color = Treatment)) +
  scale_color_manual(values = col_pal, guide = F) +
  scale_fill_manual(values = col_pal, guide = F) +
  xlab("Average daily relative humidity (%)") +
  ylab(expression(paste("Proportion ", italic("Microstegium"), " leaf area brown in September", sep = ""))) +
  temp_theme 


## Separate presence Ev May ##
may_ev_dat_sim <- may_ev_dat %>%
  select(Treatment, fungicide, Mv_present) %>%
  unique() %>%
  mutate(site = NA_character_,
         exp_plot = NA_character_,
         background = case_when(Mv_present == 0 ~ "no",
                                TRUE ~ "yes")) %>%
  mutate(plant_severity_adjusted = predict(may_ev_sepp_mod2, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(may_ev_sepp_mod2, newdata = ., re.form = NA, se.fit = T, type = "response")$se.fit)

may_ev_dat2 <- may_ev_dat %>%
  mutate(background = case_when(Mv_present == 0 ~ "no",
                                TRUE ~ "yes"))

may_ev_fig <- ggplot(may_ev_dat_sim, aes(x = background, y = plant_severity_adjusted)) +
  geom_point(data = may_ev_dat2, alpha = 0.3, aes(color = Treatment)) +
  geom_errorbar(width = 0.01, aes(group = Treatment, ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se)) +
  geom_point(size = 3, shape = 21, aes(fill = Treatment)) +
  xlab(expression(paste(italic("Microstegiun"), " in background", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Elymus"), " leaf area brown in May", sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ggtitle("May") +
  temp_theme


## Separate biomass Ev June-late August ##
jun_lau_ev_biomass_sim_dat <- foc_dat2 %>%
  filter(sp == "Ev" & !(month %in% c("may", "sep"))) %>%
  group_by(Month, Treatment, fungicide) %>%
  summarise(min_ev_biomass = min(Ev_biomass.g),
            max_ev_biomass = max(Ev_biomass.g),
            min_mv_biomass = min(Mv_biomass.g),
            max_mv_biomass = max(Mv_biomass.g)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(Month = .$Month, Treatment = .$Treatment, fungicide = .$fungicide,
            Ev_biomass.g = c(seq(.$min_ev_biomass, .$max_ev_biomass, length.out = 100), rep(0, 100)),
            Mv_biomass.g = c(rep(0, 100), seq(.$min_mv_biomass, .$max_mv_biomass, length.out = 100)))) %>%
  ungroup() %>%
  mutate(edge_severity = 0,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(plant_severity_adjusted = predict(ev_sepb_mod, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(ev_sepb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

jun_lau_ev_bio_ev_fig <- ggplot(filter(jun_lau_ev_biomass_sim_dat, Ev_biomass.g > 0), aes(x = Ev_biomass.g, y = plant_severity_adjusted)) +
  geom_point(data = jun_lau_ev_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  facet_wrap(~ Month) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " biomass (g ", m^-1, ")", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Elymus"), " leaf area brown", sep = ""))) +
  temp_theme

# look at high Ev biomass plants
filter(jun_lau_ev_dat, Ev_biomass.g > 150)
late_aug_foc_dat <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")
filter(late_aug_foc_dat, site == "D3" & plot == 10 & treatment == "fungicide")
# too few green leaves to collect from these plants
filter(ev_bio, site == "D3" & plot == 10 & treatment == "fungicide")
filter(bg_bio, site == "D3" & plot == 10 & treatment == "fungicide")
# checked emailed sheet and it seems fine, could ask Quynh for the paper copy (probably in lab though)

jun_lau_mv_bio_ev_fig <- ggplot(filter(jun_lau_ev_biomass_sim_dat, Mv_biomass.g > 0), aes(x = Mv_biomass.g, y = plant_severity_adjusted)) +
  geom_point(data = jun_lau_ev_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  facet_wrap(~ Month) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegium"), " biomass (g ", m^-1, ")", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Elymus"), " leaf area brown", sep = ""))) +
  temp_theme


## Edge severity Ev June-late August ##
jun_lau_ev_edge_sim_dat <- foc_dat2 %>%
  filter(sp == "Ev" & !(month %in% c("may", "sep")) & !(is.na(edge_severity))) %>%
  group_by(Month, Treatment, fungicide) %>%
  summarise(min_edge = min(edge_severity),
            max_edge = max(edge_severity)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(Month = .$Month, Treatment = .$Treatment, fungicide = .$fungicide,
            edge_severity = seq(.$min_edge, .$max_edge, length.out = 100))) %>%
  ungroup() %>%
  mutate(Mv_biomass.g = 0,
         Ev_biomass.g = 0,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(plant_severity_adjusted = predict(ev_sepb_mod, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(ev_sepb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

jun_lau_edge_ev_fig <- ggplot(jun_lau_ev_edge_sim_dat, aes(x = edge_severity, y = plant_severity_adjusted)) +
  geom_point(data = jun_lau_ev_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  facet_wrap(~ Month, scales = "free_x") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste("Proportion edge ", italic("Microstegium"), " leaf area brown", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Elymus"), " leaf area brown", sep = ""))) +
  temp_theme


## Total biomass Mv June-late August ##
jun_lau_mv_biomass_sim_dat <- foc_dat2 %>%
  filter(sp == "Mv" & !(month %in% c("may", "sep"))) %>%
  group_by(Month, Treatment, fungicide) %>%
  summarise(min_tot_biomass = min(total_biomass.g),
            max_tot_biomass = max(total_biomass.g)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(Month = .$Month, Treatment = .$Treatment, fungicide = .$fungicide,
            total_biomass.g = seq(.$min_tot_biomass, .$max_tot_biomass, length.out = 100))) %>%
  ungroup() %>%
  mutate(edge_severity = 0,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(plant_severity_adjusted = predict(mv_totb_mod, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(mv_totb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

jun_lau_tot_bio_mv_fig <- ggplot(jun_lau_mv_biomass_sim_dat, aes(x = total_biomass.g, y = plant_severity_adjusted)) +
  geom_point(data = jun_lau_mv_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  facet_wrap(~ Month) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste("Total biomass (g ", m^-1, ")", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Microstegium"), " leaf area brown", sep = ""))) +
  temp_theme


## Edge severity Mv June-late August ##
jun_lau_mv_edge_sim_dat <- foc_dat2 %>%
  filter(sp == "Mv" & !(month %in% c("may", "sep")) & !(is.na(edge_severity))) %>%
  group_by(Month, Treatment, fungicide) %>%
  summarise(min_edge = min(edge_severity),
            max_edge = max(edge_severity)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(Month = .$Month, Treatment = .$Treatment, fungicide = .$fungicide,
            edge_severity = seq(.$min_edge, .$max_edge, length.out = 100))) %>%
  ungroup() %>%
  mutate(total_biomass.g = 0,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(plant_severity_adjusted = predict(mv_totb_mod, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(mv_totb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

jun_lau_edge_mv_fig <- ggplot(jun_lau_mv_edge_sim_dat, aes(x = edge_severity, y = plant_severity_adjusted)) +
  geom_point(data = jun_lau_mv_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  facet_wrap(~ Month, scales = "free_x") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste("Proportion edge ", italic("Microstegium"), " leaf area brown", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Microstegium"), " leaf area brown", sep = ""))) +
  temp_theme


## Separate biomass Ev September ##
sep_ev_biomass_sim_dat <- foc_dat2 %>%
  filter(sp == "Ev" & month == "sep") %>%
  group_by(Treatment, fungicide) %>%
  summarise(min_ev_biomass = min(Ev_biomass.g),
            max_ev_biomass = max(Ev_biomass.g),
            min_mv_biomass = min(Mv_biomass.g),
            max_mv_biomass = max(Mv_biomass.g)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(Treatment = .$Treatment, fungicide = .$fungicide,
            Ev_biomass.g = c(seq(.$min_ev_biomass, .$max_ev_biomass, length.out = 100), rep(0, 100)),
            Mv_biomass.g = c(rep(0, 100), seq(.$min_mv_biomass, .$max_mv_biomass, length.out = 100)))) %>%
  ungroup() %>%
  mutate(edge_severity = 0,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(plant_severity_adjusted = predict(sep_ev_sepb_mod, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(sep_ev_sepb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

sep_ev_fig <- ggplot(filter(sep_ev_biomass_sim_dat, Mv_biomass.g > 0), aes(x = Mv_biomass.g, y = plant_severity_adjusted)) +
  geom_point(data = sep_ev_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegium"), " biomass (g ", m^-1, ")", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Elymus"), " leaf area brown", sep = ""))) +
  ggtitle("September") +
  temp_theme


## Total biomass Mv September ##
sep_mv_biomass_sim_dat <- foc_dat2 %>%
  filter(sp == "Mv" & month == "sep") %>%
  summarise(min_tot_biomass = min(total_biomass.g),
            max_tot_biomass = max(total_biomass.g)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(total_biomass.g = seq(.$min_tot_biomass, .$max_tot_biomass, length.out = 100))) %>%
  ungroup() %>%
  mutate(fungicide = 0,
         edge_severity = 0,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(plant_severity_adjusted = predict(sep_mv_totb_mod2, newdata = ., re.form = NA, type = "response"),
         plant_severity_adjusted_se = predict(sep_mv_totb_mod2, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

sep_mv_fig <- ggplot(sep_mv_biomass_sim_dat, aes(x = total_biomass.g, y = plant_severity_adjusted)) +
  geom_point(data = sep_mv_dat, aes(color = Treatment), alpha = 0.5) +
  geom_ribbon(aes(ymin = plant_severity_adjusted - plant_severity_adjusted_se, ymax = plant_severity_adjusted + plant_severity_adjusted_se), fill = "black", alpha = 0.5) +
  geom_line(color = "black") +
  scale_color_manual(values = col_pal) +
  xlab(expression(paste("Total biomass (g ", m^-1, ")", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Microstegium"), " leaf area brown", sep = ""))) +
  ggtitle("September") +
  temp_theme


#### output ####
pdf("output/focal_severity_analysis_2019_density_exp.pdf")
may_ev_fig
jun_lau_ev_bio_ev_fig
jun_lau_mv_bio_ev_fig
jun_lau_edge_ev_fig
sep_ev_fig
hum_mv_fig
jun_lau_tot_bio_mv_fig
jun_lau_edge_mv_fig
sep_mv_fig
dev.off()

mv_severity_humidity_mod <- sep_mv_totp_water_mod2
save(mv_severity_humidity_mod, file ="output/Mv_severity_humidity_model_2019_density_exp.rda")
write_csv(sep_mv_water_dat, "intermediate-data/Mv_severity_water_treatment_sep_2019_density_exp.csv")
mv_severity_biomass_edge_mod <- mv_totb_mod
save(mv_severity_biomass_edge_mod, file ="output/Mv_severity_biomass_edge_model_2019_density_exp.rda")
write_csv(jun_lau_mv_dat, "intermediate-data/Mv_severity_june_late_aug_2019_density_exp.csv")
mv_severity_september_mod <- sep_mv_totb_mod2
save(mv_severity_september_mod, file ="output/Mv_severity_september_model_2019_density_exp.rda")
ev_severity_may_mod <- may_ev_sepp_mod2
save(ev_severity_may_mod, file = "output/Ev_severity_may_model_2019_density_exp.rda")
ev_severity_biomass_edge_mod <- ev_sepb_mod
save(ev_severity_biomass_edge_mod, file ="output/Ev_severity_biomass_edge_model_2019_density_exp.rda")
write_csv(jun_lau_ev_dat, "intermediate-data/Ev_severity_june_late_aug_2019_density_exp.csv")
ev_severity_september_mod <- sep_ev_sepb_mod
save(ev_severity_september_mod, file ="output/Ev_severity_september_model_2019_density_exp.rda")
