##### info ####

# file: temp_humidity_analysis_2019_dens_exp
# author: Amy Kendig, Chris Wojan
# date last edited: 8/20/20
# goal: evaluate the effects of density on temperature and humidity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(MASS)
library(tidyverse)
library(lubridate)
library(lme4)
library(glmmTMB)
library(DHARMa) # plot residuals for glmmTMB
library(weathermetrics)

# import data
hr_dat <- read_csv("./intermediate-data/temp_humidity_hourly_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
covariates <- read_csv("intermediate-data/covariates_2018_density_exp.csv")

# figure template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.box.margin = margin(-14, -14, -14, -14),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        strip.placement = "outside",
        plot.title = element_text(size = 14, hjust = 0.5))


#### biomass and density data ####

# copied from focal_severity_analysis_2019_density_exp.R

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


#### humidity and time ####

# times that humidity extremes occur
hum_ext_dat <-  hr_dat %>%
  group_by(site, plot, day) %>%
  mutate(min_hum = min(hum_prop),
         max_hum = max(hum_prop),
         time_type = case_when(hum_prop == min_hum ~ "Minimum humidity",
                               hum_prop == max_hum ~ "Maximum humidity",
                               TRUE ~ NA_character_),
         hour = format(time, "%H") %>%
           as.numeric()) %>%
  ungroup() %>%
  filter(!is.na(time_type))

# figure
pdf("output/hours_of_humidity_extremes.pdf", width = 5, height = 3)
ggplot(hum_ext_dat, aes(x = hour)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 15, color = "red", linetype = "dashed") +
  facet_wrap(~ time_type) +
  xlab("Hour of the day") +
  ylab("Count of days separated by plot") +
  temp_theme
dev.off()


#### relative humidity duration ####

# visualize
# hr_dat %>%
#   mutate(hour = hour(time)) %>%
#   filter(day < as.Date("2019-07-29")) %>%
#   ggplot(aes(x = hour, y = rel_hum, color = as.factor(day))) +
#   geom_line() +
#   facet_grid(site ~ plot) +
#   theme_bw() +
#   theme(legend.position = "none")
# threhsold at or near 1


#### dew point, heat index, saturation time ####

# add new metrics
hr_dat2 <- hr_dat %>%
  mutate(dewpoint = humidity.to.dewpoint(rh = rel_hum, t = temp, temperature.metric = "celsius"),
         heat_ind = heat.index(t = temp, rh = rel_hum, temperature.metric = "celsius"),
         hour = format(time, "%H") %>% as.numeric(),
         day_hum = case_when(hour < 15 ~ day,
                             TRUE ~ day + 1),
         hour_hum = case_when(hour >= 15 ~ hour - 15,
                              TRUE ~ hour + 9)) 

# splines for each day and plot
# 4 was the highest df that would lead to 0 or 1 saturation time frames (rather than > 1)
# indicate whether predicted humidity switches in (1) or out (-1) of saturation
# use 1 at beginning if day starts at 100%
spline_dat <- hr_dat2 %>%
  filter(day_hum > "2019-07-03" & day_hum < "2019-10-21") %>%
  group_by(site, plot, day_hum) %>%
  mutate(pred_hum = predict(smooth.spline(hour_hum, hum_prop, df = 4), hour_hum)$y,
         pred_sat = ifelse(pred_hum >= 1, 1, 0),
         sat_switch = c(0, diff(pred_sat))) %>%
  group_by(site, plot, day_hum) %>%
  mutate(sat_switch = case_when(pred_sat == 1 & hour_hum == 0 ~ 1,
                                TRUE ~ sat_switch)) %>%
  ungroup()

# number of times predicted humidity hits 1
spline_dat %>%
  group_by(site, plot, day_hum) %>%
  summarise(start_sat = sum(sat_switch == 1)) %>%
  ungroup() %>%
  ggplot(aes(x = start_sat)) +
  geom_histogram()

spline_dat %>%
  group_by(site, plot, day_hum) %>%
  summarise(start_sat = sum(sat_switch == -1)) %>%
  ungroup() %>%
  ggplot(aes(x = start_sat)) +
  geom_histogram()

# example plot
pdf("output/example_spline_fit_humidity_time.pdf", width = 5, height = 5)
spline_dat %>%
  filter(site == "D1" & plot == 1 & day_hum == "2019-08-20") %>%
  ggplot(aes(x = hour_hum, y = hum_prop)) +
  geom_point() +
  geom_line(aes(y = pred_hum)) +
  xlab("Hour of the day") +
  ylab("Relative humidity") +
  temp_theme
dev.off()

# determine when predicted humidity hits 100%
min_hr_hum <- spline_dat %>%
  filter(sat_switch == 1) %>%
  group_by(site, plot, day_hum) %>%
  summarise(min_sat_time = min(hour_hum))

# determine when predicted humidity goes below 100%
max_hr_hum <- spline_dat %>%
  filter(sat_switch == -1) %>%
  group_by(site, plot, day_hum) %>%
  summarise(max_sat_time = max(hour_hum))

# combine humidity times
hr_hum <- full_join(min_hr_hum, max_hr_hum) %>%
  mutate(hours_sat = max_sat_time - min_sat_time)

range(hr_hum$hours_sat, na.rm = T)

# add in saturation time
hr_dat3 <- hr_dat2 %>%
  left_join(hr_hum %>%
              select(-hours_sat)) %>%
  mutate(pred_sat = ifelse(hour_hum >= min_sat_time & hour_hum < max_sat_time, 1, 0))

# range of humidities within effective saturation
hr_dat3 %>%
  filter(pred_sat == 1) %>%
  summarise(max_hum = max(hum_prop),
            min_hum = min(hum_prop),
            min_hr_hum = min(min_sat_time),
            max_hr_hum = max(max_sat_time))

# figure
hr_dat3 %>%
  filter(pred_sat == 1) %>%
  ggplot(aes(hour_hum, hum_prop, color = as.factor(day_hum))) +
  geom_line() +
  facet_grid(site ~ plot) +
  theme_bw() +
  theme(legend.position = "none")


#### summarize weather data ####

# summarize dew point by each day, using only hours at saturation
dy_dat_dew <- filter(hr_dat3, pred_sat == 1) %>%
  group_by(site, plot, treatment, day_hum) %>%
  summarise(dewp_min_hum = min(dewpoint),
            dewp_max_hum = max(dewpoint)) %>%
  ungroup() %>%
  mutate(dewp_rng_hum = dewp_max_hum - dewp_min_hum)

# daily summary for humidity data
dy_hum_dat <- hr_dat3 %>%
  group_by(site, plot, treatment, day_hum) %>%
  summarise(hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            hum_dur = sum(hum_prop == 1),
            sat_dur = unique(max_sat_time - min_sat_time) %>%
              replace_na(0)) %>%
  ungroup() %>%
  full_join(dy_dat_dew %>%
              select(site, plot, treatment, day_hum, dewp_rng_hum)) %>%
  mutate(month = case_when(day_hum < as.Date("2019-07-29") ~ "early_aug",
                           day_hum >= as.Date("2019-07-29") & day_hum < as.Date("2019-08-28") ~ "late_aug",
                           day_hum >= as.Date("2019-08-28") & day_hum < as.Date("2019-09-24") ~ "sep",
                           TRUE ~ "oct") %>%
           fct_relevel("early_aug", "late_aug", "sep", "oct"),
         dewp_rng_hum = replace_na(dewp_rng_hum, 0),
         dew_intensity = dewp_rng_hum * hum_dur,
         dew_intensity2 = dewp_rng_hum * sat_dur)

# daily summary for temperature data
dy_temp_dat <- hr_dat3 %>%
  group_by(site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hind_max = max(heat_ind)) %>%
  ungroup() %>%
  mutate(month = case_when(day < as.Date("2019-07-29") ~ "early_aug",
                           day >= as.Date("2019-07-29") & day < as.Date("2019-08-28") ~ "late_aug",
                           day >= as.Date("2019-08-28") & day < as.Date("2019-09-24") ~ "sep",
                           TRUE ~ "oct") %>%
           fct_relevel("early_aug", "late_aug", "sep", "oct"))

# summarize by month
mo_dat <- dy_hum_dat %>%
  group_by(site, plot, treatment, month) %>%
  summarise(hum_avg_se = sd(hum_avg)/sqrt(length(hum_avg)),
            hum_min_se = sd(hum_min)/sqrt(length(hum_min)),
            hum_max_se = sd(hum_max)/sqrt(length(hum_max)),
            hum_dur_se = sd(hum_dur)/sqrt(length(hum_dur)),
            dew_int_se = sd(dew_intensity)/sqrt(length(dew_intensity)),
            dew_int2_se = sd(dew_intensity2)/sqrt(length(dew_intensity2)),
            hum_avg = mean(hum_avg),
            hum_min = mean(hum_min),
            hum_max = mean(hum_max),
            hum_dur = mean(hum_dur),
            dew_intensity = mean(dew_intensity),
            dew_intensity2 = mean(dew_intensity2),
            dew_days = as.numeric(hum_max == 1)) %>%
  ungroup() %>%
  full_join(dy_temp_dat %>%
              group_by(site, plot, treatment, month) %>%
              summarise(temp_avg_se = sd(temp_avg)/sqrt(length(temp_avg)),
                        temp_min_se = sd(temp_min)/sqrt(length(temp_min)),
                        temp_max_se = sd(temp_max)/sqrt(length(temp_max)),
                        hind_max_se = sd(hind_max)/sqrt(length(hind_max)),
                        temp_avg = mean(temp_avg),
                        temp_min = mean(temp_min),
                        temp_max = mean(temp_max),
                        hind_max = mean(hind_max)) %>%
              ungroup())


#### weather metric correlations ####

# correlation matrix
filter(mo_dat) %>%
  select(c(hum_avg:dew_days, temp_avg:hind_max)) %>%
  GGally::ggpairs()

# hum_avg: hum_min and hum_dur (0.8) and other hum metrics
# hum_min: hum_dur (0.5) and temp_min (0.6)
# hum_max: hum_dur (0.7), both dew_intensities (0.6)
# hum_dur: both dew_intensities (0.8)
# dew_intensity: dew_intensity2 (1), dew_days (0.5)
# dew_intensity2: dew_days (0.5)
# temp_avg: temp_min (1), temp_max (0.8), hind_max (0.9)
# temp_min: temp_max (0.7) and hind_max (0.8)
# temp_max: hind_max (0.9)

# metric suites
# both dew_intensities, dew_days, hum_dur, hum_max, hum_avg
# temp_avg, temp_min, temp_max, hind_max


#### add non-weather info ####

dat <- left_join(mo_dat, bio) %>%
  left_join(dens_dat) %>%
  left_join(covariates) %>%
  mutate(month_name = recode(month, "early_aug" = "Early August", "late_aug" = "Late August", "sep" = "September", "oct" = "October"),
         Mv_pc_biomass = Mv_biomass.g / Mv_density,
         Ev_pc_biomass = Ev_biomass.g / Ev_density,
         total_pc_biomass = total_biomass.g / total_density)


#### canopy cover correlation ####

# Chris did most already in canopy_cover_temp_humidity_relationships

# new dew intensity metric
cor.test(dat$dew_intensity2, dat$canopy_cover.prop)
# not correlated


#### average temperature stats ####

# models
tavg_totb_mod <- glmmTMB(temp_avg ~ total_biomass.g * month_name + (1|site), data = dat, family = "gaussian")
tavg_sepb_mod <- glmmTMB(temp_avg ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "gaussian")
tavg_totd_mod <- glmmTMB(temp_avg ~ total_density * month_name + (1|site), data = dat, family = "gaussian")
tavg_sepd_mod <- glmmTMB(temp_avg ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "gaussian")
tavg_totp_mod <- glmmTMB(temp_avg ~ bg_present * month_name + (1|site), data = dat, family = "gaussian")
tavg_sepp_mod <- glmmTMB(temp_avg ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "gaussian")

# model comparison
AIC(tavg_totb_mod, tavg_sepb_mod, tavg_totd_mod, tavg_sepd_mod, tavg_totp_mod, tavg_sepp_mod)
# sep present
summary(tavg_sepp_mod)
# Mv reduces temperature
stepAIC(tavg_sepp_mod)
# keep full model
plot(simulateResiduals(tavg_sepp_mod))
# adding month to fixed effects instead of random improved this


#### min temperature stats ####

# models
tmin_totb_mod <- glmmTMB(temp_min ~ total_biomass.g * month_name + (1|site), data = dat, family = "gaussian")
tmin_sepb_mod <- glmmTMB(temp_min ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "gaussian")
tmin_totd_mod <- glmmTMB(temp_min ~ total_density * month_name + (1|site), data = dat, family = "gaussian")
tmin_sepd_mod <- glmmTMB(temp_min ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "gaussian")
tmin_totp_mod <- glmmTMB(temp_min ~ bg_present * month_name + (1|site), data = dat, family = "gaussian")
tmin_sepp_mod <- glmmTMB(temp_min ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "gaussian")

# model comparison
AIC(tmin_totb_mod, tmin_sepb_mod, tmin_totd_mod, tmin_sepd_mod, tmin_totp_mod, tmin_sepp_mod)
# total density
summary(tmin_totd_mod)
# not sig
stepAIC(tmin_totd_mod)
# keep full model
plot(simulateResiduals(tmin_totd_mod))
# sig deviation


#### max temperature stats ####

# models
tmax_totb_mod <- glmmTMB(temp_max ~ total_biomass.g * month_name + (1|site), data = dat, family = "gaussian")
tmax_sepb_mod <- glmmTMB(temp_max ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "gaussian")
tmax_totd_mod <- glmmTMB(temp_max ~ total_density * month_name + (1|site), data = dat, family = "gaussian")
tmax_sepd_mod <- glmmTMB(temp_max ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "gaussian")
tmax_totp_mod <- glmmTMB(temp_max ~ bg_present * month_name + (1|site), data = dat, family = "gaussian")
tmax_sepp_mod <- glmmTMB(temp_max ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "gaussian")

# model comparison
AIC(tmax_totb_mod, tmax_sepb_mod, tmax_totd_mod, tmax_sepd_mod, tmax_totp_mod, tmax_sepp_mod)
# sep biomass
summary(tmax_sepb_mod)
# not sig
stepAIC(tmax_sepb_mod)
# keep full model
plot(simulateResiduals(tmax_sepb_mod))


#### max temperature stats ####

# models
hind_totb_mod <- glmmTMB(hind_max ~ total_biomass.g * month_name + (1|site), data = dat, family = "gaussian")
hind_sepb_mod <- glmmTMB(hind_max ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "gaussian")
hind_totd_mod <- glmmTMB(hind_max ~ total_density * month_name + (1|site), data = dat, family = "gaussian")
hind_sepd_mod <- glmmTMB(hind_max ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "gaussian")
hind_totp_mod <- glmmTMB(hind_max ~ bg_present * month_name + (1|site), data = dat, family = "gaussian")
hind_sepp_mod <- glmmTMB(hind_max ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "gaussian")

# model comparison
AIC(hind_totb_mod, hind_sepb_mod, hind_totd_mod, hind_sepd_mod, hind_totp_mod, hind_sepp_mod)
# sep biomass
summary(hind_sepb_mod)
# not sig
stepAIC(hind_sepb_mod)
# keep full model
plot(simulateResiduals(hind_sepb_mod))


#### average humidity stats ####

# models
havg_totb_mod <- glmmTMB(hum_avg ~ total_biomass.g * month_name + (1|site), data = dat, family = "beta_family")
# can ignore warning if model converges
havg_sepb_mod <- glmmTMB(hum_avg ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "beta_family")
havg_totd_mod <- glmmTMB(hum_avg ~ total_density * month_name + (1|site), data = dat, family = "beta_family")
havg_sepd_mod <- glmmTMB(hum_avg ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "beta_family")
havg_totp_mod <- glmmTMB(hum_avg ~ bg_present * month_name + (1|site), data = dat, family = "beta_family")
havg_sepp_mod <- glmmTMB(hum_avg ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "beta_family")

# model comparison
AIC(havg_totb_mod, havg_sepb_mod, havg_totd_mod, havg_sepd_mod, havg_totp_mod, havg_sepp_mod)
# total biomass
summary(havg_totb_mod)
# increase
stepAIC(havg_totb_mod)
# keep full model
plot(simulateResiduals(havg_totb_mod))


#### min humidity stats ####

# models
hmin_totb_mod <- glmmTMB(hum_min ~ total_biomass.g * month_name + (1|site), data = dat, family = "beta_family")
hmin_sepb_mod <- glmmTMB(hum_min ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "beta_family")
hmin_totd_mod <- glmmTMB(hum_min ~ total_density * month_name + (1|site), data = dat, family = "beta_family")
hmin_sepd_mod <- glmmTMB(hum_min ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "beta_family")
hmin_totp_mod <- glmmTMB(hum_min ~ bg_present * month_name + (1|site), data = dat, family = "beta_family")
hmin_sepp_mod <- glmmTMB(hum_min ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "beta_family")

# model comparison
AIC(hmin_totb_mod, hmin_sepb_mod, hmin_totd_mod, hmin_sepd_mod, hmin_totp_mod, hmin_sepp_mod)
# separate biomass
summary(hmin_sepb_mod)
# not sig
stepAIC(hmin_sepb_mod)
# keep full model
plot(simulateResiduals(hmin_sepb_mod))


#### max humidity stats ####

# reduce humidity by a small amount (values of 1 not allowed)
dat_max <- dat %>%
  mutate(hum_max = hum_max - 1e-7)

# models
hmax_totb_mod <- glmmTMB(hum_max ~ total_biomass.g * month_name + (1|site), data = dat_max, family = "beta_family")
hmax_sepb_mod <- glmmTMB(hum_max ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat_max, family = "beta_family")
hmax_totd_mod <- glmmTMB(hum_max ~ total_density * month_name + (1|site), data = dat_max, family = "beta_family")
hmax_sepd_mod <- glmmTMB(hum_max ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat_max, family = "beta_family")
hmax_totp_mod <- glmmTMB(hum_max ~ bg_present * month_name + (1|site), data = dat_max, family = "beta_family")
hmax_sepp_mod <- glmmTMB(hum_max ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat_max, family = "beta_family")

# model comparison
AIC(hmax_totb_mod, hmax_sepb_mod, hmax_totd_mod, hmax_sepd_mod, hmax_totp_mod, hmax_sepp_mod)
# total biomass
summary(hmax_totb_mod)
# not sig
stepAIC(hmax_totb_mod)
# keep full model
plot(simulateResiduals(hmax_totb_mod))


#### humidity duration stats ####

# models
hdur_totb_mod <- glmmTMB(hum_dur ~ total_biomass.g * month_name + (1|site), data = dat, family = "gaussian")
hdur_sepb_mod <- glmmTMB(hum_dur ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "gaussian")
hdur_totd_mod <- glmmTMB(hum_dur ~ total_density * month_name + (1|site), data = dat, family = "gaussian")
hdur_sepd_mod <- glmmTMB(hum_dur ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "gaussian")
hdur_totp_mod <- glmmTMB(hum_dur ~ bg_present * month_name + (1|site), data = dat, family = "gaussian")
hdur_sepp_mod <- glmmTMB(hum_dur ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "gaussian")

# model comparison
AIC(hdur_totb_mod, hdur_sepb_mod, hdur_totd_mod, hdur_sepd_mod, hdur_totp_mod, hdur_sepp_mod)
# total biomass
summary(hdur_totb_mod)
# total biomass has a positive effect in September
stepAIC(hdur_totb_mod)
# keep both in model
plot(simulateResiduals(hdur_sepb_mod))


#### dew intensity stats ####

# models
dewi_totb_mod <- glmmTMB(dew_intensity2 ~ total_biomass.g * month_name + (1|site), data = dat, family = "gaussian")
dewi_sepb_mod <- glmmTMB(dew_intensity2 ~ (Ev_biomass.g + Mv_biomass.g) * month_name + (1|site), data = dat, family = "gaussian")
dewi_totd_mod <- glmmTMB(dew_intensity2 ~ total_density * month_name + (1|site), data = dat, family = "gaussian")
dewi_sepd_mod <- glmmTMB(dew_intensity2 ~ (Ev_density + Mv_density) * month_name + (1|site), data = dat, family = "gaussian")
dewi_totp_mod <- glmmTMB(dew_intensity2 ~ bg_present * month_name + (1|site), data = dat, family = "gaussian")
dewi_sepp_mod <- glmmTMB(dew_intensity2 ~ (Ev_present + Mv_present) * month_name + (1|site), data = dat, family = "gaussian")

# model comparison
AIC(dewi_totb_mod, dewi_sepb_mod, dewi_totd_mod, dewi_sepd_mod, dewi_totp_mod, dewi_sepp_mod)
# total biomass
summary(dewi_totb_mod)
# marginal total biomass has a positive effect in September
stepAIC(dewi_totb_mod)
# keep both in model
plot(simulateResiduals(dewi_sepb_mod))


#### per capita biomass stats ####

cor.test(~ hum_avg + total_pc_biomass, data = dat)$p.value

# correlations by plot treatment
dat %>%
  group_by(month_name, plot) %>%
  summarise(correlation = cor.test(~ hum_avg + total_pc_biomass)$estimate,
            p_value = cor.test(~ hum_avg + total_pc_biomass)$p.value) %>%
  filter(p_value < 0.05)

# figures
dat %>%
  filter(plot == 6) %>%
  ggplot(aes(hum_avg, total_pc_biomass, color = site)) +
  geom_point() +
  facet_wrap(~ month_name)
# per capita biomass decreases with increasing humidity

dat %>%
  filter(plot == 9) %>%
  ggplot(aes(hum_avg, total_pc_biomass, color = site)) +
  geom_point() +
  facet_wrap(~ month_name)
# per capita biomass increases with increasing humidity


#### figures ####

# presence/absence data
pres_sim_dat <- expand_grid(month_name = c("Early August", "Late August", "September", "October"), site = c("D1", "D2", "D3", "D4"), Ev_present = c(0, 1), Mv_present = c(0, 1)) %>%
  filter(!(Ev_present == 1 & Mv_present == 1)) %>%
  mutate(month_name = fct_relevel(month_name, "Early August", "Late August", "September", "October"),
         pred = predict(tavg_sepp_mod, newdata = .),
         pred_se = predict(tavg_sepp_mod, newdata = ., se.fit = T)$se.fit,
         background = case_when(Ev_present == 0 & Mv_present == 0 ~ "none",
                                Ev_present == 1 & Mv_present == 0 ~ "Ev",
                                Ev_present == 0 & Mv_present == 1 ~ "Mv") %>%
           fct_relevel("none", "Ev", "Mv"))

# edit raw data
dat_pres <- dat %>%
  mutate(background = case_when(Ev_present == 0 & Mv_present == 0 ~ "none",
                                Ev_present == 1 & Mv_present == 0 ~ "Ev",
                                Ev_present == 0 & Mv_present == 1 ~ "Mv") %>%
           fct_relevel("none", "Ev", "Mv"))

# average temperature
tavg_plot <- ggplot(pres_sim_dat, aes(x = background, y = pred)) +
  geom_point(data = dat_pres, alpha = 0.3, aes(y = temp_avg, color = site)) +
  geom_errorbar(width = 0.05, aes(ymin = pred - pred_se, ymax = pred + pred_se, group = site)) +
  geom_point(size = 3, shape = 21, aes(fill = site)) +
  facet_wrap(~month_name, scales = "free") +
  temp_theme +
  ylab("Average daily temperature (°C)") +
  xlab("Background plants")

# # total density data
# totd_sim_dat <- expand_grid(month = c("early_aug", "late_aug", "sep", "oct"), site = c("D1", "D2", "D3", "D4"), total_density = 7:71) %>%
#   mutate(month = fct_relevel(month, "early_aug", "late_aug", "sep", "oct"),
#          pred = predict(tmin_totd_mod, newdata = .),
#          pred_se = predict(tmin_totd_mod, newdata = ., se.fit = T)$se.fit,
#          month_name = recode(month, "early_aug" = "Early August", "late_aug" = "Late August", "sep" = "September", "oct" = "October"))
# 
# # minimum temperature
# tmin_plot <- ggplot(totd_sim_dat, aes(x = total_density, y = pred)) +
#   geom_point(data = dat, alpha = 0.3, aes(y = temp_min, color = site)) +
#   geom_ribbon(alpha = 0.5, aes(ymin = pred - pred_se, ymax = pred + pred_se, fill = site)) +
#   geom_line(aes(color = site)) +
#   facet_wrap(~month_name, scales = "free") +
#   temp_theme +
#   ylab("Minimum daily temperature (°C)") +
#   xlab(expression(paste("Total plant density (plants ", m^-1, ")", sep = "")))

# # separate biomass
# sepb_sim_dat <- dat %>%
#   group_by(month, month_name, site) %>%
#   summarise(min_Ev_biomass.g = min(Ev_biomass.g),
#             max_Ev_biomass.g = max(Ev_biomass.g),
#             min_Mv_biomass.g = min(Mv_biomass.g),
#             max_Mv_biomass.g = max(Mv_biomass.g))  %>%
#   ungroup() %>%
#   rowwise() %>%
#   do(tibble(month = .$month, month_name = .$month_name, site = .$site,
#                 Ev_biomass.g = c(seq(.$min_Ev_biomass.g, .$max_Ev_biomass.g, length.out = 100), rep(0, 100)),
#                 Mv_biomass.g = c(rep(0, 100), seq(.$min_Mv_biomass.g, .$max_Mv_biomass.g, length.out = 100)))) %>%
#   ungroup() %>%
#   mutate(tmax_pred = predict(tmax_sepb_mod, newdata = .),
#          tmax_pred_se = predict(tmax_sepb_mod, newdata = ., se.fit = T)$se.fit,
#          hmin_pred = predict(hmin_sepb_mod, newdata = ., type = "response"),
#          hmin_pred_se = predict(hmin_sepb_mod, newdata = ., type = "response", se.fit = T)$se.fit,
#          biomass.g = Ev_biomass.g + Mv_biomass.g,
#          sp = case_when(Ev_biomass.g > 0 ~ "Ev",
#                         Mv_biomass.g > 0 ~ "Mv"))

# # maximum temperature
# tmax_plot <- ggplot(filter(sepb_sim_dat, sp == "Ev"), aes(x = biomass.g, y = tmax_pred)) +
#   geom_point(data = dat, alpha = 0.3, aes(x = Ev_biomass.g, y = temp_max, color = site)) +
#   geom_ribbon(alpha = 0.5, aes(ymin = tmax_pred - tmax_pred_se, ymax = tmax_pred + tmax_pred_se, fill = site)) +
#   geom_line(aes(color = site)) +
#   facet_wrap(~month_name, scales = "free") +
#   temp_theme +
#   ylab("Maximum daily temperature (°C)") +
#   xlab(expression(paste(italic("Elymus"), " biomass (g ", m^-1, ")", sep = "")))

# rounding function
scale_fun <- function(x){
  y = 100 * x
  sprintf("%.0f", y)
}
# 
# # minimum humidity
# hmin_plot <- ggplot(filter(sepb_sim_dat, sp == "Ev"), aes(x = biomass.g, y = hmin_pred)) +
#   geom_point(data = dat, alpha = 0.3, aes(x = Ev_biomass.g, y = hum_min, color = site)) +
#   geom_ribbon(alpha = 0.5, aes(ymin = hmin_pred - hmin_pred_se, ymax = hmin_pred + hmin_pred_se, fill = site)) +
#   geom_line(aes(color = site)) +
#   facet_wrap(~month_name, scales = "free") +
#   temp_theme +
#   ylab("Minimum daily relative humidity (%)") +
#   xlab(expression(paste(italic("Elymus"), " biomass (g ", m^-1, ")", sep = ""))) +
#   scale_y_continuous(labels = scale_fun)

# total biomass
totb_sim_dat <- dat %>%
  group_by(month_name, site) %>%
  summarise(min_biomass.g = min(total_biomass.g),
            max_biomass.g = max(total_biomass.g))  %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(month_name = .$month_name, site = .$site,
            total_biomass.g = seq(.$min_biomass.g, .$max_biomass.g, length.out = 100))) %>%
  ungroup() %>%
  mutate(havg_pred = predict(havg_totb_mod, newdata = ., type = "response"),
         havg_pred_se = predict(havg_totb_mod, newdata = ., type = "response", se.fit = T)$se.fit,
         hdur_pred = predict(hdur_totb_mod, newdata = ., type = "response"),
         hdur_pred_se = predict(hdur_totb_mod, newdata = ., type = "response", se.fit = T)$se.fit)

# average relative humidity
havg_plot <- ggplot(totb_sim_dat, aes(x = total_biomass.g, y = havg_pred)) +
  geom_point(data = dat, alpha = 0.3, aes(x = total_biomass.g, y = hum_avg, color = site)) +
  geom_ribbon(alpha = 0.5, aes(ymin = havg_pred - havg_pred_se, ymax = havg_pred + havg_pred_se, fill = site)) +
  geom_line(aes(color = site)) +
  facet_wrap(~month_name, scales = "free") +
  temp_theme +
  ylab("Average daily relative humidity (%)") +
  xlab(expression(paste("Total biomass (g ", m^-1, ")", sep = ""))) +
  scale_y_continuous(labels = scale_fun)

# humidity duration
hdur_plot <- ggplot(filter(totb_sim_dat, month_name == "September"), aes(x = total_biomass.g, y = hdur_pred)) +
  geom_point(data = dat, alpha = 0.3, aes(x = total_biomass.g, y = hum_dur, color = site)) +
  geom_ribbon(alpha = 0.5, aes(ymin = hdur_pred - hdur_pred_se, ymax = hdur_pred + hdur_pred_se, fill = site)) +
  geom_line(aes(color = site)) +
  facet_wrap(~month_name) +
  temp_theme +
  ylab(expression(paste("100% relative humidity duration (hours ", day^-1, ")", sep = ""))) +
  xlab(expression(paste("Total biomass (g ", m^-1, ")", sep = "")))

#### save figures ####
pdf("output/temp_humidity_analysis_2019_density_exp.pdf")
tavg_plot
havg_plot
hdur_plot
dev.off()


#### save daily data ####

day_dat <- hr_dat %>%
  group_by(site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            hum_dur = sum(hum_prop == 1)) %>%
  ungroup() 

write_csv(day_dat, "./intermediate-data/temp_humidity_daily_2019_density_exp.csv")


#### presentation figure ####

# colors
col_pal = c("#C0A76D", "#55A48B")

# add background info
dat2 <- dat %>%
  mutate(mv_background = case_when(plot > 1 & plot <= 4 ~ "Mv",
                                   TRUE ~ "not Mv"))

# biomass values by background
dat2 %>%
  filter(month_name == "September")  %>%
  group_by(mv_background) %>%
  summarise(min_bio = min(total_biomass.g),
            max_bio = max(total_biomass.g))

# total biomass
totb_fix_sim_dat <- dat %>%
  filter(month_name == "September") %>%
  summarise(min_biomass.g = min(total_biomass.g),
            max_biomass.g = max(total_biomass.g)) %>%
  rowwise() %>%
  do(tibble(month_name = "September", 
            site = NA,
            total_biomass.g = seq(.$min_biomass.g, .$max_biomass.g, length.out = 100))) %>%
  ungroup() %>%
  mutate(havg_pred = predict(havg_totb_mod, newdata = ., re.form = NA, type = "response"),
         havg_pred_se = predict(havg_totb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

# average relative humidity
havg_pres_plot <- ggplot(totb_fix_sim_dat, aes(x = total_biomass.g, y = havg_pred)) +
  geom_point(data = filter(dat2, month_name == "September"), alpha = 0.3, color = col_pal[1], aes(x = total_biomass.g, y = hum_avg, shape = mv_background)) +
  geom_ribbon(alpha = 0.5, fill = col_pal[1], aes(ymin = havg_pred - havg_pred_se, ymax = havg_pred + havg_pred_se, fill = site)) +
  geom_line(color = col_pal[1]) +
  scale_shape_manual(values = c(15, 16), guide = F) +
  ylab("Average daily relative humidity (%)") +
  xlab(expression(paste("Total biomass (g ", m^-1, ")", sep = ""))) +
  scale_y_continuous(labels = scale_fun) +
  temp_theme

# save figure
pdf("output/average_humidity_total_biomass_figure_2019_density_exp.pdf", width = 3.5, height = 3.5)
havg_pres_plot
dev.off()


#### output ####
save(havg_totb_mod, file = "output/average_humidity_total_biomass_model_2019_density_exp.rda")
write_csv(dat, "intermediate-data/average_humidity_total_biomass_2019_density_exp.csv")