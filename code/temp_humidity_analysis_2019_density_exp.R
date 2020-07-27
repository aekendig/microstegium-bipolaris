##### info ####

# file: temp_humidity_analysis_2019_dens_exp
# author: Amy Kendig
# date last edited: 7/20/20
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

# import data
hr_dat <- read_csv("./intermediate-data/temp_humidity_hourly_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


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


#### monthly summaries ####

# from hourly data
mo_hr_dat <- hr_dat %>%
  mutate(month = case_when(day < as.Date("2019-07-29") ~ "early_aug",
                           day >= as.Date("2019-07-29") & day < as.Date("2019-08-28") ~ "late_aug",
                           day >= as.Date("2019-08-28") & day < as.Date("2019-09-24") ~ "sep",
                           TRUE ~ "oct") %>%
           fct_relevel("early_aug", "late_aug", "sep", "oct")) %>%
  group_by(site, plot, treatment, month) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            temp_avg_se = sd(temp)/sqrt(length(temp)),
            hum_avg_se = sd(hum_prop)/sqrt(length(hum_prop)),
            hum_dur = sum(hum_prop == 1)) %>%
  ungroup() %>%
  mutate(time = "hourly")

# from daily calculations
mo_dy_dat <- hr_dat %>%
  group_by(site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            hum_dur = sum(hum_prop == 1)) %>%
  ungroup() %>%
  mutate(month = case_when(day < as.Date("2019-07-29") ~ "early_aug",
                           day >= as.Date("2019-07-29") & day < as.Date("2019-08-28") ~ "late_aug",
                           day >= as.Date("2019-08-28") & day < as.Date("2019-09-24") ~ "sep",
                           TRUE ~ "oct") %>%
           fct_relevel("early_aug", "late_aug", "sep", "oct")) %>%
  group_by(site, plot, treatment, month) %>%
  summarise(temp_avg_se = sd(temp_avg)/sqrt(length(temp_avg)),
            temp_min_se = sd(temp_min)/sqrt(length(temp_min)),
            temp_max_se = sd(temp_max)/sqrt(length(temp_max)),
            hum_avg_se = sd(hum_avg)/sqrt(length(hum_avg)),
            hum_min_se = sd(hum_min)/sqrt(length(hum_min)),
            hum_max_se = sd(hum_max)/sqrt(length(hum_max)),
            temp_avg = mean(temp_avg),
            temp_min = mean(temp_min),
            temp_max = mean(temp_max),
            hum_avg = mean(hum_avg),
            hum_min = mean(hum_min),
            hum_max = mean(hum_max),
            hum_dur = mean(hum_dur)) %>%
  ungroup() %>%
  mutate(time = "daily")

# combine dataframes
mo_dat <- full_join(mo_hr_dat, mo_dy_dat) %>%
  left_join(plots)


#### visualizations to compare monthly and daily ####

# save theme
plot_theme <- theme_bw() +
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

# average temperature
avgt_plot <- ggplot(mo_dat, aes(x = background_density, 
                                y = temp_avg,
                                ymin = temp_avg - temp_avg_se,
                                ymax = temp_avg + temp_avg_se, 
                                fill = time)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.3)) +
  geom_errorbar(width = 0.1, position = position_dodge(0.3)) +
  facet_grid(month ~ background, scales = "free_x", switch = "both") +
  scale_fill_manual(values = c("black", "white")) +
  plot_theme

avgt_plot
# equivalent

# min temperature
avgt_plot %+%
  aes(y = temp_min,
      ymin = temp_min - temp_min_se,
      ymax = temp_min + temp_min_se)
# hourly is lower

# humidity duration
ggplot(mo_dat, aes(x = background_density, 
                   y = hum_dur,
                   fill = time)) +
  geom_point(size = 2, shape = 21, position = position_dodge(0.3)) +
  facet_grid(month ~ background, scales = "free_x", switch = "both") +
  scale_fill_manual(values = c("black", "white")) +
  plot_theme
# daily makes more sense


#### monthly averages ####

mo_dy_dat %>%
  filter(plot == 1) %>%
  group_by(month) %>%
  summarise(average_temp = mean(temp_avg),
            min_temp = mean(temp_min),
            max_temp = mean(temp_max),
            hum_dur = mean(hum_dur))


#### edit daily data ####

dat <- left_join(mo_dy_dat, bio) %>%
  left_join(dens_dat) %>%
  mutate(month_name = recode(month, "early_aug" = "Early August", "late_aug" = "Late August", "sep" = "September", "oct" = "October"),
         Mv_pc_biomass = Mv_biomass.g / Mv_density,
         Ev_pc_biomass = Ev_biomass.g / Ev_density,
         total_pc_biomass = total_biomass.g / total_density)


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
  plot_theme +
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
#   plot_theme +
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
#   plot_theme +
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
#   plot_theme +
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
  plot_theme +
  ylab("Average daily relative humidity (%)") +
  xlab(expression(paste("Total biomass (g ", m^-1, ")", sep = ""))) +
  scale_y_continuous(labels = scale_fun)

# humidity duration
hdur_plot <- ggplot(filter(totb_sim_dat, month_name == "September"), aes(x = total_biomass.g, y = hdur_pred)) +
  geom_point(data = dat, alpha = 0.3, aes(x = total_biomass.g, y = hum_dur, color = site)) +
  geom_ribbon(alpha = 0.5, aes(ymin = hdur_pred - hdur_pred_se, ymax = hdur_pred + hdur_pred_se, fill = site)) +
  geom_line(aes(color = site)) +
  facet_wrap(~month_name) +
  plot_theme +
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
  plot_theme

# save figure
pdf("output/average_humidity_total_biomass_figure_2019_density_exp.pdf", width = 3.5, height = 3.5)
havg_pres_plot
dev.off()


#### output ####
save(havg_totb_mod, file = "output/average_humidity_total_biomass_model_2019_density_exp.rda")
write_csv(dat, "intermediate-data/average_humidity_total_biomass_2019_density_exp.csv")