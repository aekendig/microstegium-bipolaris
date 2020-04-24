##### info ####

# file: temp_humidity_analysis_2019_dens_exp
# author: Amy Kendig
# date last edited: 4/23/20
# goal: evaluate the effects of density on temperature and humidity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest) # p-values for lme4
library(glmmTMB)
library(DHARMa) # plot residuals for glmmTMB

# import data
hr_dat <- read_csv("./intermediate-data/temp_humidity_hourly_2019_density_exp.csv")
dy_dat <- read_csv("./intermediate-data/temp_humidity_daily_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### monthly summaries ####

# from hourly data
mo_hr_dat <- hr_dat %>%
  mutate(month = month(time)) %>%
  group_by(site, plot, treatment, month) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            temp_avg_se = sd(temp)/sqrt(length(temp)),
            hum_avg_se = sd(hum_prop)/sqrt(length(hum_prop))) %>%
  mutate(time = "hourly")

# from daily data
mo_dy_dat <- dy_dat %>%
  mutate(month = month(day)) %>%
  group_by(site, plot, treatment, month) %>%
  summarise(temp_avg = mean(temp_avg),
            temp_min = mean(temp_min),
            temp_max = mean(temp_max),
            hum_avg = mean(hum_avg),
            hum_min = mean(hum_min),
            hum_max = mean(hum_max),
            temp_avg_se = sd(temp_avg)/sqrt(length(temp_avg)),
            temp_min_se = sd(temp_min)/sqrt(length(temp_min)),
            temp_max_se = sd(temp_max)/sqrt(length(temp_max)),
            hum_avg_se = sd(hum_avg)/sqrt(length(hum_avg)),
            hum_min_se = sd(hum_min)/sqrt(length(hum_min)),
            hum_max_se = sd(hum_max)/sqrt(length(hum_max))) %>%
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

# use daily stats
dat <- left_join(mo_dy_dat, plots) %>%
  mutate(month_name = month(month, label = T))


#### visualize ####

# avg temperature
davgt_plot <- ggplot(dat, aes(x = background_density, y = temp_avg)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_grid(month ~ background, scales = "free_x", switch = "both") +
  plot_theme
davgt_plot
# nothing obvious

# min temp
dmint_plot <- davgt_plot %+%
  aes(x = background_density, y = temp_min)
dmint_plot
# maybe positive effect of seedlings in October

# max temp
dmaxt_plot <- davgt_plot %+%
  aes(x = background_density, y = temp_max)
dmaxt_plot
# slight cooling with Mv in July and August

# avg humidity
davgh_plot <- davgt_plot %+%
  aes(x = background_density, y = hum_avg)
davgh_plot
# positive effect of seedlings

# min humidity
dminh_plot <- davgt_plot %+%
  aes(x = background_density, y = hum_min)
dminh_plot
# positive effect of Mv

# max humidity
dmaxh_plot <- davgt_plot %+%
  aes(x = background_density, y = hum_max)
dmaxh_plot
# positive effect of plants


#### average temperature stats ####

# Ev adult
tavg_ea_mod <- lmer(temp_avg ~ background_density * month + (1|site),
                  data = filter(dat, background == "Ev adult"))
# singular fit
summary(tavg_ea_mod)
# no variance explained by site
tavg_ea_mod <- lm(temp_avg ~ background_density * month,
                    data = filter(dat, background == "Ev adult"))
summary(tavg_ea_mod)
# month sig
tavg_ea_mod2 <- update(tavg_ea_mod, .~. -background_density:month)
summary(tavg_ea_mod2)
# month sig
# plot(tavg_ea_mod2)

# Ev seedling
tavg_es_mod <- lmer(temp_avg ~ background_density * month + (1|site),
                  data = filter(dat, background == "Ev seedling"))
# singular fit
summary(tavg_es_mod)
# no variance explained by site
tavg_es_mod <- lm(temp_avg ~ background_density * month ,
                    data = filter(dat, background == "Ev seedling"))
summary(tavg_es_mod)
# month sig
tavg_es_mod2 <- update(tavg_es_mod, .~. -background_density:month)
summary(tavg_es_mod2)
# month sig
# plot(tavg_es_mod2)

# Mv seedling
tavg_ms_mod <- lmer(temp_avg ~ background_density * month + (1|site),
                  data = filter(dat, background == "Mv seedling"))
# singular fit
summary(tavg_ms_mod)
# very little variance explained by site
tavg_ms_mod <- lm(temp_avg ~ background_density * month,
                    data = filter(dat, background == "Mv seedling"))
summary(tavg_ms_mod)
# month sig
tavg_ms_mod2 <- update(tavg_ms_mod, .~. -background_density:month)
summary(tavg_ms_mod2)
# month sig
# plot(tavg_ms_mod2)


#### min temperature stats ####

# Ev adult
tmin_ea_mod <- lmer(temp_min ~ background_density * month + (1|site),
                  data = filter(dat, background == "Ev adult"))
# singular fit
summary(tmin_ea_mod)
# very little variance explained by site
tmin_ea_mod <- lm(temp_min ~ background_density * month,
                    data = filter(dat, background == "Ev adult"))
summary(tmin_ea_mod)
# month sig
tmin_ea_mod2 <- update(tmin_ea_mod, .~. -background_density:month)
summary(tmin_ea_mod2)
# month sig
# plot(tmin_ea_mod2)

# Ev seedling
tmin_es_mod <- lmer(temp_min ~ background_density * month + (1|site),
                  data = filter(dat, background == "Ev seedling"))
# singular fit
summary(tmin_es_mod)
# very little variance explained by site
tmin_es_mod <- lm(temp_min ~ background_density * month,
                    data = filter(dat, background == "Ev seedling"))
summary(tmin_es_mod)
# month sig
tmin_es_mod2 <- update(tmin_es_mod, .~. -background_density:month)
summary(tmin_es_mod2)
# month sig
# plot(tmin_es_mod2)

# Mv seedling
tmin_ms_mod <- lmer(temp_min ~ background_density * month + (1|site),
                  data = filter(dat, background == "Mv seedling"))
# singular fit
summary(tmin_ms_mod)
# very little variance explained by site
tmin_ms_mod <- lm(temp_min ~ background_density * month,
                    data = filter(dat, background == "Mv seedling"))
summary(tmin_ms_mod)
# month sig
tmin_ms_mod2 <- update(tmin_ms_mod, .~. -background_density:month)
summary(tmin_ms_mod2)
# month sig
# plot(tmin_ms_mod2)


#### max temperature stats ####

# Ev adult
tmax_ea_mod <- lmer(temp_max ~ background_density * month + (1|site),
                  data = filter(dat, background == "Ev adult"))
summary(tmax_ea_mod)
# month sig
tmax_ea_mod2 <- update(tmax_ea_mod, .~. -background_density:month)
summary(tmax_ea_mod2)
# month sig
# plot(tmax_ea_mod2)

# Ev seedling
tmax_es_mod <- lmer(temp_max ~ background_density * month + (1|site),
                  data = filter(dat, background == "Ev seedling"))
summary(tmax_es_mod)
# month sig
tmax_es_mod2 <- update(tmax_es_mod, .~. -background_density:month)
summary(tmax_es_mod2)
# month sig
# plot(tmax_es_mod2)

# Mv seedling
tmax_ms_mod <- lmer(temp_max ~ background_density * month + (1|site),
                  data = filter(dat, background == "Mv seedling"))
summary(tmax_ms_mod)
# month sig
tmax_ms_mod2 <- update(tmax_ms_mod, .~. -background_density:month)
summary(tmax_ms_mod2)
# month sig
# plot(tmax_ms_mod2)


#### average humidity stats ####

# Ev adult
havg_ea_mod <- glmmTMB(hum_avg ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Ev adult"),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(havg_ea_mod)
# month sig
havg_ea_mod2 <- update(havg_ea_mod, .~. -background_density:month)
summary(havg_ea_mod2)
# month and density sig
# plot(simulateResiduals(havg_ea_mod2))
# high residuals for predicted

# Ev seedling
havg_es_mod <- glmmTMB(hum_avg ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Ev seedling"),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(havg_es_mod)
# month sig
havg_es_mod2 <- update(havg_es_mod, .~. -background_density:month)
summary(havg_es_mod2)
# month and density sig
# plot(simulateResiduals(havg_es_mod2))
# high residuals for predicted

# Mv seedling
havg_ms_mod <- glmmTMB(hum_avg ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Mv seedling"),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(havg_ms_mod)
# month sig
havg_ms_mod2 <- update(havg_ms_mod, .~. -background_density:month)
# this warning can be ignored as long as the model converges
summary(havg_ms_mod2)
# month and density sig
# plot(simulateResiduals(havg_ms_mod2))
# high residuals for predicted


#### min humidity stats ####

# Ev adult
hmin_ea_mod <- glmmTMB(hum_min ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Ev adult"),
                  family = beta_family)
summary(hmin_ea_mod)
# month sig
hmin_ea_mod2 <- update(hmin_ea_mod, .~. -background_density:month)
# this warning can be ignored as long as the model converges
summary(hmin_ea_mod2)
# month sig
# plot(simulateResiduals(hmin_ea_mod2))
# high residuals for predicted

# Ev seedling
hmin_es_mod <- glmmTMB(hum_min ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Ev seedling"),
                  family = beta_family)
summary(hmin_es_mod)
# month sig
hmin_es_mod2 <- update(hmin_es_mod, .~. -background_density:month)
summary(hmin_es_mod2)
# month sig
# plot(simulateResiduals(hmin_es_mod2))
# high residuals for predicted

# Mv seedling
hmin_ms_mod <- glmmTMB(hum_min ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Mv seedling"),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(hmin_ms_mod)
# month sig
hmin_ms_mod2 <- update(hmin_ms_mod, .~. -background_density:month)
summary(hmin_ms_mod2)
# month sig
# background marginally sig
# plot(simulateResiduals(hmin_ms_mod2))
# high residuals for predicted


#### max humidity stats ####

# Ev adult
hmax_ea_mod <- glmmTMB(hum_max ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Ev adult") %>%
                    mutate(hum_max = ifelse(hum_max == 1, 0.9999, hum_max)),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(hmax_ea_mod)
# month sig
hmax_ea_mod2 <- update(hmax_ea_mod, .~. -background_density:month)
summary(hmax_ea_mod2)
# month sig
# background marginally sig
# plot(simulateResiduals(hmax_ea_mod2))
# slight increase in residuals with predicted values

# Ev seedling
hmax_es_mod <- glmmTMB(hum_max ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Ev seedling") %>%
                    mutate(hum_max = ifelse(hum_max == 1, 0.9999, hum_max)),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(hmax_es_mod)
# month and density sig
hmax_es_mod2 <- update(hmax_es_mod, .~. -background_density:month)
# this warning can be ignored as long as the model converges
summary(hmax_es_mod2)
# month and density sig
# plot(simulateResiduals(hmax_es_mod2))
# slight increase in residuals with predicted values

# Mv seedling
hmax_ms_mod <- glmmTMB(hum_max ~ background_density *  month + (1|site),
                  data = filter(dat, background == "Mv seedling") %>%
                    mutate(hum_max = ifelse(hum_max == 1, 0.9999, hum_max)),
                  family = beta_family)
# this warning can be ignored as long as the model converges
summary(hmax_ms_mod)
# month sig
hmax_ms_mod2 <- update(hmax_ms_mod, .~. -background_density:month)
# this warning can be ignored as long as the model converges
summary(hmax_ms_mod2)
# month sig
# plot(simulateResiduals(hmax_ms_mod2))
# slight decrease in residuals with predicted values


#### figures ####

# simulate data
dat_sim <- tibble(background = rep(c("Ev adult", "Ev seedling", "Mv seedling"), each = 400),
                  month = rep(rep(7:10, each = 100), 3),
                  background_density = c(rep(seq(0, 8, length.out = 100), 4),
                                         rep(seq(0, 16, length.out = 100), 4),
                                         rep(seq(0, 64, length.out = 100), 4))) %>%
  mutate(site = NA,
         month_name = month(month, label = T))

# add humdity predictions to data
dat_hum <- dat_sim %>%
  mutate(hum_avg = case_when(background == "Ev adult" ~ predict(havg_ea_mod2, newdata = ., re.form = NA, type = "response"),
                          background == "Ev seedling" ~ predict(havg_es_mod2, newdata = ., re.form = NA, type = "response"),
                          background == "Mv seedling" ~ predict(havg_ms_mod2, newdata = ., re.form = NA, type = "response")),
         hum_avg_se = case_when(background == "Ev adult" ~ predict(havg_ea_mod2, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                          background == "Ev seedling" ~ predict(havg_es_mod2, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                          background == "Mv seedling" ~ predict(havg_ms_mod2, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit))

# average humidity
hum_plot <- ggplot(dat_hum, aes(x = background_density, y = hum_avg, color = month_name, fill = month_name)) +
  geom_line() +
  geom_ribbon(aes(ymin = hum_avg - hum_avg_se, ymax = hum_avg + hum_avg_se), alpha = 0.5, color = NA) +
  stat_summary(data = dat, geom = "errorbar", fun.data = "mean_se", width = 0.1, color = "black") +
  stat_summary(data = dat, geom = "point", fun = "mean", size = 3, shape = 21, color = "black") +
  facet_grid(month_name ~ background, scales = "free_x", switch = "both") +
  plot_theme +
  xlab("Background density") +
  ylab("Monthly average relative humidity")
hum_plot
# it would be nice to have curves that can saturate
