##### info ####

# file: severity_temporal_trends_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 10/28/21
# goal: temporal trends in disease severity and environmental variables


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lubridate)
library(zoo) # rollmean
library(cowplot)

# import data
# d1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
# didn't end up using above -- too few data points (3)
d2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R
# envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") 
# temp_humidity_data_processing_2019_density_exp
# didn't end up using above -- too few data points (3)
weatherD2Dat <- read_csv("data/BONWR_weather_2019.csv")


#### edit data ####

# summarize weather data by day
weatherD2Dat2 <- weatherD2Dat %>%
  filter(!is.na(TAIRGZ)) %>% # remove daily summaries
  mutate(date_time = mdy_hm(utc_valid),
         date = as_date(date_time)) %>%
  group_by(date) %>%
  summarize(temp_avg = mean(TAIRGZ),
            hum_avg = mean(XRIRGZ))

# rolling average, make long
weatherD2Dat3 <- weatherD2Dat2 %>%
  arrange(date) %>%
  mutate(temp = rollmean(temp_avg, k = 7, fill = NA),
         hum = rollmean(hum_avg, k = 7, fill = NA)) %>%
  select(-c(temp_avg, hum_avg)) %>%
  pivot_longer(cols = c(temp, hum),
               names_to = "var",
               values_to = "val") %>%
  mutate(var = fct_recode(var,
                          "humidity (%)" = "hum",
                          "temperature (F)" = "temp")) %>%
  filter(!is.na(val))
  

# format edge severity
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix) %>%
  select(month, site, plot, treatment, edge_severity) %>%
  group_by(site, plot, treatment) %>%
  mutate(months = n_distinct(month)) %>%
  ungroup() %>%
  mutate(max_months = max(months)) %>%
  filter(months == max_months) %>% # plot must have all months
  mutate(date = fct_recode(month,
                           "2019-05-07" = "may",
                           "2019-06-04" = "jun",
                           "2019-07-02" = "jul",
                           "2019-07-29" = "early_aug",
                           "2019-08-28" = "late_aug",
                           "2019-09-24" = "sep") %>%
           as_date(),
         plot_ID = paste(site, plot, treatment, sep = "_"))

# format focal severity
d2Dat2 <- d2Dat %>%
  filter(treatment != "fungicide") %>%
  group_by(site, plot, treatment) %>%
  mutate(months = n_distinct(month)) %>%
  ungroup() %>%
  mutate(max_months = max(months)) %>%
  filter(months == max_months) %>% # plot must have all months
  mutate(date = fct_recode(month,
                           "2019-05-07" = "may",
                           "2019-06-04" = "jun",
                           "2019-07-02" = "jul",
                           "2019-07-29" = "early_aug",
                           "2019-08-28" = "late_aug") %>%
           as_date(),
         plot_ID = paste(site, plot, treatment, sep = "_"),
         plant_group = paste(sp, age))

n_distinct(d2Dat2$plot_ID) # 35 plots

#### visualize ####
# 
# ggplot(edgeSevD2Dat2, aes(x = date, y = edge_severity)) +
#   stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
#   stat_summary(geom = "line", fun = "mean") +
#   stat_summary(geom = "point", size = 2, fun = "mean")

# settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, 'cm'))

# transform data
weatherFig <- ggplot(filter(weatherD2Dat3, date >= "2019-05-07" & date <= "2019-08-28"), 
                     aes(x = date, y = val, color = var)) +
  geom_line() +
  labs(y = "Environmental metric") +
  scale_color_manual(values = c("blue", "black")) +
  fig_theme +
  theme(legend.position = c(0.85, 0.2),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

severityFig <- ggplot(d2Dat2, aes(x = date, y = severity, color = plant_group)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  labs(x = "Date", y = "Disease severity") +
  fig_theme +
  theme(legend.position = c(0.13, 0.73))

pdf("output/severity_environment_time_series.pdf", width = 4, height = 3)
plot_grid(weatherFig, severityFig, nrow = 2)
dev.off()
