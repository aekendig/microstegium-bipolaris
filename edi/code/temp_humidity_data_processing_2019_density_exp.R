##### outputs ####

# temp_humidity_monthly_2019_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(caTools)
library(weathermetrics)

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


#### import data ####

# make a list of the file names
file_names = list.files(path = "data/temp-humidity-density-exp-20191025", 
                        pattern = "[.]csv")

# Import files, all in one list
dat = lapply(paste("data/temp-humidity-density-exp-20191025/", file_names, sep = ""), read.csv, header=T)


#### edit data ####

# assign plot to files
for(i in 1:length(file_names)){
  dat[[i]]$id <- sub(".csv", "", file_names[i])
}

# make list into data table
# add and remove columns
# remove data before placement and after collection
# sensors were placed by 3:00 pm on 7/2/19 and collected after 10:00 am on 10/22/19
dat2 <- do.call(rbind.data.frame, dat) %>%
  as_tibble() %>%
  rename(temp = "Temp...C..c.1", rel_hum = "RH.....c.1.2") %>%
  mutate(site = substr(id, 1, 2),
         plot = substr(id, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.factor(),
         treatment = "water",
         time = case_when(site == "D3" & plot == 6 ~ as.POSIXlt(Date.Time, format = "%m/%d/%y %H:%M"),
                          TRUE ~ as.POSIXlt(Date.Time, format = "%m/%d/%y %H:%M:%S")),
         day = case_when(site == "D3" & plot == 6 ~ as.Date(Date.Time, format = "%m/%d/%y %H:%M"),
                         TRUE ~ as.Date(Date.Time, format = "%m/%d/%y %H:%M:%S")),
         hum_prop = rel_hum / 100) %>%
#  select(-c(Date.Time)) %>%
  filter(time >= "2019-07-03 00:00:00" & time <= "2019-10-21 23:59:59")

# check days and times
dat2 %>%
  group_by(day) %>%
  summarise(min_time = min(time),
            max_time = max(time))

# check for missing data
filter(dat2, is.na(temp) | is.na(hum_prop)) 
filter(dat2, site == "D3" & plot == 6 & day == as.Date("2019-07-23")) %>% data.frame()
filter(dat2, site == "D3" & plot == 6 & day == as.Date("2019-09-12")) %>% data.frame()
# one hour missing from each day

# remove missing data
dat3 <- dat2 %>%
  filter(!is.na(temp) & !is.na(hum_prop))


#### daily summaries ####

day_dat <- dat3 %>%
  group_by(id, site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            hrs_10_35 = sum(temp > 10 & temp < 35))

ggplot(day_dat, aes(x = day, y = temp_avg, color = site)) +
  geom_point() +
  facet_wrap(~ plot) +
  temp_theme
  
ggplot(day_dat, aes(x = day, y = hum_avg, color = site)) +
  geom_point() +
  facet_wrap(~ plot) +
  temp_theme


#### humidity and time ####

# times that humidity extremes occur
hum_ext_dat <-  dat3 %>%
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
ggplot(hum_ext_dat, aes(x = hour)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 15, color = "red", linetype = "dashed") +
  facet_wrap(~ time_type) +
  xlab("Hour of the day") +
  ylab("Count of days separated by plot") +
  temp_theme


#### dew point, heat index, saturation time ####

# add new metrics
dat4 <- dat3 %>%
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
spline_dat <- dat4 %>%
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
spline_dat %>%
  filter(site == "D1" & plot == 1 & day_hum == "2019-08-20") %>%
  ggplot(aes(x = hour_hum, y = hum_prop)) +
  geom_point() +
  geom_line(aes(y = pred_hum)) +
  xlab("Hour of the day") +
  ylab("Relative humidity") +
  temp_theme

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
dat5 <- dat4 %>%
  left_join(hr_hum %>%
              select(-hours_sat)) %>%
  mutate(pred_sat = ifelse(hour_hum >= min_sat_time & hour_hum < max_sat_time, 1, 0))

# range of humidities within effective saturation
dat5 %>%
  filter(pred_sat == 1) %>%
  summarise(max_hum = max(hum_prop),
            min_hum = min(hum_prop),
            min_hr_hum = min(min_sat_time),
            max_hr_hum = max(max_sat_time))

# figure
dat5 %>%
  filter(pred_sat == 1) %>%
  ggplot(aes(hour_hum, hum_prop, color = as.factor(day_hum))) +
  geom_line() +
  facet_grid(site ~ plot) +
  theme_bw() +
  theme(legend.position = "none")


#### summarize weather data ####

# summarize dew point by each day, using only hours at saturation
dy_dat_dew <- filter(dat5, pred_sat == 1) %>%
  group_by(site, plot, treatment, day_hum) %>%
  summarise(dewp_min_hum = min(dewpoint),
            dewp_max_hum = max(dewpoint)) %>%
  ungroup() %>%
  mutate(dewp_rng_hum = dewp_max_hum - dewp_min_hum)

# daily summary for humidity data
dy_hum_dat <- dat5 %>%
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
dy_temp_dat <- dat5 %>%
  group_by(site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hind_max = max(heat_ind),
            hrs_10_35 = sum(temp > 10 & temp < 35)) %>%
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
                        hrs_10_35_se = sd(hrs_10_35)/sqrt(length(hrs_10_35)),
                        temp_avg = mean(temp_avg),
                        temp_min = mean(temp_min),
                        temp_max = mean(temp_max),
                        hind_max = mean(hind_max),
                        hrs_10_35 = mean(hrs_10_35)) %>%
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



#### save data #### 
write_csv(mo_dat, "intermediate-data/temp_humidity_monthly_2019_density_exp.csv")
