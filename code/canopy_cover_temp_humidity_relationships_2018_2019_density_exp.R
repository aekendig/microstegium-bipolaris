##### info ####

# file: canopy_cover_temp_humidity_relationships_2018_2019_dens_exp
# author: Chris Wojan (with large parts copied from Amy Kendig)
# date last edited: 8/18/20
# goal: look for relationships among canopy cover and temperature/humidity measures


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(weathermetrics)

# import data
dy_dat <- read_csv("./data/temp_humidity_daily_2019_density_exp.csv")
hr_dat <- read_csv("./data/temp_humidity_hourly_2019_density_exp.csv")
covariates <- read_csv("./data/covariates_2018_density_exp.csv")


#### process ####

# calculate dewpoint from temp and relative humidity for each hour's reading
hr_dat$dewpoint <- humidity.to.dewpoint(rh=hr_dat$rel_hum, t=hr_dat$temp,
                                        temperature.metric = "celsius")

# calculate heat index from temp and relative humidity for each hour's reading
hr_dat$heat_ind <- heat.index(t=hr_dat$temp, rh=hr_dat$rel_hum,
                              temperature.metric = "celsius")

# summarize temp, relative humidity, dewpoint, and heat index for each day
dy_dat <- hr_dat %>%
  group_by(site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop),
            hum_dur = sum(hum_prop == 1),
            dewp_avg = mean(dewpoint),
            dewp_min = min(dewpoint),
            dewp_max = max(dewpoint),
            hind_avg = mean(heat_ind),
            hind_min = min(heat_ind),
            hind_max = max(heat_ind)
            ) %>%
  ungroup()

# create a subset of the hourly data for hours at saturation
hr_dat_dew <- filter(hr_dat, hum_prop==1)

# summarize humidity by each day, using only hours at saturation
dy_dat_dew <- hr_dat_dew %>%
  group_by(site, plot, treatment, day) %>%
  summarise(dewp_avg_hum = mean(dewpoint),
            dewp_min_hum = min(dewpoint),
            dewp_max_hum = max(dewpoint)) %>%
  ungroup()

# merge the humidity summary variables at saturation with the full daily summary data
dy_dat <- merge(dy_dat, dy_dat_dew, by = c("site", "plot", "treatment", "day"), all.x=T)

# add range variables, dew intensity, month variable
dy_dat <- mutate(dy_dat,
                 temp_rng = temp_max - temp_min,
                 dewp_rng = dewp_max - dewp_min,
                 hind_rng = hind_max - hind_min,
                 dewp_rng_hum = dewp_max_hum - dewp_min_hum,
                 dew_intensity = dewp_rng_hum * hum_dur,
                 saturated = hum_max,
                 month = case_when(day < as.Date("2019-07-29") ~ "early_aug",
                                   day >= as.Date("2019-07-29") & day < as.Date("2019-08-28") ~ "late_aug",
                                   day >= as.Date("2019-08-28") & day < as.Date("2019-09-24") ~ "sep",
                                   TRUE ~ "oct") %>%
                   fct_relevel("early_aug", "late_aug", "sep", "oct"))

# for days without saturation, set dew intensity and dewpoint range at saturation to 0
dy_dat[is.na(dy_dat$dew_intensity), "dew_intensity"] <- 0
dy_dat[is.na(dy_dat$dewp_rng_hum), "dewp_rng_hum"] <- 0

# mark days without saturation
dy_dat[dy_dat$saturated<1, "saturated"] <- 0

# summarize temp, humidity, dewpoint, heat index, dew measures, by month
mo_dat <- dy_dat %>%
  group_by(site, plot, treatment, month) %>%
  summarise(temp_avg_se = sd(temp_avg)/sqrt(length(temp_avg)),
            temp_min_se = sd(temp_min)/sqrt(length(temp_min)),
            temp_max_se = sd(temp_max)/sqrt(length(temp_max)),
            temp_rng_se = sd(temp_rng)/sqrt(length(temp_rng)),
            hum_avg_se = sd(hum_avg)/sqrt(length(hum_avg)),
            hum_min_se = sd(hum_min)/sqrt(length(hum_min)),
            hum_max_se = sd(hum_max)/sqrt(length(hum_max)),
            dewp_avg_se = sd(dewp_avg)/sqrt(length(dewp_avg)),
            dewp_min_se = sd(dewp_min)/sqrt(length(dewp_min)),
            dewp_max_se = sd(dewp_max)/sqrt(length(dewp_max)),
            dewp_rng_se = sd(dewp_rng)/sqrt(length(dewp_rng)),
            dew_int_se = sd(dew_intensity)/sqrt(length(dew_intensity)),
            dew_rng_hum_se = sd(dewp_rng_hum)/sqrt(length(dewp_rng_hum)),
            hind_avg_se = sd(hind_avg)/sqrt(length(hind_avg)),
            hind_min_se = sd(hind_min)/sqrt(length(hind_min)),
            hind_max_se = sd(hind_max)/sqrt(length(hind_max)),
            hind_rng_se = sd(hind_rng)/sqrt(length(hind_rng)),
            temp_avg = mean(temp_avg),
            temp_min = mean(temp_min),
            temp_max = mean(temp_max),
            temp_rng = mean(temp_rng),
            hum_avg = mean(hum_avg),
            hum_min = mean(hum_min),
            hum_max = mean(hum_max),
            hum_dur = mean(hum_dur),
            hum_dur_sum = sum(hum_dur),
            dewp_avg = mean(dewp_avg),
            dewp_min = mean(dewp_min),
            dewp_max = mean(dewp_max),
            dewp_rng = mean(dewp_rng),
            dew_intensity = mean(dew_intensity),
            dewp_rng_hum = mean(dewp_rng_hum),
            hind_avg = mean(hind_avg),
            hind_max = mean(hind_max),
            hind_min = mean(hind_min),
            hind_rng = mean(hind_rng),
            dew_days = sum(saturated)
            ) %>%
  ungroup()

# merge canopy cover data with monthly temp/humidity data
mo_cov <- merge(mo_dat, covariates, by = c("site", "plot", "treatment"))


#### analyze ####

# create a vector of variable names for looping
dep_vars <- c("hum_dur","dew_days", "dew_intensity", "temp_max", "temp_min", "hind_max", "hum_avg")

# create an empty data frame to save correlation results
cc_cors <- data.frame()

# loop through seven variables of interest, and save their correlation with canopy cover
for(i in dep_vars){
  icors <- by(mo_cov, mo_cov$month, function(mo_cov) cor.test(mo_cov$canopy_cover.prop,mo_cov[,i]))
  idat <- data.frame()
  for(j in levels(mo_cov$month)){
    jdat <- data.frame(var=i, month=j, estimate=icors[[j]]$estimate, p.value=icors[[j]]$p.value)
    idat <- rbind(idat, jdat)
  }
  cc_cors <- rbind(cc_cors, idat)
}

# round the correlation test results for graphical display
cc_cors$estimate <- round(cc_cors$estimate, digits = 3)
cc_cors$p.value <- round(cc_cors$p.value, digits = 4)


#### visualize ####

# set theme to Amy's style
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

# plot variables of interest by canopy cover, faceting by month

ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=hum_dur))+
  geom_text(data=filter(cc_cors, var=="hum_dur"), aes(x=0.6,y=12,label=paste0("cor = ",estimate,"\n p = ",p.value)),  color="blue")+
  facet_wrap(vars(month)) +
  plot_theme


ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=dew_days))+
  geom_text(data=filter(cc_cors, var=="dew_days"), aes(x=0.6,y=26,label=paste0("cor = ",estimate,"\n p = ",p.value)), color="blue")+
  facet_wrap(vars(month))+
  plot_theme

ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=dew_intensity))+
  geom_text(data=filter(cc_cors, var=="dew_intensity"), aes(x=0.6,y=50,label=paste0("cor = ",estimate,"\n p = ",p.value)), color="blue")+
  facet_wrap(vars(month))+
  plot_theme


ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=temp_min))+
  geom_text(data=filter(cc_cors, var=="temp_min"), aes(x=0.6,y=13,label=paste0("cor = ",estimate,"\n p = ",p.value)), color="blue")+
  facet_wrap(vars(month))+
  plot_theme

ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=temp_max))+
  geom_text(data=filter(cc_cors, var=="temp_max"), aes(x=0.6,y=25,label=paste0("cor = ",estimate,"\n p = ",p.value)), color="blue")+
  facet_wrap(vars(month))+
  plot_theme

ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=hind_max))+
  geom_text(data=filter(cc_cors, var=="hind_max"), aes(x=0.6,y=32,label=paste0("cor = ",estimate,"\n p = ",p.value)), color="blue")+
  facet_wrap(vars(month))+
  plot_theme

ggplot(data=mo_cov)+
  geom_point(aes(x=canopy_cover.prop, y=hum_avg))+
  geom_text(data=filter(cc_cors, var=="hum_avg"), aes(x=0.6,y=.93,label=paste0("cor = ",estimate,"\n p = ",p.value)), color="blue")+
  facet_wrap(vars(month))+
  plot_theme
