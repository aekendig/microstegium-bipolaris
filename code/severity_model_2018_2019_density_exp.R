##### info ####

# file: severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 1/21/21
# goal: evaluate effects of transmission and environment on disease


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)
library(MuMIn)

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp.R
mvEdgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp.R

# import environmental variables
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") 
# covariate_data_processing_2018_density_exp
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") 
# temp_humidity_data_processing_2019_density_exp


#### edit data ####

# severity data
sevD1Datb <- sevD1Dat %>%
  filter(!(sp == "Ev" & site == "D2" & plot == 7 & treatment == "water" & ID == "A")) %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity)) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_severity",
              values_from = severity) %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(late_aug_severity = mean(late_aug_severity, na.rm = T),
            jul_severity = mean(jul_severity, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = sp,
              names_glue = "{sp}_{.value}",
              values_from = c(jul_severity, late_aug_severity))

sevD1Datc <- sevD1Datb %>%
  select(site, plot, treatment, Mv_late_aug_severity, Ev_jul_severity, Mv_jul_severity) %>%
  rename(severity = Mv_late_aug_severity) %>%
  mutate(sp = "Mv") %>%
  full_join(sevD1Datb %>%
              select(site, plot, treatment, Ev_late_aug_severity, Ev_jul_severity, Mv_jul_severity) %>%
              rename(severity = Ev_late_aug_severity) %>%
              mutate(sp = "Ev"))

sevD2Datb <- sevD2Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity)) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
  filter(ID %in% c("1", "2", "3", "A")) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_severity",
              values_from = severity) %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(early_aug_severity = mean(early_aug_severity, na.rm = T),
            jul_severity = mean(jul_severity, na.rm = T),
            jun_severity = mean(jun_severity, na.rm = T),
            may_severity = mean(may_severity, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = sp,
              names_glue = "{sp}_{.value}",
              values_from = c(may_severity, jun_severity, jul_severity, early_aug_severity))

sevD2Datc <- sevD2Datb %>%
  select(site, plot, treatment, Mv_early_aug_severity, Ev_jul_severity, Mv_jul_severity, Ev_jun_severity, Mv_jun_severity, Ev_may_severity) %>%
  rename(severity = Mv_early_aug_severity) %>%
  mutate(sp = "Mv") %>%
  full_join(sevD2Datb %>%
              select(site, plot, treatment, Ev_early_aug_severity, Ev_jul_severity, Mv_jul_severity, Ev_jun_severity, Mv_jun_severity, Ev_may_severity) %>%
              rename(severity = Ev_early_aug_severity) %>%
              mutate(sp = "Ev"))

edgeSevD2Dat <- mvEdgeSevD2Dat %>%
  filter(month != "sep") %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_edge_severity",
              values_from = edge_severity)

# environmental data
envDDat <- envD2Dat %>%
  select(site, plot, treatment, month, dew_intensity2, temp_avg) %>%
  rename(dew_intensity = dew_intensity2) %>%
  pivot_wider(values_from = c(dew_intensity, temp_avg),
              names_from = month,
              names_glue = "{month}_{.value}") %>%
  select(site, plot, treatment, early_aug_dew_intensity, late_aug_dew_intensity, early_aug_temp_avg, late_aug_temp_avg) %>%
  mutate_at(c("early_aug_dew_intensity", "late_aug_dew_intensity"), log) %>%
  full_join(envD1Dat)

# combine data
d1dat <- sevD1Datc %>%
  left_join(plotsD) %>%
  left_join(envD1Dat) %>%
  mutate(severity = asin(sqrt(severity)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  select(site, plot, sp, severity, treatment, fungicide, soil_moisture_jun.prop, canopy_cover.prop, mv_inf_jul.prop, Ev_jul_severity, Mv_jul_severity) %>%
  drop_na()

d2dat <- sevD2Datc %>%
  left_join(plotsD) %>%
  left_join(envDDat) %>%
  left_join(edgeSevD2Dat) %>%
  mutate(severity = asin(sqrt(severity)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  select(site, plot, sp, severity, treatment, fungicide, Ev_may_severity, Ev_jun_severity, Mv_jun_severity, Ev_jul_severity, Mv_jul_severity, soil_moisture_oct.prop, canopy_cover.prop, jul_edge_severity) %>%
  drop_na()

d2datw <- sevD2Datc %>%
  left_join(plotsD) %>%
  left_join(envDDat) %>%
  left_join(edgeSevD2Dat) %>%
  mutate(severity = asin(sqrt(severity)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  select(site, plot, sp, severity, treatment, fungicide, Ev_jul_severity, Mv_jul_severity, soil_moisture_oct.prop, canopy_cover.prop, jul_edge_severity, early_aug_dew_intensity, early_aug_temp_avg, soil_moisture_oct.prop) %>%
  drop_na()


#### visualize ####

# year 1
ggpairs(d1dat %>%
          filter(sp == "Mv") %>%
          select(treatment, soil_moisture_jun.prop, canopy_cover.prop, mv_inf_jul.prop, Ev_jul_severity, Mv_jul_severity))
# only sig is < 0.4

# year 2
ggpairs(d2dat %>%
          filter(sp == "Mv") %>%
          select(treatment, Ev_may_severity, Ev_jun_severity, Mv_jun_severity, Ev_jul_severity, Mv_jul_severity, soil_moisture_oct.prop, canopy_cover.prop, jul_edge_severity))
# mv and ev july severity 0.4
# mv july severity and canopy 0.5

ggpairs(d2datw %>%
          filter(sp == "Mv") %>%
          select(treatment, Ev_jul_severity, Mv_jul_severity, early_aug_dew_intensity, early_aug_temp_avg, soil_moisture_oct.prop, canopy_cover.prop, jul_edge_severity))
# mv and ev july severity 0.7
# ev july severity and temp 0.5
# ev july severity and canopy 0.4
# mv july severity and canopy 0.4
# dew and temp 0.5
# dew and soil moisture 0.7
# temp and soil moisture 0.5

# recombine data without correlated variables
d2dat2 <- sevD2Datc %>%
  left_join(plotsD) %>%
  left_join(envDDat) %>%
  left_join(edgeSevD2Dat) %>%
  mutate(severity = asin(sqrt(severity)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  select(site, plot, sp, severity, treatment, fungicide, Ev_may_severity, Ev_jun_severity, Mv_jun_severity, Ev_jul_severity, soil_moisture_oct.prop, canopy_cover.prop, jul_edge_severity) %>%
  drop_na()
# same sample size

d2datw2 <- sevD2Datc %>%
  left_join(plotsD) %>%
  left_join(envDDat) %>%
  left_join(edgeSevD2Dat) %>%
  mutate(severity = asin(sqrt(severity)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  select(site, plot, sp, severity, treatment, fungicide, canopy_cover.prop, jul_edge_severity, early_aug_dew_intensity, soil_moisture_oct.prop) %>%
  drop_na()
# same sample size

# recombine data without early months
d2dat3 <- sevD2Datc %>%
  left_join(plotsD) %>%
  left_join(envDDat) %>%
  left_join(edgeSevD2Dat) %>%
  mutate(severity = asin(sqrt(severity)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  select(site, plot, sp, severity, treatment, fungicide, Mv_jul_severity, Ev_jul_severity, soil_moisture_oct.prop, canopy_cover.prop, jul_edge_severity) %>%
  drop_na()


#### fit year 1 model ####

# must set for dredge function, default is "na.omit"
options(na.action = "na.fail")

d1mod1 <- lm(severity ~ sp * (fungicide + soil_moisture_jun.prop + canopy_cover.prop + mv_inf_jul.prop + Ev_jul_severity + Mv_jul_severity), 
           dat = d1dat)
summary(d1mod1)
# plot(d1mod1)
d1mod2 <- model.avg(get.models(dredge(d1mod1), subset = delta < 4))
summary(d1mod2) 
# fungicide has a larger negative effect on Mv than Ev
# Mv has lower severity than Ev
# soil moisture reduced severity for both species


#### fit year 2 models ####

d2mod1 <- lm(severity ~ sp * (fungicide + Ev_may_severity + Ev_jun_severity + Mv_jun_severity + Ev_jul_severity + soil_moisture_oct.prop + jul_edge_severity + canopy_cover.prop), 
             dat = d2dat2)
summary(d2mod1)
# plot(d2mod1)
d2mod2 <- model.avg(get.models(dredge(d2mod1), subset = delta < 4))
summary(d2mod2) 
# Ev july severity decreased Ev early August, but increased Mv
# Mv had lower severity than Ev
# fungicide reduced both

# same as when months earlier than July were left out, but soil moisture reduced severity (there were more samples)

# replace Ev july with Mv - same? need to remove canopy cover
d2mod3 <- lm(severity ~ sp * (fungicide + Ev_may_severity + Ev_jun_severity + Mv_jun_severity + Mv_jul_severity + soil_moisture_oct.prop + jul_edge_severity), 
             dat = d2dat)
summary(d2mod3)
# plot(d2mod3)
d2mod4 <- model.avg(get.models(dredge(d2mod3), subset = delta < 4))
summary(d2mod4) 
# yes, same pattern (decreased Ev early August, but increased Mv)

# model with both - need to remove canopy cover
d2mod5 <- lm(severity ~ sp * (fungicide + Ev_jul_severity + Mv_jul_severity + soil_moisture_oct.prop + jul_edge_severity), 
             dat = d2dat3)
summary(d2mod5)
# plot(d2mod5)
d2mod6 <- model.avg(get.models(dredge(d2mod5), subset = delta < 4))
summary(d2mod6)
# only Ev is significant

# water-only data with dew intensity
d2modw1 <- lm(severity ~ sp * (canopy_cover.prop + jul_edge_severity + early_aug_dew_intensity), 
             dat = d2datw2)
summary(d2modw1)
# plot(d2modw1)
d2modw2 <- model.avg(get.models(dredge(d2modw1), subset = delta < 4))
summary(d2modw2) 
# canopy increased Ev, but decreased Mv
# mv had higher severity

# water-only model with temp
d2modw3 <- lm(severity ~ sp * (canopy_cover.prop + jul_edge_severity + early_aug_temp_avg), 
              dat = d2datw)
summary(d2modw3)
# plot(d2modw3)
d2modw4 <- model.avg(get.models(dredge(d2modw3), subset = delta < 4))
summary(d2modw4) 
# canopy increased Ev, but decreased Mv


#### figures ####

# did this before May and June severity were added to year 2 model
# the figure looks like a cloud of points and the slopes are super slight

# soil data
d1soil <- d1dat %>%
  group_by(sp) %>%
  summarise(fungicide = 0,
            Ev_jul_severity = mean(Ev_jul_severity),
            Mv_jul_severity = mean(Mv_jul_severity),
            canopy_cover.prop = mean(canopy_cover.prop),
            mv_inf_jul.prop = mean(mv_inf_jul.prop)) %>%
  ungroup() %>%
  merge(tibble(soil_moisture_jun.prop = seq(min(d1dat$soil_moisture_jun.prop), 
                                         max(d1dat$soil_moisture_jun.prop),
                                         length.out = 100)),
        all = T) %>%
  mutate(year = "Year 1",
         pred_sev = predict(d1mod2, newdata = .),
         pred_sev.se = predict(d1mod2, newdata = ., se.fit = T)$se.fit) %>%
  rename(soil_moisture = soil_moisture_jun.prop) %>%
  as_tibble()

d2soil <- d2dat %>%
  group_by(sp) %>%
  summarise(fungicide = 0,
            Ev_jul_severity = mean(Ev_jul_severity),
            canopy_cover.prop = mean(canopy_cover.prop),
            jul_edge_severity = mean(jul_edge_severity)) %>%
  ungroup() %>%
  merge(tibble(soil_moisture_oct.prop = seq(min(d2dat$soil_moisture_oct.prop), 
                                            max(d2dat$soil_moisture_oct.prop),
                                            length.out = 100)),
        all = T) %>%
  mutate(year = "Year 2",
         pred_sev = predict(d2mod2, newdata = .),
         pred_sev.se = predict(d2mod2, newdata = ., se.fit = T)$se.fit) %>%
  rename(soil_moisture = soil_moisture_oct.prop) %>%
  as_tibble()

soildat <- full_join(d1soil, d2soil)

# soil moisture
d1dat %>%
  mutate(year = "Year 1") %>%
  rename(soil_moisture = soil_moisture_jun.prop) %>%
  full_join(d2dat %>%
              mutate(year = "Year 2") %>%
              rename(soil_moisture = soil_moisture_oct.prop)) %>%
  ggplot(aes(x = soil_moisture, y = severity, color = sp, fill = sp)) +
  geom_point() +
  geom_line(data = soildat, aes(y = pred_sev)) +
  geom_ribbon(data = soildat, aes(y = pred_sev, ymin = pred_sev - pred_sev.se, ymax = pred_sev + pred_sev.se), color = NA, alpha = 0.5) +
  facet_wrap(~ year, scales = "free") +
  theme_bw()
  
# Ev july
d2dat3 %>%
  ggplot(aes(x = Ev_jul_severity, y = severity)) +
  geom_point() +
  facet_wrap(~sp, scale = "free") +
  theme_bw()
