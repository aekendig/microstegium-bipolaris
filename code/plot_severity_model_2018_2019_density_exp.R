##### info ####

# file: plot_severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 5/4/21
# goal: effects of within and outside of plot severity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(betareg)
library(car)

# import data
d1dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
d2dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")


#### functions ####

# transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}


#### edit data ####

# format edge severity
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(!(month %in% c("sep"))) %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity)

# one plot in one month is missing, use nearby plots
edgeSevD2Dat3 <- tibble(site = "D1", plot = 6, treatment = "fungicide", month = "early_aug") %>%
  mutate(edge_severity = edgeSevD2Dat2 %>% 
           filter(site == "D1" & plot %in% c(4, 5) & treatment == "water" & month == "early_aug") %>%
           summarise(sev = mean(edge_severity)) %>%
           pull(sev)) %>%
  full_join(edgeSevD2Dat2) %>%
  filter(!is.na(edge_severity)) # use to check for missing values

# bg severity
bgSevD1Dat <- d1dat %>%
  filter((plot %in% 2:4 & sp == "Mv") | # select background measurements (includes focal)
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-c(severity, prop_healthy))

bgSevD2Dat <- d2dat %>%
  filter((plot %in% 2:4 & sp == "Mv") |
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-c(severity, prop_healthy))

# focal dataset
focSevD1Dat <- d1dat %>%
  filter(month != "jul") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  mutate(month = dplyr::recode(month, 
                               "late_aug" = "jul", # match prior month
                               "sep" = "late_aug")) %>%
  select(-c(severity, prop_healthy))

focSevD2Dat <- d2dat %>%
  filter(month != "may") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  mutate(month = dplyr::recode(month,  # match prior month
                               "jun" = "may",
                               "jul" = "jun",
                               "early_aug" = "jul",
                               "late_aug" = "early_aug")) %>%
  select(-c(severity, prop_healthy))

# prior severity
priorSevD1Dat <- d1dat %>%
  filter(month != "sep") %>% # remove last month
  mutate(foc_prior_disease = prop_healthy * lesions) %>%
  rename(foc_sp = sp,
         foc_age = age,
         foc_prior_prop_healthy = prop_healthy,
         foc_prior_lesions = lesions) %>%
  select(-severity)

priorSevD2Dat <- d2dat %>%
  filter(month != "late_aug") %>% # remove last month
  mutate(foc_prior_disease = prop_healthy * lesions) %>%
  rename(foc_sp = sp,
         foc_age = age,
         foc_prior_prop_healthy = prop_healthy,
         foc_prior_lesions = lesions) %>%
  select(-severity)

# combine data
sevD1Dat <- focSevD1Dat %>%
  left_join(priorSevD1Dat) %>%
  left_join(bgSevD1Dat) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         bg_sp = case_when(plot == 1 ~ foc_sp, # make the focal species the bg in 1 plots
                           TRUE ~ bg_sp),
         bg_age = case_when(plot == 1 ~ foc_age,
                           TRUE ~ bg_age),
         foc_sp_age = paste(foc_sp, foc_age, sep = "_"),
         bg_sp_age = paste(bg_sp, bg_age, sep = "_"),
         bg_lesions = case_when(foc_sp_age == bg_sp_age & plot == 1 ~ foc_lesions, # use focal lesions as bg in 1 plots
                                TRUE ~ bg_lesions),
         foc_lesions_change = foc_lesions - foc_prior_lesions,
         bg_prior_disease = bg_lesions * foc_prior_prop_healthy) %>%
  group_by(foc_sp, foc_age, bg_sp, bg_age) %>%
  mutate(mean_foc_ls_chg = mean(foc_lesions_change, na.rm = T),
         sd_foc_ls_chg = sd(foc_lesions_change, na.rm = T),
         mean_foc_pri_dis = mean(foc_prior_disease, na.rm = T),
         sd_foc_pri_dis = sd(foc_prior_disease, na.rm = T),
         mean_bg_pri_dis = mean(bg_prior_disease, na.rm = T),
         sd_bg_pri_dis = sd(bg_prior_disease, na.rm = T)) %>%
  ungroup() %>%
  mutate(foc_ls_chg_s = (foc_lesions_change - mean_foc_ls_chg) / sd_foc_ls_chg,
         foc_pr_dis_s = (foc_prior_disease - mean_foc_pri_dis) / sd_foc_pri_dis,
         bg_pr_dis_s = (bg_prior_disease - mean_bg_pri_dis) / sd_bg_pri_dis)

sevD2Dat <- focSevD2Dat %>%
  left_join(priorSevD2Dat) %>%
  left_join(bgSevD2Dat) %>%
  left_join(edgeSevD2Dat3) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         bg_sp = case_when(plot == 1 ~ foc_sp, # make the focal species the bg in 1 plots
                           TRUE ~ bg_sp),
         bg_age = case_when(plot == 1 ~ foc_age,
                            TRUE ~ bg_age),
         foc_sp_age = paste(foc_sp, foc_age, sep = "_"),
         bg_sp_age = paste(bg_sp, bg_age, sep = "_"),
         bg_lesions = case_when(foc_sp_age == bg_sp_age & plot == 1 ~ foc_lesions, # use focal lesions as bg in 1 plots
                                TRUE ~ bg_lesions),
         foc_lesions_change = foc_lesions - foc_prior_lesions,
         bg_prior_disease = bg_lesions * foc_prior_prop_healthy) %>%
  group_by(foc_sp, foc_age, bg_sp, bg_age) %>%
  mutate(mean_foc_ls_chg = mean(foc_lesions_change, na.rm = T),
         sd_foc_ls_chg = sd(foc_lesions_change, na.rm = T),
         mean_foc_pri_dis = mean(foc_prior_disease, na.rm = T),
         sd_foc_pri_dis = sd(foc_prior_disease, na.rm = T),
         mean_bg_pri_dis = mean(bg_prior_disease, na.rm = T),
         sd_bg_pri_dis = sd(bg_prior_disease, na.rm = T)) %>%
  ungroup() %>%
  mutate(foc_ls_chg_s = (foc_lesions_change - mean_foc_ls_chg) / sd_foc_ls_chg,
         foc_pr_dis_s = (foc_prior_disease - mean_foc_pri_dis) / sd_foc_pri_dis,
         bg_pr_dis_s = (bg_prior_disease - mean_bg_pri_dis) / sd_bg_pri_dis)

# modified bg data for figure
bgSevDat <- sevD1Dat %>%
  select(month, site, plot, treatment, bg_sp, bg_age, bg_sp_age, bg_lesions) %>%
  unique() %>%
  mutate(yearf = "2018") %>%
  full_join(sevD2Dat %>%
              select(month, site, plot, treatment, bg_sp, bg_age, bg_sp_age, bg_lesions) %>%
              unique() %>%
              mutate(yearf = "2019")) %>%
  filter(!is.na(bg_lesions)) %>%
  group_by(yearf, bg_sp_age) %>%
  mutate(mean_bg_ls = mean(bg_lesions),
         sd_bg_ls = sd(bg_lesions)) %>%
  ungroup() %>%
  mutate(bg_ls_s = (bg_lesions - mean_bg_ls) / sd_bg_ls,
         density_level = case_when(plot == 1 ~ "none",
                                   plot %in% c(2, 5, 8) ~ "low",
                                   plot %in% c(3, 6, 9) ~ "medium",
                                   plot %in% c(4, 7, 10) ~ "high") %>%
           fct_relevel("none", "low", "medium"))


#### visualizations ####

# prior disease and lesions
ggplot(priorSevD1Dat, aes(foc_prior_lesions, foc_prior_disease)) +
  geom_point(aes(color = month)) +
  geom_smooth() +
  facet_wrap(foc_sp ~ foc_age)
# too correlated, have to pick one of these

# lesion change distributions
ggplot(sevD1Dat, aes(foc_lesions_change)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD1Dat, aes(foc_ls_chg_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# focal disease distributions
ggplot(sevD1Dat, aes(foc_prior_disease)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD1Dat, aes(foc_pr_dis_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# bg disease distributions
ggplot(sevD1Dat, aes(bg_prior_disease)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD1Dat, aes(bg_pr_dis_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# intraspecific
ggplot(sevD1Dat, aes(foc_pr_dis_s, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(foc_pr_dis_s, foc_ls_chg_s, color = treatment, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# interspecific
ggplot(sevD1Dat, aes(bg_pr_dis_s, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp_age ~ bg_sp_age, scales = "free")

ggplot(sevD2Dat, aes(bg_pr_dis_s, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp_age ~ bg_sp_age, scales = "free")

#### start here: combine two below for one panel of figure ####

# bg by density
ggplot(filter(bgSevDat, yearf == "2018"), aes(density_level, bg_ls_s, color = bg_sp_age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(0.1)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(0.1)) +
  fig_theme +
  theme(legend.position = c(0.2, 0.95)) +
  ylab(expression(paste("Scaled lesions (area ", plot^-1, ")", sep = ""))) +
  xlab(expression(paste("Background plant density (", plot^-1, ")", sep = "")))

ggplot(filter(bgSevDat, yearf == "2019"), aes(density_level, bg_ls_s, color = bg_sp_age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(0.1)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(0.1)) +
  fig_theme +
  theme(axis.title.x = element_blank()) +
  ylab(expression(paste("Scaled lesions (g ", plot^-1, ")", sep = "")))


#### divide data ####

mvMvD1Dat <- sevD1Dat %>%
  filter(foc_sp == "Mv" & bg_sp == "Mv" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease))

evSevSD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_seedling" & bg_sp_age == "Ev_seedling" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease))

evAevAD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_adult" & bg_sp_age == "Ev_adult" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease))

mvMvD2Dat <- sevD2Dat %>%
  filter(foc_sp == "Mv" & bg_sp == "Mv" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease))

evSevSD2Dat <- sevD2Dat %>%
  filter(foc_sp_age == "Ev_seedling" & bg_sp_age == "Ev_seedling" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease))

evAevAD2Dat <- sevD2Dat %>%
  filter(foc_sp_age == "Ev_adult" & bg_sp_age == "Ev_adult" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease))

mvEvSD1Dat <- sevD1Dat %>%
  filter(foc_sp == "Mv" & bg_sp_age == "Ev_seedling" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease) & !is.na(bg_prior_disease))

mvEvAD1Dat <- sevD1Dat %>%
  filter(foc_sp == "Mv" & bg_sp_age == "Ev_adult" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease) & !is.na(bg_prior_disease))

evSMvD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_seedling" & bg_sp_age == "Mv_seedling" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease) & !is.na(bg_prior_disease))

evSEvAD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_seedling" & bg_sp_age == "Ev_adult" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease) & !is.na(bg_prior_disease))

evAMvD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_adult" & bg_sp_age == "Mv_seedling" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease) & !is.na(bg_prior_disease))

evAEvSD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_adult" & bg_sp_age == "Ev_seedling" & 
           !is.na(foc_lesions_change) & !is.na(foc_prior_disease) & !is.na(bg_prior_disease))



#### intraspecific models ####

# 2018

# mv by mv
mvMvSevD1Mod1 <- brm(foc_ls_chg_s ~ fungicide * foc_pr_dis_s + (1|plotf),
                     data = mvMvD1Dat, family = gaussian,
                     prior <- c(prior(normal(0, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99)) # divergent transitions
summary(mvMvSevD1Mod1)

# evS by evS
evSevSSevD1Mod1 <- update(mvMvSevD1Mod1, newdata = evSevSD1Dat)
summary(evSevSSevD1Mod1)

# evA by evA
evAevASevD1Mod1 <- update(mvMvSevD1Mod1, newdata = evAevAD1Dat)
summary(evAevASevD1Mod1)

# 2019

# mv by mv
mvMvSevD2Mod1 <- update(mvMvSevD1Mod1, newdata = mvMvD2Dat)
summary(mvMvSevD2Mod1)

# evS by evS
evSevSSevD2Mod1 <- update(mvMvSevD1Mod1, newdata = evSevSD2Dat)
summary(evSevSSevD2Mod1)

# evA by evA
evAevASevD2Mod1 <- update(mvMvSevD1Mod1, newdata = evAevAD2Dat)
summary(evAevASevD2Mod1)


#### interspecific models ####

# 2018

# mv by evS
mvEvSSevD1Mod1 <-brm(foc_ls_chg_s ~ fungicide * (foc_pr_dis_s + bg_pr_dis_s) + (1|plotf),
                     data = mvEvSD1Dat, family = gaussian,
                     prior <- c(prior(normal(0, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                     control = list(adapt_delta = 0.99)) # divergent transitions
summary(mvEvSSevD1Mod1)

# mv by evA
mvEvASevD1Mod1 <- update(mvEvSSevD1Mod1, newdata = mvEvAD1Dat)
summary(mvEvASevD1Mod1)



#### intraspecific density regressions ####

# function
intra_dens_fun <- function(sev, SP, AGE, year, edge = "no", edge_sev = NULL){
  
  # select data
  if(year == 2018){
    dat <- d1dat
  } else{
    dat <- d2dat
  }
  
  # select group
  # tranform column
  # add density info
  dat2 <- dat %>%
    filter(sp == SP & age == AGE & treatment == "water") %>%
    mutate(severity = transform01(.data[[sev]])) %>%
    left_join(plotsD)
  
  # add edge severity if available
  # center edge severity so that intercept is at mean value
  if(edge == "yes"){
    dat3 <- dat2 %>%
      left_join(edgeSevD2Dat3 %>%
                  select(site, plot, treatment, .data[[edge_sev]]) %>%
                  mutate(edge_severity = .data[[edge_sev]] - mean(.data[[edge_sev]])))
  } else{
    dat3 <- dat2
  }
  
  # select plots
  # add focals to density
  if(SP == "Mv"){
    dat4 <- dat3 %>%
      filter(plot %in% 1:4) %>%
      mutate(planted_density = background_density + 3)
  }else if(AGE == "seedling"){
    dat4 <- dat3 %>%
      filter(plot %in% c(1, 5:7)) %>%
      mutate(planted_density = background_density + 3)
  }else{
    dat4 <- dat3 %>%
      filter(plot %in% c(1, 8:10)) %>%
      mutate(planted_density = background_density + 1)
  }
  
  # remove missing data
  dat5 <- dat4 %>%
    filter(!is.na(severity))
  
  # fit regression
  if(edge == "yes"){
    mod <- betareg(severity ~ planted_density + edge_severity, data = dat5)
  } else{
    mod <- betareg(severity ~ planted_density, data = dat5)
  }
  
  # extract values
  out<- tibble(est = summary(mod)[[1]]$mean[2, 1],
               p = summary(mod)[[1]]$mean[2, 4],
               n = nrow(dat5))
  
  return(out)
  
}

# year 1
d1IntraFits <- d1dat %>%
  pivot_longer(cols = c(jul_severity, late_aug_severity, sep_severity),
               names_to = "month_severity",
               values_to = "severity") %>%
  select(sp, age, month_severity) %>%
  unique() %>%
  rename("SP" = sp, "AGE" = age, "sev" = month_severity) %>%
  mutate(year = 2018) %>%
  mutate(sev_mod = pmap(., intra_dens_fun)) %>%
  unnest_wider(sev_mod)
# Ev adults, late August severity

# year 2
d2IntraFits <- d2dat %>%
  pivot_longer(cols = c(jun_severity, jul_severity, early_aug_severity, late_aug_severity),
               names_to = "month_severity",
               values_to = "severity") %>%
  select(sp, age, month_severity) %>%
  unique() %>%
  left_join(tibble(month_severity = c("jun_severity", "jul_severity", "early_aug_severity", "late_aug_severity"),
                   edge_sev = c("may_edge_severity", "jun_edge_severity", "jul_edge_severity", "early_aug_edge_severity"))) %>%
  rename("SP" = sp, "AGE" = age, "sev" = month_severity) %>%
  mutate(year = 2019,
         edge = "yes") %>%
  mutate(sev_mod = pmap(., intra_dens_fun)) %>%
  unnest_wider(sev_mod)
#  none are sig

# format data
intraDat <- d1dat %>%
  filter(plot %in% 1:4 & sp == "Mv") %>%
  full_join(d1dat %>% 
              filter(plot %in% c(1, 5:7) & sp == "Ev" & age == "seedling")) %>%
  full_join(d1dat %>% 
              filter(plot %in% c(1, 8:10) & sp == "Ev" & age == "adult")) %>%
  mutate(yearf = "2018") %>%
  full_join(d2dat %>%
              filter(plot %in% 1:4 & sp == "Mv") %>%
              full_join(d2dat %>% 
                          filter(plot %in% c(1, 5:7) & sp == "Ev" & age == "seedling")) %>%
              full_join(d2dat %>% 
                          filter(plot %in% c(1, 8:10) & sp == "Ev" & age == "adult")) %>%
              mutate(yearf = "2019")) %>%
  mutate(plant_group = paste(sp, age, sep = " ") %>% dplyr::recode("Mv seedling" = "Mv"),
         density_level = case_when(plot == 1 ~ "none",
                                   plot %in% c(2, 5, 8) ~ "low",
                                   plot %in% c(3, 6, 9) ~ "medium",
                                   plot %in% c(4, 7, 10) ~ "high"))

# figure
ggplot(intraDat, aes(density_level, late_aug_severity, color = plant_group, group = interaction(plant_group, yearf))) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2, aes(shape = yearf))
# no clear pattern


#### edge severity regressions ####

# function
edge_sev_fun <- function(sev, SP, AGE, edge_sev, rev = F){
  
  # select group
  # tranform column
  # add density info
  dat2 <- d2dat %>%
    filter(sp == SP & age == AGE & treatment == "water") %>%
    mutate(severity = .data[[sev]],
           severityC = severity - mean(severity, na.rm = T),
           severity01 = transform01(severity)) %>%
    left_join(edgeSevD2Dat3 %>%
                select(site, plot, treatment, .data[[edge_sev]]) %>%
                mutate(edge_severity = .data[[edge_sev]],
                       edge_severityC = edge_severity - mean(edge_severity),
                       edge_severity01 = transform01(edge_severity)))
  
  # select plots
  # add focals to density
  if(SP == "Mv"){
    dat3 <- dat2 %>%
      filter(plot %in% 1:4)
  }else if(AGE == "seedling"){
    dat3 <- dat2 %>%
      filter(plot %in% c(1, 5:7))
  }else{
    dat3 <- dat2 %>%
      filter(plot %in% c(1, 8:10))
  }
  
  # remove missing data
  dat4 <- dat3 %>%
    filter(!is.na(severity))
  
  # fit regression
  if(rev == F){
    mod <- betareg(severity01 ~ edge_severityC, data = dat4)
  } else{
    mod <- betareg(edge_severity01 ~ severityC, data = dat4)
  }
  
  # extract values
  out<- tibble(est = summary(mod)[[1]]$mean[2, 1],
               p = summary(mod)[[1]]$mean[2, 4],
               n = nrow(dat4))
  
  return(out)
  
}

# apply function
d2EdgeFits <- d2dat %>%
  pivot_longer(cols = c(jun_severity, jul_severity, early_aug_severity, late_aug_severity),
               names_to = "month_severity",
               values_to = "severity") %>%
  select(sp, age, month_severity) %>%
  unique() %>%
  left_join(tibble(month_severity = c("jun_severity", "jul_severity", "early_aug_severity", "late_aug_severity"),
                   edge_sev = c("may_edge_severity", "jun_edge_severity", "jul_edge_severity", "early_aug_edge_severity"))) %>%
  rename("SP" = sp, "AGE" = age, "sev" = month_severity) %>%
  mutate(sev_mod = pmap(., edge_sev_fun)) %>%
  unnest_wider(sev_mod)
# Mv late august sig
# largest est for others: Ev seedling july, Ev adult july

d2RevEdgeFits <- d2dat %>%
  pivot_longer(cols = c(may_severity, jun_severity, jul_severity, early_aug_severity),
               names_to = "month_severity",
               values_to = "severity") %>%
  select(sp, age, month_severity) %>%
  unique() %>%
  filter(!(age == "seedling" & month_severity == "may_severity")) %>%
  left_join(tibble(month_severity = c("may_severity", "jun_severity", "jul_severity", "early_aug_severity"),
                   edge_sev = c("jun_edge_severity", "jul_edge_severity", "early_aug_edge_severity", "late_aug_edge_severity"))) %>%
  rename("SP" = sp, "AGE" = age, "sev" = month_severity) %>%
  mutate(rev = T) %>%
  mutate(sev_mod = pmap(., edge_sev_fun)) %>%
  unnest_wider(sev_mod)
# Mv late august, Ev adult jun, jul, early august all sig
# interpretation of above: Ev adults drive disease transmission May-July, Mv drives transmission August

# format data
edgeDat <- d2dat %>%
  filter(plot %in% 1:4 & sp == "Mv") %>%
  mutate(severity = late_aug_severity) %>%
  left_join(edgeSevD2Dat3 %>%
              select(site, plot, treatment, early_aug_edge_severity) %>%
              rename(edge_severity = early_aug_edge_severity)) %>%
  full_join(d2dat %>% 
              filter(plot %in% c(1, 5:7) & sp == "Ev" & age == "seedling") %>%
              mutate(severity = jul_severity) %>%
              left_join(edgeSevD2Dat3 %>%
                          select(site, plot, treatment, jun_edge_severity) %>%
                          rename(edge_severity = jun_edge_severity))) %>%
  full_join(d2dat %>% 
              filter(plot %in% c(1, 8:10) & sp == "Ev" & age == "adult") %>%
              mutate(severity = jul_severity) %>%
              left_join(edgeSevD2Dat3 %>%
                          select(site, plot, treatment, jun_edge_severity) %>%
                          rename(edge_severity = jun_edge_severity))) %>%
  mutate(plant_group = paste(sp, age, sep = " ") %>% dplyr::recode("Mv seedling" = "Mv"))

# figure
ggplot(edgeDat, aes(edge_severity, severity, color = plant_group)) +
  geom_point()


#### focal severity regressions ####

# don't need to format focal data if not doing intragroup plots, then just use the plot-level data

# filter for focal
fSevD2Dat <- sevD2Dat %>%
  filter(focal == 1)

unique(fSevD1Dat$focal)

# check data
fSevD1Dat %>%
  filter(leaves_infec == 0 & lesion_area.pix > 0)
fSevD1Dat %>%
  filter(leaves_infec == 0 & is.na(lesion_area.pix))
fSevD2Dat %>%
  filter(leaves_infec == 0 & lesion_area.pix > 0)
fSevD2Dat %>%
  filter(leaves_infec == 0 & is.na(lesion_area.pix))

# summarize by plot
fSevDat <- fSevD1Dat %>%
  mutate(yearf = "2018") %>%
  full_join(fSevD2Dat %>%
              mutate(yearf = "2019")) %>%
  group_by(yearf, month, site, plot, treatment, sp, age) %>%
  summarize(leaves_tot = mean(leaves_tot, na.rm = T),
            leaves_infec = mean(leaves_infec, na.rm = T),
            leaf_area.pix = mean(leaf_area.pix, na.rm = T),
            lesion_area.pix = mean(lesion_area.pix, na.rm = T)) %>%
  ungroup() %>%
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity01 = transform01(severity))

# check data
fSevDat %>%
  filter(severity > 1)

# format background data
bgDat <- d1dat %>%
  mutate(yearf = "2018") %>%
  full_join(d2dat %>%
              mutate(yearf = "2019")) %>%
  mutate(plant_group = paste(sp, age, sep = "_")) %>%
  pivot_longer(cols = c(may_severity, jun_severity, jul_severity, early_aug_severity, late_aug_severity, sep_severity),
               names_to = "month",
               values_to = "bg_severity") %>%
  mutate(month = str_replace(month, "_severity", ""),
         month = dplyr::recode(month, 
                               "sep" = "rem",
                               "late_aug" = "sep",
                               "early_aug" = "late_aug",
                               "jul" = "early_aug",
                               "jun" = "jul",
                               "may" = "jun"),
         month = case_when(month == "early_aug" & yearf == "2018" ~ "late_aug", # use July to precede late Aug in 2018 (no early Aug)
                           TRUE ~ month)) %>%
  filter(!is.na(bg_severity) & month != "rem" & treatment == "water" &
           (plot == 1 | 
              (plot %in% 2:4 & plant_group == "Mv_seedling") | 
              (plot %in% 5:7 & plant_group == "Ev_seedling") | 
              (plot %in% 8:10 & plant_group == "Ev_adult"))) %>%
  left_join(plotsD %>%
              select(plot, background_density) %>%
              unique()) %>%
  mutate(planted_density = case_when(plant_group == "Mv_seedling" ~ background_density + 3,
                                     plant_group == "Ev_seedling" ~ background_density + 3,
                                     plant_group == "Ev_adult" ~ background_density + 1),
         FOI = planted_density * bg_severity) %>%
  select(-c(sp, age))
  

# function
focal_sev_fun <- function(YEAR, MONTH, SP, AGE, BG){
  
  # select focal data
  # select group
  # tranform column
  # add density info
  dat <- bgDat %>%
    filter(plant_group == BG & yearf == YEAR & month == MONTH) %>%
    left_join(fSevDat %>%
                filter(sp == SP & age == AGE)) %>%
    mutate(planted_densityC = planted_density - mean(planted_density),
           bg_severityC = bg_severity - mean(bg_severity),
           FOIC = FOI - mean(FOI))
  
  # fit regression
  mod1 <- betareg(severity01 ~ planted_densityC + bg_severityC, data = dat)
  mod2 <- betareg(severity01 ~ FOIC, data = dat)

  # extract values
  out<- tibble(dens_est = summary(mod1)[[1]]$mean[2, 1],
               dens_p = summary(mod1)[[1]]$mean[2, 4],
               sev_est = summary(mod1)[[1]]$mean[3, 1],
               sev_p = summary(mod1)[[1]]$mean[3, 4],
               FOI_est = summary(mod2)[[1]]$mean[2, 1],
               FOI_p = summary(mod2)[[1]]$mean[2, 4],
               n = nrow(dat))
  
  return(out)
  
}

# apply function
focSevFits <- bgDat %>%
  select(plant_group, yearf, month) %>%
  unique() %>%
  merge(fSevDat %>%
          select(sp, age) %>%
          unique(), all = T) %>%
  as_tibble() %>%
  filter(!(yearf == "2019" & month == "sep")) %>%
  rename("YEAR" = yearf, "MONTH" = month, "SP" = sp, "AGE" = age, "BG" = plant_group) %>%
  mutate(sev_mod = pmap(., focal_sev_fun)) %>%
  unnest_wider(sev_mod)

focSevFits %>%
  filter(dens_p < 0.05 | sev_p < 0.05 | FOI_p < 0.05)
# positive effects
# Ev adults increase Ev adult severity through density (jul-aug 2018)
# Ev adults increase Ev seedling severity through severity (early-late aug 2019)
# Ev adult increases Mv severity through severity (may-jun, jun-jul, jul-early aug 2019)
# Ev seedling increases Mv severity through severity (jul-aug 2018), density (early-late aug 2019), and FOI (jul-early aug 2019)
# Mv increases Mv severity through severity (aug-sep 2018, jun-jul 2019, and early-late aug 2019)

# negative effects
# Ev adults decrease Ev adult severity through severity (jul-aug 2018)
# Ev seedling decreases Ev seedling severity through density (aug-sep 2018)
# Ev seedling decreases Ev adult severity through severity (aug-sep 2018)
# Mv decreases Ev adult severity through density (jul-early aug 2019)
# Mv decreases Ev seedling through density and severity (sep 2018)


