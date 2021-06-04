##### info ####

# file: plot_severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/3/21
# goal: effects of within and outside of plot severity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)
library(gridExtra)

# import data
d1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
d2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R


#### edit data ####

# # format edge severity
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
bgSevD1Dat <- d1Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") | # select background measurements (includes focal)
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-severity)

bgSevD2Dat <- d2Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") |
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-severity)

# focal dataset (includes background when background is same species)
focNextSevD1Dat <- d1Dat %>%
  filter(month != "jul") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_next_lesions = lesions) %>%
  mutate(month = dplyr::recode(month, 
                               "late_aug" = "jul", # match prior month
                               "sep" = "late_aug")) %>%
  select(-severity)

focNextSevD2Dat <- d2Dat %>%
  filter(month != "may") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_next_lesions = lesions) %>%
  mutate(month = dplyr::recode(month,  # match prior month
                               "jun" = "may",
                               "jul" = "jun",
                               "early_aug" = "jul",
                               "late_aug" = "early_aug")) %>%
  select(-severity)

# prior severity
focSevD1Dat <- d1Dat %>%
  filter(month != "sep") %>% # remove last month
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  select(-severity)

focSevD2Dat <- d2Dat %>%
  filter(month != "late_aug") %>% # remove last month
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  select(-severity)

# combine data
sevD1Dat <- focNextSevD1Dat %>%
  left_join(focSevD1Dat) %>%
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
         foc_lesions_change = log(foc_next_lesions / foc_lesions),
         foc_lesions_change = case_when(foc_lesions == 0 & foc_next_lesions == 0 ~ log(1), # 4 plots had no lesions, no change
                                        TRUE ~ foc_lesions_change),
         foc_lesions_change = case_when(month == "jul" ~ foc_lesions_change / 2, # account for longer lapse between samples 
                                        TRUE ~ foc_lesions_change)) %>%
  group_by(foc_sp, foc_age, bg_sp, bg_age) %>%
  mutate(mean_bg_les = mean(bg_lesions, na.rm = T),
         sd_bg_les = sd(bg_lesions, na.rm = T)) %>%
  ungroup() %>%
  mutate(bg_les_s = (bg_lesions - mean_bg_les) / sd_bg_les,
         edge_severity = 0) %>%
  filter(!is.na(foc_next_lesions) & !is.na(foc_lesions) & !is.na(bg_lesions))

# check for missing data
filter(sevD1Dat, is.na(foc_next_lesions))
filter(sevD1Dat, is.na(bg_lesions)) # 18 missing
filter(sevD1Dat, is.na(foc_lesions)) # 7 missing
filter(sevD1Dat, is.na(foc_lesions_change) & !is.na(foc_next_lesions) & !is.na(foc_lesions))  
filter(sevD1Dat, foc_lesions_change %in% c(Inf, -Inf))


sevD2Dat <- focNextSevD2Dat %>%
  left_join(focSevD2Dat) %>%
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
         foc_lesions_change = log(foc_next_lesions / foc_lesions),
         foc_lesions_change = case_when(foc_lesions == 0 & foc_next_lesions == 0 ~ log(1), # 4 plots had no lesions, no change
                                        TRUE ~ foc_lesions_change),
         foc_lesions_change = case_when(foc_lesions == 0 & foc_next_lesions > 0 ~ log(foc_next_lesions / 1e-6),
                                        TRUE ~ foc_lesions_change),  # used close to smallest lesion amount (6.4e-6) as denominator
         foc_lesions_change = case_when(foc_lesions > 0 & foc_next_lesions == 0 ~ log(1e-6 / foc_lesions),
                                        TRUE ~ foc_lesions_change)) %>%  # used close to smallest lesion amount (6.4e-6) as numerator
  group_by(foc_sp, foc_age, bg_sp, bg_age) %>%
  mutate(mean_bg_les = mean(bg_lesions, na.rm = T),
         sd_bg_les = sd(bg_lesions, na.rm = T)) %>%
  ungroup() %>%
  mutate(bg_les_s = (bg_lesions - mean_bg_les) / sd_bg_les) %>%
  filter(!is.na(foc_next_lesions) & !is.na(foc_lesions) & !is.na(bg_lesions))

# check for missing data
filter(sevD2Dat, is.na(foc_next_lesions))
filter(sevD2Dat, is.na(bg_lesions)) # 189 missing
filter(sevD2Dat, is.na(foc_lesions)) # 185 missing
filter(sevD2Dat, is.na(foc_lesions_change) & !is.na(foc_next_lesions) & !is.na(foc_lesions)) 
filter(sevD2Dat, foc_lesions_change %in% c(Inf, -Inf)) # 21 cases


#### visualizations ####

# lesions change distributions
ggplot(sevD1Dat, aes(foc_lesions_change)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(foc_lesions_change)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# bg lesions distributions
ggplot(sevD1Dat, aes(bg_lesions)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD1Dat, aes(bg_les_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(bg_lesions)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(bg_les_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# pairwise combinations
ggplot(sevD1Dat, aes(bg_les_s, foc_lesions_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(foc_sp_age ~ bg_sp_age, scales = "free")

ggplot(sevD2Dat, aes(bg_les_s, foc_lesions_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(foc_sp_age ~ bg_sp_age, scales = "free")
# looked at high Ev adult plot (bg_les_s > 7) - no Ev scans for late August

# edge
ggplot(sevD2Dat, aes(edge_severity, foc_lesions_change, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(~ foc_sp_age, scales = "free")


#### divide data ####

mvD1Dat <- sevD1Dat %>%
  filter(foc_sp == "Mv") %>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Mv_seedling"))

evSD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_seedling") %>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Ev_seedling"))

evAD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_adult")

mvD2Dat <- sevD2Dat %>%
  filter(foc_sp == "Mv") %>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Mv_seedling"))

evSD2Dat <- sevD2Dat %>%
  filter(foc_sp_age == "Ev_seedling")%>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Ev_seedling"))

evAD2Dat <- sevD2Dat %>%
  filter(foc_sp_age == "Ev_adult")


#### 2018 models ####

# Mv
mvSevD1Mod1 <- brm(foc_lesions_change ~ fungicide * bg_les_s * bg_sp_age + (1|plotf),
                   data = mvD1Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99))
summary(mvSevD1Mod1)

# Ev seedling
evSSevD1Mod1 <- update(mvSevD1Mod1, newdata = evSD1Dat)
summary(evSSevD1Mod1)

# Ev adult
evASevD1Mod1 <- update(mvSevD1Mod1, newdata = evAD1Dat)
summary(evASevD1Mod1)

# save
save(mvSevD1Mod1, file = "output/mv_focal_plot_transmission_2018_density_exp.rda")
save(evSSevD1Mod1, file = "output/evS_focal_plot_transmission_2018_density_exp.rda")
save(evASevD1Mod1, file = "output/evA_focal_plot_transmission_2018_density_exp.rda")


#### 2019 models ####

# Mv
mvSevD2Mod1 <- brm(foc_lesions_change ~ fungicide * (bg_les_s * bg_sp_age + edge_severity) + (1|plotf),
                   data = mvD2Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
summary(mvSevD2Mod1)

# Ev seedling
evSSevD2Mod1 <- update(mvSevD2Mod1, newdata = evSD2Dat)
summary(evSSevD2Mod1)

# Ev adult
evASevD2Mod1 <- update(mvSevD2Mod1, newdata = evAD2Dat)
summary(evASevD2Mod1)

# save
save(mvSevD2Mod1, file = "output/mv_focal_plot_transmission_2019_density_exp.rda")
save(evSSevD2Mod1, file = "output/evS_focal_plot_transmission_2019_density_exp.rda")
save(evASevD2Mod1, file = "output/evA_focal_plot_transmission_2019_density_exp.rda")


#### pairwise figure ####

# prediction function
pair_pred_fun <- function(year, foc, bg, trt){
  
  # data
  if(year == "2018"){
    dat <- sevD1Dat
  }else{
    dat <- sevD2Dat
  }
  
  # model
  if(year == "2018" & foc == "Mv_seedling"){
    mod <- mvSevD1Mod1
  }else if(year == "2018" & foc == "Ev_seedling"){
    mod <- evSSevD1Mod1
  }else if(year == "2018" & foc == "Ev_adult"){
    mod <- evASevD1Mod1
  }else if(year == "2019" & foc == "Mv_seedling"){
    mod <- mvSevD2Mod1
  }else if(year == "2019" & foc == "Ev_seedling"){
    mod <- evSSevD2Mod1
  }else if(year == "2019" & foc == "Ev_adult"){
    mod <- evASevD2Mod1
  }
  
  # subset for species pair
  dat2 <- dat %>%
    filter(foc_sp_age == foc & bg_sp_age == bg & treatment == trt)
  
  # sequence of lesion values
  # other variables
  # predicted values
  simDat <- tibble(bg_les_s = seq(min(dat2$bg_les_s), max(dat2$bg_les_s), length.out = 100)) %>%
    mutate(edge_severity = mean(dat$edge_severity),
           plotf = "A",
           fungicide = case_when(trt == "water" ~ 0,
                                 trt == "fungicide" ~ 1),
           bg_sp_age = bg) %>%
    mutate(foc_lesions_change = fitted(mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(mod, newdata = ., allow_new_levels = T)[, "Q97.5"]) %>%
    select(-bg_sp_age)
  
  # output
  return(simDat)
  
}

# predicted data
predDat <- tibble(foc_sp_age = c("Ev_adult", "Ev_seedling", "Mv_seedling")) %>%
  expand_grid(tibble(bg_sp_age = c("Ev_adult", "Ev_seedling", "Mv_seedling"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  expand_grid(tibble(year = c("2018", "2019"))) %>%
  mutate(pred = pmap(list(year, foc_sp_age, bg_sp_age, treatment), pair_pred_fun)) %>%
  unnest(pred)

# coefficients function
pair_coef_fun <- function(year, foc){
  
  # model
  if(year == "2018" & foc == "Mv_seedling"){
    mod <- mvSevD1Mod1
  }else if(year == "2018" & foc == "Ev_seedling"){
    mod <- evSSevD1Mod1
  }else if(year == "2018" & foc == "Ev_adult"){
    mod <- evASevD1Mod1
  }else if(year == "2019" & foc == "Mv_seedling"){
    mod <- mvSevD2Mod1
  }else if(year == "2019" & foc == "Ev_seedling"){
    mod <- evSSevD2Mod1
  }else if(year == "2019" & foc == "Ev_adult"){
    mod <- evASevD2Mod1
  }
  
  # rename posterior samples
  if(foc == "Mv_seedling"){
    out <- posterior_samples(mod) %>%
      rename("water_Mv_seedling" = "b_bg_les_s",
             "b_Ev_adult" = "b_bg_les_s:bg_sp_ageEv_adult",
             "b_Ev_seedling" = "b_bg_les_s:bg_sp_ageEv_seedling",
             "b_Mv_seedling_fung" = "b_fungicide:bg_les_s",
             "b_Ev_adult_fung" = "b_fungicide:bg_les_s:bg_sp_ageEv_adult",
             "b_Ev_seedling_fung" = "b_fungicide:bg_les_s:bg_sp_ageEv_seedling",) %>%
      mutate(water_Ev_adult = water_Mv_seedling + b_Ev_adult,
             water_Ev_seedling = water_Mv_seedling + b_Ev_seedling,
             fungicide_Mv_seedling = water_Mv_seedling + b_Mv_seedling_fung,
             fungicide_Ev_adult = water_Ev_adult + b_Mv_seedling_fung + b_Ev_adult_fung,
             fungicide_Ev_seedling = water_Ev_seedling + b_Mv_seedling_fung + b_Ev_seedling_fung)
  } else if(foc == "Ev_seedling"){
    out <- posterior_samples(mod) %>%
      rename("water_Ev_seedling" = "b_bg_les_s",
             "b_Ev_adult" = "b_bg_les_s:bg_sp_ageEv_adult",
             "b_Mv_seedling" = "b_bg_les_s:bg_sp_ageMv_seedling",
             "b_Ev_seedling_fung" = "b_fungicide:bg_les_s",
             "b_Ev_adult_fung" = "b_fungicide:bg_les_s:bg_sp_ageEv_adult",
             "b_Mv_seedling_fung" = "b_fungicide:bg_les_s:bg_sp_ageMv_seedling",) %>%
      mutate(water_Ev_adult = water_Ev_seedling + b_Ev_adult,
             water_Mv_seedling = water_Ev_seedling + b_Mv_seedling,
             fungicide_Ev_seedling = water_Ev_seedling + b_Ev_seedling_fung,
             fungicide_Ev_adult = water_Ev_adult + b_Ev_seedling_fung + b_Ev_adult_fung,
             fungicide_Mv_seedling = water_Mv_seedling + b_Ev_seedling_fung + b_Mv_seedling_fung)
  } else{
    out <- posterior_samples(mod) %>%
      rename("water_Ev_adult" = "b_bg_les_s",
             "b_Ev_seedling" = "b_bg_les_s:bg_sp_ageEv_seedling",
             "b_Mv_seedling" = "b_bg_les_s:bg_sp_ageMv_seedling",
             "b_Ev_adult_fung" = "b_fungicide:bg_les_s",
             "b_Ev_seedling_fung" = "b_fungicide:bg_les_s:bg_sp_ageEv_seedling",
             "b_Mv_seedling_fung" = "b_fungicide:bg_les_s:bg_sp_ageMv_seedling",) %>%
      mutate(water_Ev_seedling = water_Ev_adult + b_Ev_seedling,
             water_Mv_seedling = water_Ev_adult + b_Mv_seedling,
             fungicide_Ev_adult = water_Ev_adult + b_Ev_adult_fung,
             fungicide_Ev_seedling = water_Ev_seedling + b_Ev_adult_fung + b_Ev_seedling_fung,
             fungicide_Mv_seedling = water_Mv_seedling + b_Ev_adult_fung + b_Mv_seedling_fung)
  }
  
  # re-organize and summarize
  out2 <- out %>%
    select(water_Mv_seedling, water_Ev_adult, water_Ev_seedling,
           fungicide_Mv_seedling, fungicide_Ev_adult, fungicide_Ev_seedling) %>%
    pivot_longer(cols = everything(),
                 names_to = c("treatment", ".value"),
                 names_pattern = "(.*)_(.*_.*)") %>%
    pivot_longer(cols = -treatment,
                 names_to = "bg_sp_age",
                 values_to = "beta") %>%
    group_by(treatment, bg_sp_age) %>%
    median_hdi(beta)
    
    return(out2)
  
}

# beta values
betaDat <- tibble(foc_sp_age = c("Ev_adult", "Ev_seedling", "Mv_seedling")) %>%
  expand_grid(tibble(year = c("2018", "2019"))) %>%
  mutate(coef = pmap(list(year, foc_sp_age), pair_coef_fun)) %>%
  unnest(coef) %>%
  mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"),
         foc_sp_age = str_replace(foc_sp_age, "_", " "),
         foc_sp_age = fct_recode(foc_sp_age, Mv = "Mv seedling"),
         bg_sp_age = str_replace(bg_sp_age, "_", " "),
         bg_sp_age = fct_recode(bg_sp_age, Mv = "Mv seedling"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# add to predicted dataset
predDat2 <- predDat %>%
  mutate(foc_sp_age = str_replace(foc_sp_age, "_", " "),
         foc_sp_age = fct_recode(foc_sp_age, Mv = "Mv seedling"),
         bg_sp_age = str_replace(bg_sp_age, "_", " "),
         bg_sp_age = fct_recode(bg_sp_age, Mv = "Mv seedling"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)")) %>%
  left_join(betaDat %>%
              select(year, foc_sp_age, bg_sp_age, treatment, sig))

# raw data
sevDat <- sevD1Dat %>%
  mutate(year = "2018") %>%
  full_join(sevD2Dat %>%
              mutate(year = "2019")) %>%
  mutate(foc_sp_age = str_replace(foc_sp_age, "_", " "),
         foc_sp_age = fct_recode(foc_sp_age, Mv = "Mv seedling"),
         bg_sp_age = str_replace(bg_sp_age, "_", " "),
         bg_sp_age = fct_recode(bg_sp_age, Mv = "Mv seedling"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# sig beta values
betaDat2 <- betaDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, bg_sp_age) %>%
              summarise(bg_les_s = max(bg_les_s))) %>%
  left_join(predDat2 %>%
              group_by(year, foc_sp_age) %>%
              summarise(upper = max(upper)) %>%
              ungroup() %>%
              full_join(sevDat %>%
                          group_by(year, foc_sp_age) %>%
                          summarise(foc = max(foc_lesions_change),
                                    min = min(foc_lesions_change))) %>%
              rowwise() %>%
              mutate(foc_lesions_change = max(c(foc, upper))) %>%
              ungroup() %>%
              select(-c(foc, upper))) %>%
  rename(tran = beta) %>%
  mutate(tran = round(tran, 2),
         .lower = round(.lower, 2),
         .upper = round(.upper, 2),
         lower = case_when(tran > 0 ~ .upper, # switch because they were negative
                           TRUE ~ .lower),
         upper = case_when(tran > 0 ~ .lower,
                           TRUE ~ .upper),
         foc_lesions_change = case_when(bg_sp_age == "Ev adult" & foc_sp_age == "Mv" & year == "2018" ~ min + 0.1,
                                        bg_sp_age == "Mv" & foc_sp_age == "Ev adult" & year == "2018" ~ min + 1.7,
                                        TRUE ~ foc_lesions_change)) %>%
  select(-c(.lower, .upper))
              

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.1, 0, 0.2, 2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

yearText <- tibble(year = c("2018", "2019"),
                   bg_les_s = c(-1, -1),
                   foc_lesions_change = c(3.8, 13.5),
                   bg_sp_age = "Ev adult",
                   foc_sp_age = "Ev adult",
                   treatment = "fungicide")

textSize = 3

# water figure
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = bg_les_s, y = foc_lesions_change)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(sevDat, year == "2018"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(betaDat2, year == "2018"), 
            aes(label = paste("beta", " == ", tran, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(foc_sp_age),
             cols = vars(bg_sp_age),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab("Initial infected tissue") +
  ylab("Change in infected tissue") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# fungicide figure
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = bg_les_s, y = foc_lesions_change)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(sevDat, year == "2019"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(betaDat2, year == "2019"), 
            aes(label = paste("beta", " == ", tran, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, size = textSize) +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(foc_sp_age),
             cols = vars(bg_sp_age),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab("Initial infected tissue") +
  ylab("Change in infected tissue") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg <- get_legend(pairD1Fig)

# combine plots
combFig <- plot_grid(pairD1Fig + theme(legend.position = "none"), pairD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.9))

# combine
pdf("output/plot_severity_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()

# save beta values
write_csv(betaDat, "output/plot_severity_pairwise_coefficients_2018_2019_density_exp.csv")


#### extract effects ####

# bg_eff_fun <- function(mod, foc, yr){
#   
#   if(foc == "Mv"){
#     out <- posterior_samples(mod) %>%
#       rename("Mv" = "b_bg_pr_dis_s", 
#              "b_Ev_adult" = "b_bg_pr_dis_s:bg_sp_ageEv_adult",
#              "b_Ev_seedling" = "b_bg_pr_dis_s:bg_sp_ageEv_seedling") %>%
#       mutate(Ev_adult = Mv + b_Ev_adult,
#              Ev_seedling = Mv + b_Ev_seedling)
#   } else if(foc == "Ev seedling"){
#     out <- posterior_samples(mod) %>%
#       rename("Ev_seedling" = "b_bg_pr_dis_s", 
#              "b_Ev_adult" = "b_bg_pr_dis_s:bg_sp_ageEv_adult",
#              "b_Mv" = "b_bg_pr_dis_s:bg_sp_ageMv_seedling") %>%
#       mutate(Ev_adult = Ev_seedling + b_Ev_adult,
#              Mv = Ev_seedling + b_Mv)
#   } else{
#     out <- posterior_samples(mod) %>%
#       rename("Ev_adult" = "b_bg_pr_dis_s", 
#              "b_Ev_seedling" = "b_bg_pr_dis_s:bg_sp_ageEv_seedling",
#              "b_Mv" = "b_bg_pr_dis_s:bg_sp_ageMv_seedling") %>%
#       mutate(Ev_seedling = Ev_adult + b_Ev_seedling,
#              Mv = Ev_adult + b_Mv)
#   }
# 
#   
#   out2 <- out %>%
#     select(Mv, Ev_adult, Ev_seedling) %>%
#     pivot_longer(cols = everything(),
#                  names_to = "bg_sp_age",
#                  values_to = "bg_eff") %>%
#     group_by(bg_sp_age) %>%
#     median_hdi(bg_eff) %>%
#     mutate(focal = foc,
#            background = str_replace(bg_sp_age, "_", " "),
#            yearf = yr)
#   
#   return(out2)
#   
# }
# 
# bgEff <- bg_eff_fun(mvSevD1Mod1, "Mv", "2018") %>%
#   full_join(bg_eff_fun(evSSevD1Mod1, "Ev seedling", "2018")) %>%
#   full_join(bg_eff_fun(evASevD1Mod1, "Ev adult", "2018")) %>%
#   full_join(bg_eff_fun(mvSevD2Mod1, "Mv", "2019")) %>%
#   full_join(bg_eff_fun(evSSevD2Mod1, "Ev seedling", "2019")) %>%
#   full_join(bg_eff_fun(evASevD2Mod1, "Ev adult", "2019")) %>%
#   mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
#                          TRUE ~ "includes 0"))


#### effects figure ####

# fig_theme <- theme_bw() +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.spacing.x = unit(0,"line"),
#         axis.text.y = element_text(size = 8, color = "black"),
#         axis.text.x = element_text(size = 8, color = "black"),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10),
#         axis.line = element_line(color = "black"),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.background = element_blank(),
#         legend.position = "none",
#         legend.key.size = unit(1, "mm"),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8),
#         strip.placement = "outside")
# 
# col_pal = c("#472F7DFF", "#39568CFF", "#35B779FF")
# shape_pal = c(23, 22, 24)

# bg_dens_figD1 <- ggplot(filter(bgSevDat, yearf == "2018"), aes(density_level, bg_ls_s, color = bg_sp_age)) +
#   stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(0.7)) +
#   stat_summary(geom = "point", size = 2, fun = "mean", fill = "white", position = position_dodge(0.7), aes(shape = bg_sp_age)) +
#   geom_text(x = 0.5, y = 1.03, label = "2018", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
#   ylab("Background lesions") +
#   xlab("Background plant density") +
#   scale_color_manual(values = col_pal, name = "Plant group") +
#   scale_shape_manual(values = shape_pal, name = "Plant group") +
#   fig_theme +
#   theme(legend.position = c(0.23, 0.75))
# 
# bg_dens_figD2 <- ggplot(filter(bgSevDat, yearf == "2019"), aes(density_level, bg_ls_s, color = bg_sp_age)) +
#   stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(0.7)) +
#   stat_summary(geom = "point", size = 2, fun = "mean", fill = "white", position = position_dodge(0.7), aes(shape = bg_sp_age)) +
#   geom_text(x = 0.5, y = 1.3, label = "2019", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
#   ylab("Background lesions") +
#   xlab("Background plant density") +
#   scale_color_manual(values = col_pal, name = "Plant group") +
#   scale_shape_manual(values = shape_pal, name = "Plant group") +
#   fig_theme
# 
# les_figD1 <- ggplot(filter(bgEff, yearf == "2018"), aes(background, bg_eff, color = focal)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, position = position_dodge(0.5)) +
#   geom_point(aes(shape = focal, fill = sig), size = 2, position = position_dodge(0.5)) +
#   geom_text(x = 0.48, y = 1.06, label = "2018", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
#   xlab("Background plant group") +
#   ylab("Change in focal lesions") +
#   scale_color_manual(values = col_pal, name = "Focal plant group") +
#   scale_shape_manual(values = shape_pal, name = "Focal plant group") +
#   scale_fill_manual(values = c("white", "black"), guide = F) +
#   fig_theme
# 
# les_figD2 <- ggplot(filter(bgEff, yearf == "2019"), aes(background, bg_eff, color = focal)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, position = position_dodge(0.5)) +
#   geom_point(aes(shape = focal, fill = sig), size = 2, position = position_dodge(0.5)) +
#   geom_text(x = 0.48, y = 1.18, label = "2019", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
#   xlab("Background plant group") +
#   ylab("Change in focal lesions") +
#   scale_color_manual(values = col_pal, name = "Focal plant group") +
#   scale_shape_manual(values = shape_pal, name = "Focal plant group") +
#   scale_fill_manual(values = c("white", "black"), guide = F) +
#   fig_theme +
#   theme(legend.position = c(0.6, 0.87))
# 
# 
# # combine figure
# pdf("output/plot_severity_effect_figure_2018_2019_density_exp.pdf", width = 4.7, height = 4.7)
# plot_grid(bg_dens_figD1, bg_dens_figD2, les_figD1, les_figD2,
#           nrow = 2,
#           labels = LETTERS[1:4])
# dev.off()
