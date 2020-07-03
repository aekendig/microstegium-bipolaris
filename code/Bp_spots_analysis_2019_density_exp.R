##### info ####

# file: Bp_spots_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/2/20
# goal: analyze occurrence of Bp spots


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(DHARMa)
library(MASS)
library(tidyverse)
library(glmmTMB)

# import all raw data files
dt_may <- read_csv("data/ev_disease_may_2019_density_exp.csv")
dt_jun <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
dt_jul <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
dt_early_aug <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
dt_late_aug <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")
dt_sep <- read_csv("data/focal_disease_sep_2019_density_exp.csv")
edge_dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
temp_hum <- read_csv("intermediate-data/temp_humidity_daily_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


#### edit data ####

# follow the same protocols used in focal_severity_analysis_2019_density_exp for temp_hum, bio, and edge dat

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

# dataset function
dat_fun <- function(dat_in, mth){
  
  dat_out <- dat_in %>%
    filter(sp == "Ev" & !is.na(Bp_spots)) %>%
    mutate(month = mth,
           fungicide = case_when(treatment == "water" ~ 0,
                                 TRUE ~ 1),
           exp_plot = paste0(plot, toupper(substr(treatment, 1, 1))),
           site_exp_plot = paste(site, exp_plot, sep = " "),
           Month = recode(month, early_aug = "Early August", jul = "July", 
                          jun = "June", late_aug = "Late August", may = "May", sep = "September") %>%
             fct_relevel("May", "June", "July", "Early August", "Late August", "September"),
           Treatment = recode(treatment, water = "control (water)")) %>%
    left_join(temp_hum2) %>%
    left_join(edge_dat2) %>%
    left_join(bio) %>%
    left_join(dens_dat)
  
  return(dat_out)
}

# datasets
may_dat <- dat_fun(dt_may, "may")
may_ev_dat <- dat_fun(filter(dt_may, ID %in% c("1", "2", "3", "A")), "may")
jun_ev_dat <- dat_fun(dt_jun, "jun")
jul_ev_dat <- dat_fun(dt_jul, "jul")
eau_ev_dat <- dat_fun(dt_early_aug, "early_aug")
lau_ev_dat <- dat_fun(dt_late_aug, "late_aug")
sep_ev_dat <- dat_fun(dt_sep, "sep")

# water analysis
ev_water_dat <- full_join(eau_ev_dat, lau_ev_dat) %>%
  filter(treatment == "water")

# June-late August analysis
jun_lau_ev_dat <- full_join(jun_ev_dat, jul_ev_dat) %>%
  full_join(eau_ev_dat) %>%
  full_join(lau_ev_dat)


#### visualize ####

ggplot(ev_water_dat, aes(x = total_biomass.g, y = Bp_spots)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~month)

ggplot(ev_dat, aes(x = total_biomass.g, y = Bp_spots, color = treatment)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~month)


#### model selection for water treatment ####

# include humidity data (only available for control plots)

# divide data
sep_ev_water_dat <- filter(sep_ev_dat, treatment == "water")

## early-late August ##
ev_totb_water_mod <- glmmTMB(Bp_spots ~ Month * (total_biomass.g + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "binomial")
ev_sepb_water_mod <- glmmTMB(Bp_spots ~ Month * (Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "binomial")
ev_totd_water_mod <- glmmTMB(Bp_spots ~ Month * (total_density + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "binomial")
ev_sepd_water_mod <- glmmTMB(Bp_spots ~ Month * (Ev_density + Mv_density + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "binomial")
ev_totp_water_mod <- glmmTMB(Bp_spots ~ Month * (bg_present + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "binomial")
ev_sepp_water_mod <- glmmTMB(Bp_spots ~ Month * (Ev_present + Mv_present + edge_severity + hum_avg) + (1|site/plot), data = ev_water_dat, family = "binomial")

# compare different plant models
AIC(ev_totb_water_mod, ev_sepb_water_mod, ev_totd_water_mod, ev_sepd_water_mod, ev_totp_water_mod, ev_sepp_water_mod)
# total presence is the lowest
summary(ev_totp_water_mod)
# month:bg (*)
stepAIC(ev_totp_water_mod)
# keep full model
plot(simulateResiduals(ev_totp_water_mod))

## September ##
sep_ev_totb_water_mod <- glmmTMB(Bp_spots ~ total_biomass.g + edge_severity + hum_avg + (1|site/plot), data = sep_ev_water_dat, family = "binomial")
# can't converge
summary(sep_ev_totb_water_mod)
sep_ev_totb_water_mod <- glmmTMB(Bp_spots ~ total_biomass.g + edge_severity + hum_avg + (1|site_exp_plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_sepb_water_mod <- glmmTMB(Bp_spots ~ Ev_biomass.g + Mv_biomass.g + edge_severity + hum_avg + (1|site_exp_plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_totd_water_mod <- glmmTMB(Bp_spots ~ total_density + edge_severity + hum_avg + (1|site_exp_plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_sepd_water_mod <- glmmTMB(Bp_spots ~ Ev_density + Mv_density + edge_severity + hum_avg + (1|site_exp_plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_totp_water_mod <- glmmTMB(Bp_spots ~ bg_present + edge_severity + hum_avg + (1|site_exp_plot), data = sep_ev_water_dat, family = "binomial")
sep_ev_sepp_water_mod <- glmmTMB(Bp_spots ~ Ev_present + Mv_present + edge_severity + hum_avg + (1|site_exp_plot), data = sep_ev_water_dat, family = "binomial")

# compare different plant models
AIC(sep_ev_totb_water_mod, sep_ev_sepb_water_mod, sep_ev_totd_water_mod, sep_ev_sepd_water_mod, sep_ev_totp_water_mod, sep_ev_sepp_water_mod)
# tot biomass is the lowest
summary(sep_ev_totb_water_mod)
# none significant
stepAIC(sep_ev_totb_water_mod)
# keep full model
plot(simulateResiduals(sep_ev_totb_water_mod))


#### model selection for both treatments ####

## May all Ev ##
may_ev_all_totb_mod <- glmmTMB(Bp_spots ~ fungicide * total_biomass.g + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g) + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_totd_mod <- glmmTMB(Bp_spots ~ fungicide * total_density + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density) + (1|site/exp_plot), data = may_dat, family = "binomial")
# convergence issue
summary(may_ev_all_sepd_mod) # remove site
may_ev_all_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density) + (1|site_exp_plot), data = may_dat, family = "binomial")
# convergence issue
filter(may_dat, Mv_present == 1 & treatment == "fungicide")
# only 0's
may_ev_all_totp_mod <- glmmTMB(Bp_spots ~ fungicide * bg_present + (1|site/exp_plot), data = may_dat, family = "binomial")
may_ev_all_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present) + (1|site/exp_plot), data = may_dat, family = "binomial")
# convergence issue

# compare different plant models
AIC(may_ev_all_totb_mod, may_ev_all_sepb_mod, may_ev_all_totd_mod, may_ev_all_totp_mod)
# separate biomass is the lowest
summary(may_ev_all_sepb_mod)
stepAIC(may_ev_all_sepb_mod)
# remove Mv_biomass and interactions
may_ev_all_sepb_mod2 <- update(may_ev_all_sepb_mod, .~. - Mv_biomass.g - fungicide:Mv_biomass.g)
summary(may_ev_all_sepb_mod2)
# none sig
plot(simulateResiduals(may_ev_all_sepb_mod2))


## May Ev ##
may_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * total_biomass.g + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
summary(may_ev_totb_mod)
may_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * total_biomass.g + (1|site_exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
may_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g) + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
may_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * total_density + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
may_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density) + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
may_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * bg_present + (1|site/exp_plot), data = may_ev_dat, family = "binomial")
# convergence issue
may_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present) + (1|site/exp_plot), data = may_ev_dat, family = "binomial")

# compare different plant models
AIC(may_ev_sepb_mod, may_ev_totd_mod)
# total density
summary(may_ev_totd_mod)
# none
stepAIC(may_ev_totd_mod)
# keep full model
plot(simulateResiduals(may_ev_totd_mod))
# sig dev


## June-late August Ev ##
jun_lau_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "binomial")
jun_lau_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "binomial")
jun_lau_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "binomial")
jun_lau_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "binomial")
jun_lau_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "binomial")
jun_lau_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) * Month + (1|site/exp_plot), data = jun_lau_ev_dat, family = "binomial")

# compare different plant models
AIC(jun_lau_ev_totb_mod, jun_lau_ev_sepb_mod, jun_lau_ev_totd_mod, jun_lau_ev_sepd_mod, jun_lau_ev_totp_mod, jun_lau_ev_sepp_mod)
# total presence
summary(jun_lau_ev_totp_mod)
# fung:edge:month
stepAIC(jun_lau_ev_totp_mod)
# keep full model
plot(simulateResiduals(jun_lau_ev_totp_mod))


## September Ev ##
sep_ev_totb_mod <- glmmTMB(Bp_spots ~ fungicide * (total_biomass.g + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_sepb_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_biomass.g + Mv_biomass.g + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_totd_mod <- glmmTMB(Bp_spots ~ fungicide * (total_density + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_sepd_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_density + Mv_density + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_totp_mod <- glmmTMB(Bp_spots ~ fungicide * (bg_present + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")
sep_ev_sepp_mod <- glmmTMB(Bp_spots ~ fungicide * (Ev_present + Mv_present + edge_severity) + (1|exp_plot), data = sep_ev_dat, family = "binomial")

# compare different plant models
AIC(sep_ev_totb_mod, sep_ev_sepb_mod, sep_ev_totd_mod, sep_ev_sepd_mod, sep_ev_totp_mod, sep_ev_sepp_mod)
# separate presence
summary(sep_ev_sepp_mod)
# none
stepAIC(sep_ev_sepp_mod)
# remove fung:Ev and edge severity
sep_ev_sepp_mod2 <- update(sep_ev_sepp_mod, .~. - fungicide:Ev_present - edge_severity - fungicide:edge_severity)
summary(sep_ev_sepp_mod2)
# fung:Mv present (*)
plot(simulateResiduals(sep_ev_sepp_mod2))


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

## Edge severity Ev June-late August ##
jun_lau_ev_edge_sim_dat <- jun_lau_ev_dat %>%
  filter(!(is.na(edge_severity))) %>%
  group_by(Month, Treatment, fungicide) %>%
  summarise(min_edge = min(edge_severity),
            max_edge = max(edge_severity)) %>%
  ungroup() %>%
  rowwise() %>%
  do(tibble(Month = .$Month, Treatment = .$Treatment, fungicide = .$fungicide,
            edge_severity = seq(.$min_edge, .$max_edge, length.out = 100))) %>%
  ungroup() %>%
  mutate(bg_present = 1,
         site = NA_character_,
         exp_plot = NA_character_) %>%
  mutate(Bp_spots = predict(jun_lau_ev_totp_mod, newdata = ., re.form = NA, type = "response"),
         Bp_spots_se = predict(jun_lau_ev_totp_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

jun_lau_ev_sum_dat <- jun_lau_ev_dat %>%
  group_by(Month) %>%
  mutate(edge_bin = cut_interval(edge_severity, n = 5) %>%
           as.character()) %>%
  ungroup() %>%
  group_by(Month, edge_bin) %>%
  mutate(min_interval = parse_number(strsplit(edge_bin, ",")[[1]])[1],
         max_interval = parse_number(strsplit(edge_bin, ",")[[1]])[2],
         edge_severity = (max_interval + min_interval) / 2) %>%
  ungroup() %>%
  group_by(Month, Treatment, edge_severity) %>%
  summarise(Bp_spots_lower = mean_cl_boot(Bp_spots)$ymin,
            Bp_spots_upper = mean_cl_boot(Bp_spots)$ymax,
            Bp_spots = mean_cl_boot(Bp_spots)$y)

jun_lau_edge_ev_fig <- ggplot(jun_lau_ev_edge_sim_dat, aes(x = edge_severity, y = Bp_spots)) +
  geom_errorbar(data = jun_lau_ev_sum_dat, width = 0, aes(group = Treatment, ymin = Bp_spots_lower, ymax = Bp_spots_upper)) +
  geom_point(data = jun_lau_ev_sum_dat, size = 3, shape = 21, aes(fill = Treatment)) +
  geom_ribbon(aes(ymin = Bp_spots - Bp_spots_se, ymax = Bp_spots + Bp_spots_se, fill = Treatment), alpha = 0.5) +
  geom_line(aes(color = Treatment)) +
  facet_wrap(~ Month, scales = "free_x") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste("Proportion edge ", italic("Microstegium"), " leaf area brown", sep = ""))) +
  ylab(expression(paste("Proportion ", italic("Elymus"), " with ", italic(Bipolaris), "-like lesions", sep = ""))) +
  temp_theme


#### output ####

pdf("output/Bp_spots_analysis_2019_density_exp.pdf")
jun_lau_edge_ev_fig
dev.off()

bp_spots_edge_mod <- jun_lau_ev_totp_mod
save(bp_spots_edge_mod, file ="output/Bp_spots_ege_model_2019_density_exp.rda")
