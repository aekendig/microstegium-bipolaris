##### info ####

# file: mv_germination_disease_analysis_2018_density_exp
# author: Amy Kendig
# date last edited: 7/14/20
# goal: analyze Mv germination and disease occurrence


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(glmmTMB)

# import data
set1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
set2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
sev <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")


#### edit data ####

# notes
unique(set1$notes_check_1) # seedlings with lesions in notes
unique(set1$notes_check_2) # may be a contaminate on plate
unique(set1$notes_check_3)
unique(set1$notes_germination_check_1) # some plates were put into fridge during one day
unique(set1$notes_germination_final) 
unique(set2$notes) # may be a contaminate on plate

# check columns
unique(set1$germination_check_1)
unique(set1$germination_final)

# condense info in set 1
set1b <- set1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final),
         seeds_dark_check_1 = rowSums(cbind(seeds_dark_check_1, seeds_red_check_1, seeds_pink_check_1, seeds_green_check_1), na.rm = T),
         seeds_dark_check_2 = rowSums(cbind(seeds_dark_check_2, seeds_red_check_2, seeds_pink_check_2, seeds_green_check_2), na.rm = T),
         seeds_dark_check_3 = rowSums(cbind(seeds_dark_check_3, seeds_red_check_3, seeds_green_check_3), na.rm = T)) %>%
  select(colnames(set2)[-c(12:13)], seeds_dark_check_3, seeds_light_check_3, date_check_3, notes_check_1, notes_check_2, notes_check_3, notes_germination_check_1, notes_germination_final)

# disease severity
sevb <- sev %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  group_by(month, site, plot, treatment, sp) %>%
  summarise(plant_severity = mean(plant_severity, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Mv")

# combine data
dat <- full_join(set1b, set2) %>%
  rowwise() %>%
  mutate(seeds_dark = max(cbind(seeds_dark_check_1, seeds_dark_check_2, seeds_dark_check_3), na.rm = T),
         seeds_light = max(cbind(seeds_light_check_1, seeds_light_check_2, seeds_light_check_3), na.rm = T)) %>%
  ungroup() %>%
  mutate(germination.prop = germination_final / seeds,
         dark.prop = seeds_dark / seeds,
         light.prop = seeds_light / seeds,
         site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           recode("F" = "fungicide", "W" = "water"),
         fungicide = recode(treatment, fungicide = "1", water = "0") %>% as.numeric(),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide")) %>%
  left_join(sevb)

# add treatments
dat_treat <- dat  %>%
  left_join(treat)
dat_plots <- dat  %>%
  left_join(plots)


#### visualizations ####

# template theme
temp_theme <- theme_bw() +
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

# colors
col_pal = c("#a6611a", "#018571")

# different check times
ggplot(dat, aes(x = seeds_dark_check_1, y = seeds_dark_check_2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)
# most have more during check 2

ggplot(dat, aes(x = seeds_dark_check_2, y = seeds_dark_check_3)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)
# a handful have more during check 3

# germination with treatments
ggplot(dat_plots, aes(x = background_density, y = germination.prop)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.01, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, aes(color = Treatment)) +
  facet_wrap(~background, scales = "free_x") +
  scale_color_manual(values = col_pal) +
  temp_theme

# infection with treatments
ggplot(dat_plots, aes(x = background_density, y = dark.prop)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.01, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, aes(color = Treatment)) +
  facet_wrap(~background, scales = "free_x") +
  scale_color_manual(values = col_pal) +
  temp_theme

ggplot(dat_plots, aes(x = background_density, y = light.prop)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.01, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, aes(color = Treatment)) +
  facet_wrap(~background, scales = "free_x") +
  scale_color_manual(values = col_pal) +
  temp_theme


#### direct disease models ####

# filter data
no_back_plant_dat <- filter(dat_treat, plot == 1)

# model
mv_no_back_germ_mod <- brm(data = no_back_plant_dat, family = binomial,
                           germination_final | trials(seeds) ~ fungicide + (1|site/plot),
                           prior <- c(prior(normal(0, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd)),
                           iter = 6000, warmup = 1000, chains = 1, 
                           control = list(adapt_delta = 0.99))

# check model and add chains
summary(mv_no_back_germ_mod)
prior_summary(mv_no_back_germ_mod)
mv_no_back_germ_mod <- update(mv_no_back_germ_mod, chains = 3, 
                              control = list(adapt_delta = 0.9999))
plot(mv_no_back_germ_mod)
pp_check(mv_no_back_germ_mod, nsamples = 100)


## high Mv model ##

# filter data
hi_mv_plant_dat <- filter(dat_treat, plot == 4)

# model
mv_hi_mv_germ_mod <- brm(data = hi_mv_plant_dat, family = binomial,
                           germination_final | trials(seeds) ~ fungicide + (1|site/plot),
                           prior <- c(prior(normal(0, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd)),
                           iter = 6000, warmup = 1000, chains = 1, 
                           control = list(adapt_delta = 0.99))

# check model and add chains
summary(mv_hi_mv_germ_mod)
mv_hi_mv_germ_mod <- update(mv_hi_mv_germ_mod, chains = 3)
plot(mv_hi_mv_germ_mod)
pp_check(mv_hi_mv_germ_mod, nsamples = 100)


# model with all data (assume density of source plot doesn't affect germination)
mv_germ_mod <- brm(data = dat_treat, family = binomial,
                           germination_final | trials(seeds) ~ fungicide + (1|site/plot),
                           prior <- c(prior(normal(0, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd)),
                           iter = 6000, warmup = 1000, chains = 1, 
                           control = list(adapt_delta = 0.99))

# check model and add chains
summary(mv_germ_mod)
prior_summary(mv_germ_mod)
mv_germ_mod <- update(mv_germ_mod, chains = 3,
                      control = list(adapt_delta = 0.9999))
plot(mv_germ_mod)
pp_check(mv_germ_mod, nsamples = 100)


#### severity and germination ####

# filter data
mv_sev_2018_dat <- filter(dat_treat, treatment == "water") %>% 
  mutate(non_germination_final = seeds - germination_final) %>%
  as.data.frame()

# check for density-severity correlations
for(i in 33:35){
  mv_sev_2018_dat$severity = mv_sev_2018_dat[,i]
  print(cor.test(mv_sev_2018_dat$background_density, mv_sev_2018_dat$severity))
}
# sep is correlated 0.2

# germination-severity relationships
mv_germ_sev_jul_18_mod <- glmmTMB(cbind(germination_final, non_germination_final) ~ severity_jul + (1|site/plot), family = binomial, data = mv_sev_2018_dat)
summary(mv_germ_sev_jul_18_mod) # not sig
mv_germ_sev_lau_18_mod <- glmmTMB(cbind(germination_final, non_germination_final) ~ severity_late_aug + (1|site/plot), family = binomial, data = mv_sev_2018_dat)
summary(mv_germ_sev_lau_18_mod) # not sig
mv_germ_sev_sep_18_mod <- glmmTMB(cbind(germination_final, non_germination_final) ~ severity_sep + (1|site/plot), family = binomial, data = mv_sev_2018_dat)
summary(mv_germ_sev_sep_18_mod) # not sig


#### output ####
save(mv_no_back_germ_mod, file = "output/mv_germination_no_background_model_2018_density_exp.rda")
save(mv_hi_mv_germ_mod, file = "output/mv_germination_high_mv_model_2018_density_exp.rda")
save(mv_germ_mod, file = "output/mv_germination_all_data_model_2018_density_exp.rda")
