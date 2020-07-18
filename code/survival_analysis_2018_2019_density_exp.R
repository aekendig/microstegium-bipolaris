##### info ####

# file: survival_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/15/20
# goal: evaluate the effects of density and fungicide on survival

#### need to add severity from year 1 ####


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(glmmTMB)

# import data
daty1_0 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
daty2_0 <- read_csv("data/all_replacement_2019_density_exp.csv")
bgbio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
sevy1_0 <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
sevy1_bgev_0 <- read_csv("intermediate-data/ev_background_leaf_scans_2018_density_exp.csv")
sevy2_0 <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")


#### edit data ####

# add plant type and fungicide treatment
daty1_1 <- daty1_0 %>%
  mutate(plant_type = paste(sp, age, sep = " "),
         Treatment = recode(treatment, water = "control (water)"),
         fungicide = recode(treatment, water = 0, fungicide = 1))

daty2_1 <- daty2_0 %>%
  mutate(plant_type = paste(sp, age, sep = " "))

# disease severity
sevy1_1 <- sevy1_0 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}")

sevy1_bgev_1 <- sevy1_bgev_0 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}")

sevy2_1 <- sevy2_0 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}")

# summer survival 2018
# remove NA's 
# remove background Microstegium (not individual plants)
sumy1 <- daty1_1 %>%
  filter(month == "September" & !is.na(survival) & !(sp == "Mv" & focal == 0)) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  left_join(sevy1_1) %>%
  left_join(sevy1_bgev_1)

# survival through winter given summer survival, 2018
# remove NA's 
# remove background Microstegium (not individual plants)
winy1 <- daty1_1 %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  filter(September == 1 & !is.na(April) & !(sp == "Mv" & focal == 0)) %>%
  select(-September) %>%
  rename(survival = April) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  left_join(sevy1_1) %>%
  left_join(sevy1_bgev_1)

# add IDs to 2019 background plants
# remove unnecessary columns
daty2_2 <- daty2_1 %>%
  group_by(site, plot, treatment, plant_type, focal, replace_date) %>%
  mutate(Bg_ID = as.character(1:n()),
         ID = case_when(ID == "Bg" ~ Bg_ID,
                        TRUE ~ ID)) %>%
  ungroup() %>%
  select(-c(sp, age, Bg_ID))

# focal plant ID's
focid_0 <- tibble(plant_type = c(rep(unique(daty2_1$plant_type)[1:2], each = 3), unique(daty2_1$plant_type)[3]),
              ID = c(1, 2, 3, 1, 2, 3, "A"),
              focal = 1)

# background plant ID's
bgid_0 <- plots %>% 
  select(plot, background, background_density) %>%
  unique() %>%
  filter(plot != 1) %>%
  mutate(start = 1,
         focal = 0) %>%
  group_by(plot, background, focal) %>%
  expand(start, background_density, ID = full_seq(start:background_density, 1)) %>%
  ungroup() %>%
  select(-start) %>%
  rename(plant_type = background) %>%
  mutate(ID = as.character(ID))

# merge ID lists with plot
focid_1 <- plots %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(focid_0)

bgid_1 <- plots %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  inner_join(bgid_0)

allid <- full_join(focid_1, bgid_1)

# check that number of plants is correct
nrow(allid)
2*4*(7*3 + 4+7 + 16+7 + 64+7 + 4+7 + 8+7 + 16+7 + 2+7 + 4+7 + 8+7)
# yes - matches 
  
# summer survival 2019
sumy2 <- daty2_2 %>%
  select(-replace_date) %>%
  unique() %>%
  mutate(survival = 0) %>%
  full_join(allid) %>%
  mutate(survival = replace_na(survival, 1)) %>%
  left_join(covar) %>%
  left_join(bgbio %>%
              rename(background_biomass = biomass.g)) %>%
  mutate(background_pc_bio = background_biomass / background_density,
         Treatment = recode(treatment, water = "control (water)"),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         sp = substr(plant_type, 1, 2)) %>%
  left_join(sevy2_1)

# process focal Ev biomass for direct disease effects (plot 1)
# check for missing data
filter(ev_bio, plot == 1) %>%
  summarise(missing = sum(is.na(weight))) # none
# add up all Ev in plot
ev_bio_1 <- ev_bio %>%
  filter(plot == 1) %>%
  group_by(site, treatment) %>%
  summarise(ev_biomass = sum(weight)) %>%
  ungroup()

# no background data
no_back_plant_dat <- sumy1 %>%
  filter(plot == 1 & background == "Mv seedling") %>%
  mutate(yearf = "year 1") %>%
  full_join(sumy2 %>%
              filter(plot == 1 & background == "Mv seedling") %>%
              mutate(yearf = "year 2"))

# Ev seedling no background winter data
winy1 %>%
  filter(plot == 1 & background == "Mv seedling" & plant_type == "Ev seedling")
# only 7 plants because most didn't survive first year

winy1 %>%
  filter(plot == 1 & background == "Mv seedling" & plant_type == "Ev adult")
# only 7 plants because most didn't survive first year


#### visualize treatment effects ####

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

# template figures
temp_fig <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free_x", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

temp_fig2 <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  geom_smooth(method = "glm", aes(color = treatment)) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  scale_color_manual(values = c("black", "red"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

temp_fig3 <- ggplot(plots, aes(x = treatment, y = plot)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  temp_theme +
  xlab("Background")

# save treatment figures
pdf("./output/survival_visualize_treatment_2018_2019_density_exp.pdf")

# focal plants, summer survival, 2018
temp_fig %+%
  filter(sumy1, focal == 1) %+%
  aes(y = survival) +
  ylab("Year 1 summer focal survival")

# all plants, summer survival, 2018
temp_fig %+%
  sumy1 %+%
  aes(y = survival) +
  ylab("Year 1 summer all survival")

# focal plants, winter survival, 2018
temp_fig %+%
  filter(winy1, focal == 1) %+%
  aes(y = survival) +
  ylab("Year 1 winter focal survival")

# all plants, winter survival, 2018
temp_fig %+%
  winy1 %+%
  aes(y = survival) +
  ylab("Year 1 winter all survival")   
  
# focal plants, summer survival, 2019
temp_fig %+%
  filter(sumy2, focal == 1) %+%
  aes(y = survival) +
  ylab("Year 2 summer focal survival")

# all plants, summer survival, 2019
temp_fig %+%
  sumy2 %+%
  aes(y = survival) +
  ylab("Year 2 summer all survival") 

# focal plants, summer survival, biomass 2019
temp_fig2 %+%
  filter(sumy2, focal == 1) %+%
  aes(x = background_biomass, y = survival) +
  ylab("Year 2 summer focal survival") +
  xlab("Background biomass (g)")

# focal plants, summer survival, per capita biomass 2019
temp_fig2 %+%
  filter(sumy2, focal == 1) %+%
  aes(x = background_pc_bio, y = survival) +
  ylab("Year 2 summer focal survival") +
  xlab("Per capita ackground biomass (g)")

# focal plants, summer survival, treatment 2018
temp_fig3 %+%
  filter(sumy1, focal == 1 & background_density > 0) +
  ylab("Year 1 summer focal survival")

# focal plants, winter survival, treatment 2018
temp_fig3 %+%
  filter(sumy1, focal == 1 & background_density > 0) +
  ylab("Year 1 winter focal survival")

# focal plants, summer survival, treatment 2019
temp_fig3 %+%
  filter(sumy2, focal == 1 & background_density > 0) +
  ylab("Year 2 summer focal survival")

dev.off()

pdf("output/focal_survival_no_background_treatment_2018_2019_density_exp.pdf", width = 6, height = 3)
ggplot(no_back_plant_dat, aes(x = Treatment, y = survival)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_wrap(~plant_type) +
  scale_fill_manual(values = col_pal, guide = F) +
  ylab("Summer survival") +
  temp_theme
dev.off()

ggplot(no_back_plant_dat, aes(x = Treatment, y = survival)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_grid(yearf~plant_type) +
  scale_fill_manual(values = col_pal, guide = F) +
  ylab("Summer survival") +
  temp_theme


#### visualize covariate effects ####

# covariate template figure
temp_fig_cov <- ggplot(sumy1, aes(x = soil_moisture_jun.prop, y = survival)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "loess") +
  facet_wrap(~ plant_type) +
  temp_theme

# save covariate figures
pdf("./output/survival_visualize_covariates_2018_2019_density_exp.pdf")

# soil moisture summer survival year 1
temp_fig_cov +
  xlab("Proportion soil moisture") +
  ylab("Year 1 summer all survival")
# include this in model

# soil moisture winter survival year 1
temp_fig_cov %+%
  winy1 +
  xlab("Proportion soil moisture") +
  ylab("Year 1 winter all survival")  

# soil moisture summer survival year 2
temp_fig_cov %+%
  sumy2 +
  xlab("Proportion soil moisture") +
  ylab("Year 2 summer all survival")

# October soil moisture summer survival year 2
temp_fig_cov %+%
  sumy2 %+%
  aes(x = soil_moisture_oct.prop) +
  xlab("Proportion soil moisture in October") +
  ylab("Year 2 summer all survival")
# include this in model

# canopy cover summer survival year 1
temp_fig_cov %+%
  aes(x = canopy_cover.prop) +
  xlab("Proportion canopy cover") +
  ylab("Year 1 summer all survival")

# canopy cover winter survival year 1
temp_fig_cov %+%
  winy1 %+%
  aes(x = canopy_cover.prop) +
  xlab("Proportion canopy cover") +
  ylab("Year 1 winter all survival")  

# canopy cover summer survival year 2
temp_fig_cov %+%
  sumy2 %+%
  aes(x = canopy_cover.prop) +
  xlab("Proportion canopy cover") +
  ylab("Year 2 summer all survival")

# September biomass summer survival year 1
temp_fig_cov %+%
  aes(x = mv_sep.g) +
  xlab("September biomass") +
  ylab("Year 1 summer all survival")

# September biomass winter survival year 1
temp_fig_cov %+%
  winy1 %+%
  aes(x = mv_sep.g) +
  xlab("September biomass") +
  ylab("Year 1 winter all survival")  

# September biomass summer survival year 2
temp_fig_cov %+%
  sumy2 %+%
  aes(x = mv_sep.g) +
  xlab("September biomass") +
  ylab("Year 2 summer all survival")

dev.off()


#### visualize severity effects ####

# July without background
temp_fig_sev <- ggplot(no_back_plant_dat, aes(x = severity_jul, y = survival)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(treatment ~ plant_type, scale = "free") +
  temp_theme
temp_fig_sev
# sdecrease with Mv seedling

# August without background
temp_fig_sev %+%
  aes(x = severity_late_aug)
# not a strong effect

# winter
temp_fig_sev %+%
  winy1 %+%
  aes(x = severity_sep)
# slight reduction with water


#### direct disease models ####

# divide data
mv_no_back_plant_dat <- filter(no_back_plant_dat, plant_type == "Mv seedling")
evs_no_back_plant_dat <- filter(no_back_plant_dat, plant_type == "Ev seedling")
eva_no_back_plant_dat <- filter(no_back_plant_dat, plant_type == "Ev adult")

# Mv model
mv_no_back_surv_mod <- brm(data = mv_no_back_plant_dat, family = bernoulli,
                           survival ~ fungicide + (1|site) + (1|yearf),
                           prior <- c(prior(normal(0, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd)),
                           iter = 6000, warmup = 1000, chains = 1, 
                           control = list(adapt_delta = 0.9999))

# check model and add chains
summary(mv_no_back_surv_mod)
prior_summary(mv_no_back_surv_mod)
mv_no_back_surv_mod <- update(mv_no_back_surv_mod, chains = 3)
plot(mv_no_back_surv_mod)
pp_check(mv_no_back_surv_mod, nsamples = 100)

# Ev seedling model
evs_no_back_surv_mod <- brm(data = evs_no_back_plant_dat, family = bernoulli,
                           survival ~ fungicide + (1|site) + (1|yearf),
                           prior <- c(prior(normal(0, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd)),
                           iter = 6000, warmup = 1000, chains = 1, 
                           control = list(adapt_delta = 0.9999))

# check model and add chains
summary(evs_no_back_surv_mod)
prior_summary(evs_no_back_surv_mod)
evs_no_back_surv_mod <- update(evs_no_back_surv_mod, chains = 3)
plot(evs_no_back_surv_mod)
pp_check(evs_no_back_surv_mod, nsamples = 100)


# Ev adult model
eva_no_back_surv_mod <- brm(data = eva_no_back_plant_dat, family = bernoulli,
                            survival ~ fungicide + (1|site) + (1|yearf),
                            prior <- c(prior(normal(0, 10), class = Intercept),
                                       prior(normal(0, 10), class = b),
                                       prior(cauchy(0, 1), class = sd)),
                            iter = 6000, warmup = 1000, chains = 1, 
                            control = list(adapt_delta = 0.9999))

# check model and add chains
summary(eva_no_back_surv_mod)
prior_summary(eva_no_back_surv_mod)
eva_no_back_surv_mod <- update(eva_no_back_surv_mod, chains = 3)
plot(eva_no_back_surv_mod)
pp_check(eva_no_back_surv_mod, nsamples = 100)


#### Mv survival models ####

# separate data
mv_mv_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Mv seedling" & background == "Mv seedling")

mv_evs_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Mv seedling" & background == "Ev seedling")

mv_eva_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Mv seedling" & background == "Ev adult")

# check for missing data
sum(is.na(mv_mv_surv_dat$survival))
sum(is.na(mv_evs_surv_dat$survival))
sum(is.na(mv_eva_surv_dat$survival))

# Mv background
mv_mv_surv_mod <- brm(data = mv_mv_surv_dat, family = bernoulli,
                      survival ~ fungicide * background_density + (1|site),
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(normal(0, 10), class = b),
                                prior(cauchy(0, 1), class = sd)),
                      iter = 6000, warmup = 1000, chains = 1)

# check model and increase chains
summary(mv_mv_surv_mod)
prior_summary(mv_mv_surv_mod)
mv_mv_surv_mod <- update(mv_mv_surv_mod, chains = 3)
plot(mv_mv_surv_mod)
pp_check(mv_mv_surv_mod, nsamples = 100)

# Ev seedling background
mv_evs_surv_mod <- update(mv_mv_surv_mod, newdata = mv_evs_surv_dat,
                          control = list(adapt_delta = 0.99))

# check model and increase chains
summary(mv_evs_surv_mod)
plot(mv_evs_surv_mod)
pp_check(mv_evs_surv_mod, nsamples = 100)

# Ev adult background
mv_eva_surv_mod <- update(mv_mv_surv_mod, newdata = mv_eva_surv_dat,
                          control = list(adapt_delta = 0.99))

# check model and increase chains
summary(mv_eva_surv_mod)
plot(mv_eva_surv_mod)
pp_check(mv_eva_surv_mod, nsamples = 100)


#### Ev seedling survival models ####

# separate data
evs_mv_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Ev seedling" & background == "Mv seedling")

evs_evs_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Ev seedling" & background == "Ev seedling")

evs_eva_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Ev seedling" & background == "Ev adult")

# check for missing data
sum(is.na(evs_mv_surv_dat$survival))
sum(is.na(evs_evs_surv_dat$survival))
sum(is.na(evs_eva_surv_dat$survival))

# Mv background
evs_mv_surv_mod <- update(mv_mv_surv_mod, newdata = evs_mv_surv_dat,
                          control = list(adapt_delta = 0.99))

# check model and increase chains
summary(evs_mv_surv_mod)
plot(evs_mv_surv_mod)
pp_check(evs_mv_surv_mod, nsamples = 100)

# Ev seedling background
evs_evs_surv_mod <- update(mv_mv_surv_mod, newdata = evs_evs_surv_dat,
                           control = list(adapt_delta = 0.99))

# check model and increase chains
summary(evs_evs_surv_mod)
plot(evs_evs_surv_mod)
pp_check(evs_evs_surv_mod, nsamples = 100)

# Ev adult background
evs_eva_surv_mod <- update(mv_mv_surv_mod, newdata = evs_eva_surv_dat,
                           control = list(adapt_delta = 0.99))

# check model and increase chains
summary(evs_eva_surv_mod)
plot(evs_eva_surv_mod)
pp_check(evs_eva_surv_mod, nsamples = 100)


#### Ev adult survival models ####

# separate data
eva_mv_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Ev adult" & background == "Mv seedling")

eva_evs_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Ev adult" & background == "Ev seedling")

eva_eva_surv_dat <- sumy2 %>%
  filter(focal == 1 & plant_type == "Ev adult" & background == "Ev adult")

# check for missing data
sum(is.na(eva_mv_surv_dat$survival))
sum(is.na(eva_evs_surv_dat$survival))
sum(is.na(eva_eva_surv_dat$survival))

# Mv background
eva_mv_surv_mod <- update(mv_mv_surv_mod, newdata = eva_mv_surv_dat,
                          control = list(adapt_delta = 0.99, max_treedepth = 15))

# check model and increase chains
summary(eva_mv_surv_mod)
plot(eva_mv_surv_mod)
pp_check(eva_mv_surv_mod, nsamples = 100)

# Ev seedling background
eva_evs_surv_mod <- update(mv_mv_surv_mod, newdata = eva_evs_surv_dat)

# check model and increase chains
summary(eva_evs_surv_mod)
plot(eva_evs_surv_mod)
pp_check(eva_evs_surv_mod, nsamples = 100)

# Ev adult background
eva_eva_surv_mod <- update(mv_mv_surv_mod, newdata = eva_eva_surv_dat,
                           control = list(adapt_delta = 0.99))

# check model and increase chains
summary(eva_eva_surv_mod)
plot(eva_eva_surv_mod)
pp_check(eva_eva_surv_mod, nsamples = 100)


#### model fit figure ####

# model fits over simulated data
surv_mv_sim_dat <- tibble(background_density = rep(seq(0, 64, length.out = 300), 2),
                          fungicide = rep(c(0, 1), each = 300)) %>%
  merge(tibble(plant_type = c("Mv seedling", "Ev seedling", "Ev adult"))) %>%
  as_tibble() %>%
  mutate(site = NA) %>%
  mutate(survival = case_when(plant_type == "Mv seedling" ~ fitted(mv_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
                              plant_type == "Ev seedling" ~ fitted(evs_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
                              plant_type == "Ev adult" ~ fitted(eva_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"]),
         survival_lower = case_when(plant_type == "Mv seedling" ~ fitted(mv_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
                                    plant_type == "Ev seedling" ~ fitted(evs_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
                                    plant_type == "Ev adult" ~ fitted(eva_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"]),
         survival_upper = case_when(plant_type == "Mv seedling" ~ fitted(mv_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
                                    plant_type == "Ev seedling" ~ fitted(evs_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
                                    plant_type == "Ev adult" ~ fitted(eva_mv_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]))

surv_evs_sim_dat <- tibble(background_density = rep(seq(0, 16, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  merge(tibble(plant_type = c("Mv seedling", "Ev seedling", "Ev adult"))) %>%
  as_tibble() %>%
  mutate(site = NA) %>%
  mutate(survival = case_when(plant_type == "Mv seedling" ~ fitted(mv_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
                              plant_type == "Ev seedling" ~ fitted(evs_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
                              plant_type == "Ev adult" ~ fitted(eva_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"]),
         survival_lower = case_when(plant_type == "Mv seedling" ~ fitted(mv_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
                                    plant_type == "Ev seedling" ~ fitted(evs_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
                                    plant_type == "Ev adult" ~ fitted(eva_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"]),
         survival_upper = case_when(plant_type == "Mv seedling" ~ fitted(mv_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
                                    plant_type == "Ev seedling" ~ fitted(evs_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
                                    plant_type == "Ev adult" ~ fitted(eva_evs_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]))

surv_eva_sim_dat <- tibble(background_density = rep(seq(0, 8, length.out = 100), 2),
                           fungicide = rep(c(0, 1), each = 100)) %>%
  merge(tibble(plant_type = c("Mv seedling", "Ev seedling", "Ev adult"))) %>%
  as_tibble() %>%
  mutate(site = NA) %>%
  mutate(survival = case_when(plant_type == "Mv seedling" ~ fitted(mv_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
                              plant_type == "Ev seedling" ~ fitted(evs_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
                              plant_type == "Ev adult" ~ fitted(eva_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Estimate"]),
         survival_lower = case_when(plant_type == "Mv seedling" ~ fitted(mv_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
                                    plant_type == "Ev seedling" ~ fitted(evs_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
                                    plant_type == "Ev adult" ~ fitted(eva_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q2.5"]),
         survival_upper = case_when(plant_type == "Mv seedling" ~ fitted(mv_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
                                    plant_type == "Ev seedling" ~ fitted(evs_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
                                    plant_type == "Ev adult" ~ fitted(eva_eva_surv_mod, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]))

# combine simulated data
surv_sim_dat <- surv_mv_sim_dat %>%
  mutate(background = "Mv seedling") %>%
  full_join(surv_evs_sim_dat %>%
              mutate(background = "Ev seedling")) %>%
  full_join(surv_eva_sim_dat %>%
              mutate(background = "Ev adult")) %>%
  mutate(Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"))

# figure
pdf("output/suvival_analysis_model_fits_2019_density_exp.pdf")
ggplot(filter(sumy2, focal == 1), aes(x = background_density, y = survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = Treatment)) +
  geom_ribbon(data = surv_sim_dat, alpha = 0.5, aes(ymin = survival_lower, ymax = survival_upper, fill = Treatment)) +
  geom_line(data = surv_sim_dat, aes(color = Treatment)) +
  facet_grid(plant_type ~ background, scales = "free_x", switch = "both") +
  xlab("Background density") +
  ylab("Summer survival") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  temp_theme
dev.off()


#### severity and survival ####

# divide data
# remove extra 1 plots
mv_sev_2018_dat <- filter(sumy1, treatment == "water" & plant_type == "Mv seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_"))
evs_sev_2018_dat <- filter(sumy1, treatment == "water" & plant_type == "Ev seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_"))
eva_sev_2018_dat <- filter(sumy1, treatment == "water" & plant_type == "Ev adult" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_"))
mv_sev_2019_dat <- filter(sumy2, treatment == "water" & plant_type == "Mv seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_"))
evs_sev_2019_dat <- filter(sumy2, treatment == "water" & plant_type == "Ev seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_"))
eva_sev_2019_dat <- filter(sumy2, treatment == "water" & plant_type == "Ev adult" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_"))

# combine to test correlations with density
dens_sev_dat <- mv_sev_2018_dat %>%
  group_by(site, plot, treatment, background_density) %>%
  summarise(mvjul18 = mean(severity_jul, na.rm = T), 
            mvlau18 = mean(severity_late_aug, na.rm = T), 
            mvsep18 = mean(severity_sep, na.rm = T)) %>%
  full_join(evs_sev_2018_dat %>%
              group_by(site, plot, treatment, background_density) %>%
              summarise(evsjul18 = mean(severity_jul, na.rm = T), 
                        evslau18 = mean(severity_late_aug, na.rm = T), 
                        evssep18 = mean(severity_sep, na.rm = T))) %>%
  full_join(eva_sev_2018_dat %>%
              group_by(site, plot, treatment, background_density) %>%
              summarise(evajul18 = mean(severity_jul, na.rm = T), 
                        evalau18 = mean(severity_late_aug, na.rm = T), 
                        evasep18 = mean(severity_sep, na.rm = T))) %>%
  full_join(mv_sev_2019_dat %>%
              group_by(site, plot, treatment, background_density) %>%
              summarise(mvjun19 = mean(severity_jun, na.rm = T), 
                        mvjul19 = mean(severity_jul, na.rm = T), 
                        mveau19 = mean(severity_early_aug, na.rm = T), 
                        mvlau19 = mean(severity_late_aug, na.rm = T), 
                        mvsep19 = mean(severity_sep, na.rm = T))) %>%
  full_join(evs_sev_2019_dat %>%
              group_by(site, plot, treatment, background_density) %>%
              summarise(evsjun19 = mean(severity_jun, na.rm = T), 
                        evsjul19 = mean(severity_jul, na.rm = T), 
                        evseau19 = mean(severity_early_aug, na.rm = T), 
                        evslau19 = mean(severity_late_aug, na.rm = T), 
                        evssep19 = mean(severity_sep, na.rm = T))) %>%
  full_join(eva_sev_2019_dat %>%
              group_by(site, plot, treatment, background_density) %>%
              summarise(evajun19 = mean(severity_jun, na.rm = T), 
                        evajul19 = mean(severity_jul, na.rm = T), 
                        evaeau19 = mean(severity_early_aug, na.rm = T), 
                        evalau19 = mean(severity_late_aug, na.rm = T), 
                        evasep19 = mean(severity_sep, na.rm = T))) %>%
  ungroup() %>%
  as.data.frame()

# check for density-severity correlations
for(i in 5:28){
  dens_sev_dat$severity = dens_sev_dat[,i]
  print(cor.test(dens_sev_dat$background_density, dens_sev_dat$severity))
}
# none significantly correlated

# survival-severity relationships
mv_surv_sev_jul_18_mod <- glmmTMB(survival ~ severity_jul + (1|site/plot), family = binomial, data = mv_sev_2018_dat)
summary(mv_surv_sev_jul_18_mod) # not sig
mv_surv_sev_lau_18_mod <- glmmTMB(survival ~ severity_late_aug + (1|site/plot), family = binomial, data = mv_sev_2018_dat)
summary(mv_surv_sev_lau_18_mod) # not sig
mv_surv_sev_sep_18_mod <- glmmTMB(survival ~ severity_sep + (1|site/plot), family = binomial, data = mv_sev_2018_dat)
summary(mv_surv_sev_sep_18_mod) # convergence issue
mv_surv_sev_sep_18_mod <- glmmTMB(survival ~ severity_sep + (1|site_plot), family = binomial, data = mv_sev_2018_dat)
mv_surv_sev_sep_18_mod <- glmmTMB(survival ~ severity_sep + (1|site), family = binomial, data = mv_sev_2018_dat)
mv_surv_sev_sep_18_mod <- glmmTMB(survival ~ severity_sep + (1|plot), family = binomial, data = mv_sev_2018_dat)
# can't get sep model to converge

evs_surv_sev_jul_18_mod <- glmmTMB(survival ~ severity_jul + (1|site/plot), family = binomial, data = evs_sev_2018_dat)
summary(evs_surv_sev_jul_18_mod) # not sig
evs_surv_sev_lau_18_mod <- glmmTMB(survival ~ severity_late_aug + (1|site/plot), family = binomial, data = evs_sev_2018_dat)
summary(evs_surv_sev_lau_18_mod) # not sig
evs_surv_sev_sep_18_mod <- glmmTMB(survival ~ severity_sep + (1|site/plot), family = binomial, data = evs_sev_2018_dat)
summary(evs_surv_sev_sep_18_mod) # also doesn't converge. it doesn't work because the plant had to survive until september for the measurement to be taken

eva_surv_sev_jul_18_mod <- glmmTMB(survival ~ severity_jul + (1|site/plot), family = binomial, data = eva_sev_2018_dat)
summary(eva_surv_sev_jul_18_mod) # convergence issue
eva_surv_sev_jul_18_mod <- glmmTMB(survival ~ severity_jul + (1|site_plot), family = binomial, data = eva_sev_2018_dat)
# not sig
eva_surv_sev_lau_18_mod <- glmmTMB(survival ~ severity_late_aug + (1|site/plot), family = binomial, data = eva_sev_2018_dat) # convergence issue
eva_surv_sev_lau_18_mod <- glmmTMB(survival ~ severity_late_aug + (1|site_plot), family = binomial, data = eva_sev_2018_dat) # convergence issue
# survival likely too high in Ev adults

mv_surv_sev_jun_19_mod <- glmmTMB(survival ~ severity_jun + (1|site/plot), family = binomial, data = mv_sev_2019_dat)
summary(mv_surv_sev_jun_19_mod) # not sig
mv_surv_sev_jul_19_mod <- glmmTMB(survival ~ severity_jul + (1|site/plot), family = binomial, data = mv_sev_2019_dat)
summary(mv_surv_sev_jul_19_mod) # not sig
mv_surv_sev_eau_19_mod <- glmmTMB(survival ~ severity_early_aug + (1|site/plot), family = binomial, data = mv_sev_2019_dat)
summary(mv_surv_sev_eau_19_mod) # not sig
mv_surv_sev_lau_19_mod <- glmmTMB(survival ~ severity_late_aug + (1|site/plot), family = binomial, data = mv_sev_2019_dat)
summary(mv_surv_sev_lau_19_mod) # not sig

evs_surv_sev_jun_19_mod <- glmmTMB(survival ~ severity_jun + (1|site/plot), family = binomial, data = evs_sev_2019_dat)
summary(evs_surv_sev_jun_19_mod) # not sig
evs_surv_sev_jul_19_mod <- glmmTMB(survival ~ severity_jul + (1|site/plot), family = binomial, data = evs_sev_2019_dat)
summary(evs_surv_sev_jul_19_mod) # sig reduction
evs_surv_sev_eau_19_mod <- glmmTMB(survival ~ severity_early_aug + (1|site/plot), family = binomial, data = evs_sev_2019_dat)
summary(evs_surv_sev_eau_19_mod) # not sig
evs_surv_sev_lau_19_mod <- glmmTMB(survival ~ severity_late_aug + (1|site/plot), family = binomial, data = evs_sev_2019_dat)
summary(evs_surv_sev_lau_19_mod) # not sig

eva_surv_sev_jun_19_mod <- glmmTMB(survival ~ severity_jun + (1|site/plot), family = binomial, data = eva_sev_2019_dat)
summary(eva_surv_sev_jun_19_mod) # not sig
eva_surv_sev_jul_19_mod <- glmmTMB(survival ~ severity_jul + (1|site/plot), family = binomial, data = eva_sev_2019_dat)
summary(eva_surv_sev_jul_19_mod) # not sig
eva_surv_sev_eau_19_mod <- glmmTMB(survival ~ severity_early_aug + (1|site/plot), family = binomial, data = eva_sev_2019_dat) # convergence issue
eva_surv_sev_eau_19_mod <- glmmTMB(survival ~ severity_early_aug + (1|site_plot), family = binomial, data = eva_sev_2019_dat)
summary(eva_surv_sev_eau_19_mod) # sig reduction
eva_surv_sev_lau_19_mod <- glmmTMB(survival ~ severity_late_aug + (1|site/plot), family = binomial, data = eva_sev_2019_dat)
summary(eva_surv_sev_lau_19_mod) # not sig

  

#### output ####

# models
save(mv_no_back_surv_mod, file = "output/mv_survival_no_background_model_2018_2019_density_exp.rda")
save(evs_no_back_surv_mod, file = "output/ev_seedling_survival_no_background_model_2018_2019_density_exp.rda")
save(eva_no_back_surv_mod, file = "output/ev_adult_survival_no_background_model_2018_2019_density_exp.rda")
save(mv_mv_surv_mod, file = "output/mv_survival_mv_background_model_2019_density_exp.rda")
save(mv_evs_surv_mod, file = "output/mv_survival_ev_seedling_background_model_2019_density_exp.rda")
save(mv_eva_surv_mod, file = "output/mv_survival_ev_adult_background_model_2019_density_exp.rda")
save(evs_mv_surv_mod, file = "output/ev_seedling_survival_mv_background_model_2019_density_exp.rda")
save(evs_evs_surv_mod, file = "output/ev_seedling_survival_ev_seedling_background_model_2019_density_exp.rda")
save(evs_eva_surv_mod, file = "output/ev_seedling_survival_ev_adult_background_model_2019_density_exp.rda")
save(eva_mv_surv_mod, file = "output/ev_adult_survival_mv_background_model_2019_density_exp.rda")
save(eva_evs_surv_mod, file = "output/ev_adult_survival_ev_seedling_background_model_2019_density_exp.rda")
save(eva_eva_surv_mod, file = "output/ev_adult_survival_ev_adult_background_model_2019_density_exp.rda")
save(evs_surv_sev_jul_19_mod, file = "output/ev_seedling_survival_severity_model_jul_2019_density_exp.rda")
write_csv(evs_sev_2019_dat, "intermediate-data/ev_seedling_survival_severity_data_2019_dens_exp.csv")
save(eva_surv_sev_eau_19_mod, file = "output/ev_seedling_survival_severity_model_early_aug_2019_density_exp.rda")
write_csv(eva_sev_2019_dat, "intermediate-data/ev_adult_survival_severity_data_2019_dens_exp.csv")
