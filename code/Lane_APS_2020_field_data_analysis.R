##### info ####

# file: Lane_APS_2020_field_data_analysis
# author: Amy Kendig
# date last edited: 3/13/20
# goals:
# 1. Do the proportion of individuals with Bipolaris-looking lesions differ between Mv and Ev? Did the fungicide treatment affect this?
# 2. Do the proportions of samples in which Ashish could identify Bipolaris gigantea differ between Mv and Ev? Did the fungicide treatment affect this?
# 3. Does the effect of fungicide (and inverse, Bipolaris infection) on seed production differ for Ev and Mv?
# notes:
# biomass was only measure for Mv in 2018
# 2018 data only because we're not done with seed counts or biomass yet
# could extend this to look at correlation between severity and seed production, but severity includes non-Bipolaris lesions, especially on Ev


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lme4)
library(car)
library(cowplot)

# import infection data - check metadata and add these to the methods doc
all_dis_aug <- read_csv("./data/all_disease_seeds_late_aug_2018_density_exp.csv")
ev_dis_sep <- read_csv("./data/ev_disease_seeds_sep_2018_density_exp.csv")
mv_dis_sep <- read_csv("./data/mv_disease_sep_2018_density_exp.csv")

# import isolation data
iso_jul <- read_csv("./data/fungal_isolation_jul_2018_density_exp.csv")
iso_laug <- read_csv("./data/fungal_isolation_late_aug_2018_density_exp.csv")
iso_sep <- read_csv("./data/fungal_isolation_sep_2018_density_exp.csv")

# import seed data
ev_seeds <- read_csv("./intermediate-data/ev_processed_seeds_2018_density_exp.csv")
mv_seeds <- read_csv("./intermediate-data/mv_processed_seeds_2018_density_exp.csv")


#### edit data ####

# create an infected column for Mv (any infection) and Ev (Bp_spots_Ev)
# select focals for consistent sampling across plots
# remove dead plants 
# add a fungicide column, and a site/plot column
all_dis_aug2 <- all_dis_aug %>%
  mutate(infected = case_when(sp == "Ev" ~ Bp_spots_Ev,
                              sp == "Mv" ~ infec),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_")) %>%
  filter(ID %in% c("1", "2", "3", "A") & no_green == 0)

# same as above for September data
# combine Ev and Mv
all_dis_sep <- full_join(ev_dis_sep, mv_dis_sep) %>%
  mutate(infected = case_when(sp == "Ev" ~ Bp_spots_Ev,
                              sp == "Mv" ~ as.numeric(leaves_infec > 0)),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_")) %>%
  filter(ID %in% c("1", "2", "3", "A") & !is.na(leaves_tot))

# look at isolation data
unique(iso_jul$symptoms)
unique(iso_laug$symptoms)
unique(iso_sep$symptoms)

# remove leaves that couldn't be observed
# create a binary variables and site_plot
iso_jul2 <- iso_jul %>%
  filter(!is.na(symptoms) & !(symptoms %in% c("no sample", "not much leaf"))) %>%
  mutate(identify = dplyr::recode(gigantea, "Yes" = 1, "No" = 0),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_"))

iso_laug2 <- iso_laug %>%
  filter(!is.na(gigantea) & !is.na(symptoms) & !(symptoms %in% c("no sample", "not much leaf"))) %>%
  mutate(identify = dplyr::recode(gigantea, "Yes" = 1, "No" = 0),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_"))

iso_sep2 <- iso_sep %>%
  filter(!is.na(gigantea) & !is.na(symptoms) & !(symptoms %in% c("no sample", "not much leaf"))) %>%
  mutate(identify = dplyr::recode(gigantea, "Yes" = 1, "No" = 0),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_"))

# combine July Mv and late Aug Ev
iso_jul_aug <- iso_jul2 %>%
  filter(sp == "Mv") %>%
  full_join(iso_laug2 %>%
              filter(sp == "Ev"))

# combine Ev seed data across sampling months
# select focal plants (same sampling intensity as Mv)
# group by plot (same levels as Mv)
ev_seeds2 <- ev_seeds %>%
  filter(ID %in% c("1", "2", "3", "A")) %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(seeds = sum(seeds))

# combine seed data
# use total seeds for Mv
seeds <- mv_seeds %>%
  mutate(sp = "Mv",
         seeds = total_seeds) %>%
  select(-c(bag_seeds:total_seeds)) %>%
  full_join(ev_seeds2) %>%
  mutate(log_seeds = log10(seeds),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_"))

# combine seed data
# use bio seeds for Mv
seeds2 <- mv_seeds %>%
  mutate(sp = "Mv",
         seeds = bio_seeds) %>%
  select(-c(bag_seeds:total_seeds)) %>%
  full_join(ev_seeds2) %>%
  mutate(log_seeds = log10(seeds),
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1),
         site_plot = paste(site, plot, sep = "_"))
  

#### visualize ####

# Q1: Proportion infected August
q1plot <- ggplot(all_dis_aug2, aes(x = treatment, y = infected)) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_grid(~ sp)
q1plot

q1plot %+%
  facet_grid(site ~ sp)

q1plot %+%
  facet_grid(plot ~ sp)

# Q1: Proportion infected September
q1plot %+%
  all_dis_sep

q1plot %+%
  all_dis_sep %+%
  facet_grid(site ~ sp)

q1plot %+%
  all_dis_sep %+%
  facet_grid(plot ~ sp)

# Q2: Identified July
q1plot %+%
  iso_jul2 %+%
  aes(y = identify)

q1plot %+%
  iso_jul2 %+%
  aes(y = identify) %+%
  facet_grid(site ~ sp)

q1plot %+%
  iso_jul2 %+%
  aes(y = identify) %+%
  facet_grid(plot ~ sp)

# Q2: Identified late August
q1plot %+%
  iso_laug2 %+%
  aes(y = identify)
# only Ev

# Q2: Identified September
q1plot %+%
  iso_sep2 %+%
  aes(y = identify)

q1plot %+%
  iso_sep2 %+%
  aes(y = identify) %+%
  facet_grid(site ~ sp)
# missing some combinations

q1plot %+%
  iso_sep2 %+%
  aes(y = identify) %+%
  facet_grid(plot ~ sp)
# missing some combinations

# Q2: Identified July/August
q1plot %+%
  iso_jul_aug %+%
  aes(y = identify)
# these have the highest values across all the datasets

q1plot %+%
  iso_jul_aug %+%
  aes(y = identify) %+%
  facet_grid(site ~ sp)

q1plot %+%
  iso_jul_aug %+%
  aes(y = identify) %+%
  facet_grid(plot ~ sp)
# missing Ev in one treatment in one site

# Q3: seeds (Mv total seeds)
q1plot %+%
  seeds %+%
  aes(y = seeds)

q1plot %+%
  seeds %+%
  aes(y = seeds) +
  scale_y_log10()

q1plot %+%
  seeds %+%
  aes(y = seeds) %+%
  facet_grid(site ~ sp) +
  scale_y_log10()

q1plot %+%
  seeds %+%
  aes(y = seeds) %+%
  facet_grid(plot ~ sp) +
  scale_y_log10()

# Q3: seeds (Mv bio seeds)
q1plot %+%
  seeds2 %+%
  aes(y = seeds) +
  scale_y_log10()


#### stats ####

# Q1: Proportion infected August
# mod_q1_aug_1 <- glmer(infected ~ fungicide * sp  + (1|site/plot), 
#                     data = all_dis_aug2,
#                     family = binomial)
# singular fit, probably because there are 100% infections of Mv at some sites

mod_q1_aug_1 <- glmer(infected ~ fungicide * sp  + (1|site_plot), 
                    data = all_dis_aug2,
                    family = binomial)

mod_q1_aug_2 <- glmer(infected ~ fungicide * sp  + (1|plot), 
                    data = all_dis_aug2,
                    family = binomial)

# compare random effect structure
summary(mod_q1_aug_1) # SD = 0.73
summary(mod_q1_aug_2) # SD = 0.36

# test significance of interaction (fungicide has a unique effect on each species)
Anova(mod_q1_aug_1, type = 3)
# not significant

# remove interaction
mod_q1_aug_3 <- update(mod_q1_aug_1, .~. - fungicide:sp)
summary(mod_q1_aug_3)

# test significance of main effects
Anova(mod_q1_aug_3, type = 3)
# Mv has significantly higher Bipolaris infection
# fungicide marginally reduces the proportion of plants infected


# Q1: Proportion infected September

# skip first model tried above - won't converge with these data

# modify random effects
# mod_q1_sep_1 <- glmer(infected ~ fungicide * sp  + (1|site_plot), 
#                       data = all_dis_sep,
#                       family = binomial)
# failed to converge

# mod_q1_sep_1 <- glmer(infected ~ fungicide * sp  + (1|plot), 
#                       data = all_dis_sep,
#                       family = binomial)
# singular fit

mod_q1_sep_1 <- glm(infected ~ fungicide * sp, 
                      data = all_dis_sep,
                      family = binomial)
summary(mod_q1_sep_1)

# test significance of interaction (fungicide has a unique effect on each species)
Anova(mod_q1_sep_1, type = 3)
# significant
# Mv has has higher infection, fungicide reduces infection, but it reduces it more for Mv than Ev
# difference is not obvious in figure


# Q2: Proportion identified July
# mod_q2_jul_1 <- glmer(identify ~ fungicide * sp  + (1|site/plot),
#                     data = iso_jul2,
#                     family = binomial)
# failed to converge

mod_q2_jul_1 <- glmer(identify ~ fungicide * sp  + (1|site_plot),
                      data = iso_jul2,
                      family = binomial)

mod_q2_jul_2 <- glmer(identify ~ fungicide * sp  + (1|plot),
                      data = iso_jul2,
                      family = binomial)

mod_q2_jul_3 <- glmer(identify ~ fungicide * sp  + (1|site),
                      data = iso_jul2,
                      family = binomial)

# compare random effect structure
summary(mod_q2_jul_1) # SD = 0.50
summary(mod_q2_jul_2) # SD = 0
summary(mod_q2_jul_3) # SD = 0.80

# test significance of interaction (fungicide has a unique effect on each species)
Anova(mod_q2_jul_3, type = 3)
# not significant

# remove interaction
mod_q2_jul_4 <- update(mod_q2_jul_3, .~. - fungicide:sp)
summary(mod_q2_jul_4)

# test significance of main effects
Anova(mod_q2_jul_4, type = 3)
# Mv has significantly higher Bipolaris infection

# Q2: Proportion identified July/late August
mod_q2_jla_1 <- glmer(identify ~ fungicide * sp  + (1|site/plot),
                    data = iso_jul_aug,
                    family = binomial)

mod_q2_jla_2 <- glmer(identify ~ fungicide * sp  + (1|plot),
                      data = iso_jul_aug,
                      family = binomial)

mod_q2_jla_3 <- glmer(identify ~ fungicide * sp  + (1|site),
                      data = iso_jul_aug,
                      family = binomial)

# compare random effect structure
summary(mod_q2_jla_1) # SD = 0.13 and 0.70
summary(mod_q2_jla_2) # SD = 0.17
summary(mod_q2_jla_3) # SD = 0.70

# test significance of interaction (fungicide has a unique effect on each species)
# Anova(mod_q2_jla_1, type = 3)
# # not significant
# 
# # remove interaction
# mod_q2_jla_4 <- update(mod_q2_jla_1, .~. - fungicide:sp)
# # singular fit

# test significance of interaction for next best model
Anova(mod_q2_jla_3, type = 3)
# not significant

# remove interaction
mod_q2_jla_4 <- update(mod_q2_jla_3, .~. - fungicide:sp)
summary(mod_q2_jla_4)

# test significance of main effects
Anova(mod_q2_jla_4, type = 3)
# Mv has significantly higher Bipolaris infection
# fungicide significantly reduces it


# Q3: Seeds (Mv total seeds)

# mod_q3_tot_1 <- lmer(seeds ~ fungicide * sp  + (1|plot),
#                  data = seeds)
# singular fit

mod_q3_tot_1 <- lmer(seeds ~ fungicide * sp  + (1|site),
                 data = seeds)
summary(mod_q3_tot_1)
plot(mod_q3_tot_1)

mod_q3_tot_2 <- lmer(log_seeds ~ fungicide * sp  + (1|site),
                 data = seeds)
summary(mod_q3_tot_2)
plot(mod_q3_tot_2)
# use log-transformed

# test significance of interaction (fungicide has a unique effect on each species)
Anova(mod_q3_tot_2, type = 3)
# not significant

# remove interaction
mod_q3_tot_3 <- update(mod_q3_tot_2, .~. - fungicide:sp)
summary(mod_q3_tot_3)

# test significance of main effects
Anova(mod_q3_tot_3, type = 3)
# Mv has significantly higher seeds
# fungicide treatment has no significant effect


# Q3: Seeds (Mv bio seeds)

# mod_q3_bio_1 <- lmer(seeds ~ fungicide * sp  + (1|plot),
#                  data = seeds2)
# singular fit

mod_q3_bio_1 <- lmer(seeds ~ fungicide * sp  + (1|site),
                 data = seeds2)
summary(mod_q3_bio_1)
plot(mod_q3_bio_1)

mod_q3_bio_2 <- lmer(log_seeds ~ fungicide * sp  + (1|site),
                 data = seeds2)
summary(mod_q3_bio_2)
plot(mod_q3_bio_2)
# use log-transformed

# test significance of interaction (fungicide has a unique effect on each species)
Anova(mod_q3_bio_2, type = 3)
# not significant

# remove interaction
mod_q3_bio_3 <- update(mod_q3_bio_2, .~. - fungicide:sp)
summary(mod_q3_bio_3)

# test significance of main effects
Anova(mod_q3_bio_3, type = 3)
# Mv has significantly higher seeds
# fungicide treatment has no significant effect


#### figures ####

# base theme
base_fig <- ggplot(all_dis_sep, aes(x = treatment, y = infected, fill = Species)) +
  stat_summary(geom = "point", fun.y = "mean", size = 5, shape = 21, position = position_dodge(0.1)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.1)) +
  scale_fill_manual(values = c("white", "black")) +
  theme() +
  ylab(expression(paste("Proportion of plants with ", italic(Bipolaris), "-like lesions", sep = ""))) +
  xlab("Treatment") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10, color="black", face = "italic"),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, 0, -10, -10),
        legend.margin = margin(c(0, 1, 1, 1)),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside")

# Q1: infection by species
dat_fig_q1 <- all_dis_sep %>%
  mutate(Species = dplyr::recode(sp, Ev = "E. virginicus", Mv = "M. vimineum"))

fig_q1 <- base_fig %+%
  dat_fig_q1 +
  ylim(0, 1)
fig_q1

# Q2: identified infection by species
dat_fig_q2 <- iso_jul_aug %>%
  mutate(Species = dplyr::recode(sp, Ev = "E. virginicus", Mv = "M. vimineum"))

fig_q2 <- base_fig %+%
 dat_fig_q2 %+%
  aes(y = identify) %+%
  ylab(expression(paste("Proportion of plants with ", italic("Bipolaris gigantea"), sep = ""))) +
  ylim(0, 1)
fig_q2  

# Q3: seeds
dat_fig_q3 <- seeds %>%
  mutate(Species = dplyr::recode(sp, Ev = "E. virginicus", Mv = "M. vimineum"))

fig_q3 <- base_fig %+%
  dat_fig_q3 %+%
  aes(y = log_seeds) %+%
  ylab(expression(paste(log[10], "(seeds)", sep = "")))
fig_q3


#### output ####

# figures
pdf("./output/Lane_APS_2020_field_data_analysis.pdf", width = 4, height = 4)
fig_q1
fig_q2
fig_q3
dev.off()

# data
write_csv(dat_fig_q1, "./intermediate-data/Lane_APS_2020_field_data_analysis/bipolaris_lesions_2018_density_exp.csv")
write_csv(dat_fig_q2, "./intermediate-data/Lane_APS_2020_field_data_analysis/bipolaris_identification_2018_density_exp.csv")
write_csv(dat_fig_q3, "./intermediate-data/Lane_APS_2020_field_data_analysis/seeds_2018_density_exp.csv")
