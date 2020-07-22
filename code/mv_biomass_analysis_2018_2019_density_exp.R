##### info ####

# file: mv_biomass_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/10/20
# goal: evaluate the effects of density treatments and environmental covariates on the biomass of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(glmmTMB)

# import data
bio18 <- read_csv("./data/mv_biomass_oct_2018_density_exp.csv")
bio19 <- read_csv("./data/mv_biomass_seeds_2019_density_exp.csv")
plots <- read_csv("./data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("./intermediate-data/covariates_2018_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
sev18 <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
sev19 <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")


#### edit data ####

# notes
unique(bio18$processing_notes)
unique(bio19$process_notes)

# disease severity
sev18b <- sev18 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Mv") %>%
  group_by(site, plot, treatment) %>%
  summarise(severity_jul = mean(severity_jul, na.rm = T),
            severity_late_aug = mean(severity_late_aug, na.rm = T),
            severity_sep = mean(severity_sep, na.rm = T))

sev19b <- sev19 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Mv") %>%
  mutate(plant = as.numeric(ID)) %>%
  select(-ID)

# 2018
bio18b <- bio18 %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment),
         log_bio.g = log(bio.g),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide")) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  left_join(sev18b)

# 2019
bio19b <- bio19 %>%
  rename(bio.g = biomass_weight.g) %>%
  mutate(log_bio.g = log(bio.g),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide")) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  left_join(sev19b)

# 2019 without plot repetition
bio19t <- bio19 %>%
  rename(bio.g = biomass_weight.g) %>%
  mutate(log_bio.g = log(bio.g),
         fungicide = recode(treatment, water = 0, fungicide = 1)) %>%
  left_join(treat) %>%
  left_join(covar) %>%
  left_join(sev19b) %>%
  mutate(background = fct_relevel(background, "none"),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide"))

# no background data (plot-level to use both years)
filter(bio19t, plot == 1 & is.na(bio.g))

no_back_plot_dat <- bio18b %>%
  filter(plot == 1 & background == "Mv seedling") %>%
  select(site, plot, treatment, Treatment, bio.g, severity_jul:severity_sep) %>%
  mutate(yearf = "year 1") %>%
  full_join(filter(bio19t, plot == 1) %>%
              group_by(site, plot, treatment, Treatment) %>%
              summarise(bio.g = sum(bio.g),
                        severity_jul = mean(severity_jul, na.rm = T),
                        severity_late_aug = mean(severity_late_aug, na.rm = T),
                        severity_sep = mean(severity_sep, na.rm = T)) %>%
              mutate(yearf = "year 2")) %>%
  mutate(adj_bio.g = case_when(treatment == "fungicide" ~ bio.g / 0.88,
                               TRUE ~ bio.g))

no_back_plant_dat <- filter(bio19t, plot == 1) %>%
  mutate(adj_bio.g = case_when(treatment == "fungicide" ~ bio.g / 0.88,
                               TRUE ~ bio.g))


#### check datasets ####

# sample sizes
samps_19 <- bio19b %>%
  group_by(site, plot, treatment) %>%
  summarise(n = sum(!is.na(bio.g))) %>%
  mutate(bio.g = max(bio19b$bio.g, na.rm = T) + 1)

samps_18 <- bio18b %>%
  group_by(site, plot, treatment) %>%
  summarise(n = sum(!is.na(bio.g))) %>%
  mutate(bio.g = max(bio18b$bio.g, na.rm = T) + 1)

# visualize
samp_plot_19 <- ggplot(bio19b, aes(x = plot, y = bio.g)) +
  geom_point(alpha = 0.5) +
  facet_grid(treatment ~ site)

samp_plot_19 +
  geom_text(data = samps_19, aes(label = n), size = 2)
# 2 plants missing
  
samp_plot_19 %+%
  bio18b +
  geom_text(data = samps_18, aes(label = n), size = 2)
# 2 plots missing (no data collected from these)


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

# settings
col_pal = c("black", "#0072B2", "#56B4E9", "#D55E00")
col_pal2 = c("#a6611a", "#018571")
dodge_value = 1

# template figure
temp_fig <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

# save treatment figures
pdf("./output/mv_biomass_visualize_treatment_2018_2019_density_exp.pdf")

# 2018
plot_grid(temp_fig %+%
            bio18b %+%
            aes(y =bio.g) +
            ylab("Year 1 biomass (g per plot)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig %+%
            bio18b %+%
            aes(y = log_bio.g) +
            ylab("Year 1 log-biomass (g per plot)"),
          nrow = 2)
# biomass saturates with Mv seedling density like plot-level Mv biomass collected in 2019

# 2019
plot_grid(temp_fig %+%
            bio19b %+%
            aes(y =bio.g) +
            ylab("Year 2 biomass (g per plant)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig %+%
            bio19b %+%
            aes(y = log_bio.g) +
            ylab("Year 2 log-biomass (g per plant)"),
          nrow = 2)
# biomass decreases with increasing density, but not when moving from plot with no background to some background (increases)

dev.off()

# group all background types together
ggplot(bio19t, aes(x = background_density, y = bio.g, color = background)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5, position = position_dodge(dodge_value)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(dodge_value)) +
  facet_wrap(~ Treatment) +
  scale_color_manual(values = col_pal, name = "Neighbors") +
  temp_theme +
  xlab("Neighbor density") +
  ylab(expression(paste(italic(Microstegium), " biomass (g ", individual^-1, ")", sep = ""))) +
  theme(strip.text = element_text(size = 12))

# site-specific trends
ggplot(bio19b, aes(x = background_density, y = bio.g, color = site)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "line", fun = "mean", position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3)) +
  facet_grid(treatment ~ background, scales = "free") +
  temp_theme

ggplot(bio19b, aes(x = background_density, y = bio.g, color = site)) +
  geom_point() +
  facet_grid(treatment ~ background, scales = "free") +
  temp_theme

# no background
pdf("output/mv_biomass_no_background_treatment_2018_2019_density_exp.pdf", width = 3, height = 3)
ggplot(no_back_plot_dat, aes(x = Treatment, y = bio.g, fill = Treatment)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  scale_fill_manual(values = col_pal2, guide = F) +
  ylab("Plot-level biomass (g)") +
  ggtitle("Mv seedling") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

ggplot(no_back_plot_dat, aes(x = Treatment, y = adj_bio.g, fill = Treatment)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  scale_fill_manual(values = col_pal2, guide = F) +
  ylab("Adjusted plot-level biomass (g)") +
  ggtitle("Mv seedling") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

ggplot(no_back_plot_dat, aes(x = Treatment, y = adj_bio.g, fill = Treatment)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_wrap(~yearf) +
  scale_fill_manual(values = col_pal2, guide = F) +
  ylab("Adjusted plot-level biomass (g)") +
  ggtitle("Mv seedling") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

ggplot(no_back_plant_dat, aes(x = Treatment, y = bio.g, fill = Treatment)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  scale_fill_manual(values = col_pal2, guide = F) +
  ylab("Plant-level biomass (g)") +
  ggtitle("Mv seedling") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

ggplot(no_back_plant_dat, aes(x = Treatment, y = adj_bio.g, fill = Treatment)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  scale_fill_manual(values = col_pal2, guide = F) +
  ylab("Adjusted plant-level biomass (g)") +
  ggtitle("Mv seedling") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))
dev.off()


#### visualize severity effects ####

plot_sev_jul_viz <- ggplot(no_back_plot_dat, aes(x = severity_jul, y = bio.g, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~yearf, scales = "free") +
  temp_theme
plot_sev_jul_viz
# not enough points

plot_sev_jul_viz %+%
  aes(x = severity_late_aug)
# not enough points

plot_sev_jul_viz %+%
  aes(x = severity_sep)
# not enough points

plant_sev_jul_viz <- ggplot(no_back_plant_dat, aes(x = severity_jul, y = bio.g, color = Treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  temp_theme
plant_sev_jul_viz

plant_sev_jul_viz %+%
  aes(x = severity_early_aug)

plant_sev_jul_viz %+%
  aes(x = severity_late_aug)
# decent number of points and spread
  

#### visualize covariate effects ####

# covariate template figure
temp_fig_cov <- ggplot(bio18b, aes(x = soil_moisture_jun.prop, y = bio.g)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "loess") +
  facet_wrap(~ background, scales = "free_x") +
  temp_theme

# save covariate figures
pdf("./output/mv_biomass_visualize_covariates_2018_2019_density_exp.pdf")

# June soil moisture 
plot_grid(temp_fig_cov +
            ylab("Year 1 biomass (g per plot)") +
            theme(axis.title.x = element_blank()),
          temp_fig_cov %+%
            bio19b +
            xlab("June soil moisture") +
            ylab("Year 2 biomass (g per plant)") +
            theme(strip.text.x = element_blank()),
          nrow = 2)

# October soil moisture
plot_grid(temp_fig_cov %+%
            aes(x = soil_moisture_oct.prop) +
            ylab("Year 1 biomass (g per plot)") +
            theme(axis.title.x = element_blank()),
          temp_fig_cov %+%
            bio19b %+%
            aes(x = soil_moisture_oct.prop) +
            xlab("October soil moisture") +
            ylab("Year 2 biomass (g per plant)") +
            theme(strip.text.x = element_blank()),
          nrow = 2)

# canopy cover
plot_grid(temp_fig_cov %+%
            aes(x = canopy_cover.prop) +
            ylab("Year 1 biomass (g per plot)") +
            theme(axis.title.x = element_blank()),
          temp_fig_cov %+%
            bio19b %+%
            aes(x = canopy_cover.prop) +
            xlab("Canopy cover") +
            ylab("Year 2 biomass (g per plant)") +
            theme(strip.text.x = element_blank()),
          nrow = 2)

# September biomass
plot_grid(temp_fig_cov %+%
            aes(x = mv_sep.g) +
            ylab("Year 1 biomass (g per plot)") +
            theme(axis.title.x = element_blank()),
          temp_fig_cov %+%
            bio19b %+%
            aes(x = mv_sep.g) +
            xlab("September plot-edge biomass") +
            ylab("Year 2 biomass (g per plant)") +
            theme(strip.text.x = element_blank()),
          nrow = 2)

dev.off()


#### direct disease models ####

# control biomass
filter(no_back_plant_dat, fungicide == 0) %>%
  summarise(bio = mean(bio.g))

# model
mv_no_back_bio_mod <- brm(data = no_back_plant_dat, family = gaussian,
                          bio.g ~ fungicide + (1|site),
                          prior <- c(prior(normal(13.4, 10), class = Intercept),
                                     prior(normal(0, 10), class = b),
                                     prior(cauchy(0, 1), class = sd),
                                     prior(cauchy(0, 1), class = sigma)),
                          iter = 6000, warmup = 1000, chains = 1,
                          control = list(adapt_delta = 0.99))

# check model and add chains
summary(mv_no_back_bio_mod)
prior_summary(mv_no_back_bio_mod)
mv_no_back_bio_mod <- update(mv_no_back_bio_mod, chains = 3)
plot(mv_no_back_bio_mod)
pp_check(mv_no_back_bio_mod, nsamples = 100)

# reduction due to fungicide
1-(21.05/23.80)
-0.12*13.4

# model with direct fungicide effects
mv_no_back_bio_fung_mod <- brm(data = no_back_plant_dat, family = gaussian,
                               bio.g ~ Treatment + fungicide + (1|site),
                               prior <- c(prior(normal(13.4, 10), class = Intercept),
                                          prior(normal(0, 10), class = b),
                                          prior(normal(-1.6, 0.001), class = b, coef = "Treatmentfungicide"),
                                          prior(cauchy(0, 1), class = sd),
                                          prior(cauchy(0, 1), class = sigma)),
                          iter = 6000, warmup = 1000, chains = 3,
                          control = list(adapt_delta = 0.99))

# check model
summary(mv_no_back_bio_fung_mod)
prior_summary(mv_no_back_bio_fung_mod)
plot(mv_no_back_bio_fung_mod)
pp_check(mv_no_back_bio_fung_mod, nsamples = 100)


#### disease and density models ####

# separate data
mv_mv_bio_dat <- bio19t %>%
  filter(background == "Mv seedling")

mv_evs_bio_dat <- bio19t %>%
  filter(background == "Ev seedling" & !is.na(bio.g))

mv_eva_bio_dat <- bio19t %>%
  filter(background == "Ev adult" & !is.na(bio.g))

# check for complete values
sum(is.na(mv_mv_bio_dat$bio.g))
sum(is.na(mv_evs_bio_dat$bio.g))
sum(is.na(mv_eva_bio_dat$bio.g))

# Beverton-Holt function
bh_fun <- function(dat_in, a, y0){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$background_density)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = bio.g)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# try values
mean(filter(mv_mv_bio_dat, density_level == "low")$bio.g)
bh_fun(mv_mv_bio_dat, 0.03, 20)
mean(mv_evs_bio_dat$bio.g)
bh_fun(mv_evs_bio_dat, 0, 17)
mean(filter(mv_eva_bio_dat, density_level == "low")$bio.g)
bh_fun(mv_eva_bio_dat, 0.1, 30)

# distributions
x <- seq(0, 10, length.out = 100)
y <- dexp(x, 1/0.1)
plot(x, y, type = "l")

# mv background model
mv_mv_bio_mod <- brm(data = mv_mv_bio_dat, family = gaussian,
                     bf(bio.g ~ y0 / (1 + alpha * background_density),
                        y0 ~ 0 + Treatment + (1|site),
                        alpha ~ 0 + Treatment,
                        nl = T),
                     prior <- c(prior(normal(22, 10), nlpar = "y0", class = "b", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                prior(cauchy(0, 1), class = "sigma")),
                     iter = 6000, warmup = 1000, chains = 1,
                     control = list(adapt_delta = 0.999))
                             
# check model and increase chains
summary(mv_mv_bio_mod)
prior_summary(mv_mv_bio_mod)
mv_mv_bio_mod <- update(mv_mv_bio_mod, chains = 3)
plot(mv_mv_bio_mod)
pp_check(mv_mv_bio_mod, nsamples = 100)

# ev seedling background model
mv_evs_bio_mod <- brm(data = mv_evs_bio_dat, family = gaussian,
                      bf(bio.g ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(22, 10), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(mv_evs_bio_mod)
mv_evs_bio_mod <- update(mv_evs_bio_mod, chains = 3)
plot(mv_evs_bio_mod)
pp_check(mv_evs_bio_mod, nsamples = 100)

# ev adult background model
mv_eva_bio_mod <- brm(data = mv_eva_bio_dat, family = gaussian,
                      bf(bio.g ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(30, 10), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(mv_eva_bio_mod)
mv_eva_bio_mod <- update(mv_eva_bio_mod, chains = 3)
plot(mv_eva_bio_mod)
pp_check(mv_eva_bio_mod, nsamples = 100)

# all data combined
mv_bio_dat <- bio19t %>%
  filter(!is.na(bio.g) & background != "none") %>%
  mutate(back_trt = paste(background, treatment, sep = "_"))
unique(mv_bio_dat$back_trt)

# initial value
filter(mv_bio_dat, background == "Ev adult" & density_level == "low" & treatment == "fungicide") %>%
  summarise(bio = mean(bio.g))

# all background model
mv_bio_mod <- brm(data = mv_bio_dat, family = gaussian,
                     bf(bio.g ~ y0 / (1 + alpha * background_density),
                        y0 ~ 0 + treatment + (1|site),
                        alpha ~ 0 + back_trt,
                        nl = T),
                     prior <- c(prior(normal(40, 10), nlpar = "y0", class = "b", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                prior(cauchy(0, 1), class = "sigma")),
                     iter = 6000, warmup = 1000, chains = 1,
                     control = list(adapt_delta = 0.999))

# check model and increase chains
summary(mv_bio_mod)
prior_summary(mv_bio_mod)
mv_bio_mod <- update(mv_bio_mod, chains = 3)
plot(mv_bio_mod)
pp_check(mv_bio_mod, nsamples = 100)


#### model fit figures ####

# model fits over simulated data
bio_mv_sim_dat <- tibble(background_density = rep(seq(0, 64, length.out = 300), 2),
                          fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(bio.g = fitted(mv_mv_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
         bio_lower = fitted(mv_mv_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         bio_upper = fitted(mv_mv_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"])

bio_evs_sim_dat <- tibble(background_density = rep(seq(0, 16, length.out = 100), 2),
                           fungicide = rep(c(0, 1), each = 100)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(bio.g = fitted(mv_evs_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
         bio_lower = fitted(mv_evs_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         bio_upper = fitted(mv_evs_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"])

bio_eva_sim_dat <- tibble(background_density = rep(seq(0, 8, length.out = 100), 2),
                           fungicide = rep(c(0, 1), each = 100)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(bio.g = fitted(mv_eva_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
         bio_lower = fitted(mv_eva_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         bio_upper = fitted(mv_eva_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"])

# combine simulated data
bio_sim_dat <- bio_mv_sim_dat %>%
  mutate(background = "Mv seedling") %>%
  full_join(bio_evs_sim_dat %>%
              mutate(background = "Ev seedling")) %>%
  full_join(bio_eva_sim_dat %>%
              mutate(background = "Ev adult"))

# figure
pdf("output/mv_biomass_analysis_model_fits_2019_density_exp.pdf")
ggplot(bio19b, aes(x = background_density, y = bio.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = Treatment)) +
  geom_ribbon(data = bio_sim_dat, alpha = 0.5, aes(ymin = bio_lower, ymax = bio_upper, fill = Treatment)) +
  geom_line(data = bio_sim_dat, aes(color = Treatment)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  xlab("Background density") +
  ylab(expression(paste(italic(Microstegium), " biomass (g ", plant^-1, ")", sep = ""))) +
  scale_color_manual(values = col_pal2) +
  scale_fill_manual(values = col_pal2) +
  temp_theme
dev.off()

# full data model
bio_all_sim_dat <- tibble(background_density = c(rep(seq(0, 64, length.out = 100), 2),
                                                 rep(seq(0, 16, length.out = 100), 2),
                                                 rep(seq(0, 8, length.out = 100), 2)),
                          fungicide = rep(rep(c(0, 1), each = 100), 3),
                          background = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 200)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"),
         treatment = recode(Treatment, "control (water)" = "water"),
         back_trt = paste(background, treatment, sep = "_")) %>%
  mutate(bio.g = fitted(mv_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
         bio_lower = fitted(mv_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         bio_upper = fitted(mv_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"])

# figure
pdf("output/mv_biomass_analysis_combined_model_fit_2019_density_exp.pdf")
ggplot(bio19b, aes(x = background_density, y = bio.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = Treatment)) +
  geom_ribbon(data = bio_all_sim_dat, alpha = 0.5, aes(ymin = bio_lower, ymax = bio_upper, fill = Treatment)) +
  geom_line(data = bio_all_sim_dat, aes(color = Treatment)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  xlab("Background density") +
  ylab(expression(paste(italic(Microstegium), " biomass (g ", plant^-1, ")", sep = ""))) +
  scale_color_manual(values = col_pal2) +
  scale_fill_manual(values = col_pal2) +
  temp_theme
dev.off()


#### severity and biomass ####

# water treatment
mv_sev_dat <- filter(bio19t, treatment == "water") %>%
  as.data.frame()

# check for density-severity correlations
for(i in 24:28){
  mv_sev_dat$severity = mv_sev_dat[,i]
  print(cor.test(mv_sev_dat$background_density, mv_sev_dat$severity))
}
# none significantly correlated

# biomass-severity relationships
mv_bio_sev_jun_19_mod <- glmmTMB(bio.g ~ severity_jun + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_bio_sev_jun_19_mod) # not sig
mv_bio_sev_jul_19_mod <- glmmTMB(bio.g ~ severity_jul + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_bio_sev_jul_19_mod) # not sig
mv_bio_sev_eau_19_mod <- glmmTMB(bio.g ~ severity_early_aug + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_bio_sev_eau_19_mod) # not sig
mv_bio_sev_lau_19_mod <- glmmTMB(bio.g ~ severity_late_aug + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_bio_sev_lau_19_mod) # not sig

#### output ####
save(mv_no_back_bio_mod, file = "output/mv_biomass_no_background_model_2019_density_exp.rda")
save(mv_no_back_bio_fung_mod, file = "output/mv_biomass_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
save(mv_mv_bio_mod, file = "output/mv_biomass_mv_background_model_2019_density_exp.rda")
save(mv_evs_bio_mod, file = "output/mv_biomass_ev_seedling_background_model_2019_density_exp.rda")
save(mv_eva_bio_mod, file = "output/mv_biomass_ev_adult_background_model_2019_density_exp.rda")
save(mv_bio_mod, file = "output/mv_biomass_combined_background_model_2019_density_exp.rda")
