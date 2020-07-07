##### info ####

# file: mv_biomass_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/6/20
# goal: evaluate the effects of density treatments and environmental covariates on the biomass of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)

# import data
bio18 <- read_csv("./data/mv_biomass_oct_2018_density_exp.csv")
bio19 <- read_csv("./data/mv_biomass_seeds_2019_density_exp.csv")
plots <- read_csv("./data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("./intermediate-data/covariates_2018_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# notes
unique(bio18$processing_notes)
unique(bio19$process_notes)

# 2018
bio18b <- bio18 %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment),
         log_bio.g = log(bio.g),
         fungicide = recode(treatment, water = 0, fungicide = 1)) %>%
  left_join(plots) %>%
  left_join(covar)

# 2019
bio19b <- bio19 %>%
  rename(bio.g = biomass_weight.g) %>%
  mutate(log_bio.g = log(bio.g),
         fungicide = recode(treatment, water = 0, fungicide = 1)) %>%
  left_join(plots) %>%
  left_join(covar)

# 2019 without plot repetition
bio19t <- bio19 %>%
  rename(bio.g = biomass_weight.g) %>%
  mutate(log_bio.g = log(bio.g),
         fungicide = recode(treatment, water = 0, fungicide = 1)) %>%
  left_join(treat) %>%
  left_join(covar) %>%
  mutate(background = fct_relevel(background, "none"),
         Treatment = recode(treatment, water = "Control (water)", fungicide = "Fungicide"))


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
            ylab("Year d biomass (g per plant)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig %+%
            bio19b %+%
            aes(y = log_bio.g) +
            ylab("Year d log-biomass (g per plant)"),
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


#### fungicide effect ####

# add fungicide column
expt2 <- expt %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0"))

# model
fung_mod <- brm(data = expt2, family = gaussian,
                weight.g ~ fungicide,
                prior <- c(prior(normal(0, 100), class = Intercept),
                           prior(normal(0, 10), class = b),
                           prior(cauchy(0, 1), class = sigma)),
                iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(fung_mod)
plot(fung_mod)

# remove 0 plots
bio19m <- filter(bio19t, density_level != "none")

# translate to field experiment values
bio19m %>%
  filter(!is.na(bio.g) & density_level =="low") %>%
  ggplot(aes(x = bio.g)) +
  geom_histogram() +
  geom_vline(aes(xintercept = mean(bio.g)), color = "blue", linetype = "dashed")

mean(filter(bio19m, !is.na(bio.g) & density_level =="low")$bio.g)  


### manual model fitting ###

# Beverton-Holt
bh_fun <- function(dat_in, a){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$background_density)
  y0 = mean(filter(bio19m, !is.na(bio.g) & density_level =="low")$bio.g)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = bio.g)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# separate background types
bio19m_mv <- filter(bio19m, background == "Mv seedling")

# try values
bh_fun(filter(bio19m_mv, fungicide == 0), 0.03)

# distributions
x <- seq(0, 10, length.out = 100)
y <- dexp(x, 1)
plot(x, y, type = "l")


### Mv background 2019 ###

mv_bio_mv_bg_2019_mod <- brm(data = bio19m_mv, family = gaussian,
                             bf(bio.g ~ (y0 + f)/ (alpha * background_density),
                                y0 ~ fungicide + (1|site),
                                f ~ fungicide,
                                alpha ~ fungicide,
                                nl = T),
                             prior <- c(prior(normal(23, 10), coef = "Intercept", nlpar = "y0"),
                                        prior(normal(0, 10), coef = "fungicide", nlpar = "y0"),
                                        prior(normal(0, 0.001), coef = "Intercept", nlpar = "f"),
                                        prior(normal(-2, 1), coef = "fungicide", nlpar = "f"),
                                        prior(normal(0, 1), coef = "fungicide", nlpar = "alpha"),
                                        prior(exponential(1), nlpar = "alpha")),
                             iter = 6000, warmup = 1000, chains = 1, cores = 1)
summary(mv_bio_mv_bg_2019_mod)
plot(mv_bio_mv_bg_2019_mod)
# 143 divergent transitions
mv_bio_mv_bg_2019_mod <- update(mv_bio_mv_bg_2019_mod, control = list(adapt_delta = 0.99))
# 5 divergent transitions, 29 exceed max treedepth
mv_bio_mv_bg_2019_mod <- update(mv_bio_mv_bg_2019_mod, control = list(adapt_delta = 0.999, max_treedepth = 15))
# 5 divergent transitions
mv_bio_mv_bg_2019_mod <- update(mv_bio_mv_bg_2019_mod, control = list(adapt_delta = 0.9999, max_treedepth = 15))
# 1 divergent transition
mv_bio_mv_bg_2019_mod <- update(mv_bio_mv_bg_2019_mod, control = list(adapt_delta = 0.99999, max_treedepth = 15))
# 8 divergent transitions
mv_bio_mv_bg_2019_mod <- update(mv_bio_mv_bg_2019_mod, control = list(adapt_delta = 0.99999999, max_treedepth = 15))
# 3 divergent transitions