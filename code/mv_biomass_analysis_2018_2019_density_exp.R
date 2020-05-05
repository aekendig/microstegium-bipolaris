##### info ####

# file: mv_biomass_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 5/4/20
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


#### check datasets ####

# sample sizes
samps_19 <- bio19b %>%
  group_by(site, plot, treatment) %>%
  count() %>%
  mutate(bio.g = max(bio19b$bio.g, na.rm = T) + 1)

samps_18 <- bio18b %>%
  group_by(site, plot, treatment) %>%
  count() %>%
  mutate(bio.g = max(bio18b$bio.g, na.rm = T) + 1)

# visualize
samp_plot_19 <- ggplot(bio19b, aes(x = plot, y = bio.g)) +
  geom_point(alpha = 0.5) +
  facet_grid(treatment ~ site)

samp_plot_19 +
  geom_text(data = samps_19, aes(label = n), size = 2)
  
samp_plot_19 %+%
  bio18b +
  geom_text(data = samps_18, aes(label = n), size = 2)


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

dev.off()


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


#### separate data ####

# treatments
bio_dat_water_19 <- filter(bio19b, treatment == "water")
bio_dat_fungicide_19 <- filter(bio19b, treatment == "fungicide")


### manual model fitting ###

# sigmoid Beverton-Holt
sbh_fun <- function(dat_in, r, a, d){
  
  # extract values
  xmin = min(dat_in$background_density)
  xmax = max(dat_in$background_density)
  y0 = filter(dat_in, fungicide == 0) %>%
    summarise(mean_bio = mean(bio.g)) %>%
    as.numeric() %>%
    round()
  print(y0)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 + x^(d-1) / (1/r + a * x^d))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = bio.g)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# separate data
bio19b_mv <- filter(bio19b, background == "Mv seedling")

# try values
sbh_fun(bio19b_mv, 2, 0.05, 2)


### Mv background 2019 ###

mv_bio_mv_bg_2019_mod <- brm(data = bio19b_mv, family = gaussian,
                          bf(bio.g ~ s0 + background_density / (k + alpha * background_density^2),
                             s0 ~ fungicide + (1|site),
                             k ~ fungicide,
                             alpha ~ fungicide,
                             nl = T),
                          prior <- c(prior(normal(13, 10), nlpar = "s0"),
                                     prior(normal(0.5, 1), nlpar = "k"),
                                     prior(normal(0.05, 1), nlpar = "alpha")),
                          iter = 6000, warmup = 1000, chains = 3, cores = 2,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15))
# 325 divergent transitions





