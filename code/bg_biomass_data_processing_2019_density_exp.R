##### info ####

# file: bg_biomass_data_processing_2019_density_exp
# author: Amy Kendig
# date last edited: 5/20/20
# goal: check background biomass and relationship with covariates


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
bio <- read_csv("data/bg_biomass_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")


#### edit data ####

# combine plots
bio2 <- bio %>%
  group_by(site, plot, treatment) %>%
  summarise(biomass.g = sum(biomass.g)) %>%
  ungroup()

# look at missing data
bio2 %>%
  group_by(site) %>%
  count()
# none are missing -- there are no 1's

# add 0 data to 1 plots
bio3 <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 2),
               treatment = rep(c("water", "fungicide"), 4)) %>%
  mutate(plot = 1,
         biomass.g = 0) %>%
  full_join(bio2)

# combine with plots and covariates
bio4 <- bio3 %>%
  full_join(plots) %>%
  full_join(covar) %>%
  mutate(density_level = fct_relevel(density_level, "none", "low", "medium"))

sum(is.na(bio3$biomass.g))


#### visualize treatment effects ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside")

# save treatment figures
pdf("./output/bg_biomass_visualize_treatment_2019_density_exp.pdf")

# biomass by density
ggplot(bio4, aes(x = background_density, y = biomass.g, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_wrap(~background, scales = "free") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density") +
  ylab("Background biomass (g)")

ggplot(bio4, aes(x = background_density, y = biomass.g/background_density, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_wrap(~background, scales = "free") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density") +
  ylab("Per capita background biomass (g)")

# biomass by site
ggplot(bio4, aes(x = site, y = biomass.g, shape = density_level, fill = treatment)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  facet_wrap(~ background, scales = "free") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  scale_shape_manual(values = 21:24, name = "Density level") +
  temp_theme +
  xlab("Background density") +
  ylab("Background biomass (g)")

# per capita biomass by site
ggplot(bio4, aes(x = site, y = biomass.g / background_density, shape = density_level, fill = treatment)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  facet_wrap(~ background, scales = "free") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  scale_shape_manual(values = 21:24, name = "Density level") +
  temp_theme +
  xlab("Background density") +
  ylab("Per capita background biomass (g)")

dev.off()


#### visualize covariate effects ####

# covariate template figure
temp_fig_cov <- ggplot(bio4, aes(x = soil_moisture_jun.prop, y = biomass.g)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "loess") +
  facet_wrap(~ background, scales = "free") +
  temp_theme 

# save covariate figures
pdf("./output/bg_biomass_visualize_covariates_2019_density_exp.pdf")

# June soil moisture
plot_grid(temp_fig_cov +
            ylab("Background biomass (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(y = biomass.g / background_density) +
            xlab("June soil moisture") +
            ylab("Per capita background biomass (g)"),
          nrow = 2)

# October soil moisture
plot_grid(temp_fig_cov %+%
            aes(x = soil_moisture_oct.prop) +
            ylab("Background biomass (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(x = soil_moisture_oct.prop, y = biomass.g / background_density) +
            xlab("October soil moisture") +
            ylab("Per capita background biomass (g)"),
          nrow = 2)

# canopy cover
plot_grid(temp_fig_cov %+%
            aes(x = canopy_cover.prop) +
            ylab("Background biomass (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(x = canopy_cover.prop, y = biomass.g / background_density) +
            xlab("Canopy cover") +
            ylab("Per capita background biomass (g)"),
          nrow = 2)

# September biomass
plot_grid(temp_fig_cov %+%
            aes(x = mv_sep.g) +
            ylab("Background biomass (g)") +
            theme(axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig_cov %+%
            aes(x = mv_sep.g, y = biomass.g / background_density) +
            xlab("September biomass") +
            ylab("Per capita background biomass (g)"),
          nrow = 2)

dev.off()


#### output ####

# save data without plots and covariates
write_csv(bio3, "intermediate-data/bg_processed_biomass_2019_density_exp.csv")
