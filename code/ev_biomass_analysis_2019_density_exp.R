##### info ####

# file: ev_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 5/6/20
# goal: evaluate the effects of density treatments and environmental covariates on the seed production of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)

# import data
bio <- read_csv("./data/ev_biomass_seeds_oct_2019_density_exp.csv")
spike <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")


#### edit data ####

# notes
unique(bio$processing_notes)

# format spikelet data
spike2 <- spike %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g))

# merge with plots and covariates
bio2 <- bio %>%
  filter(!is.na(weight)) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1)) %>%
  rename(veg_weight.g = weight) %>%
  left_join(spike2) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         tot_weight.g = veg_weight.g + spikelet_weight.g) %>%
  left_join(plots) %>%
  left_join(covar)


#### check data ####

# max value
max(bio2$veg_weight.g)

# sample sizes
samps <- bio2 %>%
  group_by(site, plot, treatment) %>%
  count() %>%
  mutate(veg_weight.g = 23)

# visualize
ggplot(bio2, aes(x = plot, y = veg_weight.g)) +
  geom_point(alpha = 0.5, aes(color = age)) +
  facet_grid(treatment ~ site) +
  geom_text(data = samps, aes(label = n), size = 2)


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

# template figure
temp_fig <- ggplot(plots, aes(x = background_density, y = plot, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

# save treatment figures
pdf("./output/ev_biomass_visualize_treatment_2019_density_exp.pdf")

plot_grid(temp_fig %+%
            bio2 %+%
            aes(y = veg_weight.g) +
            ylab("Vegetative weight (g)") +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  strip.text.x = element_blank()),
          temp_fig %+%
            bio2 %+%
            aes(y = tot_weight.g) +
            ylab("Total weight (g)"),
          nrow = 2,
          rel_heights = c(0.8, 1))

dev.off()

