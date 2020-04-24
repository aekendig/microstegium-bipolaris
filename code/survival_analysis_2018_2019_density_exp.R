##### info ####

# file: survival_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/22/20
# goal: evaluate the effects of density and fungicide on survival


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
daty1_0 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
daty2_0 <- read_csv("data/all_replacement_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")


#### edit data ####

# add plant type
daty1_1 <- daty1_0 %>%
  mutate(plant_type = paste(sp, age, sep = " "))

daty2_1 <- daty2_0 %>%
  mutate(plant_type = paste(sp, age, sep = " "))

# summer survival 2018
# remove NA's 
# remove background Microstegium (not individual plants)
sumy1 <- daty1_1 %>%
  filter(month == "September" & !is.na(survival) & !(sp == "Mv" & focal == 0)) %>%
  left_join(plots) %>%
  left_join(covar)

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
  left_join(covar)

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
  
# summer survival 2018
sumy2 <- daty2_2 %>%
  select(-replace_date) %>%
  unique() %>%
  mutate(survival = 0) %>%
  full_join(allid) %>%
  mutate(survival = replace_na(survival, 1)) %>%
  left_join(covar)


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
  stat_summary(geom = "point", fun.y = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free_x", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme +
  xlab("Background density")

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

dev.off()


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


#### model decisions ####

# linearly increasing
# saturating increase
# low density peak
# hump-shaped
# no change
# intercept ranges from 0.25 to 1

# literature search
# Beverton-Holt (fish, Shepherd and Cushing 1980)
# linear (fish, Forrester 1995; plants, Schamp and Aarssen 1991)
# ANOVA (plants, Bell et al. 2006)
# logistic glm (plants, Sletvold 2005)
# linear regression of ln(-ln(1-p)) (trees, Zhu et al. 2015)

# soil moisture June increases year 1 summer survival
# soil moisture Oct does too, but less linear

# survival ~ 