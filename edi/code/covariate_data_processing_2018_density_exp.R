##### outputs ####

# covariates_2018_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)
library(betareg)

# import data
sj <- read_csv("data/soil_moisture_jun_2018_density_exp.csv") # June soil moisture
so <- read_csv("data/soil_moisture_oct_2018_density_exp.csv") # October soil moisture
cc <- read_csv("data/canopy_cover_jul_2018_density_exp.csv") # canopy cover
bj <- read_csv("data/plot_edge_mv_weight_jul_2018_density_exp.csv") # July biomass and infection
be <- read_csv("data/plot_edge_mv_weight_early_aug_2018_density_exp.csv") # early August biomass
bl <- read_csv("data/plot_edge_mv_weight_late_aug_2018_density_exp.csv") # late August biomass
bs <- read_csv("data/plot_edge_mv_weight_sep_2018_density_exp.csv") # September biomass and infection
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# June soil moisture
sj2 <- sj %>%
  mutate(soil_moisture_jun.prop = soil_moisture.vwc/100) %>%
  select(site, plot, treatment, soil_moisture_jun.prop)

# October soil moisture
so2 <- so %>%
  rowwise() %>%
  mutate(soil_moisture_oct.prop = mean(c(soil_moisture.vwc.1, soil_moisture.vwc.2, soil_moisture.vwc.3))/100) %>%
  select(site, plot, treatment, soil_moisture_oct.prop)

# Canopy cover
unique(cc$processing_notes)
cc2 <- cc %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                     type == "c" ~ (count * 1.04) / 100)) %>%
  select(site, plot, treatment, canopy_cover.prop)

# July plot edge biomass  
unique(bj$process_notes)
hist(bj$mv.g)
filter(bj, mv.g > 8)
hist(bj$mv_inf.g)
filter(bj, mv_inf.g > 1)

bj2 <- bj %>%
  select(site, plot, treatment, mv.g, mv_inf.g) %>%
  mutate(mv_jul.g = mv.g + mv_inf.g,
         mv_inf_jul.prop = mv_inf.g / mv_jul.g,
         mv_inf_jul.g = mv_inf.g) %>%
  select(site, plot, treatment, mv_jul.g, mv_inf_jul.prop, mv_inf_jul.g)

# Early August plot edge biomass
unique(be$process_notes)
hist(be$mv.g)
filter(be, mv.g > 15)

be2 <- be %>%
  select(site, plot, treatment, mv.g) %>%
  rename(mv_eau.g = mv.g)

# late August plot edge biomass
unique(bl$process_notes)
hist(bl$mv.g)

bl2 <- bl %>%
  select(site, plot, treatment, mv.g) %>%
  rename(mv_lau.g = mv.g)

# September plot edge biomass
unique(bs$process_notes)
hist(bs$mv.g)
filter(bs, mv.g > 12)
hist(bs$mv_inf.g)
filter(bs, mv_inf.g > 2.5)

bs2 <- bs %>%
  select(site, plot, treatment, mv.g, mv_inf.g) %>%
  mutate(mv_sep.g = mv.g + mv_inf.g,
         mv_inf_sep.prop = mv_inf.g / mv_sep.g,
         mv_inf_sep.g = mv_inf.g) %>%
  select(-c(mv.g, mv_inf.g))


# combine data
co <- full_join(sj2, so2) %>%
  full_join(cc2) %>%
  full_join(bj2) %>%
  full_join(be2) %>%
  full_join(bl2) %>%
  full_join(bs2)


#### select unique covariates ####

#co %>%
#  select(-c(site, plot, treatment)) %>%
#  ggpairs()
# the two soil moisture measurements are correlated - pick one: the October one has fewer extreme data and includes 3 data points
# July biomass correlated with early and late August
# early and late August biomass are correlated
# all biomass measurements have some extreme (high) values
# infected proportion in July and September are only correlated 0.3, the September data has fewer extreme values

# examine biomass values more
sum(is.na(co$mv_jul.g)) #3
sum(is.na(co$mv_eau.g)) #2
sum(is.na(co$mv_lau.g)) #3
sum(is.na(co$mv_sep.g)) #2
# Mv early August and September are uncorrelated and have all values

filter(co, mv_eau.g > 15) %>% 
  select(site, plot, treatment, mv_eau.g)
filter(co, mv_sep.g > 15) %>% 
  select(site, plot, treatment, mv_sep.g)
# these match the hand-written entries

# select unique covariates
co2 <- co %>%
  select(site:canopy_cover.prop, mv_eau.g, mv_sep.g, mv_inf_jul.prop, mv_inf_sep.prop)

# combine with plot data for figures
co3 <- left_join(co2, plots) %>%
  mutate(fungicide = recode(treatment, "water" = 0, "fungicide" = 1))


#### relationship with treatments ####

# standard plot
base_plot <- ggplot(co3, aes(x = background_density, y = soil_moisture_jun.prop, fill = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point", size = 2, shape = 21, position = position_dodge(0.5)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  xlab("Background density") +
  theme_bw() +
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

# June soil moisture
base_plot %+%
  filter(co3, !is.na(soil_moisture_jun.prop)) +
  ylab("June soil moisture")

# October soil moisture
base_plot %+%
  filter(co3, !is.na(soil_moisture_oct.prop)) %+%
  aes(y = soil_moisture_oct.prop) +
  ylab("October soil moisture")

# Canopy cover
base_plot %+%
  filter(co3, !is.na(canopy_cover.prop)) %+%
  aes(y = canopy_cover.prop) +
  ylab("Canopy cover")

# Early August biomass
base_plot %+%
  filter(co3, !is.na(mv_eau.g)) %+%
  aes(y = mv_eau.g) +
  ylab("Early August biomass")

# September biomass
base_plot %+%
  filter(co3, !is.na(mv_sep.g)) %+%
  aes(y = mv_sep.g) +
  ylab("September biomass")

# July infection
base_plot %+%
  filter(co3, !is.na(mv_inf_jul.prop)) %+%
  aes(y = mv_inf_jul.prop) +
  ylab("July proportion infected")

# September infection
base_plot %+%
  filter(co3, !is.na(mv_inf_sep.prop)) %+%
  aes(y = mv_inf_sep.prop) +
  ylab("September proportion infected")


#### statistical tests of treatments ####

# selected tests based on trends in figure

# Mv density and June soil moisture
jsm_mv_mod <- betareg(soil_moisture_jun.prop ~ background_density * fungicide, data = filter(co3, background == "Mv seedling"))
summary(jsm_mv_mod)
# not sig

# Ev adult density and October soil moisture
osm_ea_mod <- betareg(soil_moisture_oct.prop ~ background_density * fungicide, data = filter(co3, background == "Ev adult"))
summary(osm_ea_mod)
# not sig

# Ev seedling density and October soil moisture
osm_es_mod <- betareg(soil_moisture_oct.prop ~ background_density * fungicide, data = filter(co3, background == "Ev seedling"))
summary(osm_es_mod)
# not sig

# Mv density and October soil moisture
osm_mv_mod <- betareg(soil_moisture_oct.prop ~ background_density * fungicide, data = filter(co3, background == "Mv seedling"))
summary(osm_mv_mod)
# not sig

# Mv density and canopy cover
cc_mv_mod <- betareg(canopy_cover.prop ~ background_density * fungicide, data = filter(co3, background == "Mv seedling"))
summary(cc_mv_mod)
# not sig

# Ev adult density and early August biomass
eab_ea_mod <- lm(mv_eau.g ~ background_density * fungicide, data = filter(co3, background == "Ev adult"))
summary(eab_ea_mod)
# not sig

# Ev seedling density and early August biomass
eab_es_mod <- lm(mv_eau.g ~ background_density * fungicide, data = filter(co3, background == "Ev seedling"))
summary(eab_es_mod)
# not sig

# Mv density and early August biomass
eab_mv_mod <- lm(mv_eau.g ~ background_density * fungicide, data = filter(co3, background == "Mv seedling"))
summary(eab_mv_mod)
# density sig
eab_mv_mod2 <- update(eab_mv_mod, .~. -background_density:fungicide)
summary(eab_mv_mod2)
# density sig
eab_mv_mod3 <- update(eab_mv_mod2, .~. -fungicide)
summary(eab_mv_mod3)
# density sig

# Ev adult density and September biomass
sb_ea_mod <- lm(mv_sep.g ~ background_density * fungicide, data = filter(co3, background == "Ev adult"))
summary(sb_ea_mod)
# not sig

# Ev seedling density and September biomass
sb_es_mod <- lm(mv_sep.g ~ background_density * fungicide, data = filter(co3, background == "Ev seedling"))
summary(sb_es_mod)
# not sig

# Mv density and September biomass
sb_mv_mod <- lm(mv_sep.g ~ background_density * fungicide, data = filter(co3, background == "Mv seedling"))
summary(sb_mv_mod)
# not sig

# Ev seedling density and July infected
jpi_es_mod <- betareg(mv_inf_jul.prop ~ background_density * fungicide, data = filter(co3, background == "Ev seedling" & mv_inf_jul.prop > 0))
summary(jpi_es_mod)
# not sig

# Mv density and July infected
jpi_mv_mod <- betareg(mv_inf_jul.prop ~ background_density * fungicide, data = filter(co3, background == "Mv seedling" & mv_inf_jul.prop > 0))
summary(jpi_mv_mod)
# significant interaction

# Ev seedling density and September infected
spi_es_mod <- betareg(mv_inf_sep.prop ~ background_density * fungicide, data = filter(co3, background == "Ev seedling"))
summary(spi_es_mod)
# not sig

# Mv density and September infected
spi_mv_mod <- betareg(mv_inf_sep.prop ~ background_density * fungicide, data = filter(co3, background == "Mv seedling"))
summary(spi_mv_mod)
# fungicide effect
# there is a trend of this across the three densities
spi_mv_mod2 <- update(spi_mv_mod, .~. -background_density:fungicide)
summary(spi_mv_mod2)
# not sig
spi_mv_mod3 <- update(spi_mv_mod2, .~. -background_density)
summary(spi_mv_mod3)
# not sig


#### remove covariates with treatment effects ####

co4 <- co2 %>%
  select(-mv_eau.g)
# added July disease back in because I wanted it for an analysis


#### output intermediate data ####

write_csv(co4, "intermediate-data/covariates_2018_density_exp.csv")

