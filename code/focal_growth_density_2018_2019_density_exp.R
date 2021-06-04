##### info ####

# file: focal_growth_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/3/21
# goal: analyses of plant growth as a function of density


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi
library(cowplot)

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import growth data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(density = case_when(plot %in% 2:4 ~ background_density + 3,
                             plot %in% 5:7 ~ background_density + 3,
                             plot%in% 8:10 ~ background_density + 1,
                             plot == 1 ~ NA_real_)) %>%
  select(plot, treatment, background, density)

# missing data
filter(growthD1Dat, is.na(tillers_jul)) # 0
filter(growthD1Dat, is.na(tillers_jun)) # 0
filter(growthD1Dat, tillers_jul == 0 | tillers_jun == 0) # these plants are dead
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling

# combine data
growthD1Dat2 <- growthD1Dat %>%
  left_join(plotDens) %>%
  mutate(plant_growth = log(tillers_jul/tillers_jun),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(plot, str_sub(treatment, 1, 1), sep = "")) %>%
  filter(plot > 1 & tillers_jun > 0 & tillers_jul > 0)

growthD2Dat2 <- mvBioD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  full_join(evBioD2Dat %>%
              rename(biomass_weight.g = weight)) %>%
  left_join(plotDens) %>%
  mutate(plant_growth = log(biomass_weight.g),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(plot, str_sub(treatment, 1, 1), sep = "")) %>%
  filter(!is.na(biomass_weight.g) & plot > 1)


#### models ####

# initial visualization
ggplot(growthD1Dat2, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(foc ~ bg, scales = "free")

ggplot(growthD2Dat2, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(focal ~ background, scales = "free")

# fit models
growthD1Mod <- brm(data = growthD1Dat2, family = gaussian,
                   plant_growth ~ density * foc * bg * fungicide + (1|plotf),
                   prior <- c(prior(normal(0.75, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3) 

mod_check_fun(growthD1Mod)

growthD2Mod <- update(growthD1Mod, newdata = growthD2Dat2)
mod_check_fun(growthD2Mod)

# save models
save(growthD1Mod, file = "output/focal_growth_density_model_2018_density_exp.rda")
save(growthD2Mod, file = "output/focal_growth_density_model_2019_density_exp.rda")


#### values for text ####

# common terms on both sides of = were deleted

# Mv intraspecific vs. interspecific
evS_mv_ctrl_hyp = "density:focs = 0"
evS_mv_fung_hyp = "density:focs + density:focs:fungicide = 0"
evA_mv_ctrl_hyp = "density:foca = 0"
evA_mv_fung_hyp = "density:foca + density:foca:fungicide = 0"

# EvS intra vs. inter
mv_evS_ctrl_hyp = "density:focs + density:focs:bgs  = 0"
mv_evS_fung_hyp = "density:focs + density:focs:bgs + density:focs:fungicide + density:focs:bgs:fungicide  = 0"
evA_evS_ctrl_hyp = "density:focs + density:focs:bgs  = density:foca + density:foca:bgs"
evA_evS_fung_hyp = "density:focs + density:focs:bgs + density:focs:fungicide + density:focs:bgs:fungicide  = density:foca + density:foca:bgs + density:foca:fungicide + density:foca:bgs:fungicide"

# EvA intra vs. inter
mv_evA_ctrl_hyp = "density:foca + density:foca:bga  = 0"
mv_evA_fung_hyp = "density:foca + density:foca:bga + density:foca:fungicide + density:foca:bga:fungicide  = 0"
evS_evA_ctrl_hyp = "density:foca + density:foca:bga  = density:focs + density:focs:bga"
evS_evA_fung_hyp = "density:foca + density:foca:bga + density:foca:fungicide + density:foca:bga:fungicide  = density:focs + density:focs:bga + density:focs:fungicide + density:focs:bga:fungicide"

growthD1hyps <- hypothesis(growthD1Mod, 
           c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
           mv_evS_ctrl_hyp, mv_evS_fung_hyp, evA_evS_ctrl_hyp, evA_evS_fung_hyp,
           mv_evA_ctrl_hyp, mv_evA_fung_hyp, evS_evA_ctrl_hyp, evS_evA_fung_hyp))
# none are significantly different from zero

growthD2hyps <- hypothesis(growthD2Mod, 
           c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
             mv_evS_ctrl_hyp, mv_evS_fung_hyp, evA_evS_ctrl_hyp, evA_evS_fung_hyp,
             mv_evA_ctrl_hyp, mv_evA_fung_hyp, evS_evA_ctrl_hyp, evS_evA_fung_hyp))
# none are significantly different from zero

write_csv(growthD1hyps[[1]], "output/focal_growth_intra_vs_intra_comp_2018_density_exp.csv")
write_csv(growthD2hyps[[1]], "output/focal_growth_intra_vs_intra_comp_2019_density_exp.csv")


#### figure ####

# density function
dens_fun <- function(foc, bg, treatment, year){
  
  if(year == "2018"){
    dat <- growthD1Dat2 %>% filter(foc == foc & bg == bg & treatment == treatment)
  }else{
    dat <- growthD2Dat2 %>% filter(foc == foc & bg == bg & treatment == treatment)
  }
  
  density = seq(min(growthD1Dat2$density), max(growthD1Dat2$density), length.out = 100)
  
  return(density)
}

# predicted data
predD1Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2018",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(density = pmap(list(foc, bg, treatment, year), dens_fun)) %>%
  unnest(density) %>%
  mutate(plant_growth = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predD2Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2018",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(density = pmap(list(foc, bg, treatment, year), dens_fun)) %>%
  unnest(density) %>%
  mutate(plant_growth = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predDat <- predD1Dat %>%
  full_join(predD2Dat) %>%
  mutate(focal = fct_recode(foc, Mv = "m", "Ev seedling" = "s", "Ev adult" = "a"),
         background = fct_recode(bg, Mv = "m", "Ev seedling" = "s", "Ev adult" = "a"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

#### start here ####

# Mv background
mv_mv_ctrl_alpha = "density = 0"
mv_mv_ctrl_alpha = "density + density:fungicide = 0"
evS_mv_ctrl_alpha = "density + density:focs = 0"
evS_mv_fung_alpha = "density + density:fungicide + density:focs + density:focs:fungicide = 0"
evA_mv_ctrl_alpha = "density + density:foca = 0"
evA_mv_fung_alpha = "density + density:fungicide + density:foca + density:foca:fungicide = 0"

# EvS background
evS_evS_ctrl_alpha = "density + density:focs + density:bgs + density:focs:bgs = 0"
evS_evS_ctrl_alpha = "density + density:focs + density:bgs + density:focs:bgs + density:fungicide + density:focs:fungicide + density:bgs:fungicide + density:focs:bgs:fungicide  = 0"
mv_evS_ctrl_alpha = "density + density:focs + density:focs:bgs  = 0"
mv_evS_fung_alpha = "density + density:fungicide + density:focs + density:focs:bgs + density:focs:fungicide + density:focs:bgs:fungicide  = 0"
evA_evS_ctrl_alpha = "density:focs + density:focs:bgs  = density:foca + density:foca:bgs"
evA_evS_fung_alpha = "density:focs + density:focs:bgs + density:focs:fungicide + density:focs:bgs:fungicide  = density:foca + density:foca:bgs + density:foca:fungicide + density:foca:bgs:fungicide"

# EvA intra vs. inter
mv_evA_ctrl_alpha = "density:foca + density:foca:bga  = 0"
mv_evA_fung_alpha = "density:foca + density:foca:bga + density:foca:fungicide + density:foca:bga:fungicide  = 0"
evS_evA_ctrl_alpha = "density:foca + density:foca:bga  = density:focs + density:focs:bga"
evS_evA_fung_alpha = "density:foca + density:foca:bga + density:foca:fungicide + density:foca:bga:fungicide  = density:focs + density:focs:bga + density:focs:fungicide + density:focs:bga:fungicide"


# combine alphas
alphasDat <- posterior_samples(growthD1Mod) %>%
  mutate(year = "2018") %>%
  full_join(posterior_samples(growthD2Mod) %>%
              mutate(year = "2019")) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  mutate(evA_evA_water = b_density,
         evA_evS_water = b_density + b_densitybackgroundEv_seedling,
         evA_mv_water = b_density + b_densitybackgroundMv_seedling,
         evS_evA_water = b_density + b_densityfocalEv_seedling,
         evS_evS_water = evA_evS_water + b_densityfocalEv_seedling + b_densityfocalEv_seedlingbackgroundEv_seedling,
         evS_mv_water = evA_mv_water + b_densityfocalEv_seedling + b_densityfocalEv_seedlingbackgroundMv_seedling,
         mv_evA_water = b_density + b_densityfocalMv_Seedling,
         evS_evS_water = evA_evS_water + b_densityfocalEv_seedling + b_densityfocalEv_seedlingbackgroundEv_seedling,
        ) %>%
  select(year, focal, a_water, s_water, p_water, a_fungicide, s_fungicide, p_fungicide) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(year, focal),
               names_to = c(".value", "bg_treatment"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_longer(cols = -c(year, focal, bg_treatment),
               names_to = "background",
               values_to = "alpha") %>%
  group_by(year, treatment, focal, background) %>%
  mean_hdi(alpha) %>%
  ungroup() %>%
  mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"),
         focal = fct_recode(focal, "Ev adult" = "p", "Ev seedling" = "s", "Mv" = "a"),
         background = fct_recode(background, "Ev adult" = "p", "Ev seedling" = "s", "Mv" = "a"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# add to predicted dataset
predDat2 <- predDat %>%
  inner_join(alphasDat %>%
              select(year, treatment, focal, background, sig))

# raw data
figDat <- d1dat %>%
  select(site, plot, treatment, sp, age, ID, plant_group, background, density, plant_growth) %>%
  mutate(year = "2018") %>%
  full_join(d2dat %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, plant_growth) %>%
              mutate(year = "2019"))%>%
  mutate(focal = str_replace(plant_group, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)")) %>%
  filter(background != "none" & !is.na(plant_growth) & plant_growth != Inf & plant_growth != -Inf)

# sig alpha values
alphasDat2 <- alphasDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, background) %>%
              summarise(density = max(density))) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(upper = max(upper)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(raw = max(plant_growth))) %>%
              rowwise() %>%
              mutate(plant_growth = max(c(raw, upper))) %>%
              ungroup() %>%
              select(-c(raw, upper))) %>%
  rename(param = alpha) %>%
  mutate(param = round(param, 2),
         plant_growth = case_when(treatment == "fungicide" & background == "Mv" & focal == "Ev seedling" ~ plant_growth - 0.5,
                                  treatment == "fungicide" & background == "Mv" & focal == "Mv" ~ plant_growth - 0.6,
                                  TRUE ~ plant_growth + 0.1))

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.1, 0, 0.2, 2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

yearText <- tibble(year = c("2018", "2019"),
                   density = c(3.2, 3),
                   plant_growth = c(0.87, 2.7),
                   background = "Ev adult",
                   focal = "Ev adult",
                   treatment = "fungicide")

textSize = 3

# water figure
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(alphasDat2, year == "2018"), 
            aes(label = paste("alpha", " == ", param, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# fungicide figure
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(alphasDat2, year == "2019"), 
            aes(label = paste("alpha", " == ", param, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg <- get_legend(pairD1Fig)

# combine plots
combFig <- plot_grid(pairD1Fig + theme(legend.position = "none"), pairD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.9),
                     label_x = c(0, -0.01))

# combine
pdf("output/focal_growth_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()

# save alpha data
write_csv(alphasDat, "output/focal_growth_pairwise_coefficients_2018_2019_density_exp.csv")

