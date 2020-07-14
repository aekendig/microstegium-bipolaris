##### info ####

# file: ev_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/13/20
# goal: evaluate the effects of density treatments and environmental covariates on the seed production of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(glmmTMB)

# import data
bio <- read_csv("./data/ev_biomass_seeds_oct_2019_density_exp.csv")
spike <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
covar <- read_csv("intermediate-data/covariates_2018_density_exp.csv")
sev19 <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")


#### edit data ####

# severity data
sev19b <- sev19 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Ev")

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
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(treatment, water = "control (water)")) %>%
  rename(veg_weight.g = weight) %>%
  left_join(spike2) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         tot_weight.g = veg_weight.g + spikelet_weight.g) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  left_join(sev19b)

# no background data
no_back_plant_dat <- bio2 %>%
  filter(plot == 1 & background == "Mv seedling") %>%
  mutate(adj_bio.g = case_when(treatment == "fungicide" ~ veg_weight.g / 1.03,
                               TRUE ~ veg_weight.g))


#### check data ####

# max value
max(bio2$veg_weight.g)

# sample sizes
samps <- bio2 %>%
  group_by(site, plot, treatment) %>%
  summarise(n = sum(!is.na(veg_weight.g))) %>%
  mutate(veg_weight.g = 23)

# visualize
ggplot(bio2, aes(x = plot, y = veg_weight.g)) +
  geom_point(alpha = 0.5, aes(color = age)) +
  facet_grid(treatment ~ site) +
  geom_text(data = samps, aes(label = n), size = 2)
# one missing in D2


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

# colors
col_pal = c("#a6611a", "#018571")

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

# no background plots
pdf("output/ev_biomass_no_background_treatment_2018_2019_density_exp.pdf", width = 5, height = 3)
ggplot(no_back_plant_dat, aes(x = Treatment, y = veg_weight.g)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_wrap(~age) +
  scale_fill_manual(values = col_pal, guide = F) +
  ylab("Biomass (g)") +
  ggtitle("Ev") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

ggplot(no_back_plant_dat, aes(x = Treatment, y = adj_bio.g)) +
  stat_summary(geom = "bar", fun = "mean", aes(fill = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  facet_wrap(~age) +
  scale_fill_manual(values = col_pal, guide = F) +
  ylab("Adjusted biomass (g)") +
  ggtitle("Ev") +
  temp_theme +
  theme(plot.title = element_text(size = 12, hjust = 0.5))
dev.off()


#### direct disease models ####

# divide by age
evs_no_back_plant_dat = filter(no_back_plant_dat, age == "seedling")
eva_no_back_plant_dat = filter(no_back_plant_dat, age == "adult")

# control biomass
filter(evs_no_back_plant_dat, fungicide == 0) %>%
  summarise(bio = mean(veg_weight.g))
filter(eva_no_back_plant_dat, fungicide == 0) %>%
  summarise(bio = mean(veg_weight.g))

# Ev seedling model
evs_no_back_bio_mod <- brm(data = evs_no_back_plant_dat, family = gaussian,
                          veg_weight.g ~ fungicide + (1|site),
                          prior <- c(prior(normal(0.9, 10), class = Intercept),
                                     prior(normal(0, 10), class = b),
                                     prior(cauchy(0, 1), class = sd),
                                     prior(cauchy(0, 1), class = sigma)),
                          iter = 6000, warmup = 1000, chains = 1,
                          control = list(adapt_delta = 0.99))

# check model and add chains
summary(evs_no_back_bio_mod)
prior_summary(evs_no_back_bio_mod)
evs_no_back_bio_mod <- update(evs_no_back_bio_mod, chains = 3)
plot(evs_no_back_bio_mod)
pp_check(evs_no_back_bio_mod, nsamples = 100)

# increase due to fungicide
0.03*0.91

# model with direct fungicide effects
evs_no_back_bio_fung_mod <- brm(data = evs_no_back_plant_dat, family = gaussian,
                               veg_weight.g ~ Treatment + fungicide + (1|site),
                               prior <- c(prior(normal(0.91, 10), class = Intercept),
                                          prior(normal(0, 10), class = b),
                                          prior(normal(0.03, 0.001), class = b, coef = "Treatmentfungicide"),
                                          prior(cauchy(0, 1), class = sd),
                                          prior(cauchy(0, 1), class = sigma)),
                               iter = 6000, warmup = 1000, chains = 3,
                               control = list(adapt_delta = 0.99))

# check model
summary(evs_no_back_bio_fung_mod)
prior_summary(evs_no_back_bio_fung_mod)
plot(evs_no_back_bio_fung_mod)
pp_check(evs_no_back_bio_fung_mod, nsamples = 100)

# Ev adult model
eva_no_back_bio_mod <- brm(data = eva_no_back_plant_dat, family = gaussian,
                           veg_weight.g ~ fungicide + (1|site),
                           prior <- c(prior(normal(4.5, 10), class = Intercept),
                                      prior(normal(0, 10), class = b),
                                      prior(cauchy(0, 1), class = sd),
                                      prior(cauchy(0, 1), class = sigma)),
                           iter = 6000, warmup = 1000, chains = 1,
                           control = list(adapt_delta = 0.99))

# check model and add chains
summary(eva_no_back_bio_mod)
prior_summary(eva_no_back_bio_mod)
eva_no_back_bio_mod <- update(eva_no_back_bio_mod, chains = 3)
plot(eva_no_back_bio_mod)
pp_check(eva_no_back_bio_mod, nsamples = 100)

# increase due to fungicide
0.03*4.53

# model with direct fungicide effects
eva_no_back_bio_fung_mod <- brm(data = eva_no_back_plant_dat, family = gaussian,
                                veg_weight.g ~ Treatment + fungicide + (1|site),
                                prior <- c(prior(normal(4.53, 10), class = Intercept),
                                           prior(normal(0, 10), class = b),
                                           prior(normal(0.14, 0.001), class = b, coef = "Treatmentfungicide"),
                                           prior(cauchy(0, 1), class = sd),
                                           prior(cauchy(0, 1), class = sigma)),
                                iter = 6000, warmup = 1000, chains = 3,
                                control = list(adapt_delta = 0.99))

# check model
summary(eva_no_back_bio_fung_mod)
prior_summary(eva_no_back_bio_fung_mod)
plot(eva_no_back_bio_fung_mod)
pp_check(eva_no_back_bio_fung_mod, nsamples = 100)


#### Ev seedling disease and density models ####

# separate data
evs_mv_bio_dat <- bio2 %>%
  filter(background == "Mv seedling" & plot != 1 & age == "seedling")

evs_evs_bio_dat <- bio2 %>%
  filter(background == "Ev seedling" & plot != 1  & age == "seedling")

evs_eva_bio_dat <- bio2 %>%
  filter(background == "Ev adult" & plot != 1  & age == "seedling")

# check for complete values
sum(is.na(evs_mv_bio_dat$veg_weight.g))
sum(is.na(evs_evs_bio_dat$veg_weight.g))
sum(is.na(evs_eva_bio_dat$veg_weight.g))

# Beverton-Holt function
bh_fun <- function(dat_in, a, y0){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$background_density)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = veg_weight.g)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# try values
mean(filter(evs_mv_bio_dat, density_level == "low")$veg_weight.g)
bh_fun(evs_mv_bio_dat, 0.03, 2.5)
mean(evs_evs_bio_dat$veg_weight.g)
bh_fun(evs_evs_bio_dat, 0, 1.9)
mean(filter(evs_eva_bio_dat, density_level == "low")$veg_weight.g)
bh_fun(evs_eva_bio_dat, 0.3, 3.5)

# distributions
x <- seq(0, 10, length.out = 100)
y <- dexp(x, 1/0.1)
plot(x, y, type = "l")

# mv background model
evs_mv_bio_mod <- brm(data = evs_mv_bio_dat, family = gaussian,
                     bf(veg_weight.g ~ y0 / (1 + alpha * background_density),
                        y0 ~ 0 + Treatment + (1|site),
                        alpha ~ 0 + Treatment,
                        nl = T),
                     prior <- c(prior(normal(2.5, 10), nlpar = "y0", class = "b", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                prior(cauchy(0, 1), class = "sigma")),
                     iter = 6000, warmup = 1000, chains = 1,
                     control = list(adapt_delta = 0.999))

# check model and increase chains
summary(evs_mv_bio_mod)
prior_summary(evs_mv_bio_mod)
evs_mv_bio_mod <- update(evs_mv_bio_mod, chains = 3,
                         control = list(adapt_delta = 0.999, max_treedepth = 15))
plot(evs_mv_bio_mod)
pp_check(evs_mv_bio_mod, nsamples = 100)

# ev seedling background model
evs_evs_bio_mod <- brm(data = evs_evs_bio_dat, family = gaussian,
                      bf(veg_weight.g ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(1.9, 10), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(evs_evs_bio_mod)
evs_evs_bio_mod <- update(evs_evs_bio_mod, chains = 3,
                          control = list(adapt_delta = 0.99999))
plot(evs_evs_bio_mod)
pp_check(evs_evs_bio_mod, nsamples = 100)

# ev adult background model
evs_eva_bio_mod <- brm(data = evs_eva_bio_dat, family = gaussian,
                      bf(veg_weight.g ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(3.5, 10), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(evs_eva_bio_mod)
evs_eva_bio_mod <- update(evs_eva_bio_mod, chains = 3)
plot(evs_eva_bio_mod)
pp_check(evs_eva_bio_mod, nsamples = 100)


#### Ev adult disease and density models ####

# separate data
eva_mv_bio_dat <- bio2 %>%
  filter(background == "Mv seedling" & plot != 1 & age == "adult")

eva_evs_bio_dat <- bio2 %>%
  filter(background == "Ev seedling" & plot != 1  & age == "adult")

eva_eva_bio_dat <- bio2 %>%
  filter(background == "Ev adult" & plot != 1  & age == "adult")

# check for complete values
sum(is.na(eva_mv_bio_dat$veg_weight.g))
sum(is.na(eva_evs_bio_dat$veg_weight.g))
sum(is.na(eva_eva_bio_dat$veg_weight.g))

# try values
mean(filter(eva_mv_bio_dat, density_level == "low")$veg_weight.g)
bh_fun(eva_mv_bio_dat, 0.03, 10)
mean(eva_evs_bio_dat$veg_weight.g)
bh_fun(eva_evs_bio_dat, 0, 6)
mean(filter(eva_eva_bio_dat, density_level == "low")$veg_weight.g)
bh_fun(eva_eva_bio_dat, 0.3, 13)

# mv background model
eva_mv_bio_mod <- brm(data = eva_mv_bio_dat, family = gaussian,
                      bf(veg_weight.g ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(10, 10), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(eva_mv_bio_mod)
prior_summary(eva_mv_bio_mod)
eva_mv_bio_mod <- update(eva_mv_bio_mod, chains = 3) # redo
plot(eva_mv_bio_mod)
pp_check(eva_mv_bio_mod, nsamples = 100)

# ev seedling background model
eva_evs_bio_mod <- brm(data = eva_evs_bio_dat, family = gaussian,
                       bf(veg_weight.g ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(6, 10), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 6000, warmup = 1000, chains = 1,
                       control = list(adapt_delta = 0.999))

# check model and increase chains
summary(eva_evs_bio_mod)
eva_evs_bio_mod <- update(eva_evs_bio_mod, chains = 3)
plot(eva_evs_bio_mod)
pp_check(eva_evs_bio_mod, nsamples = 100)

# ev adult background model
eva_eva_bio_mod <- brm(data = eva_eva_bio_dat, family = gaussian,
                       bf(veg_weight.g ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(13, 10), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 6000, warmup = 1000, chains = 1,
                       control = list(adapt_delta = 0.999))

# check model and increase chains
summary(eva_eva_bio_mod)
eva_eva_bio_mod <- update(eva_eva_bio_mod, chains = 3)
plot(eva_eva_bio_mod)
pp_check(eva_eva_bio_mod, nsamples = 100)


#### model fit figure ####

# model fits over simulated data
bio_mv_sim_dat <- tibble(background_density = rep(seq(0, 64, length.out = 300), 2),
                         fungicide = rep(c(0, 1), each = 300)) %>%
  merge(tibble(age = c("seedling", "adult")), all = T) %>%
  as_tibble() %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(veg_weight.g = case_when(age == "seedling" ~ fitted(evs_mv_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
                                  age == "adult" ~ fitted(eva_mv_bio_mod, newdata = ., re_formula = NA)[, "Estimate"]),
         bio_lower = case_when(age == "seedling" ~ fitted(evs_mv_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
                               age == "adult" ~ fitted(eva_mv_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"]),
         bio_upper = case_when(age == "seedling" ~ fitted(evs_mv_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"],
                               age == "adult" ~ fitted(eva_mv_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"]))

bio_evs_sim_dat <- tibble(background_density = rep(seq(0, 16, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  merge(tibble(age = c("seedling", "adult")), all = T) %>%
  as_tibble() %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(veg_weight.g = case_when(age == "seedling" ~ fitted(evs_evs_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
                                  age == "adult" ~ fitted(eva_evs_bio_mod, newdata = ., re_formula = NA)[, "Estimate"]),
         bio_lower = case_when(age == "seedling" ~ fitted(evs_evs_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
                               age == "adult" ~ fitted(eva_evs_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"]),
         bio_upper = case_when(age == "seedling" ~ fitted(evs_evs_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"],
                               age == "adult" ~ fitted(eva_evs_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"]))

bio_eva_sim_dat <- tibble(background_density = rep(seq(0, 8, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  merge(tibble(age = c("seedling", "adult")), all = T) %>%
  as_tibble() %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(veg_weight.g = case_when(age == "seedling" ~ fitted(evs_eva_bio_mod, newdata = ., re_formula = NA)[, "Estimate"],
                                  age == "adult" ~ fitted(eva_eva_bio_mod, newdata = ., re_formula = NA)[, "Estimate"]),
         bio_lower = case_when(age == "seedling" ~ fitted(evs_eva_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"],
                               age == "adult" ~ fitted(eva_eva_bio_mod, newdata = ., re_formula = NA)[, "Q2.5"]),
         bio_upper = case_when(age == "seedling" ~ fitted(evs_eva_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"],
                               age == "adult" ~ fitted(eva_eva_bio_mod, newdata = ., re_formula = NA)[, "Q97.5"]))

# combine simulated data
bio_sim_dat <- bio_mv_sim_dat %>%
  mutate(background = "Mv seedling") %>%
  full_join(bio_evs_sim_dat %>%
              mutate(background = "Ev seedling")) %>%
  full_join(bio_eva_sim_dat %>%
              mutate(background = "Ev adult"))

# figure
pdf("output/ev_biomass_analysis_model_fits_2019_density_exp.pdf")
ggplot(bio2, aes(x = background_density, y = veg_weight.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = Treatment)) +
  geom_ribbon(data = bio_sim_dat, alpha = 0.5, aes(ymin = bio_lower, ymax = bio_upper, fill = Treatment)) +
  geom_line(data = bio_sim_dat, aes(color = Treatment)) +
  facet_grid(age~background, scales = "free", switch = "both") +
  xlab("Background density") +
  ylab(expression(paste(italic(Elymus), " biomass (g ", plant^-1, ")", sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  temp_theme
dev.off()


#### severity and biomass ####

# divide data
# remove extra 1 plots
evs_sev_dat <- filter(bio2, treatment == "water" & plant_type == "Ev seedling" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_")) %>%
  as.data.frame()
eva_sev_dat <- filter(bio2, treatment == "water" & plant_type == "Ev adult" & (plot != 1 | (plot == 1 & background == "Mv seedling"))) %>%
  mutate(site_plot = paste(site, plot, sep = "_")) %>%
  as.data.frame()

# check for density-severity correlations
for(i in 27:31){
  evs_sev_dat$severity = evs_sev_dat[,i]
  print(cor.test(evs_sev_dat$background_density, evs_sev_dat$severity))
  
  eva_sev_dat$severity = eva_sev_dat[,i]
  print(cor.test(eva_sev_dat$background_density, eva_sev_dat$severity))
}
# none significantly correlated

# biomass-severity relationships
evs_bio_sev_jun_19_mod <- glmmTMB(veg_weight.g ~ severity_jun + (1|site/plot), family = gaussian, data = evs_sev_dat)
summary(evs_bio_sev_jun_19_mod) # not sig
evs_bio_sev_jul_19_mod <- glmmTMB(veg_weight.g ~ severity_jul + (1|site/plot), family = gaussian, data = evs_sev_dat)
summary(evs_bio_sev_jul_19_mod) # sig reduction
evs_bio_sev_eau_19_mod <- glmmTMB(veg_weight.g ~ severity_early_aug + (1|site/plot), family = gaussian, data = evs_sev_dat)
summary(evs_bio_sev_eau_19_mod) # not sig
evs_bio_sev_lau_19_mod <- glmmTMB(veg_weight.g ~ severity_late_aug + (1|site/plot), family = gaussian, data = evs_sev_dat)
summary(evs_bio_sev_lau_19_mod) # not sig

eva_bio_sev_jun_19_mod <- glmmTMB(veg_weight.g ~ severity_jun + (1|site/plot), family = gaussian, data = eva_sev_dat) # convergence issue
eva_bio_sev_jun_19_mod <- glmmTMB(veg_weight.g ~ severity_jun + (1|site_plot), family = gaussian, data = eva_sev_dat) # convergence issue
eva_bio_sev_jul_19_mod <- glmmTMB(veg_weight.g ~ severity_jul + (1|site/plot), family = gaussian, data = eva_sev_dat) # convergence issue
eva_bio_sev_jul_19_mod <- glmmTMB(veg_weight.g ~ severity_jul + (1|site_plot), family = gaussian, data = eva_sev_dat)
summary(eva_bio_sev_jul_19_mod) # not sig
eva_bio_sev_eau_19_mod <- glmmTMB(veg_weight.g ~ severity_early_aug + (1|site/plot), family = gaussian, data = eva_sev_dat) # convergence issue
eva_bio_sev_eau_19_mod <- glmmTMB(veg_weight.g ~ severity_early_aug + (1|site_plot), family = gaussian, data = eva_sev_dat)
summary(eva_bio_sev_eau_19_mod) # not sig
eva_bio_sev_lau_19_mod <- glmmTMB(veg_weight.g ~ severity_late_aug + (1|site/plot), family = gaussian, data = eva_sev_dat) # convergence issue
eva_bio_sev_lau_19_mod <- glmmTMB(veg_weight.g ~ severity_late_aug + (1|site_plot), family = gaussian, data = eva_sev_dat)
# convergence issue


#### output ####

save(evs_no_back_bio_mod, file = "output/ev_seedling_biomass_no_background_model_2019_density_exp.rda")
save(evs_no_back_bio_fung_mod, file = "output/ev_seedling_biomass_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
save(eva_no_back_bio_mod, file = "output/ev_adult_biomass_no_background_model_2019_density_exp.rda")
save(eva_no_back_bio_fung_mod, file = "output/ev_adult_biomass_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
save(evs_mv_bio_mod, file = "output/ev_seedling_biomass_mv_background_model_2019_density_exp.rda")
save(evs_evs_bio_mod, file = "output/ev_seedling_biomass_ev_seedling_background_model_2019_density_exp.rda")
save(evs_eva_bio_mod, file = "output/ev_seedling_biomass_ev_adult_background_model_2019_density_exp.rda")
save(eva_mv_bio_mod, file = "output/ev_adult_biomass_mv_background_model_2019_density_exp.rda")
save(eva_evs_bio_mod, file = "output/ev_adult_biomass_ev_seedling_background_model_2019_density_exp.rda")
save(eva_eva_bio_mod, file = "output/ev_adult_biomass_ev_adult_background_model_2019_density_exp.rda")
save(evs_bio_sev_jul_19_mod, file = "output/ev_seedling_biomass_severity_model_jul_2019_density_exp.rda")