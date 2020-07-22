##### info ####

# file: mv_germination_disease_analysis_2018_litter_exp
# author: Amy Kendig
# date last edited: 2/3/20
# goal: evaluate the effects of litter treatments and environmental covariates on the germination and disease of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)

# import data
dat_jun <- read_csv("./data/both_germination_disease_jun_2018_litter_exp.csv")
dat_jul <- read_csv("./data/both_germination_disease_jul_2018_litter_exp.csv")
cov_plot <- read_csv("./intermediate-data/plot_covariates_2018_litter_exp.csv")
plots <- read_csv("./data/plot_treatments_2018_litter_exp.csv")


#### edit data ####

# edit plot data so that no litter is "sterilized"
# scale litter weight
# remove unnecessary variables
plots2 <- plots %>%
  mutate(sterilized = case_when(litter == "live" ~ 0,
                                TRUE ~ 1),
         litter_weight.scaled = (litter_weight.g - mean(litter_weight.g)) / sd(litter_weight.g),
         litter_microbes = recode(sterilized, "0" = "present", "1" = "absent")) %>%
  select(-c(flag_color, justification))

# June data planted + background
dat_jun_planted <- dat_jun %>%
  filter(seeds_added == "yes") %>%
  left_join(plots2) %>%
  select(-c(date, seeds_added, ev_germ, ev_infec)) %>%
  rename(mv_germ_planted_bg_jun = mv_germ,
         mv_infec_jun = mv_infec)

# June data background
dat_jun_bg <- dat_jun %>%
  filter(seeds_added == "no") %>%
  left_join(plots2) %>%
  select(-c(date, seeds_added, mv_infec, ev_germ, ev_infec)) %>%
  rename(mv_germ_bg_jun = mv_germ)

# edit July data
dat_jul2 <- dat_jul %>%
  filter(seeds_added == "yes") %>%
  left_join(plots2) %>%
  select(-c(date, seeds_added, ev_germ, ev_infec)) %>%
  rename(mv_germ_planted_bg_jul = mv_germ,
         mv_infec_jul = mv_infec,
         mv_germ_bg_jul = mv_germ_ev)

# combine data
micro <- full_join(dat_jun_planted, select(dat_jun_bg, -plot)) %>%
  full_join(dat_jul2) %>%
  mutate(mv_germ_planted_jun = mv_germ_planted_bg_jun - mv_germ_bg_jun,
         mv_germ_planted_jul = mv_germ_planted_bg_jul - mv_germ_bg_jul,
         mv_prop_infec_jun = mv_infec_jun / mv_germ_planted_bg_jun,
         mv_prop_infec_jul = mv_infec_jul / mv_germ_planted_bg_jul) %>%
  mutate(mv_germ_planted_cor_jun = case_when(mv_germ_planted_jun < 0 ~ 0,
                                             mv_germ_planted_jun > 200 ~ 200,
                                         TRUE ~ mv_germ_planted_jun),
         mv_germ_planted_cor_jul = case_when(mv_germ_planted_jul < 0 ~ 0,
                                             mv_germ_planted_jul > 200 ~ 200,
                                         TRUE ~ mv_germ_planted_jul),
         prop_germ_jun = mv_germ_planted_cor_jun/200,
         prop_germ_jul = mv_germ_planted_cor_jul/200) %>%
  left_join(cov_plot)

# to combine background data with covariates, use the original plot number (removed above)
dat_jun_bg2 <- dat_jun_bg %>%
  left_join(cov_plot)


#### raw data figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# color palettes
col_pal_mic <- c("white", "black")
col_pal_lit <- c("white", "yellow", "orange", "red")

# base summary figure
base_fig <- ggplot(micro, aes(x = litter_weight.g)) +
  stat_summary(aes(group = litter_microbes), geom = "errorbar", fun.data = mean_cl_boot, width = 5, position = position_dodge(15)) +
  stat_summary(aes(fill = litter_microbes), geom = "point", fun.y = mean, position = position_dodge(15), shape = 21, size = 3) +
  scale_fill_manual(values = col_pal_mic, name = "Litter microbes") +
  xlab("Litter weight (g)") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none") +
  xlim(0, 205)

# site summary figure
site_fig <- ggplot(micro, aes(x = site)) +
  stat_summary(aes(group = litter_density), geom = "errorbar", fun.data = mean_cl_boot, width = 0.2) +
  stat_summary(aes(fill = litter_density), geom = "point", fun.y = mean, shape = 21, size = 2) +
  scale_fill_manual(values = col_pal_lit, name = "Litter amount") +
  xlab("Site") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "none")

# scatterplots
scat_fig <- ggplot(micro, aes(x = soil_moisture.prop, y = mv_germ_planted_bg_jun)) +
  geom_point(aes(fill = site), shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("black", "blue", "purple", "red"), name = "Site")

# Background Mv germination June
fig_bg_jn <- base_fig %+%
  aes(y = mv_germ_bg_jun) +
  ylab("Background\ngerminants in June")

fig_bg_jn_site <- site_fig %+%
  aes(y = mv_germ_bg_jun) +
  ylab("Background\ngerminants in June")

# Background Mv germination July
fig_bg_jl <- base_fig %+%
  aes(y = mv_germ_bg_jul) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Background\ngerminants in July")

fig_bg_jl_site <- site_fig %+%
  aes(y = mv_germ_bg_jul) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Background\ngerminants in July")

# Background + planted Mv germination June
fig_bp_jn <- base_fig %+%
  aes(y = mv_germ_planted_bg_jun) +
  ylab("Background and planted\ngerminants in June")

fig_bp_jn_site <- site_fig %+%
  aes(y = mv_germ_planted_bg_jun) +
  ylab("Background and planted\ngerminants in June")

# Background + planted Mv germination July
fig_bp_jl <- base_fig %+%
  aes(y = mv_germ_planted_bg_jul) +
  ylab("Background and planted\ngerminants in July")

fig_bp_jl_site <- site_fig %+%
  aes(y = mv_germ_planted_bg_jul) +
  ylab("Background and planted\ngerminants in July")

# Planted Mv germination June
fig_pl_jn <- base_fig %+%
  aes(y = mv_germ_planted_cor_jun) +
  ylab("Planted\ngerminants in June")

fig_pl_jn_site <- site_fig %+%
  aes(y = mv_germ_planted_cor_jun) +
  ylab("Planted\ngerminants in June") +
  scale_fill_manual(values = col_pal_lit[2:4])

# Planted Mv germination July
fig_pl_jl <- base_fig %+%
  aes(y = mv_germ_planted_cor_jul) +
  ylab("Planted\ngerminants in July")

fig_pl_jl_site <- site_fig %+%
  aes(y = mv_germ_planted_cor_jul) +
  ylab("Planted\ngerminants in July")

# Mv infection June
fig_pi_jn <- base_fig %+%
  aes(y = mv_prop_infec_jun) +
  ylab("Proportion germinants\nwith lesions in June")

fig_pi_jn_site <- site_fig %+%
  aes(y = mv_prop_infec_jun) +
  ylab("Proportion germinants\nwith lesions in June")

# Mv infection July
fig_pi_jl <- base_fig %+%
  aes(y = mv_prop_infec_jul) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Proportion germinants\nwith lesions in July")

fig_pi_jl_site <- site_fig %+%
  aes(y = mv_prop_infec_jul) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Proportion germinants\nwith lesions in July")

# scatter plots

# June planted/bg vs. soil moisture
fig_bp_jn_soil <- scat_fig %+%
  theme(legend.position = c(0.15, 0.7)) +
  xlab("Proportion soil moisture") +
  ylab("Background and planted\ngerminants in June")

# July planted/bg vs. soil moisture
fig_bp_jl_soil <- scat_fig %+%
  aes(y = mv_germ_planted_bg_jul) +
  xlab("Proportion soil moisture") +
  ylab("Background and planted\ngerminants in July")

# June planted/bg vs. canopy
fig_bp_jn_can <- scat_fig %+%
  aes(x = canopy_cover.prop) +
  xlab("Proportion canopy cover") +
  ylab("Background and planted\ngerminants in June")

# July planted/bg vs. canopy
fig_bp_jl_can <- scat_fig %+%
  aes(x = canopy_cover.prop, y = mv_germ_planted_bg_jul) +
  xlab("Proportion canopy cover") +
  ylab("Background and planted\ngerminants in July")

# combine plots and save
pdf("./output/mv_germination_raw_2018_litter_exp.pdf", width = 6, height = 6)
plot_grid(fig_bg_jn, fig_bg_jl, fig_bp_jn, fig_bp_jl, fig_pl_jn, fig_pl_jl,
          ncol = 2,
          labels = letters[1:6],
          label_size = lg_txt)
dev.off()

pdf("./output/mv_disease_raw_2018_litter_exp.pdf", width = 6, height = 2)
plot_grid(fig_pi_jn, fig_pi_jl,
          ncol = 2,
          labels = letters[1:2],
          label_size = lg_txt) 
dev.off()

pdf("./output/mv_germination_site_raw_2018_litter_exp.pdf", width = 6, height = 6)
plot_grid(fig_bg_jn_site, fig_bg_jl_site, fig_bp_jn_site, fig_bp_jl_site, fig_pl_jn_site, fig_pl_jl_site,
          ncol = 2,
          labels = letters[1:6],
          label_size = lg_txt)
dev.off()

pdf("./output/mv_disease_site_raw_2018_litter_exp.pdf", width = 6, height = 2)
plot_grid(fig_pi_jn_site, fig_pi_jl_site,
          ncol = 2,
          labels = letters[1:2],
          label_size = lg_txt) 
dev.off()

pdf("./output/mv_germination_covariates_raw_2018_litter_exp.pdf", width = 6, height = 4)
plot_grid(fig_bp_jn_soil, fig_bp_jl_soil, fig_bp_jn_can, fig_bp_jl_can,
          ncol = 2,
          labels = letters[1:4],
          label_size = lg_txt)
dev.off()


#### statistical model ####

# Beverton-Holt function
bh_fun <- function(dat_in, a, y0){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$litter_weight.g)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = litter_weight.g, y = mv_germ_planted_cor_jul)) +
          stat_summary(geom = "point", fun = "mean", aes(color = litter_microbes)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = litter_microbes)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# try values
mean(filter(micro, litter_weight.g == 0)$mv_germ_planted_cor_jul)
bh_fun(micro, 0.03, 109)

# model
mv_germ_litter_mod <- brm(data = micro, family = gaussian,
                      bf(mv_germ_planted_cor_jul ~ y0 / (1 + alpha * litter_weight.g),
                         y0 ~ 0 + litter_microbes + (1|site),
                         alpha ~ 0 + litter_microbes,
                         nl = T),
                      prior <- c(prior(normal(109, 100), nlpar = "y0", class = "b", lb = 0, ub = 200),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.999))

# check model and increase chains
summary(mv_germ_litter_mod)
prior_summary(mv_germ_litter_mod)
mv_germ_litter_mod <- update(mv_germ_litter_mod, chains = 3)
plot(mv_germ_litter_mod)
pp_check(mv_germ_litter_mod, nsamples = 100)


#### model fit figure ####

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

# simulate data
sim_dat <- tibble(litter_weight.g = rep(0:200, 2),
                  litter_microbes = rep(c("absent", "present"), each = 201)) %>%
  mutate(mv_germ_planted_cor_jul = fitted(mv_germ_litter_mod, newdata = ., re_formula = NA)[, "Estimate"],
         germ_lower = fitted(mv_germ_litter_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         germ_upper = fitted(mv_germ_litter_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         prop_germ_jul = mv_germ_planted_cor_jul / 200,
         prop_lower = germ_lower / 200,
         prop_upper = germ_upper / 200)

# figure

pdf("output/mv_germinant_analysis_model_fits_2018_litter_exp.pdf")
ggplot(micro, aes(x = litter_weight.g, y = mv_germ_planted_cor_jul)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = litter_microbes)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = litter_microbes)) +
  geom_ribbon(data = sim_dat, alpha = 0.5, aes(ymin = germ_lower, ymax = germ_upper, fill = litter_microbes)) +
  geom_line(data = sim_dat, aes(color = litter_microbes)) +
  xlab("Litter weight (g)") +
  ylab("Planted germinants (out of 200)") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  temp_theme

ggplot(micro, aes(x = litter_weight.g, y = prop_germ_jul)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = litter_microbes)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = litter_microbes)) +
  geom_ribbon(data = sim_dat, alpha = 0.5, aes(ymin = prop_lower, ymax = prop_upper, fill = litter_microbes)) +
  geom_line(data = sim_dat, aes(color = litter_microbes)) +
  xlab("Litter weight (g)") +
  ylab("Proportion of plants that germinated") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  temp_theme
dev.off()

#### output ####
write_csv(micro, "intermediate-data/mv_germination_covariates_2018_litter_exp.csv")
save(mv_germ_litter_mod, file = "output/mv_germination_litter_model_2018_litter_exp.rda")
