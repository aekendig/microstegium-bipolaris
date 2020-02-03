##### info ####

# file: ev_performace_disease_2018_litter_exp
# author: Amy Kendig
# date last edited: 2/3/20
# goal: evaluate the effects of litter treatments and environmental covariates on the performance and disease of Elymus adults


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)

# import data
dat_may <- read_csv("./data/ev_size_may_2018_litter_exp.csv")
dat_jun <- read_csv("./data/ev_size_disease_jun_2018_litter_exp.csv")
dat_jul <- read_csv("./data/ev_size_disease_jul_2018_litter_exp.csv")
dat_sep <- read_csv("./data/ev_disease_sep_2018_litter_exp.csv")
leaf_scans <- read_csv("./intermediate-data/ev_damage_sep_2018_litter_exp.csv")
seed_counts <- read_csv("./intermediate-data/ev_processed_seeds_2018_litter_exp.csv")
cov_plot <- read_csv("./intermediate-data/plot_covariates_2018_litter_exp.csv")
plots <- read_csv("./data/plot_treatments_2018_litter_exp.csv")


#### visualize trends ####

# combine data
# dat_grow <- full_join(dat_may, dat_jun) %>%
#   full_join(dat_jul) %>%
#   full_join(dat_sep) %>%
#   mutate(plot = as.factor(plot),
#          prop_infec = leaves_infec / leaves_tot,
#          date = as.Date(as.character(date), format = "%Y%m%d"))
# 
# # height
# ggplot(dat_grow, aes(x = date, y = height.cm)) +
#   geom_line(aes(color = plot, linetype = site))
# # mostly subtle changes, not necessarily monotonic
# 
# # basal
# ggplot(dat_grow, aes(x = date, y = basal_circ.cm)) +
#   geom_line(aes(color = plot, linetype = site))
# # most are not monotonic
# 
# # tillers
# ggplot(dat_grow, aes(x = date, y = tillers)) +
#   geom_line(aes(color = plot, linetype = site))
# # 2 time points
# 
# # infection
# ggplot(dat_grow, aes(x = date, y = prop_infec)) +
#   geom_line(aes(color = plot, linetype = site))
# # only use July and September


#### edit data ####

# edit May data
unique(dat_may$field_notes)
filter(dat_may, field_notes == "not Elymus virginicus")
unique(dat_may$date)

dat_may2 <- dat_may %>%
  rename(height_may.cm = height.cm,
         basal_circ_may.cm = basal_circ.cm) %>%
  select(-c(date, field_notes))

# edit June data
unique(dat_jun$field_notes)
filter(dat_jun, field_notes == "not Elymus virginicus")
unique(dat_jun$date)

dat_jun2 <- dat_jun %>%
  rename(height_jun.cm = height.cm,
         basal_circ_jun.cm = basal_circ.cm,
         tillers_jun = tillers) %>%
  select(-c(date, infec, leaves_tot, leaves_infec, field_notes))

# edit July data
unique(dat_jul$field_notes)
filter(dat_jul, field_notes == "not Elymus virginicus")
unique(dat_jul$date)

dat_jul2 <- dat_jul %>%
  rename(height_jul.cm = height.cm,
         basal_circ_jul.cm = basal_circ.cm,
         tillers_jul = tillers, 
         leaves_tot_jul = leaves_tot, 
         leaves_infec_jul = leaves_infec) %>%
  select(-c(date, infec, field_notes))

# edit September data
unique(dat_sep$field_notes)
filter(dat_sep, field_notes == "not Elymus virginicus")
unique(dat_sep$date)

dat_sep2 <- dat_sep %>%
  rename(leaves_tot_sep = leaves_tot, 
         leaves_infec_sep = leaves_infec) %>%
  select(-c(date, infec, field_notes))

# combine date data
dates <- tibble(date = c(unique(dat_may$date), unique(dat_jun$date), unique(dat_jul$date), unique(dat_sep$date)),
                month = c("May", "June", "July", "September")) %>%
  mutate(date = as.Date(as.character(date), format = "%Y%m%d"),
         days_may = as.numeric(date - date[1]),
         days_jun = as.numeric(date - date[2]),
         days_jul = as.numeric(date - date[3]))

# edit plot data so that no litter is "sterilized"
# scale litter weight
plots2 <- plots %>%
  mutate(sterilized = case_when(litter == "live" ~ 0,
                                TRUE ~ 1),
         litter_weight.scaled = (litter_weight.g - mean(litter_weight.g)) / sd(litter_weight.g),
         litter_microbes = recode(sterilized, "0" = "present", "1" = "absent"),
         litter_density = fct_relevel(litter_density, "none", "low", "medium"))

# edit covariates
cov_plot2 <- cov_plot %>%
  mutate(litter_density = fct_relevel(litter_density, "none", "low", "medium"))

# combine Ev data
# remove L2 plot 3: not Elymus virginicus
# calculate differences and relative growth rates
# add covariates and plot info
elymus <- full_join(dat_may2, dat_jun2) %>%
  full_join(dat_jul2) %>%
  full_join(dat_sep2) %>%
  filter(!(site == "L2" & plot == 3)) %>%
  mutate(height_change_jun = height_jun.cm - height_may.cm,
         height_change_jul = height_jul.cm - height_jun.cm,
         basal_change_jun = basal_circ_jun.cm - basal_circ_may.cm,
         basal_change_jul = basal_circ_jul.cm - basal_circ_jun.cm,
         tillers_change = tillers_jul - tillers_jun,
         height_rel_change_jun = height_change_jun / dates$days_may[2],
         height_rel_change_jul = height_change_jul / dates$days_jun[3],
         basal_rel_change_jun = basal_change_jun / dates$days_may[2],
         basal_rel_change_jul = basal_change_jul / dates$days_jun[3],
         tillers_rel_change = tillers_change / dates$days_jun[3],
         prop_infec_jul = leaves_infec_jul / leaves_tot_jul,
         prop_infec_sep = leaves_infec_sep / leaves_tot_sep,
         prop_infec_change = prop_infec_sep - prop_infec_jul,
         prop_infec_rel_change = prop_infec_change / dates$days_jul[4]) %>%
  left_join(cov_plot2) %>%
  left_join(plots2)

# edit damage data
# remove L2 plot 3: not Elymus virginicus
# calculate proportion leaf area damaged
# add covariates and plot info
damage <- leaf_scans %>%
  select(site, plot, leaf_area.pix, lesion_area.pix) %>%
  full_join(dat_sep) %>%
  filter(!(site == "L2" & plot == 3)) %>%
  mutate(leaves_infec2 = ifelse(leaves_infec == 0 & !is.na(leaf_area.pix), 1, leaves_infec),
         prop_dam = lesion_area.pix * leaves_infec2 / (leaf_area.pix * leaves_tot),
         prop_dam = case_when(is.na(leaf_area.pix) & leaves_infec == 0 ~ 0,
                              TRUE ~ prop_dam)) %>%
  left_join(cov_plot2) %>%
  left_join(plots2)

# edit seed data
unique(seed_counts$spikelet_notes)
unique(seed_counts$ID_unclear)
sum(is.na(seed_counts$spikelets))
sum(is.na(seed_counts$spikelet_weight.mg))
sum(is.na(seed_counts$seeds))
sum(is.na(seed_counts$seed_weight.mg))

# combine seed data by plant
# add plants with no seeds
# add covariates and plot info
seeds <- seed_counts %>%
  group_by(site, plot, sp, ID) %>%
  summarise(spikelets = sum(spikelets),
            spikelet_weight.mg = sum(spikelet_weight.mg),
            seeds = sum(seeds)) %>%
  full_join(elymus %>%
              select(site, plot, sp, ID)) %>%
  mutate(spikelets = replace_na(spikelets, 0),
         spikelet_weight.mg = replace_na(spikelet_weight.mg, 0),
         seeds = replace_na(seeds, 0)) %>%
  left_join(cov_plot2) %>%
  left_join(plots2)

# combine seeds and damage
seed_dam <- inner_join(seeds, damage)


#### raw data figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# color palettes
col_pal_mic <- c("white", "black")
col_pal_lit <- c("white", "yellow", "orange", "red")

# base summary figure
base_fig <- ggplot(elymus, aes(x = litter_weight.g)) +
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
        legend.position = "none")

# site summary figure
site_fig <- ggplot(elymus, aes(x = site)) +
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
scat_fig <- ggplot(seed_dam, aes(x = prop_dam, y = seeds)) +
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

# June height raw data
fig_ht_jn <- base_fig %+%
  aes(y = height_rel_change_jun) +
  theme(legend.position = c(0.8, 0.3)) +
  ylab(expression(paste("Height change\nsince May (cm ", day^{-1}, ")", sep = "")))

fig_ht_jn_site <- site_fig %+%
  aes(y = height_rel_change_jun) %+%
  theme(legend.position = c(0.8, 0.3)) +
  ylab(expression(paste("Height change\nsince May (cm ", day^{-1}, ")", sep = "")))

# July height raw data
fig_ht_jl <- base_fig %+%
  aes(y = height_rel_change_jul) +
  ylab(expression(paste("Height change\nsince June (cm ", day^{-1}, ")", sep = "")))

fig_ht_jl_site <- site_fig %+%
  aes(y = height_rel_change_jul) +
  ylab(expression(paste("Height change\nsince June (cm ", day^{-1}, ")", sep = "")))

# June basal raw data
fig_bs_jn <- base_fig %+%
  aes(y = basal_rel_change_jun) +
  ylab(expression(paste("Basal circumference\nchange since May (cm ", day^{-1}, ")", sep = "")))

fig_bs_jn_site <- site_fig %+%
  aes(y = basal_rel_change_jun) +
  ylab(expression(paste("Basal circumference\nchange since May (cm ", day^{-1}, ")", sep = "")))

# July basal raw data
fig_bs_jl <- base_fig %+%
  aes(y = basal_rel_change_jul) +
  ylab(expression(paste("Basal circumference\nchange since June (cm ", day^{-1}, ")", sep = "")))

fig_bs_jl_site <- site_fig %+%
  aes(y = basal_rel_change_jul) +
  ylab(expression(paste("Basal circumference\nchange since June (cm ", day^{-1}, ")", sep = "")))

# July tiller raw data
fig_tl <- base_fig %+%
  aes(y = tillers_rel_change) +
  ylab("Change in number of tillers\nper day since June")

fig_tl_site <- site_fig %+%
  aes(y = tillers_rel_change) +
  ylab("Change in number of tillers\nper day since June")

# seeds raw data
fig_sd <- base_fig %+% seeds %+%
  aes(y = seeds) +
  ylab("Number of seeds\nproduced per year")

fig_sd_site <- site_fig %+% seeds %+%
  aes(y = seeds) +
  ylab("Number of seeds\nproduced per year")

# prop infection July
fig_pj <- base_fig %+%
  aes(y = prop_infec_jul) +
  ylab("Proportion leaves with\nlesions in July")

fig_pj_site <- site_fig %+%
  aes(y = prop_infec_jul) +
  ylab("Proportion leaves with\nlesions in July")

# prop infection September
fig_ps <- base_fig %+%
  aes(y = prop_infec_sep) +
  ylab("Proportion leaves with\nlesions in September")

fig_ps_site <- site_fig %+%
  aes(y = prop_infec_sep) +
  ylab("Proportion leaves with\nlesions in September")

# prop infection raw data
fig_pi <- base_fig %+%
  aes(y = prop_infec_rel_change) +
  ylab("Change in proportion leaves\nwith lesions per day")

fig_pi_site <- site_fig %+%
  aes(y = prop_infec_rel_change) +
  ylab("Change in proportion leaves\nwith lesions per day")

# damage raw data
fig_dm <- base_fig %+% damage %+%
  aes(y = prop_dam) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Proportion of leaf area\nwith lesions")

fig_dm_site <- site_fig %+% damage %+%
  aes(y = prop_dam) +
  theme(legend.position = c(0.2, 0.7)) +
  ylab("Proportion of leaf area\nwith lesions")

# scatterplots

# damage and soil moisture
fig_dam_soil <- scat_fig %+%
  aes(x = soil_moisture.prop, y = prop_dam) %+%
  xlab("Proportion soil moisture") +
  ylab("Proportion of leaf area\nwith lesions")

# seeds and soil moisture
fig_seed_soil <- scat_fig %+%
  aes(x = soil_moisture.prop, y = seeds) %+%
  theme(legend.position = c(0.2, 0.7)) +
  xlab("Proportion soil moisture") +
  ylab("Number of seeds")

# damage and canopy
fig_dam_can <- scat_fig %+%
  aes(x = canopy_cover.prop, y = prop_dam) %+%
  xlab("Proportion canopy cover") +
  ylab("Proportion of leaf area\nwith lesions")

# seeds and soil moisture
fig_seed_can <- scat_fig %+%
  aes(x = canopy_cover.prop, y = seeds) %+%
  xlab("Proportion canopy cover") +
  ylab("Number of seeds")

# seeds and damage
fig_seed_dam <- scat_fig +
  xlab("Proportion of leaf area with lesions") +
  ylab("Number of seeds")

# combine plots and save
pdf("./output/ev_performance_raw_2018_litter_exp.pdf", width = 6, height = 6)
plot_grid(fig_ht_jn, fig_ht_jl, fig_bs_jn, fig_bs_jl, fig_tl, fig_sd,
          ncol = 2,
          labels = letters[1:6],
          label_size = lg_txt)
dev.off()

pdf("./output/ev_disease_raw_2018_litter_exp.pdf", width = 6, height = 4)
plot_grid(fig_pj, fig_ps, fig_pi, fig_dm,
          ncol = 2,
          labels = letters[1:4],
          label_size = lg_txt) 
dev.off()

pdf("./output/ev_performance_site_raw_2018_litter_exp.pdf", width = 6, height = 6)
plot_grid(fig_ht_jn_site, fig_ht_jl_site, fig_bs_jn_site, fig_bs_jl_site, fig_tl_site, fig_sd_site,
          ncol = 2,
          labels = letters[1:6],
          label_size = lg_txt)
dev.off()

pdf("./output/ev_disease_site_raw_2018_litter_exp.pdf", width = 6, height = 4)
plot_grid(fig_pj_site, fig_ps_site, fig_pi_site, fig_dm_site,
          ncol = 2,
          labels = letters[1:4],
          label_size = lg_txt) 
dev.off()

pdf("./output/ev_performance_disease_covariates_raw_2018_litter_exp.pdf", width = 6, height = 6)
plot_grid(fig_seed_soil, fig_dam_soil, fig_seed_can, fig_dam_can, fig_seed_dam,
          ncol = 2,
          labels = letters[1:5],
          label_size = lg_txt)
dev.off()


#### output intermediate data ####
write_csv(seed_dam, "./intermediate-data/ev_damage_seeds_covariates_2018_litter_exp.csv")


#### June height models ####

# full model
# mod_ht_jn_1 <- brm(height_change_jun ~ sterilized * litter_weight.scaled + soil_moisture.centered + canopy_cover.centered + (1 |site),
#                    data = elymus,
#                    family = gaussian,
#                    prior <- c(prior(normal(0, 10), class = Intercept),
#                               prior(normal(0, 1), class = b),
#                               prior(cauchy(0, 1), class = sd)),
#                    iter = 6000, warmup = 1000, chains = 3, cores = 2, 
#                    control = list(adapt_delta = 0.9999))
# 
# # evaluate full model
# summary(mod_ht_jn_1)
# pp_check(mod_ht_jn_1, nsamples = 100)
# 
# # simplify model
# mod_ht_jn_2 <- update(mod_ht_jn_1, formula. = height_change_jun ~ sterilized * litter_weight.scaled + soil_moisture.centered + (1 |site))
# mod_ht_jn_3 <- update(mod_ht_jn_1, formula. = height_change_jun ~ sterilized * litter_weight.scaled + canopy_cover.centered + (1 |site))
# mod_ht_jn_4 <- update(mod_ht_jn_1, formula. = height_change_jun ~ sterilized * litter_weight.scaled + (1 |site))
# 
# # compare models
# loo_ht_jn <- list(loo(mod_ht_jn_1, reloo = T), loo(mod_ht_jn_2, reloo = T), loo(mod_ht_jn_3, reloo = T), loo(mod_ht_jn_4, reloo = T))
# loo::loo_compare(loo_ht_jn)
# # all essentially equal
# 
# # evaluate simplest model
# summary(mod_ht_jn_4)
# plot(mod_ht_jn_4)
# pp_check(mod_ht_jn_4, nsamples = 100)
