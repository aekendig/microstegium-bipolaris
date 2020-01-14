##### info ####

# file: mv_germination_disease_2018_litter_exp
# author: Amy Kendig
# date last edited: 1/14/20
# goal: evaluate the effects of litter treatments and environmental covariates on the germination and disease of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

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
micro <- full_join(select(dat_jun_planted, -plot), select(dat_jun_bg, -plot)) %>%
  full_join(select(dat_jul2, -plot)) %>%
  mutate(mv_germ_planted_jun = mv_germ_planted_bg_jun - mv_germ_bg_jun,
         mv_germ_planted_jul = mv_germ_planted_bg_jul - mv_germ_bg_jul,
         mv_prop_infec_jun = mv_infec_jun / mv_germ_planted_bg_jun,
         mv_prop_infec_jul = mv_infec_jul / mv_germ_planted_bg_jul)

# to combine with covariates, each dataset needs to be used separately because of the different plot numbers
dat_jun_planted2 <- dat_jun_planted %>%
  left_join(cov_plot)

dat_jun_bg2 <- dat_jun_bg %>%
  left_join(cov_plot)

dat_jul3 <- dat_jul2 %>%
  left_join(cov_plot)


#### raw data figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# microbes color palette
col_pal_mic <- c("white", "black")

# base summary figure
base_fig <- ggplot(micro, aes(x = litter_weight.g)) +
  stat_summary(aes(fill = litter_microbes), geom = "point", fun.y = mean, position = position_dodge(15), shape = 21, size = 3) +
  stat_summary(aes(group = litter_microbes), geom = "errorbar", fun.data = mean_cl_boot, width = 5, position = position_dodge(15)) +
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

# Background Mv germination June
fig_bg_jn <- base_fig %+%
  aes(y = mv_germ_bg_jun) +
  ylab("Background\ngerminants in June")

# Background Mv germination July
fig_bg_jl <- base_fig %+%
  aes(y = mv_germ_bg_jul) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Background\ngerminants in July")

# Background + planted Mv germination June
fig_bp_jn <- base_fig %+%
  aes(y = mv_germ_planted_bg_jun) +
  ylab("Background and planted\ngerminants in June")

# Background + planted Mv germination July
fig_bp_jl <- base_fig %+%
  aes(y = mv_germ_planted_bg_jul) +
  ylab("Background and planted\ngerminants in July")

# Planted Mv germination June
fig_pl_jn <- base_fig %+%
  aes(y = mv_germ_planted_jun) +
  ylab("Planted\ngerminants in June")

# Planted Mv germination July
fig_pl_jl <- base_fig %+%
  aes(y = mv_germ_planted_jul) +
  ylab("Planted\ngerminants in July")

# Mv infection June
fig_pi_jn <- base_fig %+%
  aes(y = mv_prop_infec_jun) +
  ylab("Proportion germinants\nwith lesions in June")

# Mv infection July
fig_pi_jl <- base_fig %+%
  aes(y = mv_prop_infec_jul) +
  theme(legend.position = c(0.8, 0.7)) +
  ylab("Proportion germinants\nwith lesions in July")


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
