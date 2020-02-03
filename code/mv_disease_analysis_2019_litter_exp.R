##### info ####

# file: mv_disease_analysis_2018_litter_exp
# author: Amy Kendig
# date last edited: 1/26/20
# goal: evaluate the effects of litter treatments on the and disease of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
germ_jun <- read_csv("./data/both_germination_disease_jun_2019_litter_exp.csv")
germ_jul <- read_csv("./data/both_germination_disease_jul_2019_litter_exp.csv")
leaf_jun <- read_csv("./data/mv_leaf_disease_jun_2019_litter_exp.csv")
leaf_jul <- read_csv("./data/mv_leaf_disease_jul_2019_litter_exp.csv")
area_jun <- read_csv("./intermediate-data/mv_damage_jun_2019_litter_exp.csv")
area_jul <- read_csv("./intermediate-data/mv_damage_jul_2019_litter_exp.csv")
bip <- read_csv("./data/mv_bipolaris_id_2019_litter_exp.csv")
plots <- read_csv("./data/litter_weight_apr_2019_litter_exp.csv")
# add in July leaf scans when completed


#### edit data ####

# edit plot data
# put the removed litter into the addition litter
# scale litter weight
# remove unnecessary variables
plots2 <- plots %>%
  spread(key = treatment, value = litter_weight.lb) %>%
  mutate(addition = addition + removal,
         removal = 0) %>%
  gather(key = "treatment", value = "litter_weight.lb", -c(date, site, block)) %>%
  mutate(litter_weight.scaled = (litter_weight.lb - mean(litter_weight.lb)) / sd(litter_weight.lb),
         litter_weight.g = litter_weight.lb * 453.592,
         treatment = fct_relevel(treatment, "removal", "control"),
         plot = as.factor(paste(site, block, sep = "_"))) %>%
  select(-c(date))

# June germination data
germ_jun2 <- germ_jun %>%
  select(-c(date, flag_color, ev_germ, ev_infec)) %>%
  rename(mv_germ_jun = mv_germ,
         mv_infec_jun = mv_infec)

# July germination data
germ_jul2 <- germ_jul %>%
  select(-c(date, flag_color)) %>%
  rename(mv_germ_jul = mv_germ,
         mv_infec_jul = mv_infec)

# June leaf data
leaf_jun2 <- leaf_jun %>%
  select(-c(date, flag_color)) %>%
  rename(leaves_tot_jun = leaves_tot,
         leaves_infec_jun = leaves_infec)

# July leaf data
leaf_jul2 <- leaf_jul %>%
  select(-c(date, flag_color)) %>%
  rename(leaves_tot_jul = leaves_tot,
         leaves_infec_jul = leaves_infec)

# June area data
area_jun2 <- area_jun %>%
  mutate(prop_dam_jun = lesion_area.pix / leaf_area.pix) %>%
  select(-c(date, leaf_objects, lesion_objects, leaf_area.pix, lesion_area.pix))

# July area data
# checked that all reps have a scan and leaf counts
area_jul2 <- area_jul %>%
  full_join(filter(leaf_jul2, scan == 1) %>%
              select(-c(scan, field_notes))) %>%
  mutate(leaves_infec_jul2 = ifelse(leaves_infec_jul == 0 & !is.na(leaf_area.pix), 1, leaves_infec_jul),
         prop_dam_jul = lesion_area.pix * leaves_infec_jul2 / (leaf_area.pix * leaves_tot_jul),
         prop_dam_jul = case_when(is.na(leaf_area.pix) & leaves_infec_jul == 0 ~ 0,
                              TRUE ~ prop_dam_jul)) %>%
  select(-c(date, leaf_area.pix:leaves_infec_jul2))
  

# combine data
germ <- full_join(germ_jun2, germ_jul2) %>%
  full_join(plots2) %>%
  mutate(mv_prop_infec_jun = mv_infec_jun / mv_germ_jun,
         mv_prop_infec_jul = mv_infec_jul / mv_germ_jul,
         treatment = fct_relevel(treatment, "removal", "control"))

leaf <- full_join(leaf_jun2, leaf_jul2) %>%
  full_join(plots2) %>%
  mutate(leaves_prop_infec_jun = leaves_infec_jun / leaves_tot_jun,
         leaves_prop_infec_jul = leaves_infec_jul / leaves_tot_jul,
         treatment = fct_relevel(treatment, "removal", "control"))

damage <- full_join(area_jun2, area_jul2) %>%
  full_join(plots2) %>%
  mutate(treatment = fct_relevel(treatment, "removal", "control"))

bipolaris <- full_join(bip, plots2) %>%
  mutate(prop_eyespots =  leaves_eyespots / leaves_tot,
         prop_bip = leaves_bip / leaves_tot,
         treatment = fct_relevel(treatment, "removal", "control"))


#### raw data figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# color palette
col_pal = c("blue", "purple", "red")

# base theme
base_theme <- theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none")

# base figures
base_fig_trt <- ggplot(germ, aes(x = treatment)) +
  stat_summary(aes(group = treatment), geom = "errorbar", fun.data = mean_cl_boot, width = 0.2) +
  stat_summary(aes(fill = treatment), geom = "point", fun.y = mean, size = 3, shape = 21) +
  xlab("Litter treatment") +
  base_theme +
  scale_fill_manual(values = col_pal, name = "Treatment")

base_fig_lb <- ggplot(germ, aes(x = litter_weight.g, color = treatment)) +
  geom_point(size = 2, aes(shape = plot)) +
  xlab("Litter weight (g)") +
  base_theme +
  scale_shape_manual(values = c(19, 15, 17, 18), name = "Plot") +
  scale_color_manual(values = col_pal, name = "Treatment")

# litter treatments
fig_trt <- ggplot(plots2, aes(x = treatment, y = litter_weight.g)) +
  geom_line(aes(size = treatment, group = plot), color = "black") +
  geom_point(aes(color = treatment, shape = plot), size = 3) +
  base_theme +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_shape_manual(values = c(19, 15, 17, 18), name = "Block") +
  scale_size_manual(values = c(0.5, 0.5, 0.5), guide = F) +
  xlab("Litter treatment") +
  ylab("Litter weight (g)")

# Bipolaris
fig_bi_lb <- base_fig_lb %+% bipolaris %+%
  aes(y = prop_bip) +
  ylab("Prop. of leaves\nwith B. megaspora")

fig_bi_trt <- base_fig_trt %+% bipolaris %+%
  aes(y = prop_bip) +
  ylab("Prop. of leaves\nwith B. megaspora")

# Mv proportion infected June
fig_pi_lb_jn <- base_fig_lb %+% 
  aes(y = mv_prop_infec_jun) +
  ylab("Prop. germinants\nwith lesions in June")

fig_pi_trt_jn <- base_fig_trt %+%
  aes(y = mv_prop_infec_jun) +
  ylab("Prop. germinants\nwith lesions in June")

# Mv proportion leaves infected June
fig_li_lb_jn <- base_fig_lb %+% leaf +
  aes(y = leaves_prop_infec_jun) +
  ylab("Prop. leaves\nwith lesions in June")

fig_li_trt_jn <- base_fig_trt %+% leaf +
  aes(y = leaves_prop_infec_jun) +
  ylab("Prop. leaves\nwith lesions in June")

# Mv proportion infected July
fig_pi_lb_jl <- base_fig_lb %+%
  aes(y = mv_prop_infec_jul) +
  ylab("Prop. germinants\nwith lesions in July")

fig_pi_trt_jl <- base_fig_trt %+%
  aes(y = mv_prop_infec_jul) +
  ylab("Prop. germinants\nwith lesions in July")

# Mv proportion leaves infected July
fig_li_lb_jl <- base_fig_lb %+% leaf +
  aes(y = leaves_prop_infec_jul) +
  ylab("Prop. leaves\nwith lesions in July")

fig_li_trt_jl <- base_fig_trt %+% leaf +
  aes(y = leaves_prop_infec_jul) +
  ylab("Prop. leaves\nwith lesions in July")

# Damage in June
fig_dm_lb_jn <- base_fig_lb %+% damage +
  aes(y = prop_dam_jun) +
  ylab("Prop. leaf area\nwith lesions in June")

fig_dm_trt_jn <- base_fig_trt %+% damage +
  aes(y = prop_dam_jun) +
  ylab("Prop. leaf area\nwith lesions in June")

# Damage in July
fig_dm_lb_jl <- base_fig_lb %+% damage +
  aes(y = prop_dam_jul) +
  ylab("Prop. leaf area\nwith lesions in July")

fig_dm_trt_jl <- base_fig_trt %+% damage +
  aes(y = prop_dam_jul) +
  ylab("Prop. leaf area\nwith lesions in July")

# Block effect: Mv proportion infected July
fig_pi_bl_jl <- base_fig_trt %+%
  aes(x = plot, y = mv_prop_infec_jul) +
  ylab("Prop. germinants\nwith lesions in July") +
  xlab("Block")

# Block effect: Mv proportion infected July
fig_li_bl_jl <- base_fig_trt %+% leaf %+%
  aes(x = plot, y = leaves_prop_infec_jul) +
  ylab("Prop. leaves\nwith lesions in July") +
  xlab("Block")

# Block effect: Mv proportion leaf area July
fig_dm_bl_jl <- base_fig_trt %+% damage %+%
  aes(x = plot, y = prop_dam_jul) +
  ylab("Prop. leaf area\nwith lesions in July") +
  xlab("Block")

# legend
leg <- get_legend(fig_trt +
                 theme(legend.position = "bottom", 
                       legend.direction = "horizontal",
                       legend.box.margin = margin(0, 0, 0, 0)))

leg_block <- get_legend(fig_pi_bl_jl +
                    theme(legend.position = "bottom", 
                          legend.direction = "horizontal",
                          legend.box.margin = margin(0, 0, 0, 0)))

# combine plots and save
g_plots <- plot_grid(fig_trt, fig_bi_lb,
          fig_pi_lb_jn, fig_pi_lb_jl,
          fig_li_lb_jn, fig_li_lb_jl,
          fig_dm_lb_jn, fig_dm_lb_jl,
          ncol = 2,
          labels = letters[1:8],
          label_size = lg_txt) 

trt_plots <- plot_grid(fig_trt, fig_bi_trt,
                      fig_pi_trt_jn, fig_pi_trt_jl,
                      fig_li_trt_jn, fig_li_trt_jl,
                      fig_dm_trt_jn, fig_dm_trt_jl,
                      ncol = 2,
                      labels = letters[1:8],
                      label_size = lg_txt)

block_plots <- plot_grid(fig_pi_bl_jl, fig_li_bl_jl, fig_dm_bl_jl,
                         nrow = 1,
                         labels = letters[1:3],
                         label_size = lg_txt)

pdf("./output/mv_disease_g_raw_2019_litter_exp.pdf", width = 6, height = 6)
plot_grid(g_plots, leg, 
          nrow = 2,
          rel_heights = c(1, 0.1))
dev.off()

pdf("./output/mv_disease_trt_raw_2019_litter_exp.pdf", width = 6, height = 6)
plot_grid(trt_plots, leg, 
          nrow = 2,
          rel_heights = c(1, 0.1))
dev.off()

pdf("./output/mv_disease_block_raw_2019_litter_exp.pdf", width = 6, height = 2)
plot_grid(block_plots, leg_block, 
          nrow = 2,
          rel_heights = c(1, 0.15))
dev.off()