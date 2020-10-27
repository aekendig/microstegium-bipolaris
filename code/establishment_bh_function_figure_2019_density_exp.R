##### info ####

# file: establishment_bh_function_figure_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: figure of establishment with litter


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library("tidyverse")
library("brms")
library("tidybayes")
library("cowplot")

# number of iterations
n_samps <- 1000

# call scripts
source("code/microstegium_establishment_bh_parameter_2019_density_exp.R")
source("code/elymus_seedling_establishment_bh_parameter_2019_density_exp.R")


#### simulations ####

# initiate data frames
est_df <- tibble(litter = seq(0, 628.82, length.out = 800)) %>%
  expand_grid(tibble(iter = 1:n_samps)) %>%
  mutate(mv_ndis = NA,
         ev_ndis = NA,
         mv_ydis = NA,
         ev_ydis = NA)

est_long_df <- tibble(litter = seq(0, 10000, length.out = 8000)) %>%
  expand_grid(tibble(iter = 1:n_samps)) %>%
  mutate(mv_ndis = NA,
         ev_ndis = NA,
         mv_ydis = NA,
         ev_ydis = NA)

# add establishment estimates
est_df2 <- est_df %>%
  rowwise() %>%
  mutate(mv_ndis = E_A_bh_fun(0, litter, iter),
         ev_ndis = E_S_bh_fun(0, litter, iter),
         mv_ydis = E_A_bh_fun(1, litter, iter),
         ev_ydis = E_S_bh_fun(1, litter, iter))

est_long_df2 <- est_long_df %>%
  rowwise() %>%
  mutate(mv_ndis = E_A_bh_fun(0, litter, iter),
         ev_ndis = E_S_bh_fun(0, litter, iter),
         mv_ydis = E_A_bh_fun(1, litter, iter),
         ev_ydis = E_S_bh_fun(1, litter, iter))


#### process output ####

# make long
est_df3 <- est_df2 %>%
  pivot_longer(cols = c(mv_ndis, ev_ndis, mv_ydis, ev_ydis), names_to = "sp_dis", values_to = "establishment") %>%
  mutate(species = substring(sp_dis, 1, 2) %>%
           recode("ev" = "E. virginicus", 
                  "mv" = "M. vimineum"),
         Disease = substring(sp_dis, 4, 7) %>%
           recode("ndis" = "without disease",
                  "ydis" = "with disease"),
         range = "field")

est_long_df3 <- est_long_df2 %>%
  pivot_longer(cols = c(mv_ndis, ev_ndis, mv_ydis, ev_ydis), names_to = "sp_dis", values_to = "establishment") %>%
  mutate(species = substring(sp_dis, 1, 2) %>%
           recode("ev" = "E. virginicus", 
                  "mv" = "M. vimineum"),
         Disease = substring(sp_dis, 4, 7) %>%
           recode("ndis" = "without disease",
                  "ydis" = "with disease"),
         range = "extended")

# summarise for plot
est_sum <- est_df3 %>%
  group_by(range, species, Disease, litter) %>%
  mean_qi(establishment) %>%
  ungroup()

est_long_sum <- est_long_df3 %>%
  group_by(range, species, Disease, litter) %>%
  mean_qi(establishment) %>%
  ungroup()

plot_dat <- full_join(est_sum, est_long_sum) %>%
  mutate(species = as.factor(species) %>% fct_relevel("M. vimineum"),
         Disease = as.factor(Disease) %>% fct_relevel("without disease"))


#### visualize ####

# figure theme
fig_theme <-   theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 10),
        strip.background = element_blank())

# figure
ext_fig <- ggplot(filter(plot_dat, range == "extended"), aes(x = litter, y = establishment, color = Disease)) +
  geom_line(size = 1.5) +
  geom_line(aes(y = .lower), linetype = "dotted") +
  geom_line(aes(y = .upper), linetype = "dotted") +
  facet_wrap(~ species, scales = "free") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  ylab("Proportion established") +
  xlab(expression(paste("Litter (g ", m^-2, ")", sep = ""))) +
  fig_theme +
  theme(strip.text = element_text(size = 10, face = "italic"),
        axis.title.x = element_blank(),
        legend.position = "none")

fld_fig <- ggplot(filter(plot_dat, range == "field"), aes(x = litter, y = establishment, color = Disease)) +
  geom_line(size = 1.5) +
  geom_line(aes(y = .lower), linetype = "dotted") +
  geom_line(aes(y = .upper), linetype = "dotted") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~ species, scales = "free") +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  ylab("Proportion established") +
  xlab(expression(paste("Litter (g ", m^-2, ")", sep = ""))) +
  fig_theme +
  theme(strip.text = element_blank(),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 0, 0))

pdf("output/establishment_bh_function_figure_2019_density_exp.pdf", width = 5, height = 5)
plot_grid(ext_fig, fld_fig,
          rel_heights = c(0.9, 1),
          nrow = 2,
          label_size = 10,
          labels = LETTERS[1:2])
dev.off()