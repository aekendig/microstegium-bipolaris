##### info ####

# file: alpha_beta_correlation_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/3/21
# goal: correlations between transmission effects and competitive effects


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
alphas <- read_csv("output/focal_growth_pairwise_coefficients_2018_2019_density_exp.csv")
betas <- read_csv("output/plot_severity_pairwise_coefficients_2018_2019_density_exp.csv")


#### edit data ####

# combine data
dat <- alphas %>%
  rename(sig_alpha = sig) %>%
  select(year, treatment, focal, background, alpha, sig_alpha) %>%
  full_join(betas %>%
              rename(focal = foc_sp_age,
                     background = bg_sp_age,
                     sig_beta = sig) %>%
              select(year, treatment, focal, background, beta, sig_beta)) %>%
  mutate(Significance = case_when(sig_alpha == "omits 0" & sig_beta == "omits 0" ~ "both",
                                  sig_alpha == "omits 0" & sig_beta != "omits 0" ~ "alpha only",
                                  sig_alpha != "omits 0" & sig_beta == "omits 0" ~ "beta only",
                                  TRUE ~ "neither") %>%
           as.factor())

# initial figure
ggplot(dat, aes(beta, alpha)) +
  geom_point(aes(shape = Significance)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ treatment)

ggplot(dat, aes(beta, alpha)) +
  geom_point(aes(shape = Significance)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(background ~ treatment)


#### correlation ####

# divide data
ctrl_dat <- dat %>%
  filter(treatment == "control (water)")
fung_dat <- dat %>%
  filter(treatment == "fungicide")

# correlations
cor.test(~ alpha + beta, data = ctrl_dat)
cor.test(~ alpha + beta, data = fung_dat)

# divid by species
mv_ctrl_dat <- ctrl_dat %>%
  filter(background == "Mv")

# correlations
cor.test(~ alpha + beta, data = mv_ctrl_dat)


#### figure ####

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.2, 'cm'))

col_pal = c("black", "#238A8DFF")
shape_pal = c(16, 17, 15, 3)

ctrl_fig <- ggplot(ctrl_dat, aes(beta, alpha)) +
  geom_point(aes(shape = Significance, color = treatment)) +
  geom_text(x = -0.505, y = 0.064, label = "Control (water)", hjust = 0, check_overlap = T, size = 3) +
  scale_color_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = shape_pal) +
  ylab(expression(paste("Competition coefficient (", alpha, ")", paste = ""))) +
  xlab(expression(paste("Transmission coefficient (", beta, ")", paste = ""))) +
  fig_theme +
  theme(legend.position = c(0.8, 0.2),
        plot.margin = unit(c(3, 7, 3, 3), "pt")) +
  coord_cartesian(ylim = c(-0.16, 0.06))

fung_fig <- ggplot(fung_dat, aes(beta, alpha)) +
  geom_point(aes(shape = Significance, color = treatment)) +
  geom_text(x = -0.75, y = 0.08, label = "Fungicide", hjust = 0, check_overlap = T, size = 3) +
  scale_color_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal[-3]) +
  ylab(expression(paste("Competition coefficient (", alpha, ")", paste = ""))) +
  xlab(expression(paste("Transmission coefficient (", beta, ")", paste = ""))) +
  fig_theme +
  theme(axis.title.y = element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt"))

pdf("output/alpha_beta_correlation_2018_2019_density_exp.pdf", width = 5, height = 2.5)
plot_grid(ctrl_fig, fung_fig, 
          nrow = 1,
          rel_widths = c(1, 0.9),
          labels = LETTERS[1:2])
dev.off()


#### text ####

# We found no significant correlation between competition coefficients and transmission coefficients for control (r = 0.19, t = 0.79, df = 16, P = 0.44) or fungicide (r = -0.22, t = -0.91, df = 16, P = 0.38) conditions (Fig. 4), suggesting that the negative effects of M. vimineum under control conditions (Fig. 3A) are due to disease-altered competitive effects rather than disease transmission.