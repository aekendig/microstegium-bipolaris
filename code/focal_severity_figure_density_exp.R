##### info ####

# file: focal_severity_figure_density_exp
# author: Amy Kendig
# date last edited: 7/27/20
# goal: figure of severity model results


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(glmmTMB)

# import models
load("output/average_humidity_total_biomass_model_2019_density_exp.rda")
load("output/mv_severity_humidity_model_2019_density_exp.rda")
load("output/Mv_all_biomass_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/Ev_all_biomass_model_greenhouse_fungicide_2019_density_exp.rda")
# load("output/Mv_severity_humidity_model_2019_density_exp.rda")
# load("output/Ev_severity_biomass_edge_model_2019_density_exp.rda")

# import data
hum_bio_dat <- read_csv("intermediate-data/average_humidity_total_biomass_2019_density_exp.csv")
sep_mv_water_dat <- read_csv("intermediate-data/mv_severity_water_treatment_sep_2019_density_exp.csv")
mv_all_bio_dat <- read_csv("intermediate-data/mv_all_biomass_data_2019_density_exp.csv")
ev_all_bio_dat <- read_csv("intermediate-data/ev_all_biomass_data_2019_density_exp.csv")
# mv_bio_edge_dat <- read_csv("intermediate-data/Mv_severity_june_late_aug_2019_density_exp.csv")
# ev_bio_edge_dat <- read_csv("intermediate-data/Ev_severity_june_late_aug_2019_density_exp.csv")

# examine models
summary(havg_totb_mod)
summary(mv_severity_humidity_mod)
summary(mv_all_bio_fung_mod)
summary(ev_all_bio_fung_mod)
# summary(mv_severity_biomass_edge_mod)
# summary(ev_severity_biomass_edge_mod)


#### figure settings ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.title = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.box.margin = margin(-10, -10, -10, -10))

# colors
col_pal = c("#C0A76D", "#55A48B")

# scale function
scale_fun <- function(x){
  y = 100 * x
  sprintf("%.0f", y)
}


#### humidity total biomass ####

# subset data
hum_bio_dat2 <- hum_bio_dat %>%
  filter(month_name == "September") %>%
  mutate(Treatment = recode(treatment, "water" = "control (water)"))

# min and max values
hum_bio_range <- hum_bio_dat2 %>%
  summarise(min = min(total_biomass.g, na.rm = T),
            max = max(total_biomass.g, na.rm = T))

# prediction data
hum_bio_sim_dat <- tibble(total_biomass.g = seq(hum_bio_range$min, hum_bio_range$max, length.out = 200)) %>%
  mutate(month_name = "September",
         site = NA,
         fungicide = 0,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(havg_pred = predict(havg_totb_mod, newdata = ., re.form = NA, type = "response"),
         havg_pred_se = predict(havg_totb_mod, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)

# plot
hum_bio_plot <- ggplot(hum_bio_sim_dat, aes(x = total_biomass.g, y = havg_pred)) +
  geom_point(data = hum_bio_dat2, alpha = 0.3, aes(x = total_biomass.g, y = hum_avg, color = Treatment)) +
geom_ribbon(alpha = 0.5, aes(ymin = havg_pred - havg_pred_se, ymax = havg_pred + havg_pred_se, fill = Treatment)) +
  geom_line(aes(color = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab("Relative humidity (%)") +
  xlab(expression(paste("Total biomass (g ", m^-2, ")", sep = ""))) +
  scale_y_continuous(labels = scale_fun) +
  temp_theme +
  theme(legend.position = "none")


#### Mv severity and humidity ####

# min and max values
mv_hum_range <- sep_mv_water_dat %>%
  summarise(min = min(hum_avg, na.rm = T),
            max = max(hum_avg, na.rm = T))

# prediction data
mv_hum_sim_dat <- tibble(hum_avg = seq(mv_hum_range$min, mv_hum_range$max, length.out = 200)) %>%
  mutate(site = NA,
         plot = NA,
         fungicide = 0,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(severity = predict(mv_severity_humidity_mod, newdata = ., type = "response", re.form = NA),
         severity_se = predict(mv_severity_humidity_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)

# plot
mv_hum_plot <- ggplot(mv_hum_sim_dat, aes(hum_avg, severity)) +
  geom_point(data = sep_mv_water_dat, alpha = 0.3, aes(x = hum_avg, y = plant_severity_adjusted, color = Treatment)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Microstegium), " disease severity (%)", sep = ""))) +
  xlab("Relative humidity (%)") +
  scale_y_continuous(labels = scale_fun) +
  scale_x_continuous(labels = scale_fun) +
  temp_theme +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 6.5, color = "black", vjust = -1.5),
        axis.title.x = element_text(size = 7, color = "black"))


#### Mv biomass and density ####

# prediction data
mv_dens_sim_dat <- tibble(Mv_density = rep(seq(0, 67, length.out = 200), 2),
                          fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Mv_biomass.g = fitted(mv_all_bio_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         Mv_biomass_lower = fitted(mv_all_bio_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         Mv_biomass_upper = fitted(mv_all_bio_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) 

# figure
mv_dens_plot <- ggplot(mv_all_bio_dat, aes(x = Mv_density, y = Mv_biomass.g)) + 
  geom_ribbon(data = mv_dens_sim_dat, aes(ymin = Mv_biomass_lower, ymax = Mv_biomass_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_dens_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegium"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegium"), " biomass (g ", m^-1, ")", sep = ""))) +
  temp_theme +
  theme(legend.position = "none")

#### Ev biomass and density ####

# prediction data
ev_dens_sim_dat <- tibble(Ev_density = rep(seq(0, 20, length.out = 200), 2),
                          fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Ev_biomass.g = fitted(ev_all_bio_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         Ev_biomass_lower = fitted(ev_all_bio_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         Ev_biomass_upper = fitted(ev_all_bio_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) 

# figure
ev_dens_plot <- ggplot(ev_all_bio_dat, aes(x = Ev_density, y = Ev_biomass.g)) + 
  geom_ribbon(data = ev_dens_sim_dat, aes(ymin = Ev_biomass_lower, ymax = Ev_biomass_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = ev_dens_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Elymus"), " biomass (g ", m^-1, ")", sep = ""))) +
  temp_theme + 
  theme(legend.position = "bottom",
        legend.direction = "horizontal")


# #### Mv edge severity ####
# 
# # min and max values
# mv_edge_range <- mv_bio_edge_dat %>%
#   filter(month == "late_aug") %>%
#   summarise(min = min(edge_severity, na.rm = T),
#             max = max(edge_severity, na.rm = T),
#             avg = mean(total_biomass.g, na.rm = T))
# 
# # prediction data
# mv_edge_sim_dat <- tibble(edge_severity = rep(seq(mv_edge_range$min, mv_edge_range$max, length.out = 200), 2),
#                           fungicide = rep(c(0, 1), each = 200)) %>%
#   mutate(Month = "Late August",
#          total_biomass.g = mv_edge_range$avg,
#          site = NA,
#          exp_plot = NA,
#          Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
#   mutate(severity = predict(mv_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA) * 100,
#          severity_se = predict(mv_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit * 100,
#          lesions = edge_severity * 100)
# 
# # plot
# mv_edge_plot <- ggplot(mv_edge_sim_dat, aes(lesions, severity)) +
#   geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
#   geom_line(aes(color = Treatment)) +
#   #annotate(geom = "text", label = "Late August", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
#   scale_color_manual(values = col_pal) +
#   scale_fill_manual(values = col_pal) +
#   ylab(expression(paste(italic(Microstegium), " lesions (% leaf surface)", sep = ""))) +
#   xlab("Background disease (% leaf surface)") +
#   temp_theme
# 
# 
# #### Ev edge severity ####
# 
# # min and max values
# ev_edge_range <- ev_bio_edge_dat %>%
#   filter(month == "jul") %>%
#   summarise(min = min(edge_severity, na.rm = T),
#             max = max(edge_severity, na.rm = T),
#             avg_mv = mean(Mv_biomass.g, na.rm = T),
#             avg_ev = mean(Ev_biomass.g, na.rm = T))
# 
# # prediction data
# ev_edge_sim_dat <- tibble(edge_severity = rep(seq(ev_edge_range$min, ev_edge_range$max, length.out = 200), 2),
#                           fungicide = rep(c(0, 1), each = 200)) %>%
#   mutate(Month = "July",
#          Mv_biomass.g = ev_edge_range$avg_mv,
#          Ev_biomass.g = ev_edge_range$avg_ev,
#          site = NA,
#          exp_plot = NA,
#          Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
#   mutate(severity = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA),
#          severity_se = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)
# 
# # plot
# ev_edge_plot <- ggplot(ev_edge_sim_dat, aes(edge_severity, severity)) +
#   geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
#   geom_line(aes(color = Treatment)) +
#   #annotate(geom = "text", label = "July", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
#   scale_color_manual(values = col_pal) +
#   scale_fill_manual(values = col_pal) +
#   ylab(expression(paste(italic(Elymus), " lesions", sep = ""))) +
#   xlab("Background disease") +
#   temp_theme +
#   theme(legend.position = "none")


#### combine plots ####

# legend
leg <- get_legend(ev_dens_plot)

# plots
plot_comb <- plot_grid(hum_bio_plot, mv_hum_plot, mv_dens_plot, ev_dens_plot + theme(legend.position = "none"),
                       nrow = 2,
                       labels = c(LETTERS[1:4], ""), label_size = 10,
                       vjust = c(1.1, 1.1, 1.5, 1.5),
                       hjust = c(-0.25, -0.1, -0.5, -0.1))

# combine
pdf("output/focal_severity_figure_density_exp.pdf",
    width = 4.3, height = 4.3)
plot_grid(plot_comb, leg,
          nrow = 2,
          rel_heights = c(1, 0.06))
dev.off()