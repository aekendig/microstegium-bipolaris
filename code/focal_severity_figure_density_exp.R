##### info ####

# file: focal_severity_figure_density_exp
# author: Amy Kendig
# date last edited: 7/15/20
# goal: figure of severity model results


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(glmmTMB)

# import models
load("output/Mv_severity_biomass_edge_model_2019_density_exp.rda")
load("output/Mv_severity_humidity_model_2019_density_exp.rda")
load("output/average_humidity_total_biomass_model_2019_density_exp.rda")
load("output/Ev_severity_biomass_edge_model_2019_density_exp.rda")

# import data
mv_bio_edge_dat <- read_csv("intermediate-data/Mv_severity_june_late_aug_2019_density_exp.csv")
sep_mv_water_dat <- read_csv("intermediate-data/Mv_severity_water_treatment_sep_2019_density_exp.csv")
hum_bio_dat <- read_csv("intermediate-data/average_humidity_total_biomass_2019_density_exp.csv")
ev_bio_edge_dat <- read_csv("intermediate-data/Ev_severity_june_late_aug_2019_density_exp.csv")

# examine models
summary(mv_severity_biomass_edge_mod)
summary(mv_severity_humidity_mod)
summary(havg_totb_mod)
summary(ev_severity_biomass_edge_mod)


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
col_pal = c("#a6611a", "#018571")


#### Mv edge severity ####

# min and max values
mv_edge_range <- mv_bio_edge_dat %>%
  filter(month == "late_aug") %>%
  summarise(min = min(edge_severity, na.rm = T),
            max = max(edge_severity, na.rm = T),
            avg = mean(total_biomass.g, na.rm = T))

# prediction data
mv_edge_sim_dat <- tibble(edge_severity = rep(seq(mv_edge_range$min, mv_edge_range$max, length.out = 200), 2),
                          fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Month = "Late August",
         total_biomass.g = mv_edge_range$avg,
         site = NA,
         exp_plot = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(severity = predict(mv_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA) * 100,
         severity_se = predict(mv_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit * 100,
         lesions = edge_severity * 100)

# plot
mv_edge_plot <- ggplot(mv_edge_sim_dat, aes(lesions, severity)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "Late August", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Microstegium), " lesions (% leaf surface)", sep = ""))) +
  xlab("Background disease (% leaf surface)") +
  temp_theme


#### Mv total biomass ####

# min and max values
mv_bio_range <- mv_bio_edge_dat %>%
  filter(month == "late_aug") %>%
  summarise(min = min(total_biomass.g, na.rm = T),
            max = max(total_biomass.g, na.rm = T),
            avg = mean(edge_severity, na.rm = T))

# prediction data
mv_bio_sim_dat <- tibble(total_biomass.g = rep(seq(mv_bio_range$min, mv_bio_range$max, length.out = 200), 2),
                          fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Month = "Late August",
         edge_severity = mv_bio_range$avg,
         site = NA,
         exp_plot = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(severity = predict(mv_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA),
         severity_se = predict(mv_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)

# plot
mv_bio_plot <- ggplot(mv_bio_sim_dat, aes(total_biomass.g, severity)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "Late August", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Microstegium), " lesions", sep = ""))) +
  xlab(expression(paste("Total biomass (g ", m^-2, ")", sep = ""))) +
  temp_theme +
  theme(legend.position = "none")


#### Mv humidity ####

# min and max values
mv_hum_range <- sep_mv_water_dat %>%
  summarise(min = min(hum_avg, na.rm = T),
            max = max(hum_avg, na.rm = T))

# prediction data
mv_hum_sim_dat <- tibble(hum_avg = seq(mv_hum_range$min, mv_hum_range$max, length.out = 200)) %>%
  mutate(site = NA,
         plot = NA,
         fungicide = 0,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"),
         hum_perc = hum_avg * 100) %>%
  mutate(severity = predict(mv_severity_humidity_mod, newdata = ., type = "response", re.form = NA),
         severity_se = predict(mv_severity_humidity_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)

# plot
mv_hum_plot <- ggplot(mv_hum_sim_dat, aes(hum_perc, severity)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "September", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Microstegium), " lesions", sep = ""))) +
  xlab("Relative humidity (%)") +
  temp_theme +
  theme(legend.position = "none")


#### humidity total biomass ####

# min and max values
hum_bio_range <- hum_bio_dat %>%
  filter(month_name == "September") %>%
  summarise(min = min(total_biomass.g, na.rm = T),
            max = max(total_biomass.g, na.rm = T))

# prediction data
hum_bio_sim_dat <- tibble(total_biomass.g = seq(hum_bio_range$min, hum_bio_range$max, length.out = 200)) %>%
  mutate(month_name = "September",
         site = NA,
         fungicide = 0,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(humidity = predict(havg_totb_mod, newdata = ., type = "response", re.form = NA) * 100,
         humidity_se = predict(havg_totb_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit * 100)

# plot
hum_bio_plot <- ggplot(hum_bio_sim_dat, aes(total_biomass.g, humidity)) +
  geom_ribbon(aes(ymin = humidity - humidity_se, ymax = humidity + humidity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "September", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab("Relative humidity (%)") +
  xlab(expression(paste("Total biomass (g ", m^-2, ")", sep = ""))) +
  temp_theme +
  theme(legend.position = "none")


#### Ev edge severity ####

# min and max values
ev_edge_range <- ev_bio_edge_dat %>%
  filter(month == "jul") %>%
  summarise(min = min(edge_severity, na.rm = T),
            max = max(edge_severity, na.rm = T),
            avg_mv = mean(Mv_biomass.g, na.rm = T),
            avg_ev = mean(Ev_biomass.g, na.rm = T))

# prediction data
ev_edge_sim_dat <- tibble(edge_severity = rep(seq(ev_edge_range$min, ev_edge_range$max, length.out = 200), 2),
                          fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Month = "July",
         Mv_biomass.g = ev_edge_range$avg_mv,
         Ev_biomass.g = ev_edge_range$avg_ev,
         site = NA,
         exp_plot = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(severity = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA),
         severity_se = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)

# plot
ev_edge_plot <- ggplot(ev_edge_sim_dat, aes(edge_severity, severity)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "July", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Elymus), " lesions", sep = ""))) +
  xlab("Background disease") +
  temp_theme +
  theme(legend.position = "none")


#### Ev Ev biomass ####

# min and max values
ev_ev_range <- ev_bio_edge_dat %>%
  filter(month == "early_aug") %>%
  summarise(min = min(Ev_biomass.g, na.rm = T),
            max = max(Ev_biomass.g, na.rm = T),
            avg_mv = mean(Mv_biomass.g, na.rm = T),
            avg_ed = mean(edge_severity, na.rm = T))

# prediction data
ev_ev_sim_dat <- tibble(Ev_biomass.g = rep(seq(ev_ev_range$min, ev_ev_range$max, length.out = 200), 2),
                          fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Month = "Early August",
         Mv_biomass.g = ev_ev_range$avg_mv,
         edge_severity = ev_ev_range$avg_ed,
         site = NA,
         exp_plot = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(severity = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA),
         severity_se = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)

# plot
ev_ev_plot <- ggplot(ev_ev_sim_dat, aes(Ev_biomass.g, severity)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "Early August", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Elymus), " lesions", sep = ""))) +
  xlab(expression(paste(italic(Elymus), " biomass (g ", m^-2, ")", sep = ""))) +
  temp_theme +
  theme(legend.position = "none")


#### Ev Mv biomass ####

# min and max values
ev_mv_range <- ev_bio_edge_dat %>%
  filter(month == "late_aug") %>%
  summarise(min = min(Mv_biomass.g, na.rm = T),
            max = max(Mv_biomass.g, na.rm = T),
            avg_ev = mean(Ev_biomass.g, na.rm = T),
            avg_ed = mean(edge_severity, na.rm = T))

# prediction data
ev_mv_sim_dat <- tibble(Mv_biomass.g = rep(seq(ev_mv_range$min, ev_mv_range$max, length.out = 200), 2),
                        fungicide = rep(c(0, 1), each = 200)) %>%
  mutate(Month = "Late August",
         Ev_biomass.g = ev_mv_range$avg_ev,
         edge_severity = ev_mv_range$avg_ed,
         site = NA,
         exp_plot = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(severity = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA),
         severity_se = predict(ev_severity_biomass_edge_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)

# plot
ev_mv_plot <- ggplot(ev_mv_sim_dat, aes(Mv_biomass.g, severity)) +
  geom_ribbon(aes(ymin = severity - severity_se, ymax = severity + severity_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "Late August", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Elymus), " lesions", sep = ""))) +
  xlab(expression(paste(italic(Microstegium), " biomass (g ", m^-2, ")", sep = ""))) +
  temp_theme +
  theme(legend.position = "none")


#### combine plots ####

# legend
leg <- get_legend(mv_edge_plot)

# combine
pdf("output/focal_severity_figure_density_exp2.pdf",
    width = 7, height = 4)
plot_grid(mv_edge_plot + theme(legend.position = "none"), mv_bio_plot, mv_hum_plot, hum_bio_plot,
          ev_edge_plot, ev_ev_plot, ev_mv_plot, leg,
          nrow = 2,
          labels = c(LETTERS[1:7], ""), label_size = 10)
dev.off()