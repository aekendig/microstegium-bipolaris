##### info ####

# file: severity-by-treatment
# author: Amy Kendig
# date last edited: 8/2/19
# goal: see how Mv and Ev disease severity is affected by treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)

# run leaf scan files
source("./code/ev-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), c("eleaf")))
source("./code/mv-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), c("eleaf", "mleaf")))

# run covariate file
source("./code/covariate-data-processing-2018.R")
rm(list = setdiff(ls(), c("eleaf", "mleaf", "covar")))

# import data
trt <- read_csv("./data/plot-treatments-for-figures-2018-density-exp.csv")
trt_s <- read_csv("./data/plot-treatments-2018-density-exp.csv")
ftilj <- read_csv("./data/focal-size-disease-jul-2018-density-exp.csv")
atila <- read_csv("./data/all-disease-seeds-late-aug-2018-density-exp.csv")
etils <- read_csv("./data/ev-disease-seeds-sep-2018-density-exp.csv")
mtils <- read_csv("./data/mv-disease-sep-2018-density-exp.csv")

# plot parameters
colpal = c("#3CBB75FF", "#39568CFF")
sm_txt = 12
lg_txt = 14
an_txt = 3


#### edit data ####

# July tiller data
unique(ftilj$field_notes)
filter(ftilj, field_notes == "forgot to collect infection information") # D3 8W EvA
filter(ftilj, field_notes == "tree fell on plot")

etilj <- ftilj %>%
  filter(sp == "Ev" & (is.na(field_notes) | field_notes != "tree fell on plot") & !is.na(height.cm)) %>%
  select(site, plot, treatment, sp, ID, infec, leaves_tot, leaves_infec, Bp_spots_Ev, field_notes) %>%
  mutate(month = "July")

mtilj <- ftilj %>%
  filter(sp == "Mv" & (is.na(field_notes) | field_notes != "tree fell on plot") & !is.na(height.cm)) %>%
  select(site, plot, treatment, sp, ID, infec, leaves_tot, leaves_infec, field_notes) %>%
  mutate(month = "July")

# August tiller data
unique(atila$field_notes)

etila <- atila %>%
  filter(sp == "Ev" & ID %in% c("1", "2", "3", "A") & no_green == 0) %>%
  select(site, plot, treatment, sp, ID, infec, leaves_tot, leaves_infec, Bp_spots_Ev, field_notes) %>%
  mutate(month = "August")

mtila <- atila %>%
  filter(sp == "Mv" & no_green == 0) %>%
  select(site, plot, treatment, sp, ID, infec, leaves_tot, leaves_infec, field_notes) %>%
  mutate(month = "August")

# September tiller data
unique(etils$field_notes)
unique(mtils$field_notes)

etils2 <- etils %>%
  filter(sp == "Ev" & ID %in% c("1", "2", "3", "A") & !is.na(leaves_tot)) %>%
  select(site, plot, treatment, sp, ID, leaves_tot, leaves_infec, Bp_spots_Ev, field_notes) %>%
  mutate(month = "September")

mtils2 <- mtils %>%
  filter(sp == "Mv" & ID %in% c("1", "2", "3") & !is.na(leaves_tot)) %>%
  select(site, plot, treatment, sp, ID, leaves_tot, leaves_infec, field_notes) %>%
  mutate(month = "September")

# tillers by month
mtil <- full_join(mtilj, mtila) %>%
  full_join(mtils2)

etil <- full_join(etilj, etila) %>%
  full_join(etils2)

# Mv focal leaves
mfsev <- mleaf %>%
  filter(focal == 1 & remove == 0) %>%
  select(month, site, plot, treatment, ID, lesion_area.pix, green_area.pix, leaf_area.pix) %>%
  mutate(plot = as.numeric(plot)) %>%
  full_join(mtil) %>%
  mutate(leaves_infec2 = ifelse(leaves_infec == 0 & !is.na(leaf_area.pix), 1, leaves_infec),
         ind_severity = (lesion_area.pix - green_area.pix) * leaves_infec2 / (leaf_area.pix * leaves_tot))
# leaves with no scan or no leaf counts will have NA for the severity metric
# severity only includes infected leaves (with 2 exceptions - lesions weren't detected)

mbsev <- mleaf %>%
  filter(focal == 0) %>%
  select(month, site, plot, treatment, ID, lesion_area.pix, green_area.pix, leaf_area.pix) %>%
  mutate(plot = as.numeric(plot)) %>%
  mutate(bg_severity = (lesion_area.pix - green_area.pix) / leaf_area.pix)
# also only includes infected leaves

# Ev focal leaves
# remove D2.7.Water.Ev.A (not Ev) and the one leaf duplicated to capture Bp
efsev <- eleaf %>%
  filter(focal == 1 & remove == 0 & plant != "D2_7W_EvA") %>%
  select(month, site, plot, treatment, ID, lesion_area.pix, leaf_area.pix) %>%
  mutate(plot = as.numeric(plot)) %>%
  full_join(etil) %>%
  mutate(leaves_infec2 = ifelse(leaves_infec == 0 & !is.na(leaf_area.pix), 1, leaves_infec),
         ind_severity = lesion_area.pix * leaves_infec2 / (leaf_area.pix * leaves_tot))


#### visualize ####

# all months, Mv focal
mfsev %>%
  ggplot(aes(x = treatment, y = ind_severity, color = treatment)) +
  geom_violin() +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

mfsev %>%
  ggplot(aes(x = treatment, y = ind_severity)) +
  geom_point(alpha = 0.3, aes(color = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", size = 3, shape = 21, aes(fill = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

# all months, Mv background
mbsev %>%
  ggplot(aes(x = treatment, y = bg_severity, color = treatment)) +
  geom_violin() +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

# all months, Ev focal
efsev %>%
  ggplot(aes(x = treatment, y = ind_severity, color = treatment)) +
  geom_violin() +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

# highest monthly mean
efsev %>%
  group_by(month, treatment) %>%
  summarise(mean_sev = mean(ind_severity, na.rm = T)) # July

# Ev with Bp
etil %>%
  filter(leaves_infec > 0) %>%
  ggplot(aes(x = treatment, y = Bp_spots_Ev, color = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

etil %>%
  ggplot(aes(x = treatment, y = Bp_spots_Ev, color = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

etil %>%
  group_by(month, treatment) %>%
  summarise(bp = sum(Bp_spots_Ev, na.rm = T),
            samps = n())

# Ev infected
etil %>%
  filter(!is.na(leaves_infec)) %>%
  mutate(disease = ifelse(leaves_infec > 0, 1, 0)) %>%
  ggplot(aes(x = treatment, y = disease, color = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)

mtil %>%
  filter(!is.na(leaves_infec)) %>%
  mutate(disease = ifelse(leaves_infec > 0, 1, 0)) %>%
  ggplot(aes(x = treatment, y = disease, color = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~month)


#### final figures ####

# Mv severity in September
pdf("./output/mv-severity-by-fungicide-sep-2018.pdf", height = 3, width = 3)
mfsev %>%
  filter(month == "September") %>%
  ggplot(aes(x = treatment, y = ind_severity)) +
  geom_violin(aes(color = treatment, fill = treatment), alpha = 0.5) +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.05) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.text.x = element_text(color = "black", size = lg_txt),
        plot.title = element_text(hjust = 0.5, size = lg_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylab("Infection severity") +
  ggtitle("Invader (Microstegium)")
dev.off()

# Ev severity in July
pdf("./output/ev-severity-by-fungicide-jul-2018.pdf", height = 3, width = 3)
efsev %>%
  filter(month == "July") %>%
  ggplot(aes(x = treatment, y = ind_severity)) +
  geom_violin(aes(color = treatment, fill = treatment), alpha = 0.5) +
  stat_summary(fun.y = "mean", geom = "point", size = 2) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.05) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.text.x = element_text(color = "black", size = lg_txt),
        plot.title = element_text(hjust = 0.5, size = lg_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) +
  ylab("Infection severity") +
  ggtitle("Native (Elymus)")
dev.off()

# Proportion of infected Ev with Bipolaris
pdf("./output/ev-bp-by-fungicide-sep-2018.pdf", height = 3, width = 3)
etil %>%
  filter(leaves_infec > 0 & month == "August") %>%
  ggplot(aes(x = treatment, y = Bp_spots_Ev, fill = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  stat_summary(fun.y = "mean", geom = "point", size = 4, shape = 21) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.text.x = element_text(color = "black", size = lg_txt),
        plot.title = element_text(hjust = 0.5, size = lg_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colpal) +
  ylab("Proportion infected with Bipolaris") +
  ggtitle("Native (Elymus)")
dev.off()