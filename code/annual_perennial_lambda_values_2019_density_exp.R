##### info ####

# file: annual_perennial_lambda_values_2019_density_exp
# author: Amy Kendig
# date last edited: 10/28/20
# goal: lifetime reproduction


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)

# import data
cdat <- read_csv("output/annual_perennial_kortessis_relative_abundance_constants_2019_density_exp.csv")


#### edit data ####

cdat2 <- cdat %>%
  filter(resident == "E. virginicus") %>%
  mutate(disease = fct_relevel(disease, "without disease"),
         rel_fit = l.A / l.P)


# figure
pdf("output/annual_perennial_lambda_values_2019_density_exp.pdf", width = 3, height = 3)
ggplot(cdat2, aes(disease, log(rel_fit), color = disease, fill = disease)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  stat_halfeye(point_interval = mean_hdci, .width = c(0.66, 0.95), point_size = 3, shape = 21, point_color = "white", interval_color = "black", slab_alpha = 0.7, aes(fill = disease)) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  ylab(expression(paste("Relative fitness (", italic("M. vimineum"), "/", italic("E. virginicus"), ")", sep = ""))) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 10),
        strip.background = element_blank())
dev.off()

# test different
cdat2 %>%
  select(iteration, disease, rel_fit) %>%
  mutate(disease = recode(disease, "without disease" = "no", "with disease" = "yes")) %>%
  pivot_wider(names_from = disease, values_from = rel_fit) %>%
  mutate(dis_diff = no - yes) %>%
  mean_hdi(dis_diff)
