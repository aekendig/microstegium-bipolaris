##### info ####

# file: annual_perennial_relative_abundance_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: relative abundance of annual to perennial


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(ggforce)

# number of iterations
n_samps <- 300

# call model script
source("code/annual_perennial_simulation_model_2019_density_exp.R")

# samples from parameter distributions
samps <- n_samps

# simulation time
simtimeR <- 600
simtimeI <- 1000


#### resident E. virginicus ####

# initial conditions
A0E <- 0 # initial annual population size
S0E <- 0 # initial perennial seedling population size
P0E <- 1 # initial perennial adult population size
L0E <- 0 # initial annual litter amount

# initiate lists
abundE <- list()
paramsE <- list()
abundE2 <- list()
paramsE2 <- list()

# simulation without disease
for(iter in 1:samps){
  
  mod_estE <- sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 0, iter)
  S0Ei <- mod_estE[[1]] %>% filter(time == simtimeR & species == "Perennial seedling") %>% pull(N)
  P0Ei <- mod_estE[[1]] %>% filter(time == simtimeR & species == "Perennial adult") %>% pull(N)
  L0Ei <- mod_estE[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  mod_invE <- sim_fun(A0 = 10, S0Ei, P0Ei, L0Ei, simtimeI, disease = 0, iter)
  abundE[[iter]] <- mod_invE[[1]]
  paramsE[[iter]] <- mod_invE[[2]]
  
  mod_estE2 <- sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 1, iter)
  S0Ei2 <- mod_estE2[[1]] %>% filter(time == simtimeR & species == "Perennial seedling") %>% pull(N)
  P0Ei2 <- mod_estE2[[1]] %>% filter(time == simtimeR & species == "Perennial adult") %>% pull(N)
  L0Ei2 <- mod_estE2[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  mod_invE2 <- sim_fun(A0 = 10, S0Ei2, P0Ei2, L0Ei2, simtimeI, disease = 1, iter)
  abundE2[[iter]] <- mod_invE2[[1]]
  paramsE2[[iter]] <- mod_invE2[[2]]
  
}

# convert to dataframes
abund_datE <- do.call(rbind, abundE)
param_datE <- do.call(rbind, paramsE)
abund_datE2 <- do.call(rbind, abundE2)
param_datE2 <- do.call(rbind, paramsE2)


#### process data ####

rel_abund_dat <- abund_datE %>%
  mutate(disease = "without disease") %>%
  full_join(abund_datE2 %>%
              mutate(disease = "with disease")) %>%
  pivot_wider(names_from = species, values_from = N) %>%
  rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
  mutate(Perennial = Perennial_seedling + Perennial_adult,
         tot_abund = Annual + Perennial,
         ann_rel_abund = Annual / tot_abund,
         per_rel_abund = Perennial / tot_abund,
         iterationF = as.factor(iteration))

# relative abundance summary
rel_sum_dat <- rel_abund_dat %>%
  filter(time >= simtimeI - 100) %>%
  group_by(disease, iteration, iterationF) %>%
  summarise(ann_rel_abund = mean(ann_rel_abund),
            per_rel_abund = mean(per_rel_abund)) %>%
  ungroup() %>%
  mutate(disease = fct_relevel(disease, "without disease"))

# absolute abundance
abs_abund_dat <- rel_abund_dat %>%
  select(disease, time, iteration, iterationF, Perennial, Annual, Litter) %>%
  pivot_longer(cols = c("Perennial", "Annual", "Litter"), names_to = "species", values_to = "density")


#### figures ####

pdf("output/annual_perennial_relative_abundance_time_series_2019_density_exp.pdf")
ggplot(rel_abund_dat, aes(time, ann_rel_abund, color = iterationF)) +
  geom_line()  +
  facet_wrap(~ disease) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/annual_perennial_relative_abundance_distribution_2019_density_exp.pdf", width = 4.5, height = 3)
ggplot(rel_sum_dat, aes(x = ann_rel_abund, fill = disease)) +
  geom_histogram(binwidth = 0.05, position = position_dodge()) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  xlab(expression(paste("Relative abundance (", italic("M. vimineum"), "/", italic("E. virginicus"), ")", sep = ""))) +
  ylab("Number of draws") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))
dev.off()

# pdf("output/annual_perennial_absolute_abundance_time_series_2019_density_exp.pdf")
# for(i in 1:33){
#   print(ggplot(abs_abund_dat, aes(time, density, color = species, linetype = disease)) +
#           geom_line()  +
#           facet_wrap_paginate(~ iterationF, nrow = 3, ncol = 3, page = i, scales = "free") +
#           theme_bw() +
#           theme(legend.position = "bottom",
#                 legend.direction = "horizontal",
#                 legend.box.margin = margin(-10, 0, 0, 0),
#                 legend.margin = margin(0, 0, 0, 0))) 
# }
# dev.off()

# see annual_perennial_troubleshooting_2019_density_exp for additional code


#### output ####
write_csv(rel_abund_dat, "output/annual_perennial_relative_abundance_time_series_2019_density_exp.csv")
