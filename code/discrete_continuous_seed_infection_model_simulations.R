#### set-up ####

# note: a more complex version of the model with seed infection is in discrete_AFP_model.R and discrete_continuous_model_simulations.Rmd, but I chose to just reduce germination by the maximum observed % due to seed infection

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)
library(patchwork)
library(ggbreak)

# import parameters
parms_all <- read_csv("output/model_parameters_2018_2019_density_exp.csv")
# model_parameters_2018_2019_density_exp

# continuous models
source("code/continuous_AFP_model.R")
source("code/continuous_AF_model.R")
source("code/continuous_AP_model.R")
source("code/continuous_FP_model.R")
source("code/continuous_A_model.R")
source("code/continuous_F_model.R")
source("code/continuous_P_model.R")

# discrete model
source("code/discrete_AFP_no_seed_infection_model.R")

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(-0.1, 0, -0.1, 0, unit = "cm"),
        plot.title = element_text(size = 10, hjust = -0.1, face = "bold"))

col_pal = c("black", "#238A8DFF")


#### parameters ####

# growing season days
gs_days <- 160
gs_days_seq <- seq(0, gs_days, by = 1)

# simulation years
years <- 100

# fungicide parameters
parms_fung <- parms_all %>%
  filter(Treatment == "fungicide" | is.na(Treatment))

# extract germination
germ_fung <- parms_fung %>%
  filter(Parameter == "g_A") %>%
  pull(Estimate)

# control parameters
# reduce germination by largest amt due to seed infection
parms_ctrl <- parms_all %>%
  filter(Treatment == "control" | is.na(Treatment)) %>%
  mutate(Estimate = if_else(Parameter == "g_A", germ_fung * (1-0.44), Estimate))

# conidia per 1 g litter
conidia_init <- as.numeric(parms_all[parms_all$Parameter == "h", "Estimate"])

# initial plant sizes
bio_A_init <- as.numeric(parms_all[parms_all$Parameter == "b_A", "Estimate"])
bio_F_init <- as.numeric(parms_all[parms_all$Parameter == "b_F", "Estimate"])
bio_P_init <- as.numeric(parms_all[parms_all$Parameter == "b_P", "Estimate"])


#### annual alone ####

# continuous input values
startA_cont_ctrl <- c(LogB_A = log(bio_A_init * 10), I_A = 0, D = 0, C = conidia_init)
startA_cont_fung <- c(LogB_A = log(bio_A_init * 10), I_A = 0, D = 0, C = 0)

# apply continuous function
modA_cont_ctrl <- ode(y = startA_cont_ctrl, times = gs_days_seq, 
                      func = cont_A_mod, parms = parms_ctrl) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_A = exp(LogB_A),
         S_A = B_A - I_A) %>%
  select(-LogB_A) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

modA_cont_fung <- ode(y = startA_cont_fung, times = gs_days_seq, 
                      func = cont_A_mod, parms = parms_fung) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_A = exp(LogB_A),
         S_A = B_A - I_A) %>%
  select(-LogB_A) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

# figure of continuous results
modA_cont_ctrl %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

modA_cont_fung %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

# discrete input values
modA_disc_ctrl <- disc_AFP_mod(A0 = 10, F0 = 0, P0 = 0, L0 = 0, C0 = conidia_init, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = 0, BF0 = 0, simtime = years, gs_time = gs_days, disc_parms = parms_ctrl, con_parms = parms_ctrl)

modA_disc_fung <- disc_AFP_mod(A0 = 10, F0 = 0, P0 = 0, L0 = 0, C0 = 0, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = 0, BF0 = 0, simtime = years, gs_time = gs_days, disc_parms = parms_fung, con_parms = parms_fung)

# figure of discrete results
modA_disc_ctrl %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

modA_disc_fung %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")


#### perennial seedling alone ####

# continuous input values
startF_cont_ctrl <- c(LogB_F = log(bio_F_init * 10), I_F = 0, D = 0, C = conidia_init)
startF_cont_fung <- c(LogB_F = log(bio_F_init * 10), I_F = 0, D = 0, C = 0)

# apply continuous function
modF_cont_ctrl <- ode(y = startF_cont_ctrl, times = gs_days_seq, 
                      func = cont_F_mod, parms = parms_ctrl) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_F = exp(LogB_F),
         S_F = B_F - I_F) %>%
  select(-LogB_F) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

modF_cont_fung <- ode(y = startF_cont_fung, times = gs_days_seq, 
                      func = cont_F_mod, parms = parms_fung) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_F = exp(LogB_F),
         S_F = B_F - I_F) %>%
  select(-LogB_F) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

# figure of continuous results
modF_cont_ctrl %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

modF_cont_fung %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

# discrete input values
modF_disc_ctrl <- disc_AFP_mod(A0 = 0, F0 = 10, P0 = 0, L0 = 0, C0 = conidia_init, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = 0, BF0 = 0, simtime = years, gs_time = gs_days, disc_parms = parms_ctrl, con_parms = parms_ctrl)

# multiply by 2 for invasion
# want Perennial adult to be at equilibrium
modF_disc_fung <- disc_AFP_mod(A0 = 0, F0 = 10, P0 = 0, L0 = 0, C0 = 0, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = 0, BF0 = 0, simtime = years * 2, gs_time = gs_days, disc_parms = parms_fung, con_parms = parms_fung)

# figure of discrete results
modF_disc_ctrl %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

modF_disc_fung %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")


#### perennial adult alone ####

# continuous input values
startP_cont_ctrl <- c(LogB_P = log(bio_P_init * 10), I_P = 0, D = 0, C = conidia_init)
startP_cont_fung <- c(LogB_P = log(bio_P_init * 10), I_P = 0, D = 0, C = 0)

# apply continuous function
modP_cont_ctrl <- ode(y = startP_cont_ctrl, times = gs_days_seq, 
                      func = cont_P_mod, parms = parms_ctrl) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_P = exp(LogB_P),
         S_P = B_P - I_P) %>%
  select(-LogB_P) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

modP_cont_fung <- ode(y = startP_cont_fung, times = gs_days_seq, 
                      func = cont_P_mod, parms = parms_fung) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_P = exp(LogB_P),
         S_P = B_P - I_P) %>%
  select(-LogB_P) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

# figure of continuous results
modP_cont_ctrl %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

modP_cont_fung %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

# discrete input values
modP_disc_ctrl <- disc_AFP_mod(A0 = 0, F0 = 0, P0 = 10, L0 = 0, C0 = conidia_init, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = bio_P_init * 10, BF0 = 0, simtime = years, gs_time = gs_days, disc_parms = parms_ctrl, con_parms = parms_ctrl)

modP_disc_fung <- disc_AFP_mod(A0 = 0, F0 = 0, P0 = 10, L0 = 0, C0 = 0, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = bio_P_init * 10, BF0 = 0, simtime = years, gs_time = gs_days, disc_parms = parms_fung, con_parms = parms_fung)

# figure of discrete results
modP_disc_ctrl %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

modP_disc_fung %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")


#### invasion without disease ####

# established perennial
modF_disc_fung_fin <- modF_disc_fung %>%
  filter(time == max(time))

# add annual
# double time so that disease can be added at 100 years
modInv_disc_fung <- disc_AFP_mod(A0 = 10, 
                                 F0 = filter(modF_disc_fung_fin, species == "Perennial first-year")$N, 
                                 P0 = filter(modF_disc_fung_fin, species == "Perennial adult")$N, 
                                 L0 = filter(modF_disc_fung_fin, species == "Litter")$N, 
                                 C0 = filter(modF_disc_fung_fin, species == "Conidia")$N, 
                                 IA0 = 0, 
                                 IF0 = filter(modF_disc_fung_fin, species == "Perennial first-year infected biomass")$N, 
                                 IP0 = filter(modF_disc_fung_fin, species == "Perennial adult infected biomass")$N, 
                                 BP0 = filter(modF_disc_fung_fin, species == "Perennial adult biomass")$N, 
                                 BF0 = filter(modF_disc_fung_fin, species == "Perennial first-year biomass")$N, 
                                 simtime = years*2, gs_time = gs_days, 
                                 disc_parms = parms_fung, con_parms = parms_fung)

# figure of discrete results
modInv_disc_fung %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")


#### add disease ####

# invaded community
modInv_disc_fung_fin <- modInv_disc_fung %>%
  filter(time == max(time)/2)

# add disease
modDis_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
                                 F0 = filter(modInv_disc_fung_fin, species == "Perennial first-year")$N, 
                                 P0 = filter(modInv_disc_fung_fin, species == "Perennial adult")$N, 
                                 L0 = filter(modInv_disc_fung_fin, species == "Litter")$N, 
                                 C0 = conidia_init, 
                                 IA0 = filter(modInv_disc_fung_fin, species == "Annual infected biomass")$N, 
                                 IF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year infected biomass")$N, 
                                 IP0 = filter(modInv_disc_fung_fin, species == "Perennial adult infected biomass")$N, 
                                 BP0 = filter(modInv_disc_fung_fin, species == "Perennial adult biomass")$N, 
                                 BF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year biomass")$N, 
                                 simtime = years, gs_time = gs_days, 
                                 disc_parms = parms_ctrl, con_parms = parms_ctrl)

# figure of discrete results
modDis_disc_ctrl %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")


#### combine simulations ####

# combine model outputs
mod_inv_dis <- modF_disc_fung %>%
  full_join(modInv_disc_fung %>%
              mutate(time = time + max(modF_disc_fung$time))) %>%
  mutate(treatment = "fungicide") %>%
  full_join(modDis_disc_ctrl %>%
              mutate(time = time + max(modF_disc_fung$time) + max(modInv_disc_fung$time)/2,
                     treatment = "control"))

# figure of discrete results
mod_inv_dis %>%
  ggplot(aes(x = time, y = N, color = treatment)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

# make discontinuous transitions continuous
bio_inv_dis <- mod_inv_dis %>%
  filter(time == 300 & 
           (str_detect(species, "biomass") == T | species == "Conidia")) %>%
  mutate(treatment = "control")

# add transitions
mod_inv_dis2 <- mod_inv_dis %>%
  full_join(bio_inv_dis)

# figure of discrete results
mod_inv_dis2 %>%
  ggplot(aes(x = time, y = N, color = treatment)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

# adjust time to invasion = 0
# remove early/late dynamics
# combine adult & first-year perennial
mod_inv_dis3 <- mod_inv_dis2 %>%
  filter(species %in% c("Annual biomass",
                        "Perennial adult biomass",
                        "Perennial first-year biomass") &
           time >= 190 & time < 315) %>%
  mutate(time = time - 200,
         species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass) %>%
  select(time, treatment, Annual_biomass, Perennial_biomass) %>%
  pivot_longer(cols = c(Annual_biomass, Perennial_biomass),
               names_to = "species",
               values_to = "N") %>%
  mutate(Species = if_else(species == "Annual_biomass", "invader", "competitor") %>%
           fct_relevel("invader"))

# save
write_csv(mod_inv_dis3, "output/discrete_seed_infection_time_series.csv")
mod_inv_dis3 <- read_csv("output/discrete_seed_infection_time_series.csv") %>%
  mutate(Species = fct_relevel(Species, "invader"))

# phase labels
phase_lab <- tibble(time = c(-5, 5, 107),
                    N = 400,
                    lab = c("pre-invasion", "invasion", "disease emergence"))

# figure of discrete results
ggsave("output/discrete_seed_infection_simulation_figure.pdf", 
       device = "pdf", width = 15, height = 6, units = "cm")
mod_inv_dis3 %>%
  ggplot(aes(x = time, y = N)) +
  geom_rect(xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "gray95", color = NA, show.legend = F) +
  geom_rect(xmin = 100, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray75", color = NA, show.legend = F) +
  geom_line(aes(color = treatment, linetype = Species)) +
  geom_text(data = phase_lab, check_overlap = T, size = 2.5,
            aes(label = lab)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_x_break(c(10, 100)) +
  labs(x = "Time (years)", y = expression(paste('Biomass (g/', m^2, ')', sep = ""))) +
  fig_theme
dev.off()


#### values for text ####

# extract perennial biomass values
pre_inv_bio <- mod_inv_dis3 %>%
  filter(time == -5 & species == "Perennial_biomass") %>%
  pull(N)

inv_bio <- mod_inv_dis3 %>%
  filter(time == 10 & species == "Perennial_biomass") %>%
  pull(N)

dis_inv_bio <- mod_inv_dis3 %>%
  filter(time == 110 & species == "Perennial_biomass" &
           treatment == "control") %>%
  pull(N)

# invasion effect
100 * (inv_bio - pre_inv_bio)/pre_inv_bio

# disease effect
100 * (dis_inv_bio - inv_bio)/inv_bio
