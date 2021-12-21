#### set-up ####

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

# control parameters
parms_ctrl <- parms_all %>%
  filter(Treatment == "control" | is.na(Treatment))

# fungicide parameters
parms_fung <- parms_all %>%
  filter(Treatment == "fungicide" | is.na(Treatment))

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
write_csv(mod_inv_dis3, "output/discrete_no_seed_infection_time_series.csv")
mod_inv_dis3 <- read_csv("output/discrete_no_seed_infection_time_series.csv") %>%
  mutate(Species = fct_relevel(Species, "invader"))

# phase labels
phase_lab <- tibble(time = c(-5, 5, 107),
                    N = 400,
                    lab = c("pre-invasion", "invasion", "disease emergence"))

# figure of discrete results
time_series_fig <- mod_inv_dis3 %>%
  ggplot(aes(x = time, y = N)) +
  geom_rect(xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "gray95", color = NA, show.legend = F) +
  geom_rect(xmin = 100, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray75", color = NA, show.legend = F) +
  geom_line(aes(color = treatment, linetype = Species)) +
  geom_text(data = phase_lab, check_overlap = T, size = 2.5,
            aes(label = lab)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_x_break(c(10, 100)) +
  labs(x = "Time (years)", y = expression(paste('Biomass (g/', m^2, ')', sep = "")),
       title = "A") +
  fig_theme + 
  theme(plot.title = element_text(size = 10, hjust = 0, face = "bold"))


#### invader intraspecific comp ####

# extract alpha values
alpha_AA_fung <- parms_fung %>%
  filter(Parameter == "alpha_AA") %>%
  pull(Estimate)

alpha_AA_ctrl <- parms_ctrl %>%
  filter(Parameter == "alpha_AA") %>%
  pull(Estimate)

# alpha values
alpha_AA_vals <- c(seq(alpha_AA_fung, 100 * alpha_AA_ctrl, length.out = 19),
                   alpha_AA_ctrl) %>%
  unique() %>% sort()

# empty list
modAlphaAA <- vector(mode = "list", length = length(alpha_AA_vals))

# cycle through alpha values
for(i in 1:length(alpha_AA_vals)){
  
  # update parameters
  parms_ctrl_temp <- parms_ctrl %>%
    mutate(Estimate = if_else(Parameter == "alpha_AA", alpha_AA_vals[i], Estimate))
  
  # add disease
  modAlphaAA_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
                                       F0 = filter(modInv_disc_fung_fin, species == "Perennial first-year")$N, 
                                       P0 = filter(modInv_disc_fung_fin, species == "Perennial adult")$N, 
                                       L0 = filter(modInv_disc_fung_fin, species == "Litter")$N, 
                                       C0 = conidia_init, 
                                       IA0 = filter(modInv_disc_fung_fin, species == "Annual infected biomass")$N, 
                                       IF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year infected biomass")$N, 
                                       IP0 = filter(modInv_disc_fung_fin, species == "Perennial adult infected biomass")$N, 
                                       BP0 = filter(modInv_disc_fung_fin, species == "Perennial adult biomass")$N, 
                                       BF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year biomass")$N, 
                                       simtime = years/2, gs_time = gs_days, 
                                       disc_parms = parms_ctrl_temp, con_parms = parms_ctrl_temp)
  
  # save data table
  modAlphaAA[[i]] <- modAlphaAA_disc_ctrl %>%
    mutate(alpha_AA = alpha_AA_vals[i])
}

# collapse list
modAlphaAA2 <- modAlphaAA %>%
  bind_rows()

# save file
write_csv(modAlphaAA2, "output/discrete_no_seed_infection_alphaAA_sim.csv")
modAlphaAA2 <- read_csv("output/discrete_no_seed_infection_alphaAA_sim.csv")

# check that each simulation looks reasonable and that biomass stabilizes
for(i in 1:length(alpha_AA_vals)){
  
  print(modAlphaAA2 %>%
          filter(alpha_AA == alpha_AA_vals[i]) %>%
          ggplot(aes(x = time, y = N)) +
          geom_line() +
          theme_classic() +
          facet_wrap(~ species, scales = "free_y"))
  
}

# extract equilibrium perennial biomass
per_bio_alone <- modF_disc_fung_fin %>%
  filter(species %in% c("Perennial adult biomass",
                        "Perennial first-year biomass")) %>%
  summarize(bio = sum(N)) %>%
  pull(bio)

# per_bio_alone <- 719.323

# select last time point
# combine adult & first-year perennial
modAlphaAA3 <- modAlphaAA2 %>%
  filter(species %in% c("Annual biomass",
                        "Perennial adult biomass",
                        "Perennial first-year biomass") &
           time == max(time)) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass,
         Impact = per_bio_alone/Perennial_biomass,
         Coefficient = case_when(alpha_AA == alpha_AA_ctrl ~ "control",
                                 alpha_AA == alpha_AA_fung ~ "fungicide",
                                 TRUE ~ NA_character_)) %>%
  filter(alpha_AA < 0.09) # trim to 10 values

# figure
alphaAA_fig <- ggplot(modAlphaAA3, aes(x = alpha_AA, y = Perennial_biomass)) +
  geom_hline(yintercept = per_bio_alone, linetype = "dashed", color = "gray60") +
  geom_text(x = max(modAlphaAA3$alpha_AA), y = per_bio_alone - 10, label = "Pre-invasion biomass",
            hjust = 1, vjust = 1, size = 2, check_overlap = T) +
  geom_line(color = "gray60") +
  geom_point(color = "gray60") +
  geom_point(data = filter(modAlphaAA3, !is.na(Coefficient)),
             aes(color = Coefficient), size = 3) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  labs(x = "Intraspecific invader competition", 
       y = expression(paste('Competitor biomass (g/', m^2, ')', sep = "")),
       title = "B") +
  fig_theme +
  theme(legend.position = c(0.3, 0.5))


#### invader interspecific comp ####

# extract alpha values
alpha_PA_fung <- parms_fung %>%
  filter(Parameter == "alpha_PA") %>%
  pull(Estimate)

alpha_PA_ctrl <- parms_ctrl %>%
  filter(Parameter == "alpha_PA") %>%
  pull(Estimate)

# alpha values
alpha_PA_vals <- c(seq(alpha_PA_fung, alpha_PA_ctrl, length.out = 10)) %>%
  unique() %>% sort()

# empty list
modAlphaPA <- vector(mode = "list", length = length(alpha_PA_vals))

# cycle through alpha values
for(i in 1:length(alpha_PA_vals)){
  
  # update parameters
  parms_ctrl_temp <- parms_ctrl %>%
    mutate(Estimate = if_else(Parameter == "alpha_PA", alpha_PA_vals[i], Estimate))
  
  # add disease
  modAlphaPA_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
                                       F0 = filter(modInv_disc_fung_fin, species == "Perennial first-year")$N, 
                                       P0 = filter(modInv_disc_fung_fin, species == "Perennial adult")$N, 
                                       L0 = filter(modInv_disc_fung_fin, species == "Litter")$N, 
                                       C0 = conidia_init, 
                                       IA0 = filter(modInv_disc_fung_fin, species == "Annual infected biomass")$N, 
                                       IF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year infected biomass")$N, 
                                       IP0 = filter(modInv_disc_fung_fin, species == "Perennial adult infected biomass")$N, 
                                       BP0 = filter(modInv_disc_fung_fin, species == "Perennial adult biomass")$N, 
                                       BF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year biomass")$N, 
                                       simtime = years/2, gs_time = gs_days, 
                                       disc_parms = parms_ctrl_temp, con_parms = parms_ctrl_temp)
  
  # save data table
  modAlphaPA[[i]] <- modAlphaPA_disc_ctrl %>%
    mutate(alpha_PA = alpha_PA_vals[i])
}

# collapse list
modAlphaPA2 <- modAlphaPA %>%
  bind_rows()

# save file
write_csv(modAlphaPA2, "output/discrete_no_seed_infection_alphaPA_sim.csv")
modAlphaPA2 <- read_csv("output/discrete_no_seed_infection_alphaPA_sim.csv")

# check that each simulation looks reasonable and that biomass stabilizes
for(i in 1:length(alpha_PA_vals)){
  
  print(modAlphaPA2 %>%
          filter(alpha_PA == alpha_PA_vals[i]) %>%
          ggplot(aes(x = time, y = N)) +
          geom_line() +
          theme_classic() +
          facet_wrap(~ species, scales = "free_y"))
  
}

# select last time point
# combine adult & first-year perennial
modAlphaPA3 <- modAlphaPA2 %>%
  filter(species %in% c("Annual biomass",
                        "Perennial adult biomass",
                        "Perennial first-year biomass") &
           time == max(time)) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass,
         Impact = per_bio_alone/Perennial_biomass,
         Coefficient = case_when(alpha_PA == alpha_PA_ctrl ~ "control",
                                 alpha_PA == alpha_PA_fung ~ "fungicide",
                                 TRUE ~ NA_character_))

# figure
alphaPA_fig <- ggplot(modAlphaPA3, aes(x = alpha_PA, y = Perennial_biomass)) +
  geom_hline(yintercept = per_bio_alone, linetype = "dashed", color = "gray60") +
  geom_text(x = max(modAlphaPA3$alpha_PA), y = per_bio_alone - 10, label = "Pre-invasion biomass",
            hjust = 1, vjust = 1, size = 2, check_overlap = T) +
  geom_line(color = "gray60") +
  geom_point(color = "gray60") +
  geom_point(data = filter(modAlphaPA3, !is.na(Coefficient)),
             aes(color = Coefficient), size = 3) +
  scale_color_manual(values = col_pal) +
  labs(x = "Interspecific invader competition", 
       y = expression(paste('Competitor biomass (g/', m^2, ')', sep = "")),
       title = "C") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())


#### competitor intraspecific comp ####

# extract alpha values
alpha_PP_fung <- parms_fung %>%
  filter(Parameter == "alpha_PP") %>%
  pull(Estimate)

alpha_PP_ctrl <- parms_ctrl %>%
  filter(Parameter == "alpha_PP") %>%
  pull(Estimate)

# alpha values
alpha_PP_vals <- c(seq(alpha_PP_fung/5, alpha_PP_ctrl, length.out = 9), alpha_PP_fung) %>%
  unique() %>% sort()

# empty list
modAlphaPP <- vector(mode = "list", length = length(alpha_PP_vals))

# cycle through alpha values
for(i in 1:length(alpha_PP_vals)){
  
  # update parameters
  parms_ctrl_temp <- parms_ctrl %>%
    mutate(Estimate = if_else(Parameter %in% c("alpha_PP", "alpha_FP"), alpha_PP_vals[i], Estimate))
  
  # add disease
  modAlphaPP_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
                                       F0 = filter(modInv_disc_fung_fin, species == "Perennial first-year")$N, 
                                       P0 = filter(modInv_disc_fung_fin, species == "Perennial adult")$N, 
                                       L0 = filter(modInv_disc_fung_fin, species == "Litter")$N, 
                                       C0 = conidia_init, 
                                       IA0 = filter(modInv_disc_fung_fin, species == "Annual infected biomass")$N, 
                                       IF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year infected biomass")$N, 
                                       IP0 = filter(modInv_disc_fung_fin, species == "Perennial adult infected biomass")$N, 
                                       BP0 = filter(modInv_disc_fung_fin, species == "Perennial adult biomass")$N, 
                                       BF0 = filter(modInv_disc_fung_fin, species == "Perennial first-year biomass")$N, 
                                       simtime = years/2, gs_time = gs_days, 
                                       disc_parms = parms_ctrl_temp, con_parms = parms_ctrl_temp)
  
  # save data table
  modAlphaPP[[i]] <- modAlphaPP_disc_ctrl %>%
    mutate(alpha_PP = alpha_PP_vals[i])
}

# collapse list
modAlphaPP2 <- modAlphaPP %>%
  bind_rows()

# save file
write_csv(modAlphaPP2, "output/discrete_no_seed_infection_alphaPP_sim.csv")
modAlphaPP2 <- read_csv("output/discrete_no_seed_infection_alphaPP_sim.csv")

# check that each simulation looks reasonable and that biomass stabilizes
for(i in 1:length(alpha_PP_vals)){
  
  print(modAlphaPP2 %>%
          filter(alpha_PP == alpha_PP_vals[i]) %>%
          ggplot(aes(x = time, y = N)) +
          geom_line() +
          theme_classic() +
          facet_wrap(~ species, scales = "free_y"))
  
}

# select last time point
# combine adult & first-year perennial
modAlphaPP3 <- modAlphaPP2 %>%
  filter(species %in% c("Annual biomass",
                        "Perennial adult biomass",
                        "Perennial first-year biomass") &
           time == max(time)) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass,
         Impact = per_bio_alone/Perennial_biomass,
         Coefficient = case_when(alpha_PP == alpha_PP_ctrl ~ "control",
                                 alpha_PP == alpha_PP_fung ~ "fungicide",
                                 TRUE ~ NA_character_))

# figure
alphaPP_fig <- ggplot(modAlphaPP3, aes(x = alpha_PP, y = Perennial_biomass)) +
  geom_hline(yintercept = per_bio_alone, linetype = "dashed", color = "gray60") +
  geom_text(x = max(modAlphaPP3$alpha_PP), y = per_bio_alone + 20, label = "Pre-invasion biomass",
            hjust = 1, vjust = 0, size = 2, check_overlap = T) +
  geom_line(color = "gray60") +
  geom_point(color = "gray60") +
  geom_point(data = filter(modAlphaPP3, !is.na(Coefficient)),
             aes(color = Coefficient), size = 3) +
  scale_color_manual(values = col_pal) +
  labs(x = "Intraspecific competitor competition", 
       y = expression(paste('Competitor biomass (g/', m^2, ')', sep = "")),
       title = "D") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())


#### figure ####

ggsave("output/discrete_no_seed_infection_simulation_figure.pdf", 
       device = "pdf", width = 15, height = 12, units = "cm")
time_series_fig / (alphaAA_fig + alphaPA_fig + alphaPP_fig)
dev.off()


#### add disease to perennial alone ####

# add 100 years for comparison to disease
modF_disc_fung2 <- disc_AFP_mod(A0 = 0, F0 = 10, P0 = 0, L0 = 0, C0 = 0, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = 0, BF0 = 0, 
                                simtime = years * 3, gs_time = gs_days, disc_parms = parms_fung, con_parms = parms_fung)

# established perennial
modF_disc_fung_fin2 <- modF_disc_fung2 %>%
  filter(time == 200)

# add disease
modDisF_disc_ctrl <- disc_AFP_mod(A0 = 0, 
                                  F0 = filter(modF_disc_fung_fin2, species == "Perennial first-year")$N, 
                                  P0 = filter(modF_disc_fung_fin2, species == "Perennial adult")$N, 
                                  L0 = filter(modF_disc_fung_fin2, species == "Litter")$N, 
                                  C0 = conidia_init, 
                                  IA0 = 0, 
                                  IF0 = filter(modF_disc_fung_fin2, species == "Perennial first-year infected biomass")$N, 
                                  IP0 = filter(modF_disc_fung_fin2, species == "Perennial adult infected biomass")$N, 
                                  BP0 = filter(modF_disc_fung_fin2, species == "Perennial adult biomass")$N, 
                                  BF0 = filter(modF_disc_fung_fin2, species == "Perennial first-year biomass")$N, 
                                  simtime = years, gs_time = gs_days, 
                                  disc_parms = parms_ctrl, con_parms = parms_ctrl)

# combine model outputs
mod_F_dis <- modF_disc_fung2 %>%
  mutate(treatment = "fungicide") %>%
  full_join(modDisF_disc_ctrl %>%
              mutate(time = time + 200,
                     treatment = "control"))

mod_F_dis %>%
  ggplot(aes(x = time, y = N, color = treatment)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

# make discontinuous transitions continuous
bio_F_dis <- mod_F_dis %>%
  filter(time == 200 & 
           (str_detect(species, "biomass") == T | species == "Conidia")) %>%
  mutate(treatment = "control")

# add transitions
mod_F_dis2 <- mod_F_dis %>%
  full_join(bio_F_dis)

# figure of discrete results
mod_F_dis2 %>%
  ggplot(aes(x = time, y = N, color = treatment)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

# combine adult & first-year perennial
mod_F_dis3 <- mod_F_dis2 %>%
  filter(species %in% c("Perennial adult biomass",
                        "Perennial first-year biomass")) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass) %>%
  select(time, treatment, Perennial_biomass)

# save
write_csv(mod_F_dis3, "output/discrete_no_seed_infection_perennial_only_time_series.csv")
mod_F_dis3 <- read_csv("output/discrete_no_seed_infection_perennial_only_time_series.csv")

# phase labels
phase_lab_F <- tibble(time = c(100, 250),
                      Perennial_biomass = 400,
                      lab = c("pre-disease", "disease emergence"))

# figure of discrete results
mod_F_dis3 %>%
  ggplot(aes(x = time, y = Perennial_biomass)) +
  geom_rect(xmin = 200, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray75", color = NA, show.legend = F) +
  geom_line(aes(color = treatment), linetype = "dashed") +
  geom_text(data = phase_lab_F, check_overlap = T, size = 2.5,
            aes(label = lab)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  labs(x = "Time (years)", y = expression(paste('Biomass (g/', m^2, ')', sep = ""))) +
  fig_theme


#### values for text ####

# extract perennial biomass values from invasion simulation
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

inv_bio2 <- mod_inv_dis3 %>%
  filter(time == 10 & species == "Annual_biomass") %>%
  pull(N)

dis_inv_bio2 <- mod_inv_dis3 %>%
  filter(time == 110 & species == "Annual_biomass" &
           treatment == "control") %>%
  pull(N)

# invasion effect
100 * (inv_bio - pre_inv_bio)/pre_inv_bio

# disease effect
100 * (dis_inv_bio - inv_bio)/inv_bio
100 * (dis_inv_bio2 - inv_bio2)/inv_bio2

# extract perennial biomass values from perennial alone simulation
pre_dis_bio <- mod_F_dis3 %>%
  filter(time == 199) %>%
  pull(Perennial_biomass)

post_dis_bio <- mod_F_dis3 %>%
  filter(time == 300 & treatment == "control") %>%
  pull(Perennial_biomass)

# disease effect
100 * (post_dis_bio - pre_dis_bio)/pre_dis_bio
