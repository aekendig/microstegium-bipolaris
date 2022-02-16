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
cont_params <- read_csv("output/continuous_model_parameters_2018_2019_density_exp.csv")
disc_params <- read_csv("output/discrete_model_parameters_2018_2019_density_exp.csv")
# model_parameters_2018_2019_density_exp.R

# continuous models
source("code/continuous_AFP_model.R")
source("code/continuous_AF_model.R")
source("code/continuous_AP_model.R")
source("code/continuous_FP_model.R")
source("code/continuous_A_model.R")
source("code/continuous_F_model.R")
source("code/continuous_P_model.R")

# discrete model
source("code/discrete_model.R")

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
        plot.title = element_text(size = 10, hjust = -0.1, face = "bold"),
        plot.subtitle = element_text(size = 7))

col_pal = c("black", "#238A8D")


#### parameters ####

# growing season days
gs_days <- 160
gs_days_seq <- seq(0, gs_days, by = 1)

# simulation years
years <- 100

# conidia per 1 g litter
conidia_init <- as.numeric(disc_params[disc_params$Parameter == "h", "Estimate"])

# initial plant sizes
bio_A_init <- as.numeric(disc_params[disc_params$Parameter == "b_A", "Estimate"])
bio_F_init <- as.numeric(disc_params[disc_params$Parameter == "b_F", "Estimate"])
bio_P_init <- as.numeric(disc_params[disc_params$Parameter == "b_P", "Estimate"])

# control and fungicide continuous parameters
cont_params_ctrl <- cont_params %>%
  filter(Treatment == "control" | is.na(Treatment))

cont_params_fung <- cont_params %>%
  filter(Treatment == "fungicide" | is.na(Treatment))


#### annual alone ####

# continuous input values
startA_cont <- c(LogB_A = log(bio_A_init * 10), I_A = 0, D = 0, C = conidia_init)

# apply continuous function
modA_cont_ctrl <- ode(y = startA_cont, times = gs_days_seq, 
                      func = cont_A_mod, parms = cont_params_ctrl) %>%
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

modA_cont_fung <- ode(y = startA_cont, times = gs_days_seq, 
                      func = cont_A_mod, parms = cont_params_fung) %>%
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
modA_disc_ctrl <- disc_AFP_mod(A0 = 10, F0 = 0, P0 = 0, L0 = 0, C0 = conidia_init, BP0 = 0, BF0 = 0, BA0 = 0, simtime = years, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

modA_disc_fung <- disc_AFP_mod(A0 = 10, F0 = 0, P0 = 0, L0 = 0, C0 = conidia_init, BP0 = 0, BF0 = 0, BA0 = 0, simtime = years, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

# save simulations
write_csv(modA_disc_ctrl, "output/discrete_continuous_model_A_ctrl_sim.csv")
write_csv(modA_disc_fung, "output/discrete_continuous_model_A_fung_sim.csv")

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
startF_cont <- c(LogB_F = log(bio_F_init * 10), I_F = 0, D = 0, C = conidia_init)

# apply continuous function
modF_cont_ctrl <- ode(y = startF_cont, times = gs_days_seq, 
                      func = cont_F_mod, parms = cont_params_ctrl) %>%
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

modF_cont_fung <- ode(y = startF_cont, times = gs_days_seq, 
                      func = cont_F_mod, parms = cont_params_fung) %>%
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
modF_disc_ctrl <- disc_AFP_mod(A0 = 0, F0 = 10, P0 = 0, L0 = 0, C0 = conidia_init, BP0 = 0, BF0 = 0, BA0 = 0, simtime = years, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

modF_disc_fung <- disc_AFP_mod(A0 = 0, F0 = 10, P0 = 0, L0 = 0, C0 = conidia_init, BP0 = 0, BF0 = 0, BA0 = 0, simtime = years, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

# save simulations
write_csv(modF_disc_ctrl, "output/discrete_continuous_model_F_ctrl_sim.csv")
write_csv(modF_disc_fung, "output/discrete_continuous_model_F_fung_sim.csv")

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


#### left off here with last parameter set ####


#### perennial adult alone ####

# continuous input values
startP_cont <- c(LogB_P = log(bio_P_init * 10), I_P = 0, D = 0, C = conidia_init)

# apply continuous function
modP_cont_ctrl <- ode(y = startP_cont, times = gs_days_seq, 
                      func = cont_P_mod, parms = cont_params_ctrl) %>%
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

modP_cont_fung <- ode(y = startP_cont, times = gs_days_seq, 
                      func = cont_P_mod, parms = cont_params_fung) %>%
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
modP_disc_ctrl <- disc_AFP_mod(A0 = 0, F0 = 0, P0 = 10, L0 = 0, C0 = conidia_init, BP0 = bio_P_init * 10, BF0 = 0, BA0 = 0, simtime = years, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

modP_disc_fung <- disc_AFP_mod(A0 = 0, F0 = 0, P0 = 10, L0 = 0, C0 = conidia_init, BP0 = bio_P_init * 10, BF0 = 0, BA0 = 0, simtime = years, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

# save simulations
write_csv(modP_disc_ctrl, "output/discrete_continuous_model_P_ctrl_sim.csv")
write_csv(modP_disc_fung, "output/discrete_continuous_model_P_fung_sim.csv")

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
                                 disc_parms = disc_params, cont_parms = cont_params)
write_csv(modInv_disc_fung, "output/discrete_continuous_model_invasion_no_disease.csv")
modInv_disc_fung <- read.csv("output/discrete_continuous_model_invasion_no_disease.csv")

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
                                 disc_parms = disc_params, cont_parms = cont_params)

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
# create a relative biomass
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
           fct_relevel("invader"),
         RelN = N / max(N))

# save
write_csv(mod_inv_dis3, "output/discrete_no_seed_infection_time_series.csv")
mod_inv_dis3 <- read_csv("output/discrete_no_seed_infection_time_series.csv") %>%
  mutate(Species = fct_relevel(Species, "invader"))

# phase labels
phase_lab <- tibble(time = c(-5, 5, 107),
                    RelN = 1.1,
                    lab = c("Pre-invasion", "Invasion", "Disease emergence"))

# figure of discrete results
time_series_fig <- mod_inv_dis3 %>%
  ggplot(aes(x = time, y = RelN)) +
  geom_rect(xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "gray95", color = NA, show.legend = F) +
  geom_rect(xmin = 100, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray75", color = NA, show.legend = F) +
  geom_line(aes(color = treatment, linetype = Species)) +
  geom_text(data = phase_lab, check_overlap = T, size = 2.5,
            aes(label = lab)) +
  scale_color_manual(values = col_pal, name = "Disease",
                     labels = c("present", "absent")) +
  scale_x_break(c(10, 100)) +
  labs(x = "Time (years)", y = "Relative biomass",
       title = "A") +
  fig_theme + 
  theme(plot.title = element_text(size = 10, hjust = 0, face = "bold"))


#### invader abundance and per capita effects ####

# extract alpha values
alpha_AA_fung <- parms_fung %>%
  filter(Parameter == "alpha_AA") %>%
  pull(Estimate)

alpha_AA_ctrl <- parms_ctrl %>%
  filter(Parameter == "alpha_AA") %>%
  pull(Estimate)

alpha_PA_fung <- parms_fung %>%
  filter(Parameter == "alpha_PA") %>%
  pull(Estimate)

alpha_PA_ctrl <- parms_ctrl %>%
  filter(Parameter == "alpha_PA") %>%
  pull(Estimate)

# alpha values
alpha_AA_vals <- c(seq(alpha_AA_fung, 50 * alpha_AA_ctrl, length.out = 10),
                   alpha_AA_ctrl) %>%
  unique() %>% sort()

alpha_PA_vals <- c(seq(alpha_PA_fung, alpha_PA_ctrl, length.out = 10))

alpha_AA_PA_vals <- tibble(alpha_AA = alpha_AA_vals) %>%
  expand_grid(tibble(alpha_PA = alpha_PA_vals))

# empty list
mod_AA_PA <- vector(mode = "list", length = nrow(alpha_AA_PA_vals))

# cycle through alpha values
for(i in 1:nrow(alpha_AA_PA_vals)){
  
  # update parameters
  parms_ctrl_temp <- parms_ctrl %>%
    mutate(Estimate = case_when(Parameter == "alpha_AA" ~ alpha_AA_PA_vals$alpha_AA[i], 
                                Parameter == "alpha_PA" ~ alpha_AA_PA_vals$alpha_PA[i],
                                TRUE ~ Estimate))
  
  # add disease
  mod_AA_PA_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
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
                                       disc_parms = parms_ctrl_temp, cont_parms = parms_ctrl_temp)
  
  # save data table
  mod_AA_PA[[i]] <- mod_AA_PA_disc_ctrl %>%
    mutate(alpha_AA = alpha_AA_PA_vals$alpha_AA[i],
           alpha_PA = alpha_AA_PA_vals$alpha_PA[i])
}

# collapse list
mod_AA_PA2 <- mod_AA_PA %>%
  bind_rows()

# save file
write_csv(mod_AA_PA2, "output/discrete_no_seed_infection_alphaAA_alphaPA_sim.csv")
mod_AA_PA2 <- read_csv("output/discrete_no_seed_infection_alphaAA_alphaPA_sim.csv")

# check that each simulation looks reasonable and that biomass stabilizes
pdf("output/discrete_no_seed_infection_alphaAA_alphaPA_sim.pdf")
for(i in 1:nrow(alpha_AA_PA_vals)){
  
  print(mod_AA_PA2 %>%
          filter(alpha_AA == alpha_AA_PA_vals$alpha_AA[i] &
                   alpha_PA == alpha_AA_PA_vals$alpha_PA[i]) %>%
          ggplot(aes(x = time, y = N)) +
          geom_line() +
          theme_classic() +
          facet_wrap(~ species, scales = "free_y"))
  
}
dev.off()

# extract equilibrium perennial biomass
per_bio_alone <- modF_disc_fung_fin %>%
  filter(species %in% c("Perennial adult biomass",
                        "Perennial first-year biomass")) %>%
  summarize(bio = sum(N)) %>%
  pull(bio)

# per_bio_alone <- 719.323

# select last time point
# combine adult & first-year perennial
mod_AA_PA3 <- mod_AA_PA2 %>%
  filter(species %in% c("Annual",
                        "Annual biomass",
                        "Perennial adult biomass",
                        "Perennial first-year biomass") &
           time == max(time)) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass,
         Impact = per_bio_alone/Perennial_biomass,
         Coefficient = case_when(alpha_AA == alpha_AA_ctrl & alpha_PA == alpha_PA_ctrl ~ "control",
                                 alpha_AA == alpha_AA_fung & alpha_PA == alpha_PA_fung ~ "fungicide",
                                 TRUE ~ NA_character_),
         RelAnnual_biomass = Annual_biomass / per_bio_alone)

# figure
alpha_AA_PA_fig <- ggplot(mod_AA_PA3, aes(x = RelAnnual_biomass, y = Impact, by = as.factor(alpha_PA))) +
  geom_hline(yintercept = 1, size = 0.25, linetype = "dashed") +
  geom_line(aes(color = alpha_PA)) +
  geom_point(data = filter(mod_AA_PA3, !is.na(Coefficient)),
             aes(shape = Coefficient), size = 2, show.legend = F) +
  scale_color_viridis_c(name = "Invader\neffect",
                        direction = -1,
                        breaks = c(min(mod_AA_PA3$alpha_PA),
                                    max(mod_AA_PA3$alpha_PA)),
                         labels = c("least", "most")) +
  labs(x = "Invader abundance (relative biomass)", 
       y = expression(paste("Impact (competitor relative ", biomass^-1, ")", sep = "")),
       title = "B") +
  fig_theme +
  theme(legend.margin = margin(-0.1, 0, 0, 0, unit = "cm"),
        legend.justification = "bottom") +
  coord_cartesian(ylim = c(0, max(mod_AA_PA3$Impact))) +
  guides(color = guide_colourbar(ticks = F, barheight = 6.5))


#### invader abundance and competitor per capita effects ####

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

alpha_AA_PP_vals <- tibble(alpha_AA = alpha_AA_vals) %>%
  expand_grid(tibble(alpha_PP = alpha_PP_vals))

# empty list
mod_AA_PP <- vector(mode = "list", length = nrow(alpha_AA_PP_vals))

# cycle through alpha values
for(i in 1:nrow(alpha_AA_PP_vals)){
  
  # update parameters
  parms_ctrl_temp <- parms_ctrl %>%
    mutate(Estimate = case_when(Parameter == "alpha_AA" ~ alpha_AA_PP_vals$alpha_AA[i], 
                                Parameter %in% c("alpha_PP", "alpha_FP") ~ alpha_AA_PP_vals$alpha_PP[i],
                                TRUE ~ Estimate))
  
  # add disease
  mod_AA_PP_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
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
                                      disc_parms = parms_ctrl_temp, cont_parms = parms_ctrl_temp)
  
  # save data table
  mod_AA_PP[[i]] <- mod_AA_PP_disc_ctrl %>%
    mutate(alpha_AA = alpha_AA_PP_vals$alpha_AA[i],
           alpha_PP = alpha_AA_PP_vals$alpha_PP[i])
}

# collapse list
mod_AA_PP2 <- mod_AA_PP %>%
  bind_rows()

# save file
write_csv(mod_AA_PP2, "output/discrete_no_seed_infection_alphaAA_alphaPP_sim.csv")
mod_AA_PP2 <- read_csv("output/discrete_no_seed_infection_alphaAA_alphaPP_sim.csv")

# check that each simulation looks reasonable and that biomass stabilizes
pdf("output/discrete_no_seed_infection_alphaAA_alphaPP_sim.pdf")
for(i in 1:nrow(alpha_AA_PP_vals)){
  
  print(mod_AA_PP2 %>%
          filter(alpha_AA == alpha_AA_PP_vals$alpha_AA[i] &
                   alpha_PP == alpha_AA_PP_vals$alpha_PP[i]) %>%
          ggplot(aes(x = time, y = N)) +
          geom_line() +
          theme_classic() +
          facet_wrap(~ species, scales = "free_y"))
  
}
dev.off()

# per_bio_alone <- 719.323

# select last time point
# combine adult & first-year perennial
mod_AA_PP3 <- mod_AA_PP2 %>%
  filter(species %in% c("Annual",
                        "Annual biomass",
                        "Perennial adult biomass",
                        "Perennial first-year biomass") &
           time == max(time)) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = str_replace_all(species, "-", "_")) %>%
  pivot_wider(names_from = "species",
              values_from = "N") %>%
  mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass,
         Impact = per_bio_alone/Perennial_biomass,
         Coefficient = case_when(alpha_AA == alpha_AA_ctrl & alpha_PP == alpha_PP_ctrl ~ "control",
                                 alpha_AA == alpha_AA_fung & alpha_PP == alpha_PP_fung ~ "fungicide",
                                 TRUE ~ NA_character_),
         RelAnnual_biomass = Annual_biomass / per_bio_alone) %>%
  filter(alpha_PP < 0.001 | alpha_PP > 0.0013) # remove value that's redundant with fungicide value

# figure
alpha_AA_PP_fig <- ggplot(mod_AA_PP3, aes(x = RelAnnual_biomass, y = Impact, by = as.factor(alpha_PP))) +
  geom_hline(yintercept = 1, size = 0.25, linetype = "dashed") +
  geom_line(aes(color = alpha_PP)) +
  geom_point(data = filter(mod_AA_PP3, !is.na(Coefficient)),
             aes(shape = Coefficient), size = 2) +
  scale_shape(name = "Disease\ntreatment") +
  scale_color_viridis_c(name = "Competitor\neffect",
                        direction = -1,
                        breaks = c(min(mod_AA_PP3$alpha_PP),
                                   max(mod_AA_PP3$alpha_PP)),
                        labels = c("least", "most")) +
  labs(x = "Invader abundance (relative biomass)", 
       title = "C") +
  coord_cartesian(ylim = c(0, max(mod_AA_PA3$Impact))) +
  fig_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.margin = margin(-0.1, 0, 0, 0, unit = "cm"),
        legend.box = "horizontal",
        legend.justification = "bottom") +
  guides(color = guide_colourbar(ticks = F, barheight = 6.5))


#### figure ####

fig_out <- time_series_fig / (alpha_AA_PA_fig + alpha_AA_PP_fig)
ggsave("output/discrete_no_seed_infection_simulation_figure.pdf", 
       plot = fig_out,
       device = "pdf", width = 14, height = 11, units = "cm")


#### add disease to perennial alone ####

# add 100 years for comparison to disease
modF_disc_fung2 <- disc_AFP_mod(A0 = 0, F0 = 10, P0 = 0, L0 = 0, C0 = conidia_init, BP0 = 0, BF0 = 0, 
                                simtime = years * 3, gs_time = gs_days, disc_parms = disc_params, cont_parms = cont_params)

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
                                  disc_parms = disc_params, cont_parms = cont_params)

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


#### sensitivity analysis functions ####

sens_sim <- function(parms_vary, parms_range, parms_n, filename){
  
  # parameter dataframe
  parms_df <- tibble(Parameter = parms_vary) %>%
    expand_grid(tibble(Estimate = seq(parms_range[1], parms_range[2], 
                                      length.out = parms_n)))
  
  # empty list
  mod_out <- vector(mode = "list", length = nrow(parms_df))
  
  # open figure
  pdf(paste0("output/", filename, ".pdf"))
  
  # cycle through alpha values
  for(i in 1:nrow(parms_df)){
    
    # parameter to change and its value
    parms_sub <- parms_df[i, ]
    
    # update parameters
    parms_ctrl_temp <- parms_ctrl %>%
      filter(!(Parameter %in% parms_sub$Parameter)) %>%
      full_join(parms_sub)
    
    # add disease
    mod_out_disc_ctrl <- disc_AFP_mod(A0 = filter(modInv_disc_fung_fin, species == "Annual")$N, 
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
                                      disc_parms = parms_ctrl_temp, cont_parms = parms_ctrl_temp)
    
    # create figure
    print(ggplot(mod_out_disc_ctrl, aes(x = time, y = N)) +
            geom_line() +
            theme_classic() +
            facet_wrap(~ species, scales = "free_y") +
            labs(title = paste0(parms_sub$Parameter, "=", parms_sub$Estimate)))
    
    # save data table
    mod_out[[i]] <- mod_out_disc_ctrl %>%
      mutate(Parameter = parms_sub$Parameter,
             Estimate = parms_sub$Estimate)
  }
  
  # stop figure device
  dev.off()
  
  # collapse list
  mod_out2 <- mod_out %>%
    bind_rows()
  
  # save file
  write_csv(mod_out2, paste0("output/", filename, ".csv"))
  
}

sens_sim_dat <- function(dat_in) {
  
  # recovered perennial biomass
  per_bio_alone <- 719.323
  
  # select last date
  # combine perennial life stages
  dat_out <- dat_in %>%
    filter(species %in% c("Perennial adult biomass",
                          "Perennial first-year biomass") &
             time == max(time)) %>%
    mutate(species = str_replace_all(species, " ", "_"),
           species = str_replace_all(species, "-", "_")) %>%
    pivot_wider(names_from = "species",
                values_from = "N") %>%
    mutate(Perennial_biomass = Perennial_first_year_biomass + Perennial_adult_biomass,
           Impact = per_bio_alone/Perennial_biomass)
  
  return(dat_out)
  
}


#### sensitivity analysis: survival ####

# run simulations
sens_sim(parms_vary = c("e_A", "e_P", "l_P"),
         parms_range = c(0.5, 1), 
         parms_n = 10, 
         filename = "sensitivity_analysis_e")

# import simulations
sens_e <- read_csv("output/sensitivity_analysis_e.csv")

# process data
sens_e2 <- sens_e %>%
  sens_sim_dat() %>%
  mutate(Plant_group = case_when(Parameter == "e_A" ~ "invader",
                                 Parameter == "e_P" ~ "1st yr competitor",
                                 Parameter == "l_P" ~ "adult competitor") %>%
           fct_relevel("invader"))

# figure
sens_e_fig <- ggplot(sens_e2, aes(x = Estimate, y = Impact, color = Plant_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  labs(y = expression(paste("Impact (competitor relative ", biomass^-1, ")", sep = "")),
       title = "A",
       x = "Establishment/survival") +
scale_color_viridis_d(option = "A", end = 0.9, name = "Plant group") +
  fig_theme +
  theme(plot.title = element_text(size = 10, hjust = -0.04, face = "bold"))


#### sensitivity analysis: growth rate ####

# run simulations
sens_sim(parms_vary = c("r_A", "r_F", "r_P"),
         parms_range = c(0.01, 0.05), 
         parms_n = 10, 
         filename = "sensitivity_analysis_r")

# import simulations
sens_r <- read_csv("output/sensitivity_analysis_r.csv")

# process data
sens_r2 <- sens_r %>%
  sens_sim_dat() %>%
  mutate(Plant_group = case_when(Parameter == "r_A" ~ "invader",
                                 Parameter == "r_F" ~ "1st yr competitor",
                                 Parameter == "r_P" ~ "adult competitor") %>%
           fct_relevel("invader"))

# figure
sens_r_fig <- ggplot(sens_r2, aes(x = Estimate, y = Impact, color = Plant_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  coord_cartesian(ylim = c(0, 20)) +
  labs(title = "B",
       x = "Growth rate") +
  scale_color_viridis_d(option = "A", end = 0.9) +
  fig_theme +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10, hjust = -0.04, face = "bold"))


#### sensitivity analysis: comp. effect A ####

# run simulations
sens_sim(parms_vary = c("alpha_AA", "alpha_FA"),
         parms_range = c(0, 0.003), 
         parms_n = 10, 
         filename = "sensitivity_analysis_alphaiA")

# accidentally left this out above
sens_sim(parms_vary = c("alpha_PA"),
         parms_range = c(0, 0.003), 
         parms_n = 10, 
         filename = "sensitivity_analysis_alphaPA")

sens_alphaFAA <- read_csv("output/sensitivity_analysis_alphaiA.csv")
sens_alphaPA <- read_csv("output/sensitivity_analysis_alphaPA.csv")
sens_alphaiA <- full_join(sens_alphaFAA, sens_alphaPA)

# process data
sens_alphaiA2 <- sens_alphaiA %>%
  sens_sim_dat() %>%
  mutate(Plant_group = case_when(Parameter == "alpha_AA" ~ "invader",
                                 Parameter == "alpha_FA" ~ "1st yr competitor",
                                 Parameter == "alpha_PA" ~ "adult competitor") %>%
           fct_relevel("invader"))

# alpha_FA > 0.001 and alpha_PA = 0 causes first year susceptible biomass <0, but it's ~0 
# alpha_AA = 0 causes annual susceptible biomass to fluctuate, but it's very small, especially relative to annual infected biomass

# figure
sens_alphaiA_fig <- ggplot(sens_alphaiA2, aes(x = Estimate, y = Impact, color = Plant_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  coord_cartesian(ylim = c(0, 20)) +
  labs(title = "C",
       x = "Invader effects") +
  scale_color_viridis_d(option = "A", end = 0.9) +
  fig_theme +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10, hjust = -0.04, face = "bold"))


#### sensitivity analysis: seed production ####

# run simulations
sens_sim(parms_vary = c("c_A", "c_F", "c_P"),
         parms_range = c(6, 80), 
         parms_n = 10, 
         filename = "sensitivity_analysis_c")

# import simulations
sens_c <- read_csv("output/sensitivity_analysis_c.csv")

# process data
sens_c2 <- sens_c %>%
  sens_sim_dat() %>%
  mutate(Plant_group = case_when(Parameter == "c_A" ~ "invader",
                                 Parameter == "c_F" ~ "1st yr competitor",
                                 Parameter == "c_P" ~ "adult competitor") %>%
           fct_relevel("invader"))

# figure
sens_c_fig <- ggplot(sens_c2, aes(x = Estimate, y = Impact, color = Plant_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  labs(title = "D",
       x = "Seed production") +
  scale_color_viridis_d(option = "A", end = 0.9) +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = -0.04, face = "bold"))


#### sensitivity analysis: virulence ####

# run simulations
sens_sim(parms_vary = c("v_A", "v_F", "v_P"),
         parms_range = c(0, 0.01), 
         parms_n = 10, 
         filename = "sensitivity_analysis_v")

# import simulations
sens_v <- read_csv("output/sensitivity_analysis_v.csv")

# process data
sens_v2 <- sens_v %>%
  sens_sim_dat() %>%
  mutate(Plant_group = case_when(Parameter == "v_A" ~ "invader",
                                 Parameter == "v_F" ~ "1st yr competitor",
                                 Parameter == "v_P" ~ "adult competitor") %>%
           fct_relevel("invader"))

# figure
sens_v_fig <- ggplot(sens_v2, aes(x = Estimate, y = Impact, color = Plant_group)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  labs(title = "E",
       x = "Virulence") +
  scale_color_viridis_d(option = "A", end = 0.9) +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = -0.04, face = "bold"))


#### sensitivity analysis figure ####
pdf("output/discrete_no_seed_infection_sensitivity_analysis.pdf", width = 5.5, height = 3.5)
sens_e_fig + theme(legend.position = "none",
                   axis.title.y = element_blank()) + 
  sens_r_fig + sens_alphaiA_fig +
  sens_c_fig + sens_v_fig + 
  get_legend(sens_e_fig) +
  plot_layout(nrow = 2)
grid::grid.draw(grid::textGrob(sens_e_fig$labels$y, x = 0.01, rot = 90,
                               gp = grid::gpar(fontsize = 7)))
dev.off()