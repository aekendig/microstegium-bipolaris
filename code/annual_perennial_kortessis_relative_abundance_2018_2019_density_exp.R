##### info ####

# file: annual_perennial_kortessis_relative_abundance_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/12/20
# goal: relative abundance of annual to perennial and invasion conditions


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidybayes)

# number of iterations
n_samps <- 300

# call model script
source("code/annual_perennial_kortessis_model_2018_2019_density_exp.R")

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
  
  mod_estE <- sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 0, iter, d)
  S0Ei <- mod_estE[[1]] %>% filter(time == simtimeR & species == "Perennial seedling") %>% pull(N)
  P0Ei <- mod_estE[[1]] %>% filter(time == simtimeR & species == "Perennial adult") %>% pull(N)
  L0Ei <- mod_estE[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  mod_invE <- sim_fun(A0 = 10, S0Ei, P0Ei, L0Ei, simtimeI, disease = 0, iter, d)
  abundE[[iter]] <- mod_invE[[1]]
  paramsE[[iter]] <- mod_invE[[2]]
  
  mod_estE2 <- sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 1, iter, d)
  S0Ei2 <- mod_estE2[[1]] %>% filter(time == simtimeR & species == "Perennial seedling") %>% pull(N)
  P0Ei2 <- mod_estE2[[1]] %>% filter(time == simtimeR & species == "Perennial adult") %>% pull(N)
  L0Ei2 <- mod_estE2[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  mod_invE2 <- sim_fun(A0 = 10, S0Ei2, P0Ei2, L0Ei2, simtimeI, disease = 1, iter, d)
  abundE2[[iter]] <- mod_invE2[[1]]
  paramsE2[[iter]] <- mod_invE2[[2]]
  
}

# convert to dataframes
abund_datE <- do.call(rbind, abundE)
param_datE <- do.call(rbind, paramsE)
abund_datE2 <- do.call(rbind, abundE2)
param_datE2 <- do.call(rbind, paramsE2)


#### resident M. vimineum ####

# initial conditions
A0M <- 1 # initial annual population size
S0M <- 0 # initial perennial seedling population size
P0M <- 0 # initial perennial adult population size
L0M <- 0 # initial annual litter amount

# initiate lists
abundM <- list()
paramsM <- list()
abundM2 <- list()
paramsM2 <- list()

# simulation without disease
for(iter in 1:samps){
  
  mod_estM <- sim_fun(A0M, S0M, P0M, L0M, simtimeR, disease = 0, iter, d)
  A0Mi <- mod_estM[[1]] %>% filter(time == simtimeR & species == "Annual") %>% pull(N)
  L0Mi <- mod_estM[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  mod_invM <- sim_fun(A0Mi, S0 = 5, P0 = 5, L0Mi, simtimeI, disease = 0, iter, d)
  abundM[[iter]] <- mod_invM[[1]]
  paramsM[[iter]] <- mod_invM[[2]]
  
  mod_estM2 <- sim_fun(A0M, S0M, P0M, L0M, simtimeR, disease = 1, iter, d)
  A0Mi2 <- mod_estM2[[1]] %>% filter(time == simtimeR & species == "Annual") %>% pull(N)
  L0Mi2 <- mod_estM2[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  mod_invM2 <- sim_fun(A0Mi, S0 = 5, P0 = 5, L0Mi2, simtimeI, disease = 1, iter, d)
  abundM2[[iter]] <- mod_invM2[[1]]
  paramsM2[[iter]] <- mod_invM2[[2]]
  
}

# convert to dataframes
abund_datM <- do.call(rbind, abundM)
param_datM <- do.call(rbind, paramsM)
abund_datM2 <- do.call(rbind, abundM2)
param_datM2 <- do.call(rbind, paramsM2)


#### process data ####

# combine abundance data
adat <- abund_datE %>%
  mutate(disease = "without disease",
         resident = "E. virginicus") %>%
  full_join(abund_datE2 %>%
              mutate(disease = "with disease",
                     resident = "E. virginicus")) %>%
  full_join(abund_datM %>%
              mutate(disease = "without disease",
                     resident = "M. vimineum")) %>%
  full_join(abund_datM2 %>%
              mutate(disease = "with disease",
                     resident = "M. vimineum")) %>%
  pivot_wider(names_from = species, values_from = N) %>%
  rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
  mutate(Perennial = Perennial_seedling + Perennial_adult,
         tot_abund = Annual + Perennial,
         ann_rel_abund = Annual / tot_abund,
         per_rel_abund = Perennial / tot_abund,
         iterationF = as.factor(iteration))

# time series
tdat <- adat %>%
  pivot_longer(cols = c(Perennial_seedling, Perennial_adult, Annual, Litter, Perennial, tot_abund, ann_rel_abund, per_rel_abund),
               names_to = "species",
               values_to = "density")

# combine parameter data
pdat <- param_datE %>%
  mutate(disease = "without disease",
         resident = "E. virginicus") %>%
  full_join(param_datE2 %>%
              mutate(disease = "with disease",
                     resident = "E. virginicus")) %>%
  full_join(param_datM %>%
              mutate(disease = "without disease",
                     resident = "M. vimineum")) %>%
  full_join(param_datM2 %>%
              mutate(disease = "with disease",
                     resident = "M. vimineum")) %>%
  pivot_wider(names_from = parameter, values_from = value) 

# invasion criteria
idat <- pdat %>%
  filter(time == 1) %>%
  mutate(ann_invades = ifelse(l.A / l.P > (1 + beta.A * L.P) / (1 + beta.P * L.P), 1, 0),
         per_invades = ifelse(l.P / l.A > (1 + beta.P * L.A) / (1 + beta.A * L.A), 1, 0),
         per_viable = ifelse(l.P > 1 + beta.P * L.P, 1, 0),
         outcome = case_when(ann_invades == 1 & per_invades == 0 ~ "M. vimineum only",
                             ann_invades == 1 & per_invades == 1  ~ "coexistence",
                             ann_invades == 0 & per_invades == 0  ~ "priority effects",
                             ann_invades == 0 & per_invades == 1 ~ "E. virginicus only",
                             TRUE ~ NA_character_),
         litter_producer = case_when(L.A > L.P ~ "Producer: M. vimineum",
                                     L.P > L.A ~ "Producer: E. virginicus",
                                     L.P == L.A ~ "Equal producers"),
         litter_tolerator = case_when(beta.A > beta.P ~ "Tolerant: M. vimineum",
                                      beta.P > beta.A ~ "Tolerant: E. virginicus",
                                      beta.P == beta.A ~ "Equal tolerators"),
         outcome = case_when(is.na(outcome) & l.P < 1 ~ "M. vimineum only",
                             is.na(outcome) & l.A < 1 ~ "E. virginicus only",
                             TRUE ~ outcome),
         litter_producer = case_when(is.na(litter_producer) & l.P < 1 ~ "Producer: M. vimineum",
                                     is.na(litter_producer) & l.A < 1 ~ "Producer: E. virginicus",
                                     TRUE ~ litter_producer)) %>%
  select(-time)

# relative abundance summary
rdat <- adat %>%
  filter(time >= simtimeI - 100 & tot_abund > 0) %>%
  group_by(resident, disease, iteration, iterationF) %>%
  summarise(ann_rel_abund = mean(ann_rel_abund),
            per_rel_abund = mean(per_rel_abund)) %>%
  ungroup()

# initial litter value
ldat <- adat %>%
  filter(time == 1) %>%
  select(resident, disease, iteration, Litter) %>%
  rename("eq_litter" = "Litter")

# combine data
cdat <- full_join(idat, rdat) %>%
  full_join(ldat) %>%
  #mutate(outcome = fct_relevel(outcome, "E. virginicus only"))
  mutate(outcome = fct_relevel(outcome, "E. virginicus only", "coexistence", "priority effects"),
         disease = fct_relevel(disease, "without disease"),
         Resident = ifelse(resident == "E. virginicus", "Resident: E. virginicus", "Resident: M. vimineum") %>%
           fct_relevel("Resident: M. vimineum"))

# check that relative abundance makes sense
filter(cdat, outcome == "E. virginicus only" & ann_rel_abund > 0.01) %>%
  select(resident, disease, ann_rel_abund)
filter(cdat, outcome == "E. virginicus only" & ann_rel_abund > 0.01 & resident == "E. virginicus") %>%
  select(resident, disease, ann_rel_abund)
# all Microstegium residents
filter(cdat, outcome == "M. vimineum only" & per_rel_abund > 0.01) %>%
  select(resident, disease, per_rel_abund) %>%
  data.frame()
# all Elymus residents

filter(tdat, iteration == 205 & disease == "with disease" & resident == "M. vimineum" & species %in% c("Perennial_seedling", "Perennial_adult", "Annual", "Litter")) %>%
  ggplot(aes(time, density, color = species)) +
  geom_line()
# slow Microstegium declines

filter(tdat, iteration == 173 & disease == "without disease" & resident == "E. virginicus" & species %in% c("Perennial_seedling", "Perennial_adult", "Annual", "Litter")) %>%
  ggplot(aes(time, density, color = species)) +
  geom_line()
# slow Elymus decline


#### figures ####

pdf("output/annual_perennial_kortessis_outcomes_rel_abund_2018_2019_density_exp.pdf", width = 4.5, height = 3)
cdat %>%
  filter(resident == "E. virginicus") %>%
  group_by(disease, outcome) %>%
  summarise(prop = n()/samps) %>%
  mutate(prop.round = ifelse(prop < 0.01, round(prop, 3), round(prop, 2))) %>%
  ggplot(aes(outcome, prop, fill = disease)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = prop.round), position = position_dodge(0.9), vjust = -0.25, size = 2.5) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  ylab("Probability of outcome") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.position = c(0.5, 0.85),
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        strip.text = element_text(size = 10),
        strip.background = element_blank())
dev.off()

# pdf("output/annual_perennial_kortessis_coexistence_rel_abund_2018_2019_density_exp.pdf", width = 5, height = 3)
# cdat %>% 
#   filter(outcome == "coexistence") %>%
#   ggplot(aes(disease, ann_rel_abund, fill = disease)) +
#   stat_halfeye(point_interval = mean_hdci, .width = c(0.66, 0.95), point_size = 3, shape = 21, point_color = "white", slab_alpha = 0.7, aes(fill = disease)) +
#   scale_fill_viridis_d(direction = -1, end = 0.6) +
#   facet_wrap(~ Resident) +
#   ylab(expression(paste("Relative abundance (", italic("M. vimineum"), "/", italic("E. virginicus"), ")", sep = ""))) +
#   theme_bw() +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 9, color = "black"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 10),
#         legend.position = "none",
#         strip.text = element_text(size = 10),
#         strip.background = element_blank())
# dev.off()

cdat %>% 
  filter(outcome == "priority effects") %>%
  ggplot(aes(disease, ann_rel_abund, color = disease)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  facet_wrap(~ Resident)

cdat %>% 
  filter(outcome == "E. virginicus only") %>%
  ggplot(aes(disease, ann_rel_abund, color = disease)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  facet_wrap(~ Resident)

cdat %>% 
  filter(outcome == "M. vimineum only") %>%
  ggplot(aes(disease, ann_rel_abund, color = disease)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  facet_wrap(~ Resident)

# pdf("output/annual_perennial_kortessis_elymus_wins_example.pdf", width = 3, height = 3)
# tdat %>%
#   filter(iteration == 145 & disease == "without disease" & resident == "M. vimineum" & species %in% c("Annual", "Perennial")) %>%
#   ggplot(aes(time, density, color = species)) +
#   geom_line(size = 1.5) +
#   xlab("Time") +
#   ylab("Density") +
#   scale_color_viridis_d(begin = 0.3, end = 0.8, name = "Species", labels = c("M. vimineum", "E. virginicus")) +
#   theme_bw() +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.text = element_text(size = 9, face = "italic"),
#         legend.title = element_text(size = 10),
#         axis.text = element_text(size = 10, color = "black"),
#         axis.title = element_text(size = 11),
#         legend.position = c(0.73, 0.8))
# dev.off()

# pdf("output/annual_perennial_kortessis_coexistence_example.pdf", width = 3, height = 3)
# tdat %>%
#   filter(iteration == 197 & disease == "with disease" & resident == "E. virginicus" & species %in% c("Annual", "Perennial")) %>%
#   ggplot(aes(time, density, color = species)) +
#   geom_line(size = 1.5) +
#   xlab("Time") +
#   ylab("Density") +
#   scale_color_viridis_d(begin = 0.3, end = 0.8, name = "Species", labels = c("M. vimineum", "E. virginicus")) +
#   theme_bw() +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.text = element_text(size = 9, face = "italic"),
#         legend.title = element_text(size = 10),
#         axis.text = element_text(size = 10, color = "black"),
#         axis.title = element_text(size = 11),
#         legend.position = c(0.5, 0.2))
# dev.off()

# pdf("output/annual_perennial_kortessis_microstegium_wins_example.pdf", width = 3, height = 3)
# tdat %>%
#   filter(iteration == 102 & disease == "without disease" & resident == "E. virginicus" & species %in% c("Annual", "Perennial")) %>%
#   ggplot(aes(time, density, color = species)) +
#   geom_line(size = 1.5) +
#   xlab("Time") +
#   ylab("Density") +
#   scale_color_viridis_d(begin = 0.3, end = 0.8, name = "Species", labels = c("M. vimineum", "E. virginicus")) +
#   theme_bw() +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.text = element_text(size = 9, face = "italic"),
#         legend.title = element_text(size = 10),
#         axis.text = element_text(size = 10, color = "black"),
#         axis.title = element_text(size = 11),
#         legend.position = c(0.3, 0.8))
# dev.off()

pdf("output/annual_perennial_kortessis_rel_abund_2018_2019_density_exp.pdf", width = 4.5, height = 3)
cdat %>%
  filter(resident == "E. virginicus") %>%
  ggplot(aes(ann_rel_abund, fill = disease)) +
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
        legend.position = c(0.5, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))
dev.off()


#### output ####
write_csv(cdat, "output/annual_perennial_kortessis_relative_abundance_constants_2018_2019_density_exp.csv")
write_csv(tdat, "output/annual_perennial_kortessis_relative_abundance_time_series_2018_2019_density_exp.csv")
write_csv(pdat, "output/annual_perennial_kortessis_relative_abundance_parameters_2018_2019_density_exp.csv")
