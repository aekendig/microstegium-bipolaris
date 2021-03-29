##### info ####

# file: annual_perennial_relative_abundance_2019_density_exp
# author: Amy Kendig
# date last edited: 3/29/20
# goal: relative abundance of annual to perennial


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(ggforce)

# call model script
source("code/annual_perennial_simulation_model_2019_density_exp.R")

# samples from parameter distributions
samps <- sample(1:15000, 5000)

# simulation time
simtimeR <- 600
simtimeI <- 200

# function to simulate invasions
inv_fun <- function(post_draw, g_S_DE, b_A_DE, b_S_DE, alpha_SA_DE){
  
  # establish residents
  mod_est <- sim_fun(A0 = A0E, S0 = S0E, P0 = P0E, L0 = L0E, simtime = simtimeR, iter = post_draw, 
                     g.S.DE = g_S_DE, b.A.DE = b_A_DE, b.S.DE = b_S_DE, alpha.SA.DE = alpha_SA_DE)
  
  # extract final abundances
  S0Ei <- mod_est[[1]] %>% filter(time == simtimeR & species == "Perennial seedling") %>% pull(N)
  P0Ei <- mod_est[[1]] %>% filter(time == simtimeR & species == "Perennial adult") %>% pull(N)
  L0Ei <- mod_est[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  
  # invasion
  mod_inv <- sim_fun(A0 = 10, S0 = S0Ei, P0 = P0Ei, L0 = L0Ei, simtime = simtimeI, iter = post_draw, 
                     g.S.DE = g_S_DE, b.A.DE = b_A_DE, b.S.DE = b_S_DE, alpha.SA.DE = alpha_SA_DE)
  
  # output
  return(mod_inv)
}

# initial conditions
A0E <- 0 # initial annual population size
S0E <- 0 # initial perennial seedling population size
P0E <- 1 # initial perennial adult population size
L0E <- 0 # initial annual litter amount

# initiate lists
abund_no <- list()
params_no <- list()
abund_bo <- list()
params_bo <- list()
abund_boA <- list()
params_boA <- list()
abund_pr <- list()
params_pr <- list()
abund_prA <- list()
params_prA <- list()
abund_an <- list()
params_an <- list()
abund_bo2 <- list()
params_bo2 <- list()
abund_boA2 <- list()
params_boA2 <- list()
abund_pr2 <- list()
params_pr2 <- list()
abund_prA2 <- list()
params_prA2 <- list()
abund_an2 <- list()
params_an2 <- list()


#### simulations ####

for(i in 1:length(samps)){
  
  print(i)
  
  # no disease impacts
  mod_no <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = 1, b_S_DE = 1, alpha_SA_DE = 1)
  abund_no[[i]] <- mod_no[[1]]
  params_no[[i]] <- mod_no[[2]]
  
  # disease impacts on both
  mod_bo <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = mvBioDE, b_S_DE = evBioDE, alpha_SA_DE = alphaSADE)
  abund_bo[[i]] <- mod_bo[[1]]
  params_bo[[i]] <- mod_bo[[2]]
  
  # disease impacts on both, includes positive effects
  mod_boA <- inv_fun(post_draw = samps[i], g_S_DE = evGermDE, b_A_DE = mvBioDE, b_S_DE = evBioDE, alpha_SA_DE = alphaSADE)
  abund_boA[[i]] <- mod_boA[[1]]
  params_boA[[i]] <- mod_boA[[2]]
  
  # disease impacts on perennial only
  mod_pr <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = 1, b_S_DE = evBioDE, alpha_SA_DE = alphaSADE)
  abund_pr[[i]] <- mod_pr[[1]]
  params_pr[[i]] <- mod_pr[[2]]
  
  # disease impacts on perennial only, includes positive effects
  mod_prA <- inv_fun(post_draw = samps[i], g_S_DE = evGermDE, b_A_DE = 1, b_S_DE = evBioDE, alpha_SA_DE = alphaSADE)
  abund_prA[[i]] <- mod_prA[[1]]
  params_prA[[i]] <- mod_prA[[2]]
  
  # disease impacts on annual only
  mod_an <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = mvBioDE, b_S_DE = 1, alpha_SA_DE = 1)
  abund_an[[i]] <- mod_an[[1]]
  params_an[[i]] <- mod_an[[2]]
  
  # disease impacts on both x 2
  mod_bo2 <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = mvBioDE/2, b_S_DE = evBioDE/2, alpha_SA_DE = alphaSADE*2)
  abund_bo2[[i]] <- mod_bo2[[1]]
  params_bo2[[i]] <- mod_bo2[[2]]
  
  # disease impacts on both x 2, includes positive effects
  mod_boA2 <- inv_fun(post_draw = samps[i], g_S_DE = evGermDE*2, b_A_DE = mvBioDE/2, b_S_DE = evBioDE/2, alpha_SA_DE = alphaSADE*2)
  abund_boA2[[i]] <- mod_boA2[[1]]
  params_boA2[[i]] <- mod_boA2[[2]]
  
  # disease impacts on perennial only x 2
  mod_pr2 <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = 1, b_S_DE = evBioDE/2, alpha_SA_DE = alphaSADE*2)
  abund_pr2[[i]] <- mod_pr2[[1]]
  params_pr2[[i]] <- mod_pr2[[2]]
  
  # disease impacts on perennial only x 2, includes positive effects
  mod_prA2 <- inv_fun(post_draw = samps[i], g_S_DE = evGermDE*2, b_A_DE = 1, b_S_DE = evBioDE/2, alpha_SA_DE = alphaSADE*2)
  abund_prA2[[i]] <- mod_prA2[[1]]
  params_prA2[[i]] <- mod_prA2[[2]]
  
  # disease impacts on annual only x 2
  mod_an2 <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = mvBioDE/2, b_S_DE = 1, alpha_SA_DE = 1)
  abund_an2[[i]] <- mod_an2[[1]]
  params_an2[[i]] <- mod_an2[[2]]
  
}


# convert to dataframes
abund_dat_no <- do.call(rbind, abund_no)
param_dat_no <- do.call(rbind, params_no)
abund_dat_bo <- do.call(rbind, abund_bo)
param_dat_bo <- do.call(rbind, params_bo)
abund_dat_pr <- do.call(rbind, abund_pr)
param_dat_pr <- do.call(rbind, params_pr)
abund_dat_prA <- do.call(rbind, abund_prA)
param_dat_prA <- do.call(rbind, params_prA)
abund_dat_an <- do.call(rbind, abund_an)
param_dat_an <- do.call(rbind, params_an)
abund_dat_bo2 <- do.call(rbind, abund_bo2)
param_dat_bo2 <- do.call(rbind, params_bo2)
abund_dat_pr2 <- do.call(rbind, abund_pr2)
param_dat_pr2 <- do.call(rbind, params_pr2)
abund_dat_prA2 <- do.call(rbind, abund_prA2)
param_dat_prA2 <- do.call(rbind, params_prA2)
abund_dat_an2 <- do.call(rbind, abund_an2)
param_dat_an2 <- do.call(rbind, params_an2)


#### process data ####

# combine scenarios
rel_abund_dat <- abund_dat_no %>%
  mutate(disease = "observed",
         impacts = "none") %>%
  full_join(abund_dat_no %>%
              mutate(disease = "observed x 2",
                     impacts = "none")) %>%
  full_join(abund_dat_bo %>%
              mutate(disease = "observed",
                     impacts = "both (negative only)")) %>%
  full_join(abund_dat_boA %>%
              mutate(disease = "observed",
                     impacts = "both")) %>%
  full_join(abund_dat_pr %>%
              mutate(disease = "observed",
                     impacts = "perennial (negative only)")) %>%
  full_join(abund_dat_prA %>%
              mutate(disease = "observed",
                     impacts = "perennial")) %>%
  full_join(abund_dat_an %>%
              mutate(disease = "observed",
                     impacts = "annual")) %>%
  full_join(abund_dat_bo2 %>%
              mutate(disease = "observed x 2",
                     impacts = "both (negative only)")) %>%
  full_join(abund_dat_boA2 %>%
              mutate(disease = "observed x 2",
                     impacts = "both")) %>%
  full_join(abund_dat_pr2 %>%
              mutate(disease = "observed x 2",
                     impacts = "perennial (negative only)")) %>%
  full_join(abund_dat_prA2 %>%
              mutate(disease = "observed x 2",
                     impacts = "perennial")) %>%
  full_join(abund_dat_an2 %>%
              mutate(disease = "observed x 2",
                     impacts = "annual")) %>%
  pivot_wider(names_from = species, values_from = N) %>%
  rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
  mutate(Perennial = Perennial_seedling + Perennial_adult,
         tot_abund = Annual + Perennial,
         ann_rel_abund = Annual / tot_abund,
         per_rel_abund = Perennial / tot_abund,
         iterationF = as.factor(iteration))

# relative abundance summary
rel_sum_dat <- rel_abund_dat %>%
  filter(time >= simtimeI - 50) %>%
  group_by(disease, impacts, iteration, iterationF) %>%
  summarise(ann_rel_abund = mean(ann_rel_abund),
            per_rel_abund = mean(per_rel_abund)) %>%
  ungroup() %>%
  mutate(outcome = case_when(ann_rel_abund == 1 ~ "annual only",
                             per_rel_abund == 1 ~ "perennial only",
                             TRUE ~ "mixed"),
         impacts = fct_relevel(impacts, "none", "annual", "perennial (negative only)", "both (negative only)", "perennial", "both"))

# remove positive effects
rel_sum_dat2 <- rel_sum_dat %>%
  filter(!(impacts %in% c("perennial", "both"))) %>%
  mutate(impacts = recode(impacts, "perennial (negative only)" = "perennial",
                          "both (negative only)" = "both"))

# relative abundance hdi
rel_hdi_dat <- rel_sum_dat %>%
  group_by(disease, impacts) %>%
  mean_hdci(ann_rel_abund, .width = 0.95)

# absolute abundance
abs_abund_dat <- rel_abund_dat %>%
  select(disease, impacts, time, iteration, iterationF, Perennial, Annual, Litter) %>%
  pivot_longer(cols = c("Perennial", "Annual", "Litter"), names_to = "species", values_to = "density")


#### figures ####

# pdf("output/annual_perennial_relative_abundance_time_series_2019_density_exp.pdf")
# ggplot(rel_abund_dat, aes(time, ann_rel_abund, color = iterationF)) +
#   geom_line()  +
#   facet_grid(impacts ~ disease) +
#   theme_bw() +
#   theme(legend.position = "none")
# dev.off()

# pdf("output/annual_perennial_relative_abundance_distribution_2019_density_exp.pdf", width = 4.5, height = 3)

ggplot(rel_sum_dat, aes(x = outcome, fill = impacts)) +
  geom_bar(position = "dodge") +
  ylab("Number of draws") +
  facet_wrap(~ disease, nrow = 1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.position = c(0.35, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))

ggplot(rel_sum_dat2, aes(x = outcome, fill = impacts)) +
  geom_bar(position = "dodge") +
  ylab("Number of draws") +
  facet_wrap(~ disease, nrow = 1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.position = c(0.35, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))

# dev.off()

ggplot(rel_hdi_dat, aes(disease, ann_rel_abund, color = impacts)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, position = position_dodge(0.1)) +
  geom_point(size = 2, position = position_dodge(0.1))
  
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
