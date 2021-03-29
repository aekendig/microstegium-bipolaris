##### info ####

# file: annual_perennial_disease_effects_2019_density_exp
# author: Amy Kendig
# date last edited: 3/29/20
# goal: hypothetical disease effects on relative abundance


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(ggforce)

# call model script
source("code/annual_perennial_simulation_model_2019_density_exp.R")

# samples from parameter distributions
n_samps <- 2000
samps <- sample(1:15000, n_samps)

# disease effects
DE <- c(1, 5, 10, 100, 1000)

# simulation time
simtimeR <- 500
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
P0E <- 10 # initial perennial adult population size
L0E <- 0 # initial annual litter amount

# initiate lists
abund_no <- list()
abund_bo <- list()
abund_pr <- list()
abund_an <- list()


#### simulations ####

for(i in 1:length(samps)){
  
  print(i)
  
  # no disease impacts
  mod_no <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = 1, b_S_DE = 1, alpha_SA_DE = 1)
  abund_no[[i]] <- mod_no[[1]] %>% filter(time > (simtimeI - 50))
  
  # initiate lists
  abund_bo_j <- list()
  abund_pr_j <- list()
  abund_an_j <- list()
  
  for(j in 1:length(DE)){
    
    # disease impacts on both
    mod_bo <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = mvBioDE/DE[j], b_S_DE = evBioDE/DE[j], alpha_SA_DE = alphaSADE*DE[j])
    abund_bo_j[[j]] <- mod_bo[[1]] %>% filter(time > (simtimeI - 50)) %>% mutate(DE = DE[j])
    
    # disease impacts on perennial only
    mod_pr <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = 1, b_S_DE = evBioDE/DE[j], alpha_SA_DE = alphaSADE*DE[j])
    abund_pr_j[[j]] <- mod_pr[[1]] %>% filter(time > (simtimeI - 50)) %>% mutate(DE = DE[j])
    
    # disease impacts on annual only
    mod_an <- inv_fun(post_draw = samps[i], g_S_DE = 1, b_A_DE = mvBioDE/DE[j], b_S_DE = 1, alpha_SA_DE = 1)
    abund_an_j[[j]] <- mod_an[[1]] %>% filter(time > (simtimeI - 50)) %>% mutate(DE = DE[j])
    
  }
  
  abund_bo[[i]] <- do.call(rbind, abund_bo_j)
  abund_pr[[i]] <- do.call(rbind, abund_pr_j)
  abund_an[[i]] <- do.call(rbind, abund_an_j)
  
}  

# convert to dataframes
abund_dat_no <- do.call(rbind, abund_no) %>% mutate(DE = 0)
abund_dat_bo <- do.call(rbind, abund_bo)
abund_dat_pr <- do.call(rbind, abund_pr)
abund_dat_an <- do.call(rbind, abund_an)


#### process data ####

# combine scenarios
rel_abund_dat <- abund_dat_no %>%
  mutate(impacts = "none") %>%
  full_join(abund_dat_bo %>%
              mutate(impacts = "both")) %>%
  full_join(abund_dat_pr %>%
              mutate(impacts = "perennial")) %>%
  full_join(abund_dat_an %>%
              mutate(impacts = "annual")) %>%
  pivot_wider(names_from = species, values_from = N) %>%
  rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
  mutate(Perennial = Perennial_seedling + Perennial_adult,
         tot_abund = Annual + Perennial,
         ann_rel_abund = Annual / tot_abund,
         per_rel_abund = Perennial / tot_abund,
         iterationF = as.factor(iteration))

# relative abundance summary
rel_sum_dat <- rel_abund_dat %>%
  group_by(impacts, DE, iteration, iterationF) %>%
  summarise(ann_rel_abund = mean(ann_rel_abund),
            per_rel_abund = mean(per_rel_abund)) %>%
  ungroup() %>%
  mutate(outcome = case_when(ann_rel_abund == 1 ~ "annual only",
                             per_rel_abund == 1 ~ "perennial only",
                             is.na(ann_rel_abund) & is.na(per_rel_abund) ~ "neither",
                             ann_rel_abund < 1 & per_rel_abund < 1 ~ "mixed") %>%
           fct_relevel("perennial only"),
         impacts = fct_relevel(impacts, "none", "annual", "perennial", "both"),
         community = case_when(outcome == "annual only" ~ "invasive only",
                               outcome == "perennial only" ~ "native only",
                               TRUE ~ as.character(outcome)) %>%
           fct_relevel("native only", "invasive only"))



#### figures ####


ggplot(rel_sum_dat, aes(x = outcome, fill = impacts)) +
  geom_bar(position = "dodge") +
  ylab("Number of draws") +
  facet_wrap(~ DE, nrow = 1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.position = c(0.1, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))

pdf("output/annual_perennial_disease_effects_2019_density_exp.pdf", width = 2.75, height = 2.75)
rel_sum_dat %>%
  filter(impacts %in% c("none", "both")) %>%
  group_by(DE, community) %>%
  summarise(draws = n()/n_samps) %>%
  ggplot(aes(DE, draws, color = community)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_log10() +
  ylab("Probability of community type") +
  xlab("Disease intensity (x observed)") +
  scale_color_manual(values = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4"),
                     name = "Community type") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = c(0.65, 0.55),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
dev.off()

pdf("output/example_posterior_distribution_germination.pdf", width = 2.75, height = 2.75)
ggplot(g.A.samps, aes(x = g.A)) +
  geom_density(size = 2, color = "darkgoldenrod2")+
  theme_bw() +
  xlab("Germination") +
  ylab("Posterior density") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = c(0.65, 0.55),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
dev.off()

pdf("output/example_posterior_distribution_establishment.pdf", width = 2.75, height = 2.75)
ggplot(e.A.samps, aes(x = e.A)) +
  geom_density(size = 2, color = "palegreen4")+
  theme_bw() +
  xlab("Establishment") +
  ylab("Posterior density") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = c(0.65, 0.55),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
dev.off()

pdf("output/example_posterior_distribution_litter.pdf", width = 2.75, height = 2.75)
ggplot(beta.A.samps, aes(x = beta.A)) +
  geom_density(size = 2, color = "dodgerblue1")+
  theme_bw() +
  xlab("Litter effects on establishment") +
  ylab("Posterior density") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        legend.position = c(0.65, 0.55),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
dev.off()



#### output ####
write_csv(rel_abund_dat, "output/annual_perennial_disease_effects_2019_density_exp.csv")
