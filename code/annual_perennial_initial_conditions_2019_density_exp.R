##### info ####

# file: annual_perennial_initial_conditions_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: annual and perennial dominance across a range of initial conditions


#### haven't edited script from relative abundance analysis yet ####


#### set-up ####

# clear all existing data
rm(list=ls())

# number of iterations
n_samps <- 200

# call model script
source("code/annual_perennial_simulation_model_2019_density_exp.R")


#### conditions ####

# simulation time
simtime <- 600

# initial conditions
A0 <- 1 # initial annual population size
S0 <- 0 # initial perennial seedling population size
P0 <- 1 # initial perennial adult population size
L0 <- 0 # initial annual litter amount

# samples from parameter distributions
samps <- n_samps


#### no disease simulations ####

# initiate list
abund <- list()
params <- list()

# save parameter values from each sample
for(iter in 1:samps){
  
  mod <- sim_fun(A0, S0, P0, L0, simtime, disease = 0, iter)
  abund[[iter]] <- mod[[1]]
  params[[iter]] <- mod[[2]]
  
}

# convert to dataframes
abund_dat <- do.call(rbind, abund)
param_dat <- do.call(rbind, params)


#### disease simulations ####

# initiate list
abund2 <- list()
params2 <- list()

# save parameter values from each sample
for(iter in 1:samps){
  
  mod2 <- sim_fun(A0, S0, P0, L0, simtime, disease = 1, iter)
  abund2[[iter]] <- mod2[[1]]
  params2[[iter]] <- mod2[[2]]
  
}

# convert to dataframes
abund_dat2 <- do.call(rbind, abund2)
param_dat2 <- do.call(rbind, params2)


#### relative abundance ####

# relative abundance
rel_abund_dat <- abund_dat %>%
  mutate(disease = "without disease") %>%
  full_join(abund_dat2 %>%
              mutate(disease = "with disease")) %>%
  pivot_wider(names_from = species, values_from = N) %>%
  rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
  mutate(Perennial = Perennial_seedling + Perennial_adult,
         tot_abund = Annual + Perennial,
         ann_rel_abund = ifelse(tot_abund > 0, Annual / tot_abund, 0),
         per_rel_abund = ifelse(tot_abund > 0, Perennial / tot_abund, 0),
         iterationF = as.factor(iteration))

# figure
ggplot(rel_abund_dat, aes(time, ann_rel_abund, color = iterationF)) +
  geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~ disease)

# relative abundance summary
rel_sum_dat <- rel_abund_dat %>%
  filter(time >= simtime - 100 & tot_abund > 0) %>%
  group_by(disease, iteration, iterationF) %>%
  summarise(ann_rel_abund = mean(ann_rel_abund),
            per_rel_abund = mean(per_rel_abund))

# figure
pdf("output/annual_perennial_relative_abundance_2019_density_exp.pdf")
ggplot(rel_sum_dat, aes(x = ann_rel_abund, fill = disease)) +
  geom_histogram(binwidth = 0.01, position = position_dodge()) +
  xlab("Annual relative abundance") +
  ylab("Number of iterations") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8))
dev.off()