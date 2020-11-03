##### info ####

# file: annual_perennial_kortessis_initial_conditions_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: annual and perennial dominance across a range of initial conditions

#### set-up ####

# clear all existing data
rm(list=ls())

# number of iterations
n_samps <- 200

# call model script
source("code/annual_perennial_kortessis_model_2019_density_exp.R")


#### conditions ####

# simulation time
simtime <- 1000

# samples from parameter distributions
samps <- n_samps


#### simulations ####

# initiate dataframe
init_df <- tibble(A0 = 10^(1:5)) %>%
  expand_grid(tibble(S0 = 10^(1:5)/2,
                     P0 = 10^(1:5)/2)) %>%
  expand_grid(tibble(disease = c(0, 1))) %>%
  expand_grid(tibble(iter = 1:samps)) %>%
  mutate(L0 = 0)

# save relative abundance
init_df2 <- init_df %>%
  rowwise() %>%
  mutate(rel_abund = sim_fun(A0, S0, P0, L0, simtime, disease, iter, d)[[1]] %>%
           filter(time >= simtime - 200) %>%
           pivot_wider(names_from = species, values_from = N) %>%
           rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
           mutate(rel_abund = Annual / (Annual + Perennial_seedling + Perennial_adult)) %>%
           pull(rel_abund) %>%
           mean())


#### process data ####

# add columns
init_df3 <- init_df2 %>%
  ungroup() %>%
  mutate(Disease = ifelse(disease == 0, "without disease", "with disease") %>% 
           fct_relevel("without disease"),
         P_init = P0 + S0,
         A_init = A0) %>%
  group_by(A_init, P_init, Disease) %>%
  summarise(rel_abund = mean(rel_abund, na.rm = T))

# figure
pdf("output/annual_perennial_kortessis_initial_conditions_2019_density_exp.pdf", width = 6, height = 3)
ggplot(init_df3, aes(x = log10(A_init), y = log10(P_init), fill = rel_abund)) +
  geom_point(shape = 22, size = 20) +
  xlab(expression(paste("Initial ", italic("M. vimineum"), " density (", log[10], "(plants ", m^-2, "))", sep = ""))) +
  ylab(expression(paste("Initial ", italic("E. virginicus"), " density (", log[10], "(plants ", m^-2, "))", sep = ""))) +
  facet_wrap(~ Disease) +
  scale_fill_viridis_c(name = "M. vimineum\nrelative\nabundance", direction = -1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 10),
        strip.background = element_blank())
dev.off()