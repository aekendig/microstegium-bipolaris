##### info ####

# file: annual_perennial_equilibrium_litter_2019_density_exp
# author: Amy Kendig
# date last edited: 10/28/20
# goal: single-species equilibrium litter across decomposition gradient


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidybayes)

# number of iterations
n_samps <- 300

# call model script
source("code/annual_perennial_kortessis_model_2019_density_exp.R")


#### conditions ####

# simulation time
simtime <- 2

# initial conditions
A0 <- 0 # initial annual population size
S0 <- 0 # initial perennial seedling population size
P0 <- 0 # initial perennial adult population size
L0 <- 0 # initial annual litter amount

# samples from parameter distributions
samps <- n_samps


#### simulations ####

# initiate data frame
lit_df <- tibble(d_in = seq(0.05, 1, by = 0.05)) %>%
  expand_grid(tibble(iter = 1:n_samps)) %>%
  mutate(A0 = A0,
         S0 = S0,
         P0 = P0,
         L0 = L0,
         simtime = simtime,
         mv_ndis = NA,
         ev_ndis = NA,
         mv_ydis = NA,
         ev_ydis = NA)

# add litter estimates
lit_df2 <- lit_df %>%
  rowwise() %>%
  mutate(mv_ndis = sim_fun(A0, S0, P0, L0, simtime, 0, iter, d_in)[[2]] %>% filter(parameter == "L.A") %>% pull(value) %>% unique(),
         ev_ndis = sim_fun(A0, S0, P0, L0, simtime, 0, iter, d_in)[[2]] %>% filter(parameter == "L.P") %>% pull(value) %>% unique(),
         mv_ydis = sim_fun(A0, S0, P0, L0, simtime, 1, iter, d_in)[[2]] %>% filter(parameter == "L.A") %>% pull(value) %>% unique(),
         ev_ydis = sim_fun(A0, S0, P0, L0, simtime, 1, iter, d_in)[[2]] %>% filter(parameter == "L.P") %>% pull(value) %>% unique())


#### process output ####

# make long
lit_df3 <- lit_df2 %>%
  pivot_longer(cols = c(mv_ndis, ev_ndis, mv_ydis, ev_ydis), names_to = "sp_dis", values_to = "litter") %>%
  mutate(species = substring(sp_dis, 1, 2) %>%
           recode("ev" = "E. virginicus", 
                  "mv" = "M. vimineum"),
         Disease = substring(sp_dis, 4, 7) %>%
           recode("ndis" = "without disease",
                  "ydis" = "with disease"))

# summarise for plot
lit_sum <- lit_df3 %>%
  group_by(species, Disease, d_in) %>%
  mean_qi(litter, na.rm = T) %>%
  ungroup() %>%
  mutate(species = as.factor(species) %>% fct_relevel("M. vimineum"),
         Disease = as.factor(Disease) %>% fct_relevel("without disease"))


#### visualize ####

# figure
pdf("output/ annual_perennial_equilibrium_litter_figure_2019_density_exp.pdf", width = 4.5, height = 3)
ggplot(lit_sum, aes(x = d_in, y = litter, color = Disease)) +
  geom_line(size = 1.5) +
  geom_line(aes(y = .lower), linetype = "dotted") +
  geom_line(aes(y = .upper), linetype = "dotted") +
  facet_wrap(~ species) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  ylab("Equilibrium litter production") +
  xlab(expression(paste("Decomposition (", year^-1, ")", sep = ""))) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        strip.text = element_text(size = 10, face = "italic"),
        strip.background = element_blank(),
        legend.position = c(0.82, 0.83))
dev.off()


#### output ####
write_csv(lit_df2, "output/annual_perennial_equilibrium_litter_2019_density_exp.csv")
