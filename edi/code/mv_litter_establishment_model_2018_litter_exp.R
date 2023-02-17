##### outputs ####

# mv_litter_establishment_data_2018_litter_exp.csv
# mv_litter_establishment_model_2018_litter_exp.rda


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
estL1Dat <- read_csv("data/both_germination_disease_jul_2018_litter_exp.csv")
plotsL <- read_csv("data/plot_treatments_2018_litter_exp.csv")


#### edit data ####

# edit variables
# remove unnecessary variables
plotsL2 <- plotsL %>%
  mutate(sterilized = case_when(litter == "live" ~ 0,
                                TRUE ~ 1),
         sterilizedF = ifelse(sterilized == 0, "live", "sterilized") %>%
           fct_relevel("sterilized"),
         litter.g.m2 = litter_weight.g,
         litter.g.cm2 = litter.g.m2/10000) %>%
  select(-c(flag_color, justification, litter_weight.g))

# select plots with seeds added only
mvEstL1Dat <- estL1Dat %>%
  filter(seeds_added == "yes") %>%
  left_join(plotsL2) %>%
  select(-c(date, seeds_added, ev_germ, ev_infec)) %>%
  rename(mv_germ_planted_bg = mv_germ,
         mv_germ_bg = mv_germ_ev) %>%
  mutate(mv_germ_planted = mv_germ_planted_bg - mv_germ_bg,
         mv_planted = 200 + mv_germ_bg,
         mv_germ_planted_cor = case_when(mv_germ_planted < 0 ~ 0,
                                         mv_germ_planted > 200 ~ 200,
                                         TRUE ~ mv_germ_planted),
         mv_planted_cor = case_when(mv_planted < mv_germ_planted_bg ~ mv_germ_planted_bg,
                                    TRUE ~ mv_planted),
         prop_germ_num_cor = mv_germ_planted_cor/200,
         prop_germ_den_cor = mv_germ_planted_bg/mv_planted_cor)


#### model ####

# initial fit
mvEstL1Mod1 <- brm(mv_germ_planted_bg | trials(mv_planted_cor) ~ litter.g.cm2 * sterilized + (1|site),
                   data = mvEstL1Dat, family = binomial,
                   prior = c(prior(normal(0, 10), class = Intercept),
                             prior(normal(0, 10), class = b)),
                   iter = 6000, warmup = 1000, chains = 1)
# 28 divergent transitions
summary(mvEstL1Mod1)

# increase chains
mvEstL1Mod2 <- update(mvEstL1Mod1, chains = 3,
                      control = list(adapt_delta = 0.9999,
                                     max_treedepth = 15))
summary(mvEstL1Mod2)
plot(mvEstL1Mod2)
pp_check(mvEstL1Mod2, ndraws = 50)


#### output ####

save(mvEstL1Mod2, file = "output/mv_litter_establishment_model_2018_litter_exp.rda")
write_csv(mvEstL1Dat, "intermediate-data/mv_litter_establishment_data_2018_litter_exp.csv")
