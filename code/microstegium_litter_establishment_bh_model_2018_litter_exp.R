##### info ####

# file: microstegium_litter_establishment_bh_model_2018_litter_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: estimate the effect of litter on Microstegium establishment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
estL1Dat <- read_csv("./data/both_germination_disease_jul_2018_litter_exp.csv")
plotsL <- read_csv("./data/plot_treatments_2018_litter_exp.csv")


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

# add a row for no litter with the live treatment (will repeat the same data twice)
plotsL3 <- plotsL2 %>%
  full_join(tibble(plot = 13, litter = "none", litter_density = "none", seeds_added = "yes", sterilized = 0, sterilizedF = "live", litter.g.m2 = 0, litter.g.cm2 = 0))

# select plots with seeds added only
mvEstL1BhDat <- estL1Dat %>%
  filter(seeds_added == "yes") %>%
  left_join(plotsL3) %>%
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

# check repeated plots with no litter
filter(mvEstL1BhDat, litter.g.m2 == 0) %>% data.frame()


#### model ####

# initial fit
mvEstL1BhMod1 <- brm(bf(prop_germ_den_cor ~ maxEst / (1 + betaL * litter.g.cm2),
                        maxEst ~ 1 + (1|site),
                        betaL ~ 0 + sterilizedF,
                        nl = T),
                   data = mvEstL1BhDat, family = gaussian,
                   prior = c(prior(normal(1, 1), nlpar = "maxEst", class = "b", ub = 1),
                             prior(exponential(0.5), nlpar = "betaL", lb = 0),
                             prior(cauchy(0, 1), class = "sigma"),
                             prior(cauchy(0, 1), nlpar = "maxEst", class = "sd")),
                   iter = 6000, warmup = 1000, chains = 1)
# 6 divergent transitions
summary(mvEstL1BhMod1)

# increase chains
mvEstL1BhMod2 <- update(mvEstL1BhMod1, chains = 3,
                        control = list(adapt_delta = 0.9999))
summary(mvEstL1BhMod2)
plot(mvEstL1BhMod2)
pp_check(mvEstL1BhMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(sterilizedF = c("live", "sterilized") %>%
                   fct_relevel("sterilized")) %>%
  expand_grid(tibble(litter.g.cm2 = seq(0, 200/10000, length.out = 100))) %>%
  mutate(Est = fitted(mvEstL1BhMod2, newdata = ., type = "response", re_formula = NA)[, "Estimate"],
         Est_lower = fitted(mvEstL1BhMod2, newdata = ., type = "response", re_formula = NA)[, "Q2.5"],
         Est_upper = fitted(mvEstL1BhMod2, newdata = ., type = "response", re_formula = NA)[, "Q97.5"],
         litter.g.m2 = litter.g.cm2 * 10000)

# figure
pdf("output/microstegium_litter_establishment_bh_figure_2018_litter_exp.pdf", width = 3.5, height = 3.5)
ggplot(mvEstL1BhDat, aes(x = litter.g.m2, y = prop_germ_den_cor, fill = sterilizedF, color = sterilizedF)) +
  geom_ribbon(data = fitDat, aes(y = Est, ymin = Est_lower, ymax = Est_upper), alpha = 0.6, color = NA) +
  geom_line(data = fitDat, aes(y = Est)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(5)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(5)) +
  scale_color_viridis_d(name = "Litter treatment", direction = -1, end = 0.6) +
  scale_fill_viridis_d(name = "Litter treatment", direction = -1, end = 0.6) +
  ggtitle(expression(italic("M. vimineum"))) +
  ylab("Proportion established") +
  xlab(expression(paste("Litter (g ", m^-2, ")", sep = ""))) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = c(0.5, 0.18))
dev.off()


#### output ####

save(mvEstL1BhMod2, file = "output/microstegium_litter_establishment_bh_model_2018_litter_exp.rda")
