##### info ####

# file: microstegium_litter_establishment_model_2018_litter_exp
# author: Amy Kendig
# date last edited: 10/24/20
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

# select plots with seeds added only
mvEstL1Dat <- estL1Dat %>%
  filter(seeds_added == "yes" & litter != "none") %>%
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
                             prior(normal(0, 10), class = b),
                             prior(cauchy(0, 1), class = sd)),
                   iter = 6000, warmup = 1000, chains = 1)
# 10 divergent transitions
summary(mvEstL1Mod1)

# increase chains
mvEstL1Mod2 <- update(mvEstL1Mod1, chains = 3,
                      control = list(adapt_delta = 0.9999,
                                     max_treedepth = 15))
summary(mvEstL1Mod2)
plot(mvEstL1Mod2)
pp_check(mvEstL1Mod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(litter.g.cm2 = rep(seq(50/10000, 200/10000, length.out = 100), 2),
                 sterilized = rep(c(0, 1), each = 100)) %>%
  mutate(mv_planted_cor = 747) %>%
  mutate(Est = fitted(mvEstL1Mod2, newdata = ., type = "response", re_formula = NA)[, "Estimate"]/mv_planted_cor,
         Est_lower = fitted(mvEstL1Mod2, newdata = ., type = "response", re_formula = NA)[, "Q2.5"]/mv_planted_cor,
         Est_upper = fitted(mvEstL1Mod2, newdata = ., type = "response", re_formula = NA)[, "Q97.5"]/mv_planted_cor,
         litter.g.m2 = litter.g.cm2 * 10000,
         sterilizedF = ifelse(sterilized == 0, "live", "sterilized") %>%
           fct_relevel("sterilized"))

# figure
pdf("output/microstegium_litter_establishment_figure_2018_litter_exp.pdf", width = 3.5, height = 3.5)
ggplot(mvEstL1Dat, aes(x = litter.g.m2, y = prop_germ_den_cor, fill = sterilizedF, color = sterilizedF)) +
  geom_ribbon(data = fitDat, aes(y = Est, ymin = Est_lower, ymax = Est_upper), alpha = 0.6, color = NA) +
  geom_line(data = fitDat, aes(y = Est)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(5)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(5)) +
  scale_color_viridis_d(name = "Litter treatment", direction = -1, end = 0.6) +
  scale_fill_viridis_d(name = "Litter treatment", direction = -1, end = 0.6) +
  ggtitle(expression(italic("M. vimineum"))) +
  ylab("Proportion established") +
  xlab(expression(paste("Litter (g ", m^-2, ")", sep = ""))) +
  coord_cartesian(ylim = c(0.2, 1)) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = c(0.8, 0.83))
dev.off()


#### output ####

save(mvEstL1Mod2, file = "output/microstegium_litter_establishment_model_2018_litter_exp.rda")
