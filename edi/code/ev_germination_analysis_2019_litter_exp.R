##### info ####

# file: ev_germination_analysis_2018_litter_exp
# author: Amy Kendig
# date last edited: 12/1/21
# goal: evaluate the effects of litter treatments on the germination and disease of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)

# import data
germ_jun <- read_csv("./data/both_germination_disease_jun_2019_litter_exp.csv")
plots <- read_csv("./data/litter_weight_apr_2019_litter_exp.csv")
# add in leaf scans when completed


#### edit data ####

# edit plot data
# put the removed litter into the addition litter
# scale litter weight
# remove unnecessary variables
plots2 <- plots %>%
  spread(key = treatment, value = litter_weight.lb) %>%
  mutate(addition = addition + removal,
         removal = 0) %>%
  gather(key = "treatment", value = "litter_weight.lb", -c(date, site, block)) %>%
  mutate(litter_weight.scaled = (litter_weight.lb - mean(litter_weight.lb)) / sd(litter_weight.lb),
         litter_weight.g = litter_weight.lb * 453.592,
         litter.g.m2 = litter_weight.g / 4, # 4 because the plots are 2m^2
         litter.g.cm2 = litter.g.m2 / 10000,
         treatment = fct_relevel(treatment, "removal", "control"),
         plot = as.factor(paste(site, block, sep = "_"))) %>%
  select(-c(date))

# Ev planting data
plant <- tibble(treatment = c("removal", "control", "addition"),
                ev_tot = c(50, 26, 26)) 

# June germination data
# none of the germinants were infected
germ <- germ_jun %>%
  select(-c(date, flag_color, mv_germ, mv_infec)) %>%
  full_join(plots2) %>%
  full_join(plant) %>%
  mutate(ev_prop_germ = ev_germ / ev_tot,
         treatment = fct_relevel(treatment, "removal", "control"))


#### raw data figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# color palette
col_pal = c("blue", "purple", "red")

# base theme
base_theme <- theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none")

# litter treatments
fig_trt <- ggplot(plots2, aes(x = treatment, y = litter_weight.g)) +
  geom_line(aes(size = treatment, group = plot), color = "black") +
  geom_point(aes(color = treatment, shape = plot), size = 3) +
  base_theme +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_shape_manual(values = c(19, 15, 17, 18), name = "Plot") +
  scale_size_manual(values = c(0.5, 0.5, 0.5), guide = F) +
  xlab("Litter treatment") +
  ylab("Litter weight (g)")

# germination
fig_germ_trt <- ggplot(germ, aes(x = treatment, y = ev_prop_germ)) +
  stat_summary(aes(group = treatment), geom = "errorbar", fun.data = mean_cl_boot, width = 0.2) +
  stat_summary(aes(fill = treatment), geom = "point", fun.y = mean, size = 3, shape = 21) +
  xlab("Litter treatment") +
  base_theme +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  ylab("Proportion seeds\nemerged in June")

fig_germ_lb <- ggplot(germ, aes(x = litter_weight.g, y = ev_prop_germ, color = treatment)) +
  geom_point(size = 2, aes(shape = plot)) +
  xlab("Litter weight (g)") +
  base_theme +
  scale_shape_manual(values = c(19, 15, 17, 18), name = "Plot") +
  scale_color_manual(values = col_pal, name = "Treatment") +
  ylab("Proportion seeds\nemerged in June")

# legend
leg <- get_legend(fig_trt +
                 theme(legend.position = "bottom", 
                       legend.direction = "horizontal",
                       legend.box.margin = margin(0, 0, 0, 0)))

# combine plots and save
g_plots <- plot_grid(fig_trt, fig_germ_lb, fig_germ_trt,
          nrow = 1,
          labels = letters[1:3],
          label_size = lg_txt) 

pdf("./output/ev_germination_raw_2019_litter_exp.pdf", width = 6, height = 2)
plot_grid(g_plots, leg, 
          nrow = 2,
          rel_heights = c(1, 0.15))
dev.off()


#### numbers for text ####
germ %>% group_by(treatment) %>%
  summarise(mean = mean_cl_boot(ev_prop_germ)[ ,1],
            lower = mean_cl_boot(ev_prop_germ)[ ,2],
            upper = mean_cl_boot(ev_prop_germ)[ ,3])


#### germination model ####

# priors
germ %>%
  filter(litter.g.m2 == 0) %>%
  summarize(mean = mean(ev_prop_germ))

x <- seq(0, 0.4, length.out = 20)
y <- 0.035/(1 + 10 * x)
plot(x, y, type = "l")

val <- seq(0, 20, length.out = 50)
dens <- dexp(val, 0.1)
plot(val, dens, type = "l")

# initial fit
evEstL2BhMod1 <- brm(bf(ev_prop_germ ~ maxEst / (1 + betaL * litter.g.cm2),
                        maxEst ~ 1,
                        betaL ~ 1,
                        nl = T),
                     data = germ, family = gaussian,
                     prior = c(prior(normal(0.04, 1), nlpar = "maxEst", class = "b", lb = 0),
                               prior(exponential(0.1), nlpar = "betaL", lb = 0)),
                     iter = 6000, warmup = 1000, chains = 1)
summary(evEstL2BhMod1)

# increase chains
evEstL2BhMod2 <- update(evEstL2BhMod1, chains = 3,
                        control = list(adapt_delta = 0.99))
summary(evEstL2BhMod2)
plot(evEstL2BhMod2)
pp_check(evEstL2BhMod1, ndraws = 50)


#### visualize ####

evEstPredDat <- tibble(litter.g.cm2 = seq(min(germ$litter.g.cm2), max(germ$litter.g.cm2), 
                                          length.out = 100)) %>%
  mutate(Est = fitted(evEstL2BhMod2, newdata = ., type = "response", re_formula = NA)[, "Estimate"],
         Est_lower = fitted(evEstL2BhMod2, newdata = ., type = "response", re_formula = NA)[, "Q2.5"],
         Est_upper = fitted(evEstL2BhMod2, newdata = ., type = "response", re_formula = NA)[, "Q97.5"],
         litter.g.m2 = litter.g.cm2 * 10000)

ggplot(germ, aes(x = litter.g.m2, y = ev_prop_germ)) +
  geom_ribbon(data = evEstPredDat, aes(y = Est, ymin = Est_lower, ymax = Est_upper), alpha = 0.3, color = NA) +
  geom_line(data = evEstPredDat, aes(y = Est)) +
  geom_point()


#### output ####
save(evEstL2BhMod2, file = "output/elymus_litter_establishment_bh_model_2019_litter_exp.rda")
write_csv(germ, "intermediate-data/elymus_litter_establishment_data_2018_litter_exp.csv")
