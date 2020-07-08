##### info ####

# file: mv_seed_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/7/20
# goal: analyze Mv seeds


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)

# import data
plant_dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
plot_dat <- read_csv("intermediate-data/mv_plot_level_seeds_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# add density to data for visualization
plot_dat2 <- left_join(plot_dat, plots)

# data for stats
plant_dat2 <- left_join(plant_dat, treat) %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0") %>% as.numeric(),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide"),
         exp_plot = paste(plot, toupper(substr(treatment, 1, 1)), sep = ""),
         mv_density = case_when(plot %in% c(1, 5:10) ~ 3,
                                TRUE ~ background_density + 3),
         ev_seedling_density = case_when(plot %in% c(5:7) ~ background_density + 3,
                                         TRUE ~ 3),
         ev_adult_density = case_when(plot %in% c(8:10) ~ background_density + 1,
                                         TRUE ~ 1))

# split by background
mv_dat <- plant_dat2 %>%
  filter(background %in% c("none", "Mv seedling"))

ev_seedling_dat <- plant_dat2 %>%
  filter(background %in% c("none", "Ev seedling"))

ev_adult_dat <- plant_dat2 %>%
  filter(background %in% c("none", "Ev adult"))


#### visualize ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside")

# treatment figure
seed_viz <- ggplot(plot_dat2, aes(x = background_density, y = seeds, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  xlab("Background density") +
  temp_theme
seed_viz
# seeds drop with Microstegium density, but are pretty unaffected by Elymus density

# flower seeds
seed_viz %+%
  aes(y = flower_seeds)
# similar pattern

# stem seeds
seed_viz %+%
  aes(y = stem_seeds)
# similar pattern

# individual seeds
ggplot(mv_dat, aes(x = mv_density, y = seeds, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme


#### Mv background model ####

# Beverton-Holt simulation
bh_fun <- function(dat_in, a, y0){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$mv_density)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = mv_density, y = seeds)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# try values of a
bh_fun(mv_dat, 0.05, 1400)

# distributions
x = 0:10
y = dexp(x, 0.5)
plot(x, y, type = "l")

# greenhouse experiment: fungicide reduced seed heads
30 / 37
1 - 0.81

# model with greenhouse info
mv_seed_fung_mod <- brm(data = mv_dat, family = gaussian,
                        bf(seeds ~ (s0 * fung) / (1 + alpha * mv_density),
                           s0 ~ fungicide + (1|site/exp_plot),
                           fung ~ fungicide,
                           alpha ~ 0 + Treatment,
                           nl = T),
                        prior <- c(prior(normal(1400, 100), nlpar = "s0", coef = "Intercept"), 
                                   prior(normal(0, 100), nlpar = "s0", coef = "fungicide"),
                                   prior(cauchy(0, 1), nlpar = "s0", class = sd),
                                   prior(normal(1, 0.001), nlpar = "fung", coef = "Intercept"),
                                   prior(normal(-0.19, 0.001), nlpar = "fung", coef = "fungicide"),
                                   prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                   prior(cauchy(0, 1), class = sigma)),
                        iter = 6000, warmup = 1000, chains = 1, cores = 1,
                        control = list(adapt_delta = 0.99))
                   
# check model and increase chains
summary(mv_seed_fung_mod)
prior_summary(mv_seed_fung_mod)
mv_seed_fung_mod <- update(mv_seed_fung_mod, chains = 3)
plot(mv_seed_fung_mod)
pp_check(mv_seed_fung_mod, nsamples = 100)                           

# model without greenhouse info
mv_seed_mod <- brm(data = mv_dat, family = gaussian,
                        bf(seeds ~ s0 / (1 + alpha * mv_density),
                           s0 ~ fungicide + (1|site/exp_plot),
                           alpha ~ 0 + Treatment,
                           nl = T),
                        prior <- c(prior(normal(1400, 100), nlpar = "s0", coef = "Intercept"), 
                                   prior(normal(0, 100), nlpar = "s0", coef = "fungicide"),
                                   prior(cauchy(0, 1), nlpar = "s0", class = sd),
                                   prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                   prior(cauchy(0, 1), class = sigma)),
                        iter = 6000, warmup = 1000, chains = 1, cores = 1,
                        control = list(adapt_delta = 0.99))

# check model and increase chains
summary(mv_seed_mod)
prior_summary(mv_seed_mod)
mv_seed_mod <- update(mv_seed_mod, chains = 3)
plot(mv_seed_mod)
pp_check(mv_seed_mod, nsamples = 100)    


#### Ev seedling background model ####

# average seeds
mean(ev_seedling_dat$seeds)

# try values of a (want flat line)
bh_fun(mv_dat, 0, 1400)

# distributions
x = 0:10
y = dexp(x, 100) # reciprocal of the stan rate
plot(x, y, type = "l")

# model without greenhouse info
mv_seed_ev_seedling_mod <- brm(data = ev_seedling_dat, family = gaussian,
                               bf(seeds ~ s0 / (1 + alpha * ev_seedling_density),
                                  s0 ~ fungicide + (1|site/exp_plot),
                                  alpha ~ 0 + Treatment,
                                  nl = T),
                               prior <- c(prior(normal(1037, 100), nlpar = "s0", coef = "Intercept"), 
                                          prior(normal(0, 100), nlpar = "s0", coef = "fungicide"),
                                          prior(cauchy(0, 1), nlpar = "s0", class = sd),
                                          prior(exponential(0.01), nlpar = "alpha", lb = 0),
                                          prior(cauchy(0, 1), class = sigma)),
                               iter = 6000, warmup = 1000, chains = 1, cores = 1,
                               control = list(adapt_delta = 0.99))
                                    
# check model and increase chains
summary(mv_seed_ev_seedling_mod)
prior_summary(mv_seed_ev_seedling_mod)
mv_seed_ev_seedling_mod <- update(mv_seed_ev_seedling_mod, chains = 3)
plot(mv_seed_ev_seedling_mod)
pp_check(mv_seed_ev_seedling_mod, nsamples = 100)

# model with greenhouse info
mv_seed_ev_seedling_fung_mod <- brm(data = ev_seedling_dat, family = gaussian,
                                    bf(seeds ~ (s0 * fung) / (1 + alpha * ev_seedling_density),
                                       s0 ~ fungicide + (1|site/exp_plot),
                                       fung ~ fungicide,
                                       alpha ~ 0 + Treatment,
                                       nl = T),
                                    prior <- c(prior(normal(1037, 100), nlpar = "s0", coef = "Intercept"), 
                                               prior(normal(0, 100), nlpar = "s0", coef = "fungicide"),
                                               prior(cauchy(0, 1), nlpar = "s0", class = sd),
                                               prior(normal(1, 0.001), nlpar = "fung", coef = "Intercept"),
                                               prior(normal(-0.19, 0.001), nlpar = "fung", coef = "fungicide"),
                                               prior(exponential(0.01), nlpar = "alpha", lb = 0),
                                               prior(cauchy(0, 1), class = sigma)),
                                    iter = 6000, warmup = 1000, chains = 1, cores = 1,
                                    control = list(adapt_delta = 0.99))
                        

# check model and increase chains
summary(mv_seed_ev_seedling_fung_mod)
prior_summary(mv_seed_ev_seedling_fung_mod)
mv_seed_ev_seedling_fung_mod <- update(mv_seed_ev_seedling_fung_mod, chains = 3)
plot(mv_seed_ev_seedling_fung_mod)
pp_check(mv_seed_ev_seedling_fung_mod, nsamples = 100)


#### Ev adult background model ####

# average seeds
mean(ev_adult_dat$seeds, na.rm = T)

# model with greenhouse info
mv_seed_ev_adult_fung_mod <- brm(data = ev_adult_dat, family = gaussian,
                                    bf(seeds ~ (s0 * fung) / (1 + alpha * ev_adult_density),
                                       s0 ~ fungicide + (1|site/exp_plot),
                                       fung ~ fungicide,
                                       alpha ~ 0 + Treatment,
                                       nl = T),
                                    prior <- c(prior(normal(1034, 100), nlpar = "s0", coef = "Intercept"), 
                                               prior(normal(0, 100), nlpar = "s0", coef = "fungicide"),
                                               prior(cauchy(0, 1), nlpar = "s0", class = sd),
                                               prior(normal(1, 0.001), nlpar = "fung", coef = "Intercept"),
                                               prior(normal(-0.19, 0.001), nlpar = "fung", coef = "fungicide"),
                                               prior(exponential(0.01), nlpar = "alpha", lb = 0),
                                               prior(cauchy(0, 1), class = sigma)),
                                    iter = 6000, warmup = 1000, chains = 1, cores = 1,
                                    control = list(adapt_delta = 0.999))


# check model and increase chains
summary(mv_seed_ev_adult_fung_mod)
mv_seed_ev_adult_fung_mod <- update(mv_seed_ev_adult_fung_mod, chains = 3)
# model has high R-hat values and lots of divergent transitions

# seeds in control treatment
mean(filter(ev_adult_dat, treatment == "water")$seeds, na.rm = T)

# linear model for Ev adult without fungicide effect
# need to do for fungicide effece model
mv_seed_ev_adult_mod <- brm(data = ev_adult_dat, family = gaussian,
                            seeds ~ ev_adult_density * fungicide + (1|site/exp_plot),
                            prior <- c(prior(normal(1236, 100), class = "Intercept"), 
                                       prior(normal(0, 10), class = "b"),
                                       prior(cauchy(0, 1), class = "sigma")),
                            iter = 6000, warmup = 1000, chains = 1, cores = 1)

# check model and increase chains
summary(mv_seed_ev_adult_mod)
mv_seed_ev_adult_mod <- update(mv_seed_ev_adult_mod, chains = 3)
plot(mv_seed_ev_adult_mod)
pp_check(mv_seed_ev_adult_mod, nsamples = 100)  

# fungicide effect
1234 * 0.19                                 

# linear model for Ev adult with fungicide effect
mv_seed_ev_adult_fung_mod <- brm(data = ev_adult_dat, family = gaussian,
                            seeds ~ fungicide + ev_adult_density * Treatment + (1|site/exp_plot),
                            prior <- c(prior(normal(1234, 100), class = "Intercept"), 
                                       prior(normal(0, 10), class = "b"),
                                       prior(normal(-234, 0.001), coef = "fungicide"),
                                       prior(cauchy(0, 1), class = "sigma")),
                            iter = 6000, warmup = 1000, chains = 1, cores = 1)

# check model and increase chains
summary(mv_seed_ev_adult_fung_mod)
mv_seed_ev_adult_fung_mod <- update(mv_seed_ev_adult_fung_mod, chains = 3)
plot(mv_seed_ev_adult_fung_mod)
pp_check(mv_seed_ev_adult_fung_mod, nsamples = 100)  



#### figures ####

# colors
col_pal = c("#a6611a", "#018571")

# simulate data
mv_sim_dat <- tibble(mv_density = rep(seq(0, 67, length.out = 300), 2),
                     fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(mv_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(mv_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"])

mv_fung_sim_dat <- tibble(mv_density = rep(seq(0, 67, length.out = 300), 2),
                          fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_seed_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(mv_seed_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(mv_seed_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"])

mv_ev_seedling_sim_dat <- tibble(ev_seedling_density = rep(seq(0, 19, length.out = 300), 2),
                                      fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_seed_ev_seedling_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(mv_seed_ev_seedling_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(mv_seed_ev_seedling_mod, newdata = ., re_formula = NA)[, "Q97.5"])

mv_ev_seedling_fung_sim_dat <- tibble(ev_seedling_density = rep(seq(0, 19, length.out = 300), 2),
                          fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_seed_ev_seedling_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(mv_seed_ev_seedling_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(mv_seed_ev_seedling_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"])

mv_ev_adult_sim_dat <- tibble(ev_adult_density = rep(seq(0, 9, length.out = 300), 2),
                                   fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_seed_ev_adult_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(mv_seed_ev_adult_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(mv_seed_ev_adult_mod, newdata = ., re_formula = NA)[, "Q97.5"])

mv_ev_adult_fung_sim_dat <- tibble(ev_adult_density = rep(seq(0, 9, length.out = 300), 2),
                                      fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(Treatment = recode(fungicide, "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_seed_ev_adult_fung_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(mv_seed_ev_adult_fung_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(mv_seed_ev_adult_fung_mod, newdata = ., re_formula = NA)[, "Q97.5"])

# Mv plot
mv_fig <- ggplot(mv_dat, aes(x = mv_density, y = seeds)) + 
  geom_ribbon(data = mv_sim_dat, aes(ymin = seeds_lower, ymax = seeds_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegiun"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " seeds (", plant^-1, ")", sep = ""))) +
  temp_theme

mv_fung_fig <- ggplot(mv_dat, aes(x = mv_density, y = seeds)) + 
  geom_ribbon(data = mv_fung_sim_dat, aes(ymin = seeds_lower, ymax = seeds_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_fung_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Microstegiun"), " density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " seeds (", plant^-1, ")", sep = ""))) +
  temp_theme

# Ev seedling plot
ev_seedling_fig <- ggplot(ev_seedling_dat, aes(x = ev_seedling_density, y = seeds)) + 
  geom_ribbon(data = mv_ev_seedling_sim_dat, aes(ymin = seeds_lower, ymax = seeds_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_ev_seedling_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " seedling density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " seeds (", plant^-1, ")", sep = ""))) +
  temp_theme

ev_seedling_fung_fig <- ggplot(ev_seedling_dat, aes(x = ev_seedling_density, y = seeds)) + 
  geom_ribbon(data = mv_ev_seedling_fung_sim_dat, aes(ymin = seeds_lower, ymax = seeds_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_ev_seedling_fung_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " seedling density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " seeds (", plant^-1, ")", sep = ""))) +
  temp_theme

# Ev adult plot
ev_adult_fig <- ggplot(ev_adult_dat, aes(x = ev_adult_density, y = seeds)) + 
  geom_ribbon(data = mv_ev_adult_sim_dat, aes(ymin = seeds_lower, ymax = seeds_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_ev_adult_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " adult density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " seeds (", plant^-1, ")", sep = ""))) +
  temp_theme

ev_adult_fung_fig <- ggplot(ev_adult_dat, aes(x = ev_adult_density, y = seeds)) + 
  geom_ribbon(data = mv_ev_adult_fung_sim_dat, aes(ymin = seeds_lower, ymax = seeds_upper, fill = Treatment), alpha = 0.3) +
  geom_line(data = mv_ev_adult_fung_sim_dat, aes(color = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5, aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, aes(fill = Treatment)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste(italic("Elymus"), " adult density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste(italic("Microstegiun"), " seeds (", plant^-1, ")", sep = ""))) +
  temp_theme


#### output ####
pdf("output/mv_seed_analysis_2019_density_exp.pdf")
mv_fig
ev_seedling_fig
ev_adult_fig
dev.off()

pdf("output/mv_seed_analysis_greenhouse_fungicide_2019_density_exp.pdf")
mv_fung_fig
ev_seedling_fung_fig
ev_adult_fung_fig
dev.off()

save(mv_seed_mod, file ="output/Mv_seed_model_mv_background_2019_density_exp.rda")
save(mv_seed_ev_seedling_mod, file ="output/Mv_seed_model_ev_seedling_background_2019_density_exp.rda")
save(mv_seed_ev_adult_mod, file ="output/Mv_seed_model_ev_adult_background_2019_density_exp.rda")
save(mv_seed_fung_mod, file ="output/Mv_seed_model_mv_background_greenhouse_fungicide_2019_density_exp.rda")
save(mv_seed_ev_seedling_fung_mod, file ="output/Mv_seed_model_ev_seedling_background_greenhouse_fungicide_2019_density_exp.rda")
save(mv_seed_ev_adult_fung_mod, file ="output/Mv_seed_model_ev_adult_background_greenhouse_fungicide_2019_density_exp.rda")
