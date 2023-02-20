#### outputs ####

# litter_reu_mv_establishment_model.rda
# litter_reu_ev_establishment_model.rda

#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(brms)

# import data
MvEstDat <- read_csv("data/litter_reu_mv_establishment_data.csv")
EvEstDat <- read_csv("data/litter_reu_ev_establishment_data.csv")


#### edit data ####

litterDat <- MvEstDat %>%
  select(Litter.g) %>%
  unique() %>%
  arrange() %>%
  mutate(Litter.g.m2 = c(0, 50, 100, 200))

# remove competition treatments
# rename columns
# add litter column
MvEstDat2 <- MvEstDat %>%
  filter(Competition == 0) %>%
  rename(PropEst = PropEstMvDenCor) %>%
  left_join(litterDat)

EvEstDat2 <- EvEstDat %>%
  filter(Competition == 0) %>%
  rename(PropEst = PropEstEv) %>%
  left_join(litterDat)

# initial visualizations
ggplot(MvEstDat2, aes(x = Litter.g.m2, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()

ggplot(EvEstDat2, aes(x = Litter.g.m2, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()


#### fit models ####

# beta prior
x <- 0:200
y <- 0.85/(1 + 0.001 * x)
plot(x, y, type = "l")

val <- seq(0, 1, length.out = 50)
dens <- dexp(val, 1)
plot(val, dens, type = "l")

# Mv
mv_bh_mod <- brm(data = MvEstDat2, family = gaussian,
                 bf(PropEst ~ e0/(1 + beta * Litter.g.m2),
                    e0 ~ 1, 
                    beta ~ 1, 
                    nl = T),
                 prior <- c(prior(normal(0.85, 1), nlpar = "e0", lb = 0),
                            prior(exponential(1), nlpar = "beta", lb = 0),
                            prior(cauchy(0, 1), class = sigma)),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)

prior_summary(mv_bh_mod)
summary(mv_bh_mod, digits = 5)
mv_bh_mod <- update(mv_bh_mod, chains = 3)
summary(mv_bh_mod)
plot(mv_bh_mod)
fixef(mv_bh_mod)

# Ev
ev_bh_mod <- update(mv_bh_mod, newdata = EvEstDat2)

summary(ev_bh_mod)
plot(ev_bh_mod)
fixef(ev_bh_mod)

# save
save(mv_bh_mod, file = "output/litter_reu_mv_establishment_model.rda")
save(ev_bh_mod, file = "output/litter_reu_ev_establishment_model.rda")
