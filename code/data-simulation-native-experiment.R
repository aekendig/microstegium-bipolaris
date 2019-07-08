##### info ####

# file: data-simulation-native-experiment
# author: Amy Kendig, Vida Svahnstrom
# date last edited: 7/1/19
# goal: create fake data and work through analyses


#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)


#### create data ####

# dataframe with treatments
dat <- data.frame(
  inoculation = rep(c(0, 1), each = 60),
  sp = rep(rep(c("Ev", "Es", "Pc"), each = 20), 2),
  density = rep(rep(c(0, 2, 10, 50, 100), each = 4), 6),
  rep = rep(1:4, 30)
)

# examine dataframe
head(dat)
str(dat)
tail(dat)

# save base dataframe
write.csv(dat, "./output/data-simulation-native-experiment.csv", row.names = F)

# simulate biomasss data
set.seed(42)
(spbio <- rnorm(3, 20, 2))

set.seed(42)
dat$biomass <- with(dat, 
                   spbio[as.numeric(sp)] 
                   - abs(rnorm(1, sd = 0.1)) * inoculation 
                   - abs(rnorm(1)) * density/20
                   - rnorm(120)) 

# make a categorical inoculation column
dat <- dat %>%
  mutate(infection = ifelse(inoculation == 0, "no", "yes"))

# check
head(dat)
str(dat)
aggregate(biomass ~ inoculation * sp * density, data = dat, FUN = mean)


#### visualize ####

ggplot(dat, aes(x = density, y = biomass, color = infection)) +
  geom_point(alpha = 0.5) + 
  facet_wrap(~sp) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 2) +
  geom_smooth(method = "lm")

ggplot(dat, aes(x = biomass)) +
  geom_histogram()

filter(dat, biomass < 15)


#### stats ####

mod <- lm(biomass ~ inoculation * sp * density, data = dat)
summary(mod)
plot(mod)


