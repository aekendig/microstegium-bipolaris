##### info ####

# file: mv-seed-sterilization
# author: Amy Kendig, Penny Reif
# date last edited: 4/25/19
# goal: choose a sterilization protocol to isolate specific fungi


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(rethinking)

# import data
s <- read_csv("./data/mv-seed-sterilization-trial-2-density-exp.csv")
g <- read_csv("./data/mv-seed-sterilization-germination-trial-2-density-exp.csv")


#### edit data ####

s <- s %>%
  mutate(
    seeds_infected.num = recode(seeds_infected, L = 2, M = 8, H = 16)
  )


#### visualize ####

# fungal growth
s %>%
  filter(!is.na(cover.prop)) %>%
  ggplot(aes(x = time.min, y = cover.prop, colour = fungi)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2)) + 
  facet_wrap(~as.factor(concentration.perc))

# fungal growth without dark fungi
s %>%
  filter(!is.na(cover.prop) & fungi != "dark") %>%
  ggplot(aes(x = time.min, y = cover.prop, colour = fungi)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2)) + 
  facet_wrap(~as.factor(concentration.perc), scales = "free")

# fungal growth on seeds
s %>%
  filter(!is.na(seeds_infected.num)) %>%
  ggplot(aes(x = time.min, y = seeds_infected.num, colour = fungi)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2)) + 
  facet_wrap(~as.factor(concentration.perc), scales = "free")

# germination rate
g %>%
  ggplot(aes(x = time.min, y = germination, colour = as.factor(concentration.perc))) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2))


#### save figure ####

# fungal growth
pdf("./output/mv-seed-sterilization-cover.pdf", width = 5, height = 5)
s %>%
  filter(!is.na(cover.prop)) %>%
  ggplot(aes(x = time.min, y = cover.prop, colour = fungi)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2)) + 
  facet_wrap(~as.factor(concentration.perc)) +
  xlab("Time (min)") +
  ylab("Proportion of the plate covered") +
  ggtitle("Bleach concentration (percent)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = c(0.15, 0.8),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("./output/mv-seed-sterilization-seeds.pdf", width = 6.5, height = 5)
s %>%
  ggplot(aes(x = time.min, y = seeds_infected.num, colour = fungi)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2)) + 
  facet_wrap(~as.factor(concentration.perc)) +
  xlab("Time (min)") +
  ylab("Approximate number of seeds infected") +
  ggtitle("Bleach concentration (percent)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = c(0.2, 0.8),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("./output/mv-seed-sterilization-germination.pdf", width = 6.5, height = 5)
g %>%
  ggplot(aes(x = time.min, y = germination)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.2)) + 
  facet_wrap(~as.factor(concentration.perc)) +
  xlab("Time (min)") +
  ylab("Number of seeds germinated") +
  ggtitle("Bleach concentration (percent)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 12, hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


#### fit models ####

# edit germination dataframe
germ <- data.frame(g) %>%
  mutate(
    conc = scale(concentration.perc),
    time = scale(time.min)
  )

# germination models
m.g <- map2stan(
  alist(
    germination ~ dbinom(30, p),
    logit(p) <- a + bc * conc + bt * time,
    a ~ dnorm(0, 10),
    bc ~ dnorm(0, 10),
    bt ~ dnorm(0, 10)),
  data = germ, chains = 2, iter = 2500, warmup = 500)

precis(m.g) # looks like the effects of time and concentration are week

m.g.c <- map2stan(
  alist(
    germination ~ dbinom(30, p),
    logit(p) <- a + bc * conc,
    a ~ dnorm(0, 10),
    bc ~ dnorm(0, 10),
    bt ~ dnorm(0, 10)),
  data = germ, chains = 2, iter = 2500, warmup = 500)

m.g.t <- map2stan(
  alist(
    germination ~ dbinom(30, p),
    logit(p) <- a + bt * time,
    a ~ dnorm(0, 10),
    bc ~ dnorm(0, 10),
    bt ~ dnorm(0, 10)),
  data = germ, chains = 2, iter = 2500, warmup = 500)

m.g.a <- map2stan(
  alist(
    germination ~ dbinom(30, p),
    logit(p) <- a,
    a ~ dnorm(0, 10)),
  data = germ, chains = 2, iter = 2500, warmup = 500)

# compare germination models
compare(m.g.a, m.g.c, m.g.t, m.g) # intercept model is the best

# make plate number
germ$plate = 1:nrow(germ)

# edit seed infection dataframe
seed <- data.frame(s) %>%
  mutate(cover = cover.prop * 100,
         d = case_when(fungi == "dark" ~ 1,
                              TRUE ~ 0),
         w = case_when(fungi == "white" ~ 1,
                       TRUE ~ 0)) %>%
  full_join(select(germ, -germination))

# number of seeds infected models
m.s <- map2stan(
  alist(
    seeds_infected.num ~ dbinom(30, p),
    logit(p) <- a +
      a_plate[plate] +
      bc * conc +
      bt * time +
      bd * d +
      bw * w +
      bcd * conc * d +
      btd * time * d + 
      bcw * conc * w + 
      btw * time * w,
    a ~ dnorm(0, 10),
    a_plate[plate] ~ dnorm(0, sigma_plate),
    bc ~ dnorm(0, 10),
    bt ~ dnorm(0, 10),
    bd ~ dnorm(0, 10),
    bw ~ dnorm(0, 10),
    bcd ~ dnorm(0, 10),
    btd ~ dnorm(0, 10),
    bcw ~ dnorm(0, 10),
    btw ~ dnorm(0, 10),
    sigma_plate ~ dcauchy(0, 1)),
  data = seed, chains = 2, iter = 2500, warmup = 500)

precis(m.s) # looks like concentration affects dark fungi, but not the others

m.s.d <- map2stan(
  alist(
    seeds_infected.num ~ dbinom(30, p),
    logit(p) <- a +
      a_plate[plate] +
      bc * conc +
      bt * time +
      bd * d +
      bw * w +
      btd * time * d + 
      bcw * conc * w + 
      btw * time * w,
    a ~ dnorm(0, 10),
    a_plate[plate] ~ dnorm(0, sigma_plate),
    bc ~ dnorm(0, 10),
    bt ~ dnorm(0, 10),
    bd ~ dnorm(0, 10),
    bw ~ dnorm(0, 10),
    btd ~ dnorm(0, 10),
    bcw ~ dnorm(0, 10),
    btw ~ dnorm(0, 10),
    sigma_plate ~ dcauchy(0, 1)),
  data = seed, chains = 2, iter = 2500, warmup = 500)

# compare germination models
compare(m.s, m.s.d) # slightly higher weight given to the one with the interaction

# edit data for cover model
seed2 <- seed %>%
  filter(!is.na(cover)) %>%
  mutate(plate2 = paste(concentration.perc, time.min, replicate, sep = "") %>% as.factor() %>% as.numeric()) 

# infected cover models
m.c <- map2stan(
  alist(
    cover ~ dnorm(mu, sigma),
    mu <- a +
      a_plate[plate2] +
      bc * conc +
      bt * time +
      bd * d +
      bw * w +
      bcd * conc * d +
      btd * time * d + 
      bcw * conc * w + 
      btw * time * w,
    a ~ dnorm(0, 10),
    a_plate[plate2] ~ dnorm(0, sigma_plate),
    bc ~ dnorm(0, 10),
    bt ~ dnorm(0, 10),
    bd ~ dnorm(0, 10),
    bw ~ dnorm(0, 10),
    bcd ~ dnorm(0, 10),
    btd ~ dnorm(0, 10),
    bcw ~ dnorm(0, 10),
    btw ~ dnorm(0, 10),
    sigma ~ dunif(0, 100),
    sigma_plate ~ dcauchy(0, 1)),
  data = seed2, chains = 2, iter = 2500, warmup = 500)

precis(m.c) # dark fungi is higher
plot(m.c)
plot(precis(m.c))

# plot output
pdf("./output/mv-seed-sterilization-models.pdf", width = 8, height = 5)
par(mfrow = c(1, 3))
plot(precis(m.g), main = "germinated seeds")
plot(precis(m.s), main = "infected seeds")
plot(precis(m.c), main = "area covered")
dev.off()
