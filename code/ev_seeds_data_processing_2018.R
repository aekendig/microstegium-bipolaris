##### info ####

# file: ev_seeds_data_processing_2018
# author: Amy Kendig, Chris Wojan
# date last edited: 1/20/21
# goal: edit Ev seed data and check for errors
# background: spikelet counts for all samples, seed counts for a subset


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import all raw data files
di <- read_csv("./data/ev_spikelets_2018_density_litter_exp.csv")
di
de <- read_csv("./data/ev_seed_subset_2018_density_exp.csv")
de
dea <- read_csv("./data/ev_seeds_early_aug_2018_density_exp.csv")
dea
dla <- read_csv("./data/all_disease_seeds_late_aug_2018_density_exp.csv")
dla
ds <- read_csv("./data/ev_disease_seeds_sep_2018_density_exp.csv")
ds
dl <- read_csv("./data/ev_seeds_2018_litter_exp.csv")
dl


#### edit data ####

# check columns
unique(di$site)
unique(di$plot)
unique(di$plant)
filter(di, is.na(plant))$site
unique(di$collect_date)

unique(de$site)
unique(de$plot)
unique(de$treatment)
unique(de$age)
unique(de$ID)
filter(de, is.na(ID))
unique(de$collect_date)

# add and modify columns in di (spikelet counts)
di <- di %>%
  mutate(
    experiment = str_extract(site, "[aA-zZ]+") %>% recode("D" = "density", "L" = "litter"),
    treatment = str_extract(plot, "[aA-zZ]+") %>% recode("W" = "water", "F" = "fungicide"),
    plot = str_extract(plot, "[0-9]+") %>% as.numeric(),
    month = case_when(
      collect_date == 20180711 ~ "July",
      collect_date > 20180711 & collect_date < 20180827 ~ "early August",
      collect_date > 20180803 & collect_date < 20180924 ~ "late August",
      collect_date > 20180830 & collect_date < 20181022 ~ "September",
      collect_date == 20181022 ~ "October"),
    plant = ifelse(is.na(plant), "A", plant),
    ID = str_remove(plant, "EV"),
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
    age = case_when(ID == "A" ~ "adult",
                    focal == 0 & plot > 7 ~ "adult",
                    TRUE ~ "seedling"),
    ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID)
  )

# check that modifications worked
unique(di$experiment)
select(di, collect_date, month) %>% unique()
select(di, plant, age, ID, focal) %>% unique() %>% data.frame
filter(di, plant == "?")
filter(di, focal == 0) %>% select(plot, age) %>% unique()

# add and modify columns in de (seed subset count)
de <- de %>%
  mutate(
    treatment = recode(treatment, "W" = "water", "F" = "fungicide"),
    age = recode(age, "A" = "adult", "S" = "seedling"),
    ID = ifelse(is.na(ID), "A", ID),
    month = case_when(
      collect_date == 20180711 ~ "July",
      collect_date > 20180711 & collect_date < 20180827 ~ "early August",
      collect_date > 20180803 & collect_date < 20180924 ~ "late August",
      collect_date > 20180830 & collect_date < 20181022 ~ "September",
      collect_date == 20181022 ~ "October"),
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
    ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID)
  ) %>%
  mutate(collect_date = ifelse(site == "D3" & plot == 7 & treatment == "water" & month == "late August" & ID == "R2" & age == "seedling", 20180828, collect_date)) # correct to match spikelet data

# merge
d <- full_join(di, de)
d
filter(d, is.na(spikelets)) %>% data.frame()

# notes
unique(d$spikelet_notes)

# check some samples
d %>% filter(spikelet_notes == "unknown plant")
d %>% filter(spikelet_notes %in% c("1 OF 2", "2 OF 2")) %>% data.frame()

# add a removal column, sp column, edit notes, and check for duplicates
d <- d %>%
  mutate(
    ID_unclear = case_when(spikelet_notes == "unknown plant" ~ 1,
                       TRUE ~ 0), 
    spikelet_notes = ifelse(spikelet_notes %in% c("1 OF 2", "2 OF 2"), "2 samples", spikelet_notes),
    sp = "Ev"
  ) %>%
  group_by(site, plot, treatment, ID, age, month, focal) %>%
  mutate(reps = n()) %>%
  ungroup()

# add seed count from notes (Keith processed)
d2 <- d %>%
  filter(spikelet_notes == "actual seed count: 8 ; actual seed weight: 19mg") %>%
  mutate(seeds =  8,
         seed_weight.mg = 19,
         seed_date = 20190201,
         spikelet_notes = "processed by Keith Clay") %>%
  full_join(filter(d, spikelet_notes != "actual seed count: 8 ; actual seed weight: 19mg" | is.na(spikelet_notes)))

# assign month to field data
dea$month <- "early August"
dla$month <- "late August"
ds$month <- "September"

# make litter data long
dl <- dl %>%
  gather(key = month, value = seeds_collected, -c(site, plot, field_notes)) %>%
  mutate(
    month = case_when(month == "seeds_early_aug" ~ "early August",
                      month == "seeds_late_aug" ~ "late August",
                      month == "seeds_sep" ~ "September"),
    ID = "A",
    sp = "Ev"
  ) # note that late August data may be unreliable (in metadata)

# combine field data
df <- full_join(dl, dea) %>%
  full_join(dla) %>%
  full_join(ds)

# edit field data columns
unique(df$ID)

df <- df %>%
  mutate(
    experiment = str_extract(site, "[aA-zZ]+") %>% recode("D" = "density", "L" = "litter"),
    treatment = recode(treatment, "Water" = "water", "Fungicide" = "fungicide"),
    ID = case_when(ID == "A" ~ "A",
                   TRUE ~ str_remove(ID, "A")),
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
    age = case_when(ID == "A" ~ "adult",
                    focal == 0 & plot > 7 ~ "adult",
                    TRUE ~ "seedling"),
  )


#### check data ####

# collection record
dc <- full_join(d2, df)

# records in which seeds were collected and recorded
dc %>%
  filter((!is.na(spikelets) | !is.na(seeds)) & seeds_collected == 1) %>%
  dim() 
# 397

# records in which they weren't recorded
dcc <- dc %>%
  filter((!is.na(spikelets) | !is.na(seeds)) & (is.na(seeds_collected) | seeds_collected == 0) & !(month %in% c("July", "October")))
dim(dcc)
# 12
data.frame(dcc)
# recording mistakes or mislabelled seed bags

# recorded as collected, but missing
dcm <- dc %>%
  filter(is.na(spikelets) & seeds_collected == 1)
dim(dcm)
# 9
data.frame(dcm)
# one is okay

# look for overlap with not recorded and missing
dcm %>% select(site, plot, treatment, age, ID, month, spikelet_notes)
dcc %>% select(site, plot, treatment, age, ID, month, spikelet_notes)

# indicate which ones were probably mislabelled bags
dcm <- dcm %>%
  mutate(
    probably_collected = case_when(plot %in% c(10, 4) | (plot == 7 & site == "D3") ~ 1))
dcc <- dcc %>%
  mutate(
    probably_mislabelled = case_when((plot == 10 & site == "D2") | plot %in% c(4, 7) ~ 1)
  )

# look for overlaps between missing and duplicates
dcm %>% filter(is.na(probably_collected)) %>% select(site, plot, treatment, age, ID, month, spikelet_notes)
d2 %>% filter(reps > 1) %>% select(site, plot, treatment, age, ID, month, collect_date, spikelet_notes)
# some of the duplicates may have been mislabelled seed bags (D2 6W samples, D4, 4F EvA, D1 7F, D1 7W - all of them except the ones noted to have 2 samples)

# combine separated samples
colnames(d2)
group_cols <- d2 %>% select(site:spikelet_date, spikelet_notes:seed_date, ID_unclear:sp) %>% names

d3 <- d2 %>%
  filter(spikelet_notes == "2 samples") %>%
  group_by(.dots = group_cols) %>%
  summarise(spikelets = sum(spikelets),
            spikelet_weight.mg = sum(spikelet_weight.mg),
            seeds = sum(seeds),
            seed_weight.mg = sum(seed_weight.mg),
            reps = 1
  )  %>%
  ungroup() %>%
  full_join(filter(d2, spikelet_notes != "2 samples" | is.na(spikelet_notes)))

# remove plants that were probably mislabelled
dcc %>% filter(probably_mislabelled == 1) %>% data.frame()

d3 <- d3 %>%
  ungroup() %>%
  mutate(
    ID_unclear = case_when(reps > 1 |
                         (site == "D2" & plot == 10 & treatment == "water" & plant == "EVAR4" & month == "early August") |
                         (site == "D2" & plot == 10 & treatment == "water" & plant == "EVAR1" & month == "late August") |
                         (site == "D4" & plot == 7 & treatment == "water" & plant == "EV1" & month == "late August") |
                         (site == "D4" & plot == 7 & treatment == "water" & plant == "EV3" & month == "late August") |
                         (site == "D2" & plot == 4 & treatment == "water" & plant == "EV2" & month == "late August") ~ 1,
                       TRUE ~ ID_unclear)
  )


#### visualize ####

# spikelet number
d3 %>%
  ggplot(aes(x = spikelets)) +
  geom_histogram()

d3 %>%
  filter(spikelets > 200) %>%
  select(age, month)
# adults during late August and September

# spikelet weight
d3 %>%
  ggplot(aes(x = spikelet_weight.mg)) +
  geom_histogram()

d3 %>%
  filter(spikelet_weight.mg > 3500) %>%
  select(age, month)
# adults in late August

d3 %>%
  ggplot(aes(x = spikelets, y = spikelet_weight.mg, colour = age, shape = month)) +
  geom_point()

# large weight:number
d3 %>% filter(spikelet_weight.mg > 4000)
# Chris re-weighed these

# seed count
d3 %>%
  ggplot(aes(x = seeds)) +
  geom_histogram()

d3 %>%
  filter(seeds > 150) %>%
  select(age, month)
# adults in late August/September

# seed weight
d3 %>%
  ggplot(aes(x = seed_weight.mg)) +
  geom_histogram()

d3 %>%
  filter(seed_weight.mg > 600)%>%
  select(age, month)
# same as above

# seed and spikelet
d3 %>%
  ggplot(aes(x = spikelets, y = seeds)) +
  geom_point(aes(colour = age, shape = treatment)) +
  geom_smooth(method = lm)

d3 %>%
  ggplot(aes(x = spikelets, y = seeds)) +
  geom_smooth(method = lm)

# seed and spikelet weight
d3 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seed_weight.mg))  +
  geom_point(aes(colour = age, shape = treatment)) +
  geom_smooth(method = lm)

# seed count and spikelet weight
d3  %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_smooth(method = lm, size = 0.5, colour = "black")

# seed count and spikelet weight again with site
d3 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = site, shape = treatment)) +
  geom_smooth(method = lm)

# seed count and spikelet weight by treatment
d3 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = site, shape = age)) +
  facet_wrap(~treatment) +
  geom_smooth(method = lm)
# look very similar

# ratio over time
d3 %>%
  ggplot(aes(x = collect_date, y = seeds/spikelet_weight.mg)) + 
  geom_point(aes(colour = age, shape = treatment), size = 2) + 
  geom_smooth(method = "lm")
# not a strong trend

# ratio over spikelet count
d3 %>%
  ggplot(aes(x = spikelets, y = seeds/spikelet_weight.mg)) + 
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_smooth(method = "lm")
# hard tell because sample size gets smaller with more spikelets, but not a strong trend


#### stats ####

# regression fit 1: seeds and spikelet counts
m1 <- lm(seeds ~ spikelets, data = d3)
summary(m1) # R2 is 0.82

m1b <- lm(seeds ~ 0 + spikelets, data = d3)
summary(m1b) # R2 is 0.90
AIC(m1, m1b) # b is slightly lower

# regression fit 2: seed and spikelet weight
m2 <- lm(seed_weight.mg ~ spikelet_weight.mg, data = d3)
summary(m2) # R2 is 0.93

m2b <- lm(seed_weight.mg ~ 0 + spikelet_weight.mg, data = d3)
summary(m2b) # R2 is 0.95
AIC(m2, m2b)  # b is slightly lower

# regression fit 3: seed count and spikelet weight
m3 <- lm(seeds ~ spikelet_weight.mg, data = d3)
summary(m3) # R2 is 0.94

m3b <- lm(seeds ~ 0 + spikelet_weight.mg, data = d3)
summary(m3b) # R2 is 0.97
AIC(m3, m3b)  # b is slightly lower

# refit m3b without the plant that is missing spikelet counts
m3c <- lm(seeds ~ 0 + spikelet_weight.mg, data = filter(d3, !is.na(spikelets)))

# compare regressions for estimating seeds
AIC(m1b, m3c)

# modify m3 to include more info
m4 <- lm(seeds ~ spikelet_weight.mg * treatment, data = d2)
summary(m4)
m5 <- lm(seeds ~ spikelet_weight.mg + treatment, data = d2)
summary(m5)
AIC(m3, m4, m5)
# m3 is the best fit model and adding treatment hardly changes the slope of the relationship.

# visualize zero intercept model
m3b_plot <- d3  %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_abline(intercept = 0, slope = coef(m3b)[1], linetype = "dashed") +
  theme_bw() +
  xlab("spikelet weight (mg)")
m3b_plot

# residuals
d3 %>%
  filter(!is.na(seeds)) %>%
  mutate(resid = as.numeric(resid(m3b))) %>%
  ggplot(aes(x = spikelet_weight.mg, y = resid)) +
  geom_point(aes(colour = site, shape = age)) + 
  facet_wrap(~treatment)
# may underestimate water more than fungicide

# add seed estimate
d4 <- d3 %>%
  mutate(
    seeds = ifelse(is.na(seeds), spikelet_weight.mg * coef(m3b)[1], seeds)
  )

# format data for saving
d4d <- d4 %>%
  filter(experiment == "density") %>%
  select(month, collect_date, spikelet_date, seed_date, site, plot, treatment, plant, sp, ID, age, focal, spikelets, spikelet_weight.mg, seeds, seed_weight.mg, reps, ID_unclear, spikelet_notes)

d4l <- d4 %>%
  filter(experiment == "litter") %>%
  select(month, collect_date, spikelet_date, seed_date, site, plot, treatment, plant, sp, ID, age, focal, spikelets, spikelet_weight.mg, seeds, seed_weight.mg, reps, ID_unclear, spikelet_notes)


#### save output ####

# seeds prediction relationship
pdf("./output/ev_seeds_data_processing_2018_seeds_spikelet_weight.pdf", width = 5, height = 4)
m3b_plot
dev.off()

# density and litter experiment
write_csv(d4d, "./intermediate-data/ev_processed_seeds_2018_density_exp.csv")
write_csv(d4l, "./intermediate-data/ev_processed_seeds_2018_litter_exp.csv")
write_csv(d3, "./intermediate-data/ev_spikelets_seeds_2018.csv")
