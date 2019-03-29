##### info ####

# file: ev-seeds-data-processing
# author: Amy Kendig, Chris Wojan
# date last edited: 3/18/19
# goal: edit Ev seed data and check for errors
# background: spikelet counts for all samples, seed counts for a subset


#### set up ####

# load packages
#library(tidyverse)

# import all raw data files
di <- read_csv("./data/ev-spikelets-2018-density-litter-exp.csv")
di
de <- read_csv("./data/ev-seed-subset-2018-density-exp.csv")
de
dj <- read_csv("./data/ev-seeds-survival-jul-2018-density-exp.csv")
dj
dea <- read_csv("./data/ev-seeds-early-aug-2018-density-exp.csv")
dea
dla <- read_csv("./data/focal-disease-seeds-late-aug-2018-density-exp.csv")
dla
ds <- read_csv("./data/ev-disease-seeds-sep-2018-density-exp.csv")
ds
dl <- read_csv("./data/ev-seeds-2018-litter-exp.csv")
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

# add and modify columns in di
di <- di %>%
  mutate(
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
select(di, collect_date, month) %>% unique()
select(di, plant, age, ID, focal) %>% unique() %>% data.frame()
filter(di, focal == 0) %>% select(plot, age) %>% unique()

# add and modify columns in de
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
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0)
  ) %>%
  mutate(collect_date = ifelse(site == "D3" & plot == 7 & treatment == "water" & month == "late August" & ID == "R2" & age == "seedling", 20180828, collect_date)) # correct to match spikelet data

# merge
d <- full_join(di, de)
d

# notes
unique(d$spikelet_notes)

# check some samples
d %>% filter(spikelet_notes == "unknown plant")
d %>% filter(spikelet_notes %in% c("1 OF 2", "2 OF 2")) %>% data.frame()

# add a removal column, sp column, edit notes, and check for duplicates
d <- d %>%
  mutate(
    remove = case_when(spikelet_notes == "unknown plant" ~ 1,
                       TRUE ~ 0), 
    spikelet_notes = ifelse(spikelet_notes %in% c("1 OF 2", "2 OF 2"), "2 samples", spikelet_notes),
    sp = "Ev"
  ) %>%
  group_by(site, plot, treatment, ID, age, month, focal) %>%
  mutate(reps = n())

# assign month to field data
dj$month <- "July"
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
  )

# combine field data
df <- full_join(dj, dea) %>%
  full_join(dla) %>%
  full_join(ds) %>%
  full_join(dl)

# edit field data columns
unique(df$ID)

df <- df %>%
  mutate(
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
dc <- full_join(d, df)

# records in which seeds were collected and recorded
dc %>%
  filter(!is.na(spikelets) & seeds_collected == 1) %>%
  dim() 
# 404

# records in which they weren't recorded
dcc <- dc %>%
  filter(!is.na(spikelets) & (is.na(seeds_collected) | seeds_collected == 0))
# 16
data.frame(dcc)
# D2 10W EvAR1 didn't have data collected in the field
# October plants - didn't write down seed collection
# others = recording mistakes?

# recorded as collected, but missing
dcm <- dc %>%
  filter(is.na(spikelets) & seeds_collected == 1)
# 14
data.frame(dcm)

# look for overlaps between missing and duplicates
dcm %>% select(site, plot, treatment, age, ID, month, spikelet_notes)
d %>% filter(reps > 1) %>% select(site, plot, treatment, age, ID, month, collect_date, spikelet_notes)
# 5 of them may be duplicates - don't include duplicates unless noted

# combine separated samples
group_cols <- d %>% select(site:spikelet_date, spikelet_notes:seed_date, remove:sp) %>% names

d2 <- d %>%
  filter(spikelet_notes == "2 samples") %>%
  group_by(.dots = group_cols) %>%
  summarise(spikelets = sum(spikelets),
            spikelet_weight.mg = sum(spikelet_weight.mg),
            seeds = sum(seeds),
            seed_weight.mg = sum(seed_weight.mg),
            reps = 1
  )  %>%
  full_join(filter(d, spikelet_notes != "2 samples" | is.na(spikelet_notes)))

# remove duplicates
d2 <- d2 %>%
  ungroup() %>%
  mutate(
    remove = ifelse(reps > 1, 1, remove)
  )


# spikelet number
d2 %>%
  ggplot(aes(x = spikelets)) +
  geom_histogram()

d2 %>%
  filter(spikelets > 200)
# adults during late August and September

# spikelet weight
d2 %>%
  ggplot(aes(x = spikelet_weight.mg)) +
  geom_histogram()

d2 %>%
  filter(spikelet_weight.mg > 3500)
# adults in late August

d2 %>%
  ggplot(aes(x = spikelets, y = spikelet_weight.mg, colour = age, shape = month)) +
  geom_point()

# large weight:number
d2 %>% filter(spikelet_weight.mg > 4000)
# Chris will re-weigh these

# seed count
d2 %>%
  ggplot(aes(x = seeds)) +
  geom_histogram()

d2 %>%
  filter(seeds > 200)
# adult in September

# seed weight
d2 %>%
  ggplot(aes(x = seed_weight.mg)) +
  geom_histogram()

d2 %>%
  filter(seed_weight.mg > 600)
# same as above

# seed and spikelet
d2 %>%
  ggplot(aes(x = spikelets, y = seeds)) +
  geom_point(aes(colour = age, shape = treatment)) +
  geom_smooth(method = lm)

# regression fit
m1 <- lm(seeds ~ spikelets, data = d2)
summary(m1) # R2 is 0.87

# seed and spikelet weight
d2 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seed_weight.mg))  +
  geom_point(aes(colour = age, shape = treatment)) +
  geom_smooth(method = lm)

# regression fit
m2 <- lm(seed_weight.mg ~ spikelet_weight.mg, data = d2)
summary(m2) # R2 is 0.91

# seed count and spikelet weight
d2  %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_smooth(method = lm, size = 0.5, colour = "black")
# Chris is processing a few more in the gap

# regression fit
m3 <- lm(seeds ~ spikelet_weight.mg, data = d2)
summary(m3) # R2 is 0.94

# seed count and spikelet weight again with site
d2 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = site, shape = treatment)) +
  geom_smooth(method = lm)

# residuals
d2 %>%
  filter(!is.na(seeds)) %>%
  mutate(resid = as.numeric(resid(m3))) %>%
  ggplot(aes(x = spikelet_weight.mg, y = resid)) +
  geom_point(aes(colour = site, shape = age)) + 
  facet_wrap(~treatment)
# may underestimate water more than fungicide

# seed count and spikelet weight by treatment
d2 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = site, shape = age)) +
  facet_wrap(~treatment) +
  geom_smooth(method = lm)

# regression fit
m4 <- lm(seeds ~ spikelet_weight.mg * treatment, data = d2)
summary(m4)
m5 <- lm(seeds ~ spikelet_weight.mg + treatment, data = d2)
summary(m5)
AIC(m3, m4, m5)
# m3 is the best fit model and adding treatment hardly changes the slope of the relationship. it does decrease the estimate with water, but this effect seems minimal when looking at the figures. also, the simple model underestimates water.

# add ratios and seed estimate
d3 <- d2 %>%
  mutate(
    count_ratio = seeds / spikelets,
    weight_ratio = seed_weight.mg / spikelet_weight.mg, 
    cw_ratio = seeds / spikelet_weight.mg,
    collect_date = as.Date(as.character(collect_date), "%Y%m%d")
  ) %>%
  mutate(
    seeds = ifelse(is.na(seeds), coef(m3)[1] + spikelet_weight.mg * coef(m3)[2], seeds)
  )

# ratio over time
d3 %>%
  ggplot(aes(x = collect_date, y = cw_ratio)) + 
  geom_point(aes(colour = age, shape = treatment), size = 2) + 
  geom_smooth(method = "lm")
# not a strong trend

# ratio over spikelet count
d3 %>%
  ggplot(aes(x = spikelets, y = cw_ratio)) + 
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_smooth(method = "lm")
# hard tell because sample size gets smaller with more spikelets


## Save data
eseeds = d3
