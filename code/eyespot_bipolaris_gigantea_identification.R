##### info ####

# file: eyespot_bipolaris_gigantea_identification
# author: Amy Kendig
# date last edited: 11/24/20
# goal: proportion of leaves with eyespots and positive infection


#### set-up ####

# late August data cleaning
source("code/fungal_isolation_data_processing_late_aug_2018_density_exp.R")
# clears existing data
# loads tidyverse

# remove extra values
aug18 <- samps2
rm(list = setdiff(ls(), "aug18"))

# import data
jul18 <- read_csv("data/fungal_isolation_jul_2018_density_exp.csv")
sep18 <- read_csv("data/fungal_isolation_sep_2018_density_exp.csv")
lit19 <- read_csv("data/mv_bipolaris_id_2019_litter_exp.csv")
jul19 <- read_csv("data/fungal_isolation_jul_2019_density_exp.csv")
# mv_aug19 <- read_csv("data/mv_fungal_isolation_early_aug_2019_density_exp.csv") # no eyespot information
# missing two sites (asked Brett)
ev_aug19 <- read_csv("data/ev_fungal_isolation_early_aug_2019_density_exp.csv")
# missing half of D2 and D3 (ask Ashish)
# missing late August (ask Brett and Ashish)


#### edit data ####

# eyespot info
unique(jul18$symptoms)
unique(aug18$symptoms)
unique(sep18$symptoms)

# gigantea info
jul18 %>%
  select(observation, gigantea) %>%
  unique() %>%
  data.frame()

aug18 %>%
  select(observation, gigantea) %>%
  unique() %>%
  data.frame()

sep18 %>%
  select(observation, gigantea) %>%
  unique() %>%
  data.frame()

# remove missing symptoms (no symptoms or no leaves)
# create eyespots column (don't count negative statements)
dat18 <- jul18 %>%
  filter(!is.na(symptoms)) %>%
  mutate(eyespots = as.numeric(str_detect(symptoms, "eye"))) %>%
  select(month, site, plot, treatment, sp, eyespots, gigantea) %>%
  full_join(aug18 %>%
              filter(!is.na(symptoms)) %>%
              mutate(eyespots = as.numeric(str_detect(symptoms, "eye"))) %>%
              select(month, site, plot, treatment, sp, eyespots, gigantea)) %>%
  full_join(sep18 %>%
              filter(!is.na(symptoms)) %>%
              mutate(eyespots = as.numeric(str_detect(symptoms, "eye")),
                     gigantea = case_when(observation == "conidiophores" ~ "Yes",
                                          observation == "no conidiophores" ~ "No",
                                          TRUE ~ gigantea)) %>%
              select(month, site, plot, treatment, sp, eyespots, gigantea)) %>%
  mutate(experiment = "density",
         year = 2018)

# remove samples with zero leaves
# create eyespots column based on leaves_with_eyespots
dat19 <- lit19 %>%
  filter(leaves_tot > 0) %>%
  mutate(eyespots = as.numeric(leaves_eyespots > 0),
         gigantea = ifelse(leaves_bip > 0, "Yes", "No"),
         experiment = "litter",
         sp = "Mv",
         plot = block) %>%
  select(month, site, plot, treatment, sp, eyespots, gigantea, experiment) %>%
  full_join(jul19 %>%
              filter(!is.na(symptoms)) %>%
              mutate(eyespots = as.numeric(leaves_with_eyespots > 0),
                     gigantea = ifelse(leaves_with_signs > 0, "Yes", "No"),
                     experiment = "density") %>%
              select(month, site, plot, treatment, sp, eyespots, gigantea, experiment)) %>%
  full_join(ev_aug19 %>%
              filter(leaves > 0) %>%
              mutate(eyespots = as.numeric(leaves_with_eyespots > 0),
                     gigantea = ifelse(leaves_with_signs > 0, "Yes", "No"),
                     experiment = "density") %>%
              select(month, site, plot, treatment, sp, eyespots, gigantea, experiment)) %>%
    mutate(year = 2019)

# combine data
dat <- full_join(dat18, dat19) %>%
  filter(!is.na(gigantea)) %>%
  mutate(giganteaN = ifelse(gigantea == "Yes", 1, 0))


#### analysis ####

# Aug 2018 missing data on gigantea: could mean it was there and not recorded or was not there

# proportion with eyespots that had gigantea
sum(dat$giganteaN)/sum(dat$eyespots) # 66%

# analyze by species
dat %>%
  group_by(sp, eyespots) %>%
  summarise(gig = sum(giganteaN),
            n = n()) %>%
  mutate(prop = gig/n)
