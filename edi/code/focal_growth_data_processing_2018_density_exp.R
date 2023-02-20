##### outputs ####

# focal_processed_growth_2018_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lubridate)

# import size data
# start with June because individuals were replaced in May
fjn <- read_csv("data/focal_size_disease_jun_2018_density_exp.csv") 
fjl <- read_csv("data/focal_size_disease_jul_2018_density_exp.csv")


#### edit data ####

# June
fjn
unique(fjn$field_notes)
unique(subset(fjn, field_notes == "dead")$height.cm)
unique(subset(fjn, field_notes == "dead")$tillers)
unique(subset(fjn, field_notes == "dead")$basal_circ.cm)
unique(subset(fjn, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)

# July
fjl
unique(fjl$field_notes)
unique(subset(fjl, field_notes == "dead")$height.cm)
unique(subset(fjl, field_notes == "dead")$tillers)
unique(subset(fjl, field_notes == "dead")$basal_circ.cm)
unique(subset(fjl, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)

# combine data
dat <- fjn %>%
  filter(field_notes != "tree fell on plot" | is.na(field_notes)) %>%
  mutate(month = "jun") %>%
  full_join(fjl %>%
              filter(field_notes != "tree fell on plot" | is.na(field_notes)) %>%
              mutate(month = "jul")) %>%
  mutate(date = as.Date(as.character(date), format = "%Y%m%d"),
         height.cm = case_when(field_notes == "dead" ~ 0,
                               TRUE ~ height.cm),
         tillers = case_when(field_notes == "dead" ~ 0,
                             TRUE ~ tillers),
         basal_circ.cm = case_when(field_notes == "dead" ~ 0,
                                   TRUE ~ basal_circ.cm)) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  mutate(days = max(date) - min(date),
         vals = n()) %>%
  ungroup()

# look at dead plants
dat %>%
  filter(field_notes == "dead") %>%
  select(site:ID, height.cm:basal_circ.cm, month) %>%
  pivot_wider(names_from = month,
              values_from = c(height.cm, tillers, basal_circ.cm),
              names_glue = "{.value}_{month}") %>%
  data.frame()
# a lot have NA values

filter(dat, (is.na(height.cm) | is.na(tillers) | (is.na(basal_circ.cm) & sp == "Ev")) & field_notes != "dead") %>%
  data.frame()
# three plants with measures forgotten

# fix missing data
dat2 <- dat %>%
  mutate(height.cm = replace_na(height.cm, 0),
         tillers = replace_na(tillers, 0),
         basal_circ.cm = ifelse(is.na(basal_circ.cm) & 
                                  (field_notes != "forgot to measure basal circumference" | is.na(field_notes)),
                                0, basal_circ.cm))

# make wide
datw <- dat2 %>%
  pivot_wider(-c(date, infec, leaves_tot, leaves_infec, field_notes, Bp_spots_Ev, seeds_observed, vals),
              names_from = month,
              values_from = c(height.cm, tillers, basal_circ.cm),
              names_glue = "{.value}_{month}") %>%
  mutate(height_growth = (height.cm_jul - height.cm_jun)/height.cm_jun,
         height_growth = ifelse(height.cm_jun == 0, 0, height_growth),
         tiller_growth = (tillers_jul - tillers_jun)/tillers_jun,
         tiller_growth = ifelse(tillers_jun == 0, 0, tiller_growth),
         basal_growth = (basal_circ.cm_jul - basal_circ.cm_jun)/basal_circ.cm_jun,
         basal_growth = ifelse(basal_circ.cm_jun == 0, 0, basal_growth))


#### visualize ####

datw %>%
  ggplot(aes(x = height_growth)) +
  geom_histogram()

filter(datw, height_growth > 7.5) %>% data.frame()
filter(dat, site == "D3" & plot == 5 & treatment == "water" & sp == "Ev" & ID == "A") %>% data.frame()
# this seems like a very unrealistic change in growth
# checked scanned sheets and these are correct
# I suspect the June value was recorded incorrectly

datw %>%
  ggplot(aes(x = tiller_growth)) +
  geom_histogram()

datw %>%
  ggplot(aes(x = basal_growth)) +
  geom_histogram()
filter(datw, basal_growth > 2.5) %>% data.frame()
# different plant than above

# update data
datw2 <- datw %>%
  mutate(height_growth = case_when( site == "D3" & plot == 5 & treatment == "water" & sp == "Ev" & ID == "A" ~ NA_real_,
                                    TRUE ~ height_growth))

# save data
write_csv(datw2, "intermediate-data/focal_processed_growth_2018_density_exp.csv")
