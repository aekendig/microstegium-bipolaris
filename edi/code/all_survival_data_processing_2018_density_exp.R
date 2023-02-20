##### outputs ####

# all_processed_survival_2018_density_exp.csv

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import seeds data
eseeds <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")

# import survival data
fjn <- read_csv("data/focal_size_disease_jun_2018_density_exp.csv") # focal June, height and notes
fjl <- read_csv("data/focal_size_disease_jul_2018_density_exp.csv") # focal July, height and notes
fea <- read_csv("data/focal_status_early_aug_2018_density_exp.csv") # focal early August, status
ala <- read_csv("data/all_disease_seeds_late_aug_2018_density_exp.csv") # all late August, no green and leaves tot
es <- read_csv("data/ev_disease_seeds_sep_2018_density_exp.csv") # all Ev Sep, leaves tot
ms <- read_csv("data/mv_disease_sep_2018_density_exp.csv") # all Mv Sep, leaves tot
ew <- read_csv("data/ev_winter_survival_apr_2019_density_exp.csv") # all Ev over winter


#### edit data ####

# June
unique(fjn$field_notes)
unique(subset(fjn, field_notes == "dead")$height.cm)
unique(subset(fjn, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)
fjn <- fjn %>%
  mutate(survival = ifelse(is.na(height.cm), 0, 1),
         age = case_when(sp == "Ev" & ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1)

# July
unique(fjl$field_notes)
unique(subset(fjl, field_notes == "dead")$height.cm)
unique(subset(fjl, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)
fjl <- fjl %>%
  mutate(survival = ifelse(is.na(height.cm), 0, 1),
         age = case_when(sp == "Ev" & ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1)

# early August
unique(fea$status)
fea <- fea %>%
  mutate(survival = ifelse(status %in% c("dead", "tree fell on plot"), 0, 1),
         age = case_when(sp == "Ev" & ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1)

# late August
unique(ala$no_green)
unique(subset(ala, no_green == 1)$leaves_tot)
unique(subset(ala, is.na(no_green))$leaves_tot)
head(subset(ala, is.na(no_green))) %>% data.frame() # one plant is not Ev (remove from all: D2 7W Ev A), the other lost its tag and is Mv  (D3 2W Mv 1) 
unique(subset(ala, seeds_collected == 1)$no_green) # 1 or 0, will correct with actual seeds
ala <- ala %>%
  mutate(survival = ifelse(no_green == 0, 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# Ev September
unique(subset(es, is.na(leaves_tot))$field_notes) # three are unsure if plant survived
es <- es %>%
  mutate(survival = ifelse(!is.na(leaves_tot), 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# Mv September
unique(subset(ms, is.na(leaves_tot))$field_notes)
ms <- ms %>%
  mutate(survival = ifelse(!is.na(leaves_tot), 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# over winter
ew <- ew %>%
  mutate(focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# select data and add month
fjn <- fjn %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "June")

fjl <- fjl %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "July")

fea <- fea %>%
  select(site, plot, treatment, sp, age, ID, focal, survival) %>%
  mutate(month = "early August")

ala <- ala %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "late August")

es <- es %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "September")

ms <- ms %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "September")

ew <- ew %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "April")

# seeds
filter(eseeds, is.na(seeds) | seeds == 0) # all included survived

# survival by seed production
eseeds2 <- eseeds %>%
  filter(ID_unclear == 0) %>%
  mutate(survival_seeds = 1) %>%
  select(site, plot, treatment, sp, age, ID, focal, survival_seeds, month) %>%
  unique()

# combine all data
d_dens <- full_join(fjn, fjl) %>%
  full_join(fea) %>%
  full_join(ala) %>%
  full_join(es) %>%
  full_join(ms) %>%
  full_join(ew) %>%
  full_join(eseeds2)

# examine problematic data
filter(d_dens, site == "D3" & plot == 2 & treatment == "water" & ID == "1" & sp == "Mv")

# remove compromised data
d_dens <- d_dens %>%
  filter(!(site == "D4" & plot == 8 & treatment == "fungicide") &
           !(site == "D4" & plot == 10 & treatment == "fungicide") &
           !(site == "D2" & plot == 7 & treatment == "water" & ID == "A" & sp == "Ev") &
           !(site == "D3" & plot == 2 & treatment == "water" & ID == "1" & sp == "Mv" & month == "September"))

# change survival values if seeds were produced in that month
filter(d_dens, (is.na(survival) | survival == 0) & survival_seeds == 1)
d_dens <- d_dens %>%
  mutate(survival_seeds = replace_na(survival_seeds, 0),
         survival = ifelse(survival_seeds == 1, 1, survival))

# check for NA's 
sum(is.na(d_dens$survival)) #17

# change survival values if plant was noted as surviving at a later month or dying in an earlier month
# add column for if seeds were ever produced
d_dens2 <- d_dens %>%
  mutate(month = recode(month, "late August" = "lAug", "early August" = "eAug")) %>%
  select(-c(field_notes, survival_seeds)) %>%
  spread(month, survival) %>%
  mutate(October = case_when(April == 1 ~ 1, TRUE ~ October),
         September = case_when(October == 1 ~ 1, TRUE ~ September),
         lAug = case_when(September == 1 ~ 1, TRUE ~ lAug),
         eAug = case_when(lAug == 1 ~ 1, TRUE ~ eAug),
         July = case_when(eAug == 1 ~ 1, TRUE ~ July),
         June = case_when(July == 1 ~ 1, TRUE ~ June),
         July = case_when(is.na(July) & June == 0 ~ 0, TRUE ~ July),
         eAug = case_when(is.na(eAug) & July == 0 ~ 0, TRUE ~ eAug),
         lAug = case_when(is.na(lAug) & eAug == 0 ~ 0, TRUE ~ lAug),
         September = case_when(is.na(September) & lAug == 0 ~ 0, TRUE ~ September),
         October = case_when(is.na(October) & September == 0 ~ 0, TRUE ~ October),
         April = case_when(is.na(April) & October == 0 ~ 0, TRUE ~ April)) %>%
  gather("month", "survival", -c(site, plot, treatment, sp, age, ID, focal)) %>%
  mutate(month = recode(month, "lAug" = "late August", "eAug" = "early August")) %>%
  full_join(select(d_dens, site, plot, treatment, sp, age, ID, focal, month, field_notes)) %>%
  mutate(month = factor(month, levels = c("June", "July", "early August", "late August", "September", "October", "April"))) %>%
  full_join(d_dens %>%
              group_by(site, plot, treatment, sp, age, ID, focal) %>%
              summarise(seeds_produced = as.numeric(sum(survival_seeds, na.rm = T) > 0)))

# check NA values
sum(is.na(d_dens2$survival)) #968
filter(d_dens2, is.na(survival) & !(sp == "Mv" & month %in% c("October", "April"))) %>% data.frame() 
# these are plants that we didn't assess (background) at specific time points or they have notes

# make sure none come back to life
d_dens2 %>% 
  mutate(plant = paste(site, plot, treatment, sp, ID, age, focal, sep = ".")) %>%
  ggplot(aes(x = as.numeric(month), y = survival, colour = plant)) +
  geom_line() +
  theme(legend.position = "none")

# mean survival through summer
d_dens2 %>%
  filter(month == "September" & !is.na(survival)) %>%
  group_by(sp, age) %>%
  summarise(sum(survival) / length(survival))

# survival through summer and winter
d_dens2 %>%
  filter(month == "April" & !is.na(survival)) %>%
  group_by(sp, age) %>%
  summarise(sum(survival) / length(survival))

# survival through winter given summer survival
d_dens2 %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  filter(September == 1 & !is.na(April)) %>%
  group_by(sp, age) %>%
  summarise(sum(April) / length(April))


#### output intermediate data ####

write_csv(d_dens2, "intermediate-data/all_processed_survival_2018_density_exp.csv")
