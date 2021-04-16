##### info ####

# file: leaf_scans_data_processing_2019_density_exp
# author: Amy Kendig
# date last edited: 1/8/21
# goal: combine raw 2019 leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: mv_leaf_damage_severity_2019.ijm, ev_leaf_damage_severity_2019.ijm


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

#### need to update files below - should be from Jan 2021 ####

# import all raw data files
ls_ev_may <- read_csv("data/ev_leaf_scans_may_2019_density_exp.csv") # need to re-run
ls_mv_may <- read_csv("data/mv_leaf_scans_may_2019_density_exp.csv") # need to re-run
ls_ev_jun <- read_csv("data/ev_leaf_scans_jun_2019_density_exp.csv")
ls_mv_jun <- read_csv("data/mv_leaf_scans_jun_2019_density_exp.csv")
ls_ev_jul <- read_csv("data/ev_leaf_scans_jul_2019_density_exp.csv")
ls_mv_jul <- read_csv("data/mv_leaf_scans_jul_2019_density_exp.csv")
ls_ev_early_aug <- read_csv("data/ev_leaf_scans_early_aug_2019_density_exp.csv")
ls_mv_early_aug <- read_csv("data/mv_leaf_scans_early_aug_2019_density_exp.csv")
ls_ev_late_aug <- read_csv("data/ev_leaf_scans_late_aug_2019_density_exp.csv")
ls_mv_late_aug <- read_csv("data/mv_leaf_scans_late_aug_2019_density_exp.csv")
ls_ev_sep <- read_csv("data/ev_leaf_scans_sep_2019_density_exp.csv")
ls_mv_sep <- read_csv("data/mv_leaf_scans_sep_2019_density_exp.csv")

# import data files to check collection
dt_may <- read_csv("data/ev_disease_may_2019_density_exp.csv")
dt_jun <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
dt_jul <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
dt_early_aug <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
dt_late_aug <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")
dt_sep <- read_csv("data/focal_disease_sep_2019_density_exp.csv")


#### edit data ####

# function for adding columns
col_fun <- function(dat){
  
  dat2 <- dat %>%
    mutate(plant = gsub(".*:","",Slice),
           part = gsub(":.*$", "", Slice) %>% 
             recode(lesions = "lesion"),
           site = substr(plant, 1, 2),
           plot = substr(plant, 4, 5) %>% 
             gsub("[^[:digit:]]", "", .) %>% 
             as.numeric(),
           treatment = substr(plant, 5, 6) %>% 
             gsub("[^[:alpha:]]", "", .) %>% 
             as.factor() %>%
             recode("F" = "fungicide", "W" = "water"),
           ID = case_when(sp == "Ev" ~ substr(plant, 10, 11) %>% 
                            gsub(c("_"), "", .),
                          sp == "Mv" ~ ifelse(grepl("Edge", Slice, fixed = T),
                                              "Edge",
                                              substr(plant, 10, 11) %>% 
                                                gsub(c("_"), "", .))),
           focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
           age = case_when(sp == "Ev" ~ ifelse(ID == "A" | (focal == 0 & plot > 7),
                                               "adult", 
                                               "seedling"),
                           sp == "Mv" ~ "seedling")
           ) %>%
    rename(area = "Total Area")
  
  return(dat2)
  
}

# add columns

# adjust for inconsistent naming before running function
# remove scan that has no field info or uncropped image
ls_ev_may2 <- col_fun(mutate(ls_ev_may, month = "may", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D3_9W_EvR4_May" ~ "leaf:D3_9W_Ev_R4_May",
                                               Slice == "lesions:D3_9W_EvR4_May" ~ "lesions:D3_9W_Ev_R4_May",
                                               Slice == "leaf:D2_9E_Ev_R2_May" ~ "leaf:D2_9W_Ev_R2_May",
                                               Slice == "lesions:D2_9E_Ev_R2_May" ~ "lesions:D2_9W_Ev_R2_May",
                                               Slice == "leaf:D2_9E_Ev_R3_May" ~ "leaf:D2_9W_Ev_R3_May",
                                               Slice == "lesions:D2_9E_Ev_R3_May" ~ "lesions:D2_9W_Ev_R3_May",
                                               Slice == "leaf:D2_4F_Ev_A_May" ~ "leaf:D2_5F_Ev_A_May",
                                               Slice == "lesions:D2_4F_Ev_A_May" ~ "lesions:D2_5F_Ev_A_May",
                                               TRUE ~ Slice))) %>%
  mutate(ID = case_when(focal == 0 ~ substr(plant, 10, 12) %>%
                          gsub(c("_"), "", .),
                        TRUE ~ ID)) %>%
  filter(!plant %in% c("D1_9W_Ev_R8_May", "D2_8F_Ev_A_May"))

ls_mv_may2 <- col_fun(mutate(ls_mv_may, month = "may", sp = "Mv"))

# adjust for inconsistent naming before running function
ls_ev_jun2 <- col_fun(mutate(ls_ev_jun, month = "jun", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D1_8W_Ev_a_June" ~ "leaf:D1_8W_Ev_A_June",
                                               Slice == "lesions:D1_8W_Ev_a_June" ~ "lesions:D1_8W_Ev_A_June",
                                               Slice == "leaf:D1_9W_Ev_a_June" ~ "leaf:D1_9W_Ev_A_June",
                                               Slice == "lesions:D1_9W_Ev_a_June" ~ "lesions:D1_9W_Ev_A_June",
                                               TRUE ~ Slice)))

# adjust for inconsistent naming before running function
ls_mv_jun2 <- col_fun(mutate(ls_mv_jun, month = "jun", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D3_W_Mv_1_June" ~ "leaf:D3_8W_Mv_1_June",
                                               Slice == "lesions:D3_W_Mv_1_June" ~ "lesions:D3_8W_Mv_1_June",
                                               Slice == "leaf:D1_3F_edge_Mv_June" ~ "leaf:D1_3F_Edge_Mv_June",
                                               Slice == "lesions:D1_3F_edge_Mv_June" ~ "lesions:D1_3F_Edge_Mv_June",
                                               Slice == "leaf:d3_2W_Mv_2_june" ~ "leaf:D3_2W_Mv_2_june",
                                               Slice == "lesions:d3_2W_Mv_2_june" ~ "lesions:D3_2W_Mv_2_june",
                                               TRUE ~ Slice)))

# rename mis-labelled photos before running function
ls_ev_jul2 <- col_fun(mutate(ls_ev_jul, month = "jul", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D4_1W_Ev_A_July" ~ "leaf:D4_1W_Ev_3_July",
                                               Slice == "lesions:D4_1W_Ev_A_July" ~ "lesions:D4_1W_Ev_3_July",
                                               Slice == "leaf:D4_4W_Ev_A_July" ~ "leaf:D4_4W_Ev_1_July",
                                               Slice == "lesions:D4_4W_Ev_A_July" ~ "lesions:D4_4W_Ev_1_July",
                                               TRUE ~ Slice))) %>%
  filter(plant != "D2_6W_Ev_July") # scans were unexpected and missing from the same plot. The ID's may be mixed up within the plot, uncropped image that was left in folder 

# adjust for inconsistent naming before running function
ls_mv_jul2 <- col_fun(mutate(ls_mv_jul, month = "jul", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D3_9w_Mv_3_July" ~ "leaf:D3_9W_Mv_3_July",
                                               Slice == "lesions:D3_9w_Mv_3_July" ~ "lesions:D3_9W_Mv_3_July",
                                               TRUE ~ Slice))) %>%
  filter(plant != "D1_10W_Mv_July") # uncropped image that was left in folder


ls_ev_early_aug2 <- col_fun(mutate(ls_ev_early_aug, month = "early_aug", sp = "Ev"))

ls_mv_early_aug2 <- col_fun(mutate(ls_mv_early_aug, month = "early_aug", sp = "Mv"))

# rename mis-labelled photos before running function
ls_ev_late_aug2 <- col_fun(mutate(ls_ev_late_aug, month = "late_aug", sp = "Ev",
                                  Slice = case_when(Slice == "leaf:D1_1W_Ev_4_LateAug" ~ "leaf:D1_1W_Ev_A_LateAug",
                                                    Slice == "lesions:D1_1W_Ev_4_LateAug" ~ "lesions:D1_1W_Ev_A_LateAug",
                                                    Slice == "leaf:D1_4W_Ev_4_LateAug" ~ "leaf:D1_4W_Ev_A_LateAug",
                                                    Slice == "lesions:D1_4W_Ev_4_LateAug" ~ "lesions:D1_4W_Ev_A_LateAug",
                                                    Slice == "leaf:D3_2W_Ev_4_LateAug" ~ "leaf:D3_2W_Ev_A_LateAug",
                                                    Slice == "lesions:D3_2W_Ev_4_LateAug" ~ "lesions:D3_2W_Ev_A_LateAug",
                                                    Slice == "leaf:D3_1F_Ev_3_LatAug" ~ "leaf:D3_1F_Ev_A_LatAug",
                                                    Slice == "lesions:D3_1F_Ev_3_LatAug" ~ "lesions:D3_1F_Ev_A_LatAug",
                                                    Slice == "leaf:D3_1F_Ev_2_LatAug" ~ "leaf:D3_1F_Ev_3_LatAug",
                                                    Slice == "lesions:D3_1F_Ev_2_LatAug" ~ "lesions:D3_1F_Ev_3_LatAug",
                                                    Slice == "leaf:D4_2W_Ev_3_LatAug" ~ "leaf:D4_2W_Ev_A_LatAug",
                                                    Slice == "lesions:D4_2W_Ev_3_LatAug" ~ "lesions:D4_2W_Ev_A_LatAug",
                                                    Slice == "leaf:D4_2W_Ev_2_LatAug" ~ "leaf:D4_2W_Ev_3_LatAug",
                                                    Slice == "lesions:D4_2W_Ev_2_LatAug" ~ "lesions:D4_2W_Ev_3_LatAug",
                                                    TRUE ~ Slice))) %>%
  filter(!(plant %in% c("D4_2F_Ev_1_LateAug", "D3_7W_Ev_3_LateAug")))  # scans were unexpected and missing from the same plot. The ID's may be mixed up within the plot.

# remove scans that were not a leaf or mis-labelled and couldn't be reconciled
ls_mv_late_aug2 <- col_fun(mutate(ls_mv_late_aug, month = "late_aug", sp = "Mv"))

ls_ev_sep2 <- col_fun(mutate(ls_ev_sep, month = "sep", sp = "Ev")) %>%
  mutate(ID = case_when(ID == "S" ~ substr(plant, 9, 10) %>% 
                          gsub(c("_"), "", .),
                        TRUE ~ ID),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = ifelse(ID == "A" | (focal == 0 & plot > 7),
                      "adult", 
                      "seedling"))

ls_mv_sep2 <- col_fun(mutate(ls_mv_sep, month = "sep", sp = "Mv")) %>%
  mutate(ID = case_when(ID == "S" ~ substr(plant, 9, 10) %>% 
                          gsub(c("_"), "", .),
                        TRUE ~ ID),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0))

# combine with field data
dt_may2 <- full_join(dt_may, full_join(ls_ev_may2, ls_mv_may2)) %>%
  filter(ID != "Edge")
dt_jun2 <- full_join(dt_jun, full_join(ls_ev_jun2, ls_mv_jun2)) %>%
  filter(ID != "Edge")
dt_jul2 <- full_join(dt_jul, full_join(ls_ev_jul2, ls_mv_jul2)) %>%
  filter(ID != "Edge")
dt_early_aug2 <- full_join(dt_early_aug, full_join(ls_ev_early_aug2, ls_mv_early_aug2)) %>%
  filter(ID != "Edge")
dt_late_aug2 <- full_join(dt_late_aug, full_join(ls_ev_late_aug2, ls_mv_late_aug2)) %>%
  filter(ID != "Edge")
dt_sep2 <- full_join(dt_sep, full_join(ls_ev_sep2, ls_mv_sep2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_may2, is.na(scan)) %>% select(field_notes) %>% unique() # missing from field data
filter(dt_may2, is.na(scan) & !is.na(area)) # D1 9W Ev R8 doesn't have field info, didn't see uncropped scan - removed. Uncropped scan for D2 4F Ev A has the plot cut off, but I'm not sure what other plot it could be. This leaf may have come from D2 5F because the two scans were done together - switched. D2 8F Ev A is R2, this leaf was labelled twice - removed.
filter(dt_may2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_may2, scan == 1 & is.na(area)) # scan missing for D2 5F Ev A (uncropped scan indicates that bag was empty) - replaced with D2 4F

filter(dt_jun2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_jun2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_jun2, scan == 1 & is.na(area)) # scans missing for D4 6W Mv3 and D3 9F Ev1 (no uncropped scans available and I double checked OneDrive for these)

filter(dt_jul2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_jul2, is.na(scan) & is.na(field_notes)) %>% select(plant) %>% unique()
filter(dt_jul2, scan == 0 & !is.na(area)) # unexpected scans: D3 8F Ev3 (no uncropped scans available to check, but this plant had 7 leaves, so it's possible the scan = 0 is a mistake) and D4 9F Ev1 (had 7 leaves, scan = 0 may be mistake)
filter(dt_jul2, scan == 1 & is.na(area)) # missing scans: D3 4W Mv2, D1 6F EvA, D1 7F Ev1 (no uncropped scans available to check)

filter(dt_early_aug2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_early_aug2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_early_aug2, scan == 1 & is.na(area)) # missing scan: D1 5W Mv3 (checked uncropped scan and only two Mv leaves received from D1 5W)

filter(dt_late_aug2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_late_aug2, scan == 0 & !is.na(area)) # unexpected scans: D2 6W Ev3 (seems correct from uncropped scan and other information)
filter(dt_late_aug2, scan == 1 & is.na(area)) # missing scans: 
# D1 2W Ev1
# D1 3W Ev3
# D1 5W Ev1
# D1 6W Ev2
# D1 7W Ev3
# D1 9W EvA (missing the uncropped for all of the D1 W Ev scans)
# D3 7W Ev2 (not in uncropped scan, only Ev 1)

filter(dt_sep2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_sep2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_sep2, scan == 1 & is.na(area)) %>% select(site, plot, treatment, sp) %>% unique() %>% data.frame() # missing scans: the D1 and D2 scans were in the box that got lost in the mail. I'm not sure if the D3 and D4 scans were also lost in the mail or if they were not uploaded properly to OneDrive

# combine edge data
dt_edge <- full_join(ls_mv_may2, ls_mv_jun2) %>%
  full_join(ls_mv_jul2) %>%
  full_join(ls_mv_early_aug2) %>%
  full_join(ls_mv_late_aug2) %>%
  full_join(ls_mv_sep2) %>%
  filter(ID == "Edge")

# combine data
# remove unrelated variables
dat <- full_join(dt_may2, dt_jun2) %>%
  full_join(dt_jul2) %>%
  full_join(dt_early_aug2) %>%
  full_join(dt_late_aug2) %>%
  full_join(dt_sep2) %>%
  select(-c(microbiome, seeds_collected, path_ID, flower_bags))
# row number matches summed row numbers

# spread by part
datw <- dat %>%
  filter(!is.na(part)) %>% # all are missing area values
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(date:field_notes, month:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area)

datw_edge <- dt_edge %>%
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(month:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area)

# check that all plots were collected for Mv
datw_edge %>%
  group_by(month, site, treatment) %>%
  summarise(n = length(unique(plant))) %>%
  filter(n != 10)
data.frame() 
# D1 6F early Aug uncropped scan missing 
# Sep scans lost in mail: D1 8F, 9F and 10F, D1 9W and 10W, D4 9F and 10F, D4 9W and 10W


#### check values ####

# text values
unique(datw$site)
unique(datw$plot)
unique(datw$treatment)
unique(datw$sp)
unique(datw$ID)

# leaf area
datw %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# manually check extremes
datw %>%
  filter(sp == "Ev" & leaf_area.pix > 8e5) %>% data.frame()

# percent lesion area 
datw %>%
  ggplot(aes(x = lesion_area.pix/leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")
# a lot of 1's

# leaf area edge
datw_edge %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# percent lesion area edge
datw_edge %>%
  ggplot(aes(x = lesion_area.pix/leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")


#### fungicide effects ####

ggplot(datw, aes(x = treatment, y = lesion_area.pix/leaf_area.pix)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", size = 2, fun.data = "mean_cl_boot") +
  facet_grid(sp ~ month)


#### save data ####

write_csv(datw, "./intermediate-data/all_leaf_scans_2019_density_exp.csv")
write_csv(datw_edge, "./intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
