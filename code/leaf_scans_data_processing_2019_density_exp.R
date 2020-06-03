##### info ####

# file: leaf_scans_data_processing_2019_density_exp
# author: Amy Kendig
# date last edited: 6/2/20
# goal: combine raw 2019 leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: mv_leaf_damage_severity.ijm


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import all raw data files
ls_ev_may <- read_csv("data/ev_leaf_scans_may_2019_density_exp.csv")
ls_mv_may <- read_csv("data/mv_leaf_scans_may_2019_density_exp.csv")
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

# # indicate whether images were edited
# # remove edited from name
# # remove month from name
# may <- may %>%
#   mutate(edited = ifelse(grepl("edited", slice, fixed = T), 1, 0),
#          slice = gsub("_edited", "", slice) %>% gsub("_May", "", .) %>% gsub("_may", "", .),
#          month = "May")
# nrow(may) # 834
# sum(may$edited) # 153

# # create list of edited scans
# may_ed <- may %>%
#   filter(edited == 1) %>%
#   select(slice)
# 
# # check that all are contained in unedited list
# may_ed %>% filter(!(slice %in% filter(may, edited == 0)$slice))
# # one sample is not - looked through files and can't figure out why, but it's fine because we have the edited one
# 
# # remove unedited scans
# may <- may %>%
#   filter(edited == 1 | !(slice %in% may_ed$slice))
# nrow(may) # 684
# 834 - 153 + 3

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
# adjust for inconsistent naming
ls_ev_may2 <- col_fun(mutate(ls_ev_may, month = "may", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D3_9W_EvR4_May" ~ "leaf:D3_9W_Ev_R4_May",
                                               Slice == "lesions:D3_9W_EvR4_May" ~ "lesions:D3_9W_Ev_R4_May",
                                               Slice == "leaf:D2_9E_Ev_R2_May" ~ "leaf:D2_9W_Ev_R2_May",
                                               Slice == "lesions:D2_9E_Ev_R2_May" ~ "lesions:D2_9W_Ev_R2_May",
                                               Slice == "leaf:D2_9E_Ev_R3_May" ~ "leaf:D2_9W_Ev_R3_May",
                                               Slice == "lesions:D2_9E_Ev_R3_May" ~ "lesions:D2_9W_Ev_R3_May",
                                               TRUE ~ Slice))) %>%
  mutate(ID = case_when(focal == 0 ~ substr(plant, 10, 12) %>%
                          gsub(c("_"), "", .),
                        TRUE ~ ID))
ls_mv_may2 <- col_fun(mutate(ls_mv_may, month = "may", sp = "Mv"))
ls_ev_jun2 <- col_fun(mutate(ls_ev_jun, month = "jun", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D1_8W_Ev_a_June" ~ "leaf:D1_8W_Ev_A_June",
                                               Slice == "lesions:D1_8W_Ev_a_June" ~ "lesions:D1_8W_Ev_A_June",
                                               Slice == "leaf:D1_9W_Ev_a_June" ~ "leaf:D1_9W_Ev_A_June",
                                               Slice == "lesions:D1_9W_Ev_a_June" ~ "lesions:D1_9W_Ev_A_June",
                                               TRUE ~ Slice)))
ls_mv_jun2 <- col_fun(mutate(ls_mv_jun, month = "jun", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D3_W_Mv_1_June" ~ "leaf:D3_8W_Mv_1_June",
                                               Slice == "lesions:D3_W_Mv_1_June" ~ "lesions:D3_8W_Mv_1_June",
                                               Slice == "leaf:D1_3F_edge_Mv_June" ~ "leaf:D1_3F_Edge_Mv_June",
                                               Slice == "lesions:D1_3F_edge_Mv_June" ~ "lesions:D1_3F_Edge_Mv_June",
                                               Slice == "leaf:d3_2W_Mv_2_june" ~ "leaf:D3_2W_Mv_2_june",
                                               Slice == "lesions:d3_2W_Mv_2_june" ~ "lesions:D3_2W_Mv_2_june",
                                               TRUE ~ Slice)))
ls_ev_jul2 <- col_fun(mutate(ls_ev_jul, month = "jul", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D4_1W_Ev_A_July" ~ "leaf:D4_1W_Ev_3_July",
                                               Slice == "lesions:D4_1W_Ev_A_July" ~ "lesions:D4_1W_Ev_3_July",
                                               Slice == "leaf:D4_4W_Ev_A_July" ~ "leaf:D4_4W_Ev_1_July",
                                               Slice == "lesions:D4_4W_Ev_A_July" ~ "lesions:D4_4W_Ev_1_July",
                                               TRUE ~ Slice)))  # scans were unexpected and missing from the same plot, reassigned names, but the ID's may be mixed up within the plot. re-analyze in fiji for uncropped image
ls_mv_jul2 <- col_fun(mutate(ls_mv_jul, month = "jul", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D3_9w_Mv_3_July" ~ "leaf:D3_9W_Mv_3_July",
                                               Slice == "lesions:D3_9w_Mv_3_July" ~ "lesions:D3_9W_Mv_3_July",
                                               TRUE ~ Slice))) %>%
  filter(plant != "D1_10W_Mv_July") # uncropped image that was left in folder
ls_ev_early_aug2 <- col_fun(mutate(ls_ev_early_aug, month = "early_aug", sp = "Ev"))
ls_mv_early_aug2 <- col_fun(mutate(ls_mv_early_aug, month = "early_aug", sp = "Mv"))
ls_ev_late_aug2 <- col_fun(mutate(ls_ev_late_aug, month = "late_aug", sp = "Ev",
                                  Slice = case_when(Slice == "leaf:D1_1W_Ev_4_LateAug" ~ "leaf:D1_1W_Ev_A_LateAug",
                                                    Slice == "lesions:D1_1W_Ev_4_LateAug" ~ "lesions:D1_1W_Ev_A_LateAug",
                                                    Slice == "leaf:D1_4W_Ev_4_LateAug" ~ "leaf:D1_4W_Ev_A_LateAug",
                                                    Slice == "lesions:D1_4W_Ev_4_LateAug" ~ "lesions:D1_4W_Ev_A_LateAug",
                                                    Slice == "leaf:D3_2W_Ev_4_LateAug" ~ "leaf:D3_2W_Ev_A_LateAug",
                                                    Slice == "lesions:D3_2W_Ev_4_LateAug" ~ "lesions:D3_2W_Ev_A_LateAug",
                                                    TRUE ~ Slice)))
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

# check against field data
filter(dt_may2, is.na(scan)) %>% select(field_notes) %>% unique() # missing from field data
filter(dt_may2, is.na(scan) & is.na(field_notes)) # D1 9W Ev R8 doesn't have field info
filter(dt_may2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_may2, scan == 1 & is.na(area)) # scan missing for D2 5F Ev A

filter(dt_jun2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_jun2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_jun2, scan == 1 & is.na(area)) # scans missing for D4 6W Mv3 and D3 9F Ev1

filter(dt_jul2, is.na(scan)) %>% select(field_notes) %>% unique() # missing from field data
filter(dt_jul2, is.na(scan) & is.na(field_notes)) %>% select(plant) %>% unique()
filter(dt_jul2, scan == 0 & !is.na(area)) # unexpected scans: D3 8F Ev3
filter(dt_jul2, scan == 1 & is.na(area)) # missing scans: D3 4W Mv2, D1 6F EvA, D1 7F Ev1

filter(dt_early_aug2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
filter(dt_early_aug2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_early_aug2, scan == 1 & is.na(area)) # missing scan: D1 5W Mv3. Checked uncropped scand and only two Mv leaves received from D1 5W

filter(dt_late_aug2, is.na(scan)) %>% select(field_notes) %>% unique() # none missing from field data
### start here: checking below for possible manual corrections ####
filter(dt_late_aug2, scan == 0 & !is.na(area)) # unexpected scans: D2 6W Ev3 (seems correct), D3 1F Ev2, D3 7W Ev3, D4 2F Ev1, D4 2W Ev2
filter(dt_late_aug2, scan == 1 & is.na(area)) # missing scans: D1 2W Ev1, D1 3W Ev3, D1 5W Ev1, D1 6W Ev2, D1 7W Ev3, D1 9W EvA, D3 1F EvA, D3 7W Ev2, D3 8W EvA, D4 2W EvA

# spread by part
datw <- dat %>%
  select(-c(slice)) %>%
  gather(variable, value, -(edited:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value)

# indicate values to remove
datw <- datw %>%
  mutate(remove = ifelse(plant %in% c("D3_3F_Ev_A", "D3_6F_Ev_A", "D3_10W_Ev_R8") & month == "May",1, 0))


#### check values ####

# leaf area
datw %>%
  filter(remove == 0) %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# manually check extremes
datw %>%
  filter(sp == "Ev" & leaf_area.pix > 7.5e5) %>% data.frame()

# percent lesion area 
datw %>%
  filter(remove == 0) %>%
  ggplot(aes(x = (lesion_area.pix - green_area.pix)/leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# manually check extremes
datw %>%
  filter(remove == 0) %>%
  filter(sp == "Ev" & (lesion_area.pix - green_area.pix)/leaf_area.pix > 0.6) %>% data.frame()

# check that all plots were collect for Mv
datw %>%
  filter(ID == "Edge") %>%
  group_by(month, site, treatment) %>%
  summarise(n = length(unique(plant)))

# check that all Ev were collected
ev_may %>%
  filter(scan == 1) %>%
  left_join(datw) %>%
  filter(is.na(plant)) %>%
  data.frame()
# missing scans - will ask Laney about these

datw %>%
  filter(sp == "Ev") %>%
  left_join(ev_may) %>%
  filter(is.na(leaves_tot)) %>% 
  data.frame()
# shouldn't have scans - could be the ones above, but mislabelled


#### save data ####

leaf <- datw