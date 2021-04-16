##### info ####

# file: leaf_scans_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 4/13/21
# goal: combine raw 2018 leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: mv_leaf_damage_severity_2018.ijm, ev_leaf_damage_severity_2018.ijm


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import all raw data files
ls_ev_jul <- read_csv("data/ev_leaf_scans_jul_2018_density_exp.csv")
ls_mv_jul <- read_csv("data/mv_leaf_scans_jul_2018_density_exp.csv")
ls_ev_late_aug <- read_csv("data/ev_leaf_scans_late_aug_2018_density_exp.csv")
ls_mv_late_aug <- read_csv("data/mv_leaf_scans_late_aug_2018_density_exp.csv")
ls_ev_sep <- read_csv("data/ev_leaf_scans_sep_2018_density_exp.csv")
ls_mv_sep <- read_csv("data/mv_leaf_scans_sep_2018_density_exp.csv")

# import data files to check collection
dt_jul <- read_csv("data/focal_size_disease_jul_2018_density_exp.csv")
dt_late_aug <- read_csv("data/all_disease_seeds_late_aug_2018_density_exp.csv")
dt_sep_ev <- read_csv("data/ev_disease_seeds_sep_2018_density_exp.csv")
dt_sep_mv <- read_csv("data/mv_disease_sep_2018_density_exp.csv")

#### need to analyze missing scans ####
#### add code for separating transect samples from Sep 2018 and separating Bp samples ####

# not sure how I handled edited - figure that out

# list of file to check scans
fl_ev_jul <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-jul-2018-density-exp/scans/ev"))
fl_mv_jul <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-jul-2018-density-exp/scans/mv"))
fl_ev_aug <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-aug-2018-density-exp/scans/ev"))
fl_mv_aug <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-aug-2018-density-exp/scans/mv"))
fl_ev_sep <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-sep-2018-density-exp/scans/ev"))
fl_mv_sep <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-sep-2018-density-exp/scans/mv"))


#### check for missing scans ####

# function
scan_fun <- function(scans, files) {
  
  files2 <- files %>%
    mutate(plant = str_replace(file, ".tiff", ""))
  
  scans2 <- scans %>%
    transmute(plant = gsub(".*:","",Slice))
  
  mis_files <- files2 %>%
    anti_join(scans2)
  
  return(mis_files)
  
}

# apply function
scan_fun(ls_ev_jul, fl_ev_jul)


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
           treatment = substr(plant, 4, 6) %>% 
             gsub("[^[:alpha:]]", "", .) %>% 
             as.factor() %>%
             recode("F" = "fungicide", "W" = "water"),
           ID = case_when(sp == "Ev" ~ substr(plant, 9, 13) %>%
                            gsub(c("E|v|A|_"), "", .),
                          sp == "Mv" ~ substr(plant, 9, 11) %>% 
                                                gsub(c("M|v|_"), "", .)),
           ID = case_when(sp == "Ev" & ID == "" ~ "A",
                          TRUE ~ ID),
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

# July
ls_ev_jul2 <- col_fun(mutate(ls_ev_jul, month = "jul", sp = "Ev")) %>%
  mutate(treatment = recode(treatment, "FE" = "fungicide")) %>%
  filter(ID != "Bp") %>% # remove extra leaf collected to exemplify Bipolaris infection
  full_join(dt_jul %>%
              filter(sp == "Ev"))

unique(ls_ev_jul2$part)
unique(ls_ev_jul2$site)
unique(ls_ev_jul2$plot)
unique(ls_ev_jul2$treatment)
# filter(ls_ev_jul2, treatment == "FE")
unique(ls_ev_jul2$ID)
# filter(ls_ev_jul2, ID == "") %>%
#   select(plant) %>%
#   unique() %>%
#   data.frame()
# filter(ls_ev_jul2, ID == "Bp")

ls_mv_jul2 <- col_fun(mutate(ls_mv_jul, month = "jul", sp = "Mv")) %>%
  mutate(ID = case_when(ID == "B" ~ "Bg",
                        TRUE ~ ID)) %>%
  full_join(dt_jul %>%
              filter(sp == "Mv"))

unique(ls_mv_jul2$part)
unique(ls_mv_jul2$site)
unique(ls_mv_jul2$plot)
unique(ls_mv_jul2$treatment)
unique(ls_mv_jul2$ID)
# filter(ls_mv_jul2, ID == "B")


# Late August
# D2 1W EvA (0 leaves infected) should probably be D2 1F EvA (5 leaves infected)
ls_ev_late_aug2 <- col_fun(ls_ev_late_aug %>%
                             mutate(month = "late_aug",
                                    sp = "Ev",
                                    Slice = case_when(Slice == "leaf:D2_1W_EvA" ~ "leaf:D2_1F_EvA",
                                                      Slice == "lesions:D2_1W_EvA" ~ "lesions:D2_1F_EvA",
                                                      TRUE ~ Slice))) %>%
  full_join(dt_late_aug %>%
              filter(sp == "Ev") %>%
              mutate(ID = case_when(ID != "A" & substr(ID, 1, 1) == "A" ~ gsub("A", "", ID),
                                    TRUE ~ ID)))

unique(ls_ev_late_aug2$part)
unique(ls_ev_late_aug2$site)
unique(ls_ev_late_aug2$plot)
unique(ls_ev_late_aug2$treatment)
unique(ls_ev_late_aug2$ID)

# D2 10F Mv2 (0 leaves infected) should probably be D2 10F Mv3 (1 leaf infected)
ls_mv_late_aug2 <- col_fun(mutate(ls_mv_late_aug, month = "late_aug", 
                                  sp = "Mv",
                                  Slice = case_when(Slice == "leaf:D2_10F_Mv2" ~ "leaf:D2_10F_Mv3",
                                                    Slice == "lesions:D2_10F_Mv2" ~ "lesions:D2_10F_Mv3",
                                                    TRUE ~ Slice))) %>%
  mutate(ID = case_when(ID %in% c("", "g", "v") ~ "Bg",
                        TRUE ~ ID)) %>%
  full_join(dt_late_aug %>%
              filter(sp == "Mv"))

unique(ls_mv_late_aug2$part)
unique(ls_mv_late_aug2$site)
unique(ls_mv_late_aug2$plot)
unique(ls_mv_late_aug2$treatment)
unique(ls_mv_late_aug2$ID)
# filter(ls_mv_late_aug2, ID %in% c("", "g", "v")) %>% select(plant, ID) %>% data.frame()


# September
ls_ev_sep2 <- col_fun(mutate(ls_ev_sep, month = "sep", sp = "Ev")) %>%
  mutate(ID = case_when(ID == "." ~ "A",
                        TRUE ~ ID)) %>%
  full_join(dt_sep_ev %>%
              mutate(ID = case_when(ID != "A" & substr(ID, 1, 1) == "A" ~ gsub("A", "", ID),
                                    TRUE ~ ID)))

unique(ls_ev_sep2$part)
unique(ls_ev_sep2$site)
unique(ls_ev_sep2$plot)
unique(ls_ev_sep2$treatment)
unique(ls_ev_sep2$ID)
# filter(ls_ev_sep2, ID == ".") # looked at file name

ls_mv_sep2 <- col_fun(mutate(ls_mv_sep, month = "sep", sp = "Mv")) %>%
  mutate(treatment = recode(treatment, "w" = "water"),
         ID = case_when(ID == "B" ~ "Bg",
                        TRUE ~ ID)) %>%
  full_join(dt_sep_mv %>%
              filter(ID %in% c("1", "2", "3")))

unique(ls_mv_sep2$part)
unique(ls_mv_sep2$site)
unique(ls_mv_sep2$plot)
unique(ls_mv_sep2$treatment)
# filter(ls_mv_sep2, treatment == "w")
unique(ls_mv_sep2$ID)
# filter(ls_mv_sep2, ID  == "B") %>% select(plant, ID) %>% data.frame()


#### check against field data ####

# missing scans for infected plants
filter(ls_ev_jul2, is.na(area) & !is.na(leaves_infec)) %>% filter(leaves_infec > 0)
# D2 10F Ev2
# D3 10W Ev1
# D3 10W Ev3
# D4 3F Ev2
# D4 9F Ev1
# D4 4W Ev1
# unexpected scans
filter(ls_ev_jul2, !is.na(area) & (is.na(leaves_infec) | leaves_infec == 0) & focal == 1) %>% data.frame()
# none of them seem to be mislabels of above

# missing scans for infected plants
filter(ls_mv_jul2, is.na(area) & !is.na(leaves_infec)) %>% filter(leaves_infec > 0)
# D3 6W Mv1
# D3 6W Mv2
# D4 2F Mv2
# D4 6F Mv2
# D4 7F Mv1
# D4 1W Mv3
# D4 4W Mv3
# unexpected scans
filter(ls_mv_jul2, !is.na(area) & (is.na(leaves_infec) | leaves_infec == 0) & focal == 1) %>% select(plant, leaves_tot, leaves_infec) %>% unique() %>% data.frame()
# too many options to decide if any are mis-labelled
filter(ls_mv_jul2, plant == "D4_9F_Mv_3") # this plant died, so it's severity will be negated by NA leaves, but I'm not sure if it's mis-labelled

# missing scans for infected plants
filter(ls_ev_late_aug2, is.na(area) & !is.na(leaves_infec)) %>% filter(leaves_infec > 0)
# unexpected scans
filter(ls_ev_late_aug2, !is.na(area) & (is.na(leaves_infec) | leaves_infec == 0) & focal == 1) %>% data.frame() %>% select(plant, leaves_tot, leaves_infec) %>% unique()

# missing scans for infected plants
filter(ls_mv_late_aug2, is.na(area) & !is.na(leaves_infec)) %>% filter(leaves_infec > 0)
# unexpected scans
filter(ls_mv_late_aug2, !is.na(area) & (is.na(leaves_infec) | leaves_infec == 0) & focal == 1) %>% data.frame() %>% select(plant, leaves_tot, leaves_infec) %>% unique()

# missing scans for infected plants
filter(ls_ev_sep2, is.na(area) & !is.na(leaves_infec)) %>% filter(leaves_infec > 0)
# unexpected scans
filter(ls_ev_sep2, !is.na(area) & (is.na(leaves_infec) | leaves_infec == 0) & focal == 1) %>% data.frame() %>% select(plant, leaves_tot, leaves_infec) %>% unique()

# missing scans for infected plants
filter(ls_mv_sep2, is.na(area) & !is.na(leaves_infec)) %>% filter(leaves_infec > 0)
# unexpected scans
filter(ls_mv_sep2, !is.na(area) & (is.na(leaves_infec) | leaves_infec == 0) & focal == 1) %>% data.frame() %>% select(plant, leaves_tot, leaves_infec) %>% unique()


#### combine data ####

# combine background data
dt_mv_bg <- full_join(ls_mv_jul2, ls_mv_late_aug2) %>%
  full_join(ls_mv_sep2) %>%
  filter(ID == "Bg") %>%
  select(-c(date:seeds_collected))

dt_ev_bg <- full_join(ls_ev_jul2, ls_ev_late_aug2) %>%
  full_join(ls_ev_sep2) %>%
  filter(focal == 0) %>%
  select(-c(date:infec, seeds_observed))

# combine individual plant data
dat <- full_join(ls_ev_jul2, ls_mv_jul2) %>%
  full_join(ls_ev_late_aug2) %>%
  full_join(ls_mv_late_aug2) %>%
  full_join(ls_ev_sep2) %>%
  full_join(ls_mv_sep2) %>%
  filter(focal == 1) %>%
  select(-c(date:infec, seeds_observed))

# spread by part
datw <- dat %>%
  filter(!is.na(area)) %>%
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(month:seeds_collected)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area)

datw_mv_bg <- dt_mv_bg %>%
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(month:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area)

datw_ev_bg <- dt_ev_bg %>%
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(month:seeds_collected)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area)


#### check values ####

# leaf counts
filter(datw, leaf_count > 2)
filter(datw_ev_bg, leaf_count > 2)
filter(datw_mv_bg, leaf_count > 12)

# leaf area
datw %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# manually check extremes
datw %>%
  filter(sp == "Mv" & leaf_area.pix > 7e5) %>% data.frame()

# percent lesion area 
datw %>%
  ggplot(aes(x = lesion_area.pix/leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")


#### save data ####

write_csv(datw, "./intermediate-data/focal_leaf_scans_2018_density_exp.csv")
write_csv(datw_ev_bg, "./intermediate-data/ev_background_leaf_scans_2018_density_exp.csv")
write_csv(datw_mv_bg, "./intermediate-data/mv_background_leaf_scans_2018_density_exp.csv")
