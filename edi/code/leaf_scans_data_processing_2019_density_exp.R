##### info ####

# all_leaf_scans_2019_density_exp.csv
# mv_edge_leaf_scans_2019_density_exp.csv


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

# list of file to check scans (not archived with data and code, large image files)
# fl_ev_may <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-may-2019-density-exp/scans/ev"))
# fl_mv_may <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-may-2019-density-exp/scans/mv"))
# fl_ev_jun <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-jun-2019-density-exp/scans/ev"))
# fl_mv_jun <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-jun-2019-density-exp/scans/mv"))
# fl_ev_jul <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-jul-2019-density-exp/scans/ev"))
# fl_mv_jul <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-jul-2019-density-exp/scans/mv"))
# fl_ev_early_aug <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-early-aug-2019-density-exp/scans/ev"))
# fl_mv_early_aug <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-early-aug-2019-density-exp/scans/mv"))
# fl_ev_late_aug <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-late-aug-2019-density-exp/scans/ev"))
# fl_mv_late_aug <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-late-aug-2019-density-exp/scans/mv"))
# fl_ev_sep <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-sep-2019-density-exp/scans/ev"))
# fl_mv_sep <- tibble(file = list.files(path = "../leaf-scans/leaf-scans-sep-2019-density-exp/scans/mv"))


#### check for missing scans ####

# function
scan_fun <- function(scans, files) {
  
  files2 <- files %>%
    mutate(plant = str_replace(file, ".tiff|.tif|.png", ""))
  
  scans2 <- scans %>%
    transmute(plant = gsub(".*:","",Slice))
  
  mis_files <- files2 %>%
    anti_join(scans2)
  
  return(mis_files)
  
}

# apply function
# scan_fun(ls_ev_may, fl_ev_may)
# scan_fun(ls_mv_may, fl_mv_may)
# scan_fun(ls_ev_jun, fl_ev_jun)
# scan_fun(ls_mv_jun, fl_mv_jun)
# scan_fun(ls_ev_jul, fl_ev_jul)
# scan_fun(ls_mv_jul, fl_mv_jul)
# scan_fun(ls_ev_early_aug, fl_ev_early_aug)
# scan_fun(ls_mv_early_aug, fl_mv_early_aug)
# scan_fun(ls_ev_late_aug, fl_ev_late_aug)
# scan_fun(ls_mv_late_aug, fl_mv_late_aug)
# scan_fun(ls_ev_sep, fl_ev_sep)
# scan_fun(ls_mv_sep, fl_mv_sep)


#### edit data function ####

col_fun <- function(dat){
  
  dat2 <- dat %>%
    mutate(plant = gsub(".*:","",Slice),
           plant = str_replace(plant, "edited|edit", ""),
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


#### edit May data ####

# adjust for inconsistent naming before running function
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
  filter(!(plant %in% c("D1_9W_Ev_R8_May", "D2_8F_Ev_A_May"))) # remove scan that has no field info and one that was labelled wrong (correctly labelled scan is included)

ls_mv_may2 <- col_fun(mutate(ls_mv_may, month = "may", sp = "Mv"))

# combine
dt_may2 <- full_join(dt_may, full_join(ls_ev_may2, ls_mv_may2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_may2, is.na(scan)) %>% select(Slice, field_notes) %>% unique() # missing from field data
filter(dt_may2, is.na(scan) & !is.na(area)) %>% select(Slice, leaves_tot, field_notes) # D1 9W Ev R8 doesn't have field info, didn't see uncropped scan - removed. Uncropped scan for D2 4F Ev A has the plot cut off, but I'm not sure what other plot it could be. This leaf may have come from D2 5F because the two scans were done together - switched.
filter(dt_may2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_may2, scan == 1 & is.na(area)) # scan missing for D2 5F Ev A (uncropped scan indicates that bag was empty) - replaced with D2 4F


#### edit June data ####

# adjust for inconsistent naming before running function
ls_ev_jun2 <- col_fun(mutate(ls_ev_jun, month = "jun", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D1_8W_Ev_a_June_edited" ~ "leaf:D1_8W_Ev_A_June",
                                               Slice == "lesions:D1_8W_Ev_a_June_edited" ~ "lesions:D1_8W_Ev_A_June",
                                               Slice == "leaf:D1_9W_Ev_a_June" ~ "leaf:D1_9W_Ev_A_June",
                                               Slice == "lesions:D1_9W_Ev_a_June" ~ "lesions:D1_9W_Ev_A_June",
                                               Slice == "leaf:D4_9W_Ev_3orA_June_edited" ~ "leaf:D4_9W_Ev_3_June",
                                               Slice == "lesions:D4_9W_Ev_3orA_June_edited" ~ "lesions:D4_9W_Ev_3_June",
                                               TRUE ~ Slice)))

# adjust for inconsistent naming before running function
ls_mv_jun2 <- col_fun(mutate(ls_mv_jun, month = "jun", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D3_W_Mv_1_June" ~ "leaf:D3_8W_Mv_1_June",
                                               Slice == "lesions:D3_W_Mv_1_June" ~ "lesions:D3_8W_Mv_1_June",
                                               Slice == "leaf:D1_3F_edge_Mv_June" ~ "leaf:D1_3F_Edge_Mv_June",
                                               Slice == "lesions:D1_3F_edge_Mv_June" ~ "lesions:D1_3F_Edge_Mv_June",
                                               TRUE ~ Slice)))

# combine
dt_jun2 <- full_join(dt_jun, full_join(ls_ev_jun2, ls_mv_jun2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_jun2, is.na(scan)) %>% select(Slice, field_notes) %>% unique() # none missing from field data
filter(dt_jun2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_jun2, scan == 1 & is.na(area)) # scans missing for D4 6W Mv3, D3 9F Ev1 (no uncropped scans available and I double checked OneDrive for these). Scans were intentionally removed due to poor quality for D3 6F Mv1. D4 9W EvA is a possible ID of another scan.


#### edit July data ####

# rename mis-labelled photos before running function
ls_ev_jul2 <- col_fun(mutate(ls_ev_jul, month = "jul", sp = "Ev",
                             Slice = case_when(Slice == "leaf:D4_1W_Ev_A_July" ~ "leaf:D4_1W_Ev_3_July",
                                               Slice == "lesions:D4_1W_Ev_A_July" ~ "lesions:D4_1W_Ev_3_July",
                                               Slice == "leaf:D4_4W_Ev_A_July" ~ "leaf:D4_4W_Ev_1_July",
                                               Slice == "lesions:D4_4W_Ev_A_July" ~ "lesions:D4_4W_Ev_1_July",
                                               TRUE ~ Slice)))

# adjust for inconsistent naming before running function
ls_mv_jul2 <- col_fun(mutate(ls_mv_jul, month = "jul", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D3_9w_Mv_3_July" ~ "leaf:D3_9W_Mv_3_July",
                                               Slice == "lesions:D3_9w_Mv_3_July" ~ "lesions:D3_9W_Mv_3_July",
                                               Slice == "leaf:D4_7F_Mv_1or2_July" ~ "leaf:D4_7F_Mv_1_July",
                                               Slice == "lesions:D4_7F_Mv_1or2_July" ~ "lesions:D4_7F_Mv_1_July",
                                               TRUE ~ Slice)))

# combine
dt_jul2 <- full_join(dt_jul, full_join(ls_ev_jul2, ls_mv_jul2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_jul2, is.na(scan)) %>% select(Slice, leaves_tot, field_notes) %>% unique() # none missing from field data
filter(dt_jul2, is.na(scan) & is.na(field_notes)) %>% select(plant) %>% unique()
filter(dt_jul2, scan == 0 & !is.na(area)) # unexpected scans: D3 8F Ev3 (no uncropped scans available to check, but this plant had 7 leaves, so it's possible the scan = 0 is a mistake) and D4 9F Ev1 (had 7 leaves, scan = 0 may be mistake)
filter(dt_jul2, scan == 1 & is.na(area)) # missing scans: D1 6F EvA, D1 7F Ev1 (no uncropped scans available to check), D4 7F Mv2 is a possible ID of another scan. D1 1F Mv3, D1 8F EvA, D4 6W Ev1 all covered in dirt - removed.


#### edit early August data ####

ls_ev_early_aug2 <- col_fun(mutate(ls_ev_early_aug, month = "early_aug", sp = "Ev"))

ls_mv_early_aug2 <- col_fun(mutate(ls_mv_early_aug, month = "early_aug", sp = "Mv"))

# combine
dt_early_aug2 <- full_join(dt_early_aug, full_join(ls_ev_early_aug2, ls_mv_early_aug2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_early_aug2, is.na(scan)) %>% select(Slice, leaves_tot, field_notes) %>% unique() # none missing from field data
filter(dt_early_aug2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_early_aug2, scan == 1 & is.na(area)) %>% select(site, plot, treatment, sp, ID, field_notes) # missing scan: D1 5W Mv3 (checked uncropped scan and only two Mv leaves received from D1 5W)


#### edit late August data ####

# rename mis-labelled photos before running function
ls_ev_late_aug2 <- col_fun(mutate(ls_ev_late_aug, month = "late_aug", sp = "Ev",
                                  Slice = case_when(Slice == "leaf:D1_1W_Ev_4_LateAug_edited" ~ "leaf:D1_1W_Ev_A_LateAug",
                                                    Slice == "lesions:D1_1W_Ev_4_LateAug_edited" ~ "lesions:D1_1W_Ev_A_LateAug",
                                                    Slice == "leaf:D1_4W_Ev_4_LateAug_edited" ~ "leaf:D1_4W_Ev_A_LateAug",
                                                    Slice == "lesions:D1_4W_Ev_4_LateAug_edited" ~ "lesions:D1_4W_Ev_A_LateAug",
                                                    Slice == "leaf:D3_2W_Ev_4_LateAug_edited" ~ "leaf:D3_2W_Ev_A_LateAug",
                                                    Slice == "lesions:D3_2W_Ev_4_LateAug_edited" ~ "lesions:D3_2W_Ev_A_LateAug",
                                                    Slice == "leaf:D3_1F_Ev_3_LatAug" ~ "leaf:D3_1F_Ev_A_LatAug",
                                                    Slice == "lesions:D3_1F_Ev_3_LatAug" ~ "lesions:D3_1F_Ev_A_LatAug",
                                                    Slice == "leaf:D3_1F_Ev_2_LatAug" ~ "leaf:D3_1F_Ev_3_LatAug",
                                                    Slice == "lesions:D3_1F_Ev_2_LatAug" ~ "lesions:D3_1F_Ev_3_LatAug",
                                                    Slice == "leaf:D3_7W_Ev_3_LateAug" ~ "leaf:D3_7W_Ev_2_LateAug",
                                                    Slice == "lesions:D3_7W_Ev_3_LateAug" ~ "lesions:D3_7W_Ev_2_LateAug",
                                                    Slice == "leaf:D4_2W_Ev_3_LatAug_edited" ~ "leaf:D4_2W_Ev_A_LatAug",
                                                    Slice == "lesions:D4_2W_Ev_3_LatAug_edited" ~ "lesions:D4_2W_Ev_A_LatAug",
                                                    Slice == "leaf:D4_2W_Ev_2_LatAug" ~ "leaf:D4_2W_Ev_3_LatAug",
                                                    Slice == "lesions:D4_2W_Ev_2_LatAug" ~ "lesions:D4_2W_Ev_3_LatAug",
                                                    TRUE ~ Slice))) %>%
  filter(plant != "D4_2F_Ev_1_LateAug")  # scans was unexpected and of very poor quality

# remove scans that were not a leaf or mis-labelled and couldn't be reconciled
ls_mv_late_aug2 <- col_fun(mutate(ls_mv_late_aug, month = "late_aug", sp = "Mv"))

# combine
dt_late_aug2 <- full_join(dt_late_aug, full_join(ls_ev_late_aug2, ls_mv_late_aug2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_late_aug2, is.na(scan)) %>% select(Slice, leaves_tot, field_notes) %>% unique() # none missing from field data
filter(dt_late_aug2, scan == 0 & !is.na(area)) # unexpected scans: D2 6W Ev3 (seems correct from uncropped scan and other information). D3 7W Ev3 - switched with Ev2 (missing). D4 2F Ev1 (tiny, dead leaf, better not to include).
filter(dt_late_aug2, scan == 1 & is.na(area)) # missing scans: 
# D1 2W Ev1
# D1 3W Ev3
# D1 5W Ev1
# D1 6W Ev2
# D1 7W Ev3
# D1 9W EvA (missing the uncropped for all of the D1 W Ev scans)
# D3 7W Ev2 (not in uncropped scan, only Ev 1)


#### edit September data ####

ls_ev_sep2 <- col_fun(mutate(ls_ev_sep, month = "sep", sp = "Ev")) %>%
  mutate(ID = case_when(ID == "S" ~ substr(plant, 9, 10) %>% 
                          gsub(c("_"), "", .),
                        TRUE ~ ID),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = ifelse(ID == "A" | (focal == 0 & plot > 7),
                      "adult", 
                      "seedling"))

ls_mv_sep2 <- col_fun(mutate(ls_mv_sep, month = "sep", sp = "Mv",
                             Slice = case_when(Slice == "leaf:D2_8W_Mv1or2_Sept19" ~ "leaf:D2_8W_Mv1_Sept19",
                                               Slice == "lesions:D2_8W_Mv1or2_Sept19" ~ "lesions:D2_8W_Mv1_Sept19",
                                               TRUE ~ Slice))) %>%
  mutate(ID = case_when(ID == "S" ~ substr(plant, 9, 10) %>% 
                          gsub(c("_"), "", .),
                        TRUE ~ ID),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0))

# combine
dt_sep2 <- full_join(dt_sep, full_join(ls_ev_sep2, ls_mv_sep2)) %>%
  filter(ID != "Edge")

# check against field data
filter(dt_sep2, is.na(scan)) %>% select(Slice, leaves_tot, field_notes) %>% unique() # none missing from field data
filter(dt_sep2, scan == 0 & !is.na(area)) # no unexpected scans
filter(dt_sep2, scan == 1 & is.na(area)) %>% select(site, plot, treatment, sp) %>% unique() %>% data.frame() # missing scans: the D1 and D2 scans were in the box that got lost in the mail. I'm not sure if the D3 and D4 scans were also lost in the mail or if they were not uploaded properly to OneDrive


#### combine months ####

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
# D4_9W_Ev_3_June needs note that it could be A
# D4_7F_Mv_1_July needs note that it could be 2
datw <- dat %>%
  filter(!is.na(part)) %>% # all are missing area values
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(date:field_notes, month:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area) %>%
  select(date, month, site, plot, treatment, sp, ID, focal, age, leaves_tot:Bp_spots, leaf_area.pix, lesion_area.pix, leaf_count) %>%
  mutate(scan_notes = case_when(month == "jun" & site == "D4" & plot == 9 & treatment == "water" & sp == "Ev" & ID == 3 ~ "could also be EvA",
                                month == "jul" & site == "D4" & plot == 7 & treatment == "fungicide" & sp == "Mv" & ID == 1 ~ "could also be Mv2",
                                TRUE ~ NA_character_))

datw_edge <- dt_edge %>%
  rename(count = Count) %>%
  select(-c(Slice)) %>%
  gather(variable, value, -c(month:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  rename(leaf_area.pix = leaf_area,
         lesion_area.pix = lesion_area) %>%
  select(month, site, plot, treatment, sp, ID, focal, age, leaf_area.pix, lesion_area.pix, leaf_count)


#### check values ####

# check that all plots were collected for Mv
datw_edge %>%
  group_by(month, site, treatment) %>%
  summarise(n = length(unique(plot))) %>%
  filter(n != 10)
data.frame() 
# D1 6F early Aug uncropped scan missing 
# Sep scans lost in mail: D1 8F, 9F and 10F, D1 9W and 10W, D4 9F and 10F, D4 9W and 10W

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

write_csv(datw, "intermediate-data/all_leaf_scans_2019_density_exp.csv")
write_csv(datw_edge, "intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
