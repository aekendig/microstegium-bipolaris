##### info ####

# file: leaf_scans_data_processing_litter_exp_2018
# author: Amy Kendig
# date last edited: 1/26/20
# goal: process raw 2019 Microstegium leaf scan data
# background: leaf scans were analyzed using FIJI, script: mv_leaf_damage_severity.ijm


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import raw data file
fiji_jn <- read_delim("./data/leaf-scans-litter-exp-20190606/text-output/mv_leaf_scan_text_output_jun_2019_litter_exp.tsv", delim = "\t")
fiji_jl <- read_delim("./data/leaf-scans-litter-exp-20190705/text-output/mv_leaf_scan_text_output_jul_2019_litter_exp.tsv", delim = "\t")


#### edit data ####

# combine data
fiji <- fiji_jn %>%
  mutate(date = 20190606) %>%
  rbind(fiji_jl %>%
          mutate(date = 20190705))

# select count and area measurements
fiji2 <- fiji %>%
  select(Slice:'Total Area', date)

# rename columns
colnames(fiji2) <- c("slice", "objects", "area.pix", "date")

# add new columns
fiji3 <- fiji2 %>%
  mutate(
    plant = gsub(".*:","", slice) %>% gsub("Litter_", "", .),
    part = gsub(":.*$", "", slice) %>% recode(lesions = "lesion"),
    site = substr(plant, 1, 2),
    treatment = substr(plant, 4, 4) %>% recode(A = "addition", C = "control", R = "removal"),
    block = substr(plant, 5, 5),
    rep = case_when(date == 20190705 ~ parse_number(substr(plant, 7, 15)),
                    TRUE ~ 1),
    sp = "Mv"
  ) %>%
  select(-c(slice, plant))
# can safely ignore warning

# spread by part
fijiw <- fiji3 %>%
  gather(variable, value, -(date:sp)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value)


#### check values ####

# leaf area
fijiw %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram()

filter(fijiw, leaf_area.pix > 500000) %>% data.frame()
# okay

# lesion area 
fijiw %>%
  ggplot(aes(x = lesion_area.pix)) +
  geom_histogram()

filter(fijiw, lesion_area.pix > 8e5) %>% data.frame()
# okay

# leaf objects
fijiw %>%
  ggplot(aes(x = leaf_objects)) +
  geom_histogram()

filter(fijiw, leaf_objects == 12) 
# all are leaves

# lesion objects
fijiw %>%
  ggplot(aes(x = lesion_objects)) +
  geom_histogram()

# area damaged
fijiw %>%
  mutate(damage = lesion_area.pix/(leaf_area.pix + lesion_area.pix)) %>%
  ggplot(aes(x = damage)) +
  geom_histogram()


#### save data ####

# separate by month
fiji_jnw <- filter(fijiw, date == 20190606)
fiji_jlw <- filter(fijiw, date == 20190705)

# save files
write_csv(fiji_jnw, "./intermediate-data/mv_damage_jun_2019_litter_exp.csv")
write_csv(fiji_jlw, "./intermediate-data/mv_damage_jul_2019_litter_exp.csv")



