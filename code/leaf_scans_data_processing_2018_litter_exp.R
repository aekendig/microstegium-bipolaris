##### info ####

# file: leaf_scans_data_processing_litter_exp_2018
# author: Amy Kendig
# date last edited: 1/6/20
# goal: process raw 2018 Elymus leaf scan data
# background: leaf scans were analyzed using FIJI, script: ev_leaf_damage_severity.ijm


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import raw data file
fiji <- read_delim("./data/leaf-scans-litter-exp-20180927/text-output/ev_leaf_scan_text_output_sep_2018.tsv", delim = "\t")


#### edit data ####

# select count and area measurements
fiji2 <- fiji %>%
  select(Slice:'Total Area')

# rename columns
colnames(fiji2) <- c("slice", "objects", "area.pix")

# add new columns
fiji3 <- fiji2 %>%
  mutate(
    date = 20180927,
    plant = gsub(".*:","", slice),
    part = gsub(":.*$", "", slice),
  ) %>%
  mutate(
    part = recode(part, lesions = "lesion"),
    site = substr(plant, 1, 2),
    plot = substr(plant, 4, 5) %>% gsub("[^[:digit:]]", "", .),
    age = "adult",
    ID = "A"
  )

# spread by part
fijiw <- fiji3 %>%
  select(-c(slice)) %>%
  gather(variable, value, -(date:ID)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value)


#### check values ####

# leaf area
fijiw %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram()

# lesion area 
fijiw %>%
  ggplot(aes(x = lesion_area.pix)) +
  geom_histogram()

# leaf objects
fijiw %>%
  ggplot(aes(x = leaf_objects)) +
  geom_histogram()

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

write_csv(fijiw, "./intermediate-data/ev_damage_sep_2018_litter_exp.csv")



