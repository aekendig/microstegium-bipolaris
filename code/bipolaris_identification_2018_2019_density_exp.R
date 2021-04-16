##### info ####

# file: bipolaris_identification_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/14/21
# goal: proportion of leaves Bipolaris


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
jul19 <- read_csv("data/fungal_isolation_jul_2019_density_exp.csv")
mv_aug19 <- read_csv("data/mv_fungal_isolation_early_aug_2019_density_exp.csv")
# missing two sites (asked Brett)
ev_aug19 <- read_csv("data/ev_fungal_isolation_early_aug_2019_density_exp.csv")
# missing late August (ask Brett and Ashish)

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

col_pal = c("#440154FF", "#472F7DFF", "#39568CFF", "#35B779FF")
shape_pal = c(21, 23, 22, 24)


#### edit data ####

# July 2018
unique(jul18$symptoms)
unique(jul18$observation)
unique(jul18$gigantea)
jul18 %>%
  filter(is.na(gigantea))

jul18b <- jul18 %>%
  filter(!is.na(gigantea)) %>%
  mutate(bipolaris = case_when(gigantea == "Yes" ~ 1,
                               gigantea == "No" ~ 0))

# August 2018
colnames(aug18)
unique(aug18$symptoms)
unique(aug18$observation)
unique(aug18$gigantea)
unique(aug18$small_Bipolaris_isolated)
aug18 %>%
  filter(is.na(gigantea))
# use gigantea and small_Bipolaris_isolated columns

aug18b <- aug18 %>%
  mutate(bipolaris = case_when(gigantea == "Yes" ~ 1,
                               gigantea == "No" ~ 0,
                               small_Bipolaris_isolated == 1 ~ 1,
                               !is.na(symptoms) & is.na(gigantea) & is.na(small_Bipolaris_isolated) ~ 0,
                               TRUE ~ NA_real_)) %>%
  filter(!is.na(bipolaris))

# September 2018
unique(sep18$symptoms) # no sample
unique(sep18$observation)
unique(sep18$gigantea)
sep18 %>%
  filter(is.na(gigantea))

sep18b <- sep18 %>%
  mutate(bipolaris = case_when(gigantea == "Yes" ~ 1,
                               gigantea == "No" ~ 0,
                               !is.na(symptoms) & symptoms != "no sample" & is.na(gigantea) ~ 0,
                               TRUE ~ NA_real_)) %>%
  filter(!is.na(bipolaris))

# combine 2018 data
dat18 <- jul18b %>%
  full_join(aug18b) %>%
  full_join(sep18b)

#### Start here ####

# July 2019
# asked Brett and Ashish

# early Aug 2019



#### figure ####

ggplot(dat18, aes(sp, bipolaris, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(0.1)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(0.1))
# maybe don't include treatments because leaves with eyespots were specifically chosen

ggplot(dat18, aes(sp, bipolaris, color = sp)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", size = 2, fun = "mean", aes(shape = sp)) +
  scale_color_manual(values = col_pal[c(1, 4)]) +
  scale_shape_manual(values = shape_pal[c(1, 4)]) +
  ylab(expression(paste("Proportion of diseased leaves with ", italic(Bipolaris), " infection", sep = ""))) +
  fig_theme
# could put year on the x-axis and make species the color if I get 2019 data

