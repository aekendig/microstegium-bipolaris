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
jul19 <- read_csv("data/fungal_isolation_jul_2019_density_exp.csv") # Mv only
mv_aug19_1 <- read_csv("data/mv_fungal_isolation_early_aug_2019_density_exp.csv")
mv_aug19_2 <- read_csv("data/mv_fungal_isolation_2_early_aug_2019_density_exp.csv")
ev_aug19 <- read_csv("data/ev_fungal_isolation_early_aug_2019_density_exp.csv")
# late August Ev leaves or datasheet were lost

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

not_na_fun <- function(x){
  !is.na(x)
}


#### edit data ####

# July 2018
unique(jul18$symptoms) # no sample
unique(jul18$observation)
filter(jul18, is.na(observation))
unique(jul18$gigantea)
jul18 %>%
  filter(is.na(gigantea))
jul18 %>%
  filter(if_all(c(symptoms, observation, gigantea), is.na)) 
# remove these two rows
filter(jul18, plot == 8 & treatment == "fungicide" & site == "D4")
# also remove D4 8F (tree fell on plot)

jul18b <- jul18 %>%
  filter(if_any(c(symptoms, observation, gigantea), not_na_fun)) %>% # keep columns with any info
  filter(!(site == "D4" & plot == 8 & treatment == "fungicide")) %>%
  filter(symptoms != "no sample" | is.na(symptoms)) %>%
  mutate(bipolaris = case_when(gigantea == "Yes" ~ 1,
                               gigantea == "No" ~ 0))

filter(jul18b, is.na(gigantea))

# August 2018
colnames(aug18)
unique(aug18$symptoms)
filter(aug18, is.na(symptoms)) %>%
  data.frame() # some had Pyricularia isolated - check all columns for missing
unique(aug18$observation)
unique(aug18$gigantea)
unique(aug18$small_Bipolaris_isolated)
aug18 %>%
  filter(is.na(gigantea)) %>% 
  data.frame()
# use gigantea and small_Bipolaris_isolated columns

aug18b <- aug18 %>%
  filter(if_any(c(symptoms, observation, gigantea, gigantea_isolated, Pyricularia_isolated, small_Bipolaris_isolated, Curvularia_isolated), not_na_fun)) %>%  # keep columns with any info
  mutate(bipolaris = case_when(gigantea == "Yes" | small_Bipolaris_isolated == 1 ~ 1,
                               TRUE ~ 0))

filter(aug18b, is.na(gigantea) & is.na(small_Bipolaris_isolated)) # all had symptoms or Pyricularia isolations

# September 2018
unique(sep18$symptoms) # no sample
unique(sep18$observation)
unique(sep18$gigantea)
sep18 %>%
  filter(is.na(gigantea))

sep18b <- sep18 %>%
  filter(if_any(c(symptoms, observation, gigantea), not_na_fun)) %>% # keep columns with any info
  filter(symptoms != "no sample" | is.na(symptoms)) %>%
  mutate(bipolaris = case_when(gigantea == "Yes" ~ 1,
                               TRUE ~ 0))

filter(sep18b, is.na(gigantea)) %>% select(symptoms, observation) %>% unique()

# combine 2018 data
dat18 <- jul18b %>%
  full_join(aug18b) %>%
  full_join(sep18b)

# summarize by month, treatment, and species
site_sum18 <- dat18 %>%
  filter(bipolaris == 1) %>%
  group_by(month, sp, treatment) %>%
  summarize(sites = n_distinct(site)) %>%
  ungroup() %>%
  pivot_wider(names_from = sp,
              values_from = sites) %>%
  mutate(Ev = replace_na(Ev, 0)) %>%
  full_join(dat18 %>%
              group_by(month, sp, treatment) %>%
              summarize(sites = n_distinct(site)) %>%
              ungroup() %>%
              pivot_wider(names_from = sp,
                          values_from = sites,
                          names_glue = "{sp}_sites")) %>%
  mutate(month = paste0(month, " 2018"),
         E.virginicus = if_else(Ev_sites < 4, paste0(Ev, " (", Ev_sites, ")"), as.character(Ev)),
         M.vimineum = if_else(Mv_sites < 4, paste0(Mv, " (", Mv_sites, ")"), as.character(Mv)),
         treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control")) %>%
  select(month, treatment, E.virginicus, M.vimineum) %>%
  arrange(month, treatment)

# July 2019 (Mv only)
unique(jul19$notes) # no sample (empty bag)
unique(jul19$symptoms)
unique(jul19$gigantea)
filter(jul19, is.na(gigantea)) # removes missing data

jul19b <- jul19 %>%
  filter(if_any(c(symptoms, leaves_with_eyespots, gigantea, leaves_with_gigantea, gigantea_isolated, notes), not_na_fun)) %>% # keep columns with any info
  filter(notes != "no sample (empty bag)" | is.na(notes)) %>%
  mutate(bipolaris = case_when(gigantea == "Yes" ~ 1,
                               gigantea == "No" ~ 0))

filter(jul19b, is.na(gigantea))

# early Aug 2019
unique(mv_aug19_1$notes)
filter(mv_aug19_1, notes == "plate not found/lost")
unique(mv_aug19_1$gigantea)
unique(mv_aug19_1$small_Bipolaris)
filter(mv_aug19_1, is.na(gigantea)) # removes missing data
filter(mv_aug19_1, is.na(small_Bipolaris))

unique(mv_aug19_2$gigantea)
unique(mv_aug19_2$small_Bipolaris) # no missing data

unique(ev_aug19$notes) # empty bag
unique(ev_aug19$leaves_with_gigantea) # empty bag
unique(ev_aug19$symptoms)
filter(ev_aug19, is.na(leaves_with_gigantea)) %>%
  select(leaves, symptoms, notes) %>%
  unique() # leaves without eyespots weren't asessed for signs, maybe?
filter(ev_aug19, !is.na(leaves_with_gigantea)) %>%
  select(leaves, symptoms, notes) %>%
  unique() 

aug19b <- mv_aug19_1 %>%
  filter(if_any(c(Alternaria, gigantea, small_Bipolaris, Cladosporium, Fusarium, Pyricularia, Curvularia, Pyricularia_isolated, gigantea_isolated, small_Bipolaris_isolated, notes), not_na_fun)) %>% # keep columns with any info
  filter(notes != "plate not found/lost" | is.na(notes)) %>%
  full_join(mv_aug19_2 %>%
              filter(if_any(c(gigantea, small_Bipolaris, Pyricularia), not_na_fun))) %>%
  mutate(bipolaris = case_when(gigantea == "Yes" | small_Bipolaris == "Yes" ~ 1,
                               TRUE ~ 0)) %>%
  full_join(ev_aug19 %>%
              filter(if_any(c(symptoms, leaves, leaves_with_eyespots, leaves_with_gigantea, Cladosporium_isolated, notes), not_na_fun)) %>%
              filter(notes!= "empty bag" | is.na(notes)) %>%
              filter(leaves_with_gigantea != "empty bag" | is.na(leaves_with_gigantea)) %>%
              mutate(leaves_with_gigantea = as.numeric(leaves_with_gigantea),
                     bipolaris = if_else(leaves_with_gigantea > 0, 1, 0),
                     bipolaris = replace_na(bipolaris, 0)))

filter(aug19b, sp == "Mv" & (is.na(gigantea) | is.na(small_Bipolaris)))
filter(aug19b, sp == "Ev" & is.na(leaves_with_gigantea))

# late August (Ev only)
# missing

# combine 2019 data
dat19 <- jul19b %>%
  full_join(aug19b)

# summarize by month, treatment, and species
site_sum19 <- dat19 %>%
  filter(bipolaris == 1) %>%
  group_by(month, sp, treatment) %>%
  summarize(sites = n_distinct(site)) %>%
  ungroup() %>%
  pivot_wider(names_from = sp,
              values_from = sites) %>%
  mutate(Ev = replace_na(Ev, 0)) %>%
  full_join(dat19 %>%
              group_by(month, sp, treatment) %>%
              summarize(sites = n_distinct(site)) %>%
              ungroup() %>%
              pivot_wider(names_from = sp,
                          values_from = sites,
                          names_glue = "{sp}_sites") %>%
              mutate(Ev_sites = replace_na(Ev_sites, 0))) %>%
  mutate(month = paste0(month, " 2019"),
         E.virginicus = if_else(Ev_sites < 4, paste0(Ev, " (", Ev_sites, ")"), as.character(Ev)),
         M.vimineum = if_else(Mv_sites < 4, paste0(Mv, " (", Mv_sites, ")"), as.character(Mv)),
         treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control")) %>%
  select(month, treatment, E.virginicus, M.vimineum) %>%
  arrange(desc(month), treatment)

write_csv(bind_rows(site_sum18, site_sum19), "output/bipolaris_identification_sites_density_exp.csv")

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

