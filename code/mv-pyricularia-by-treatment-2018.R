##### info ####

# file: mv-pyricularia-by-treatment-2018
# author: Amy Kendig
# date last edited: 4/18/19
# goal: see if the treatment affected possible Pyricularia infection


#### set up ####

# run leaf scan file (clears workspace and loads tidyverse)
source("./code/mv-leaf-scans-data-processing-2018.R")

# clear everything except final data
rm(list = setdiff(ls(), "mleaf"))

# import data
d <- read_csv("./data/mv-leaf-scans-pyricularia-2018-density-exp.csv")



#### edit data ####

# look at identifying info
d %>%
  select(leaf) %>%
  data.frame()

# pull out identifying info
d <- d %>%
  mutate(
    site = substr(leaf, 1, 2),
    treatment = gsub("D", "", leaf) %>% str_extract("[A-Z]+") %>% recode("W" = "water", "F" = "fungicide"),
    plot = substr(leaf, 4, 5) %>% gsub("[^[:digit:]]", "", .),
    ID = substr(leaf, 7, 8) %>% gsub("[^[:digit:]]", "", .),
    pyricularia = 1
  )

# merge with full dataset

d2 <- mleaf %>%
  filter(focal == 1) %>%
  ungroup() %>%
  full_join(d) %>%
  mutate(pyricularia = replace_na(pyricularia, 0),
         pyricularia.f = as.factor(pyricularia) %>% recode("1" = "yes", "0" = "no"),
         month = factor(month, levels = c("July", "August", "September")),
         lesion.perc = (lesion_area.pix - green_area.pix) / leaf_area.pix) 


#### visualize ####

# treatment effect on pyricularia
d2 %>%
  ggplot(aes(treatment, pyricularia)) +
  stat_summary(fun.data = "mean_cl_boot", na.rm = T) +
  facet_wrap(~month)
# appears to be suppressed by fungicide

# relationship between pyricularia occurrence and total lesions
d2 %>%
  ggplot(aes(treatment, lesion.perc, colour = pyricularia.f)) +
  stat_summary(fun.data = "mean_cl_boot", na.rm = T, position = position_dodge(0.3)) +
  facet_wrap(~month)


#### save figure ####

pdf("./output/mv-pyricularia-by-treatment.pdf")
d2 %>%
  ggplot(aes(treatment, pyricularia)) +
  stat_summary(fun.data = "mean_cl_boot", na.rm = T) +
  facet_wrap(~month) +
  xlab("Treatment") +
  ylab("Proportion of leaves with possible Pyricularia") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank())
dev.off()

pdf("./output/mv-pyricularia-by-treatment-damage.pdf")
d2 %>%
  ggplot(aes(treatment, lesion.perc, fill = pyricularia.f)) +
  stat_summary(fun.data = "mean_cl_boot", na.rm = T, position = position_dodge(0.3), shape = 21) +
  facet_wrap(~month) +
  xlab("Treatment") +
  ylab("Proportion of leaf area damaged") +
  scale_fill_manual(values = c("white", "black"), name = "Pyricularia") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.1, 0.8),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
dev.off()