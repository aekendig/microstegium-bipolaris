##### info ####

# file: survival_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/5/20
# goal: evaluate the effects of density and fungicide on survival


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
daty1_0 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
daty2_0 <- read_csv("data/all_replacement_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsfig <- read_csv("data/plot_treatments_for_figures_2018_2019_density_exp.csv")


#### edit data ####

# add plant type
daty1_1 <- daty1_0 %>%
  mutate(plant_type = paste(sp, age, sep = " "))

daty2_1 <- daty2_0 %>%
  mutate(plant_type = paste(sp, age, sep = " "))

# summer survival 2018, figure
# remove NA's 
# remove background Microstegium (not individual plants)
sumfigy1 <- daty1_1 %>%
  filter(month == "September" & !is.na(survival) & !(sp == "Mv" & focal == 0)) %>%
  left_join(plotsfig)

# survival through winter given summer survival, 2018, figure
# remove NA's 
# remove background Microstegium (not individual plants)
winfigy1 <- daty1_1 %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  filter(September == 1 & !is.na(April) & !(sp == "Mv" & focal == 0)) %>%
  select(-September) %>%
  rename(survival = April) %>%
  left_join(plotsfig)

# add IDs to 2019 background plants
# remove unnecessary columns
daty2_2 <- daty2_1 %>%
  group_by(site, plot, treatment, plant_type, focal, replace_date) %>%
  mutate(Bg_ID = as.character(1:n()),
         ID = case_when(ID == "Bg" ~ Bg_ID,
                        TRUE ~ ID)) %>%
  ungroup() %>%
  select(-c(sp, age, Bg_ID))

# focal plant ID's
focid <- tibble(plant_type = c(rep(unique(daty2_1$plant_type)[1:2], each = 3), unique(daty2_1$plant_type)[3]),
              ID = c(1, 2, 3, 1, 2, 3, "A"),
              focal = 1)

# background plant ID's
bgid <- plots %>% 
  select(plot, background, background_density) %>%
  unique() %>%
  filter(background != "none") %>%
  mutate(start = 1,
         focal = 0) %>%
  group_by(plot, background, focal) %>%
  expand(start, background_density, ID = full_seq(start:background_density, 1)) %>%
  ungroup() %>%
  select(-start) %>%
  rename(plant_type = background) %>%
  mutate(ID = as.character(ID))

# merge ID lists with plot, figures
focidfig <- plotsfig %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(focid)

bgidfig <- plotsfig %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  inner_join(bgid)

allidfig <- full_join(focidfig, bgidfig)

# check that number of plants is correct
2*4*(7*3 + 4+7 + 16+7 + 64+7 + 4+7 + 8+7 + 16+7 + 2+7 + 4+7 + 8+7)
# yes
  
# summer survival 2018, figure
sumfigy2 <- daty2_2 %>%
  select(-replace_date) %>%
  unique() %>%
  mutate(survival = 0) %>%
  full_join(allidfig) %>%
  mutate(survival = replace_na(survival, 1))


#### visualize ####

# template figure
temp_fig <- ggplot(plotsfig, aes(x = background_density, y = plot, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun.y = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_grid(plant_type ~ background, scales = "free_x", switch = "both") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside") +
    xlab("Background density")

# save
pdf("./output/survival_analysis_2018_2019_density_exp.pdf")

# focal plants, summer survival, 2018
temp_fig %+%
  filter(sumfigy1, focal == 1) %+%
  aes(y = survival) +
  ylab("Year 1 summer focal survival")

# all plants, summer survival, 2018
temp_fig %+%
  sumfigy1 %+%
  aes(y = survival) +
  ylab("Year 1 summer all survival")

# focal plants, winter survival, 2018
temp_fig %+%
  filter(winfigy1, focal == 1) %+%
  aes(y = survival) +
  ylab("Year 1 winter focal survival")

# all plants, winter survival, 2018
temp_fig %+%
  winfigy1 %+%
  aes(y = survival) +
  ylab("Year 1 winter all survival")   
  
# focal plants, summer survival, 2019
temp_fig %+%
  filter(sumfigy2, focal == 1) %+%
  aes(y = survival) +
  ylab("Year 2 summer focal survival")

# all plants, summer survival, 2019
temp_fig %+%
  sumfigy2 %+%
  aes(y = survival) +
  ylab("Year 2 summer all survival") 

dev.off()

