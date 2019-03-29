##### info ####

# file: transmission-by-treatment
# author: Amy Kendig
# date last edited: 3/25/19
# goal: see how Mv and Ev disease transmission is affected by treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lme4)

# run seed data files
source("./code/ev-leaf-scans-data-processing-2018.R")
#source("./code/mv-leaf-scans-data-processing-2018.R")

# clear everything except seed data
rm(list = setdiff(ls(), c("eleaf")))

# import other data files
plots <- read_csv("./data/plot-treatments-2018-density-exp.csv")
foc <- read.csv("~/Google Drive/Microstegium Bipolaris/Data/Field Experiment Data/DensExp_FocalData_Processed_012419.csv")
bgev <- read_csv("~/Google Drive/Microstegium Bipolaris/Data/Field Experiment Data/DensExp_BgTaggedEvData_Processed_012419.csv")


#### modify data ####

# Ev tillers
et1 <- bgev %>%
  select(Month, Treatment, Site, PlotID, ID, LeavesTot, LeavesInfec, FocalSp) %>%
  rename(month = Month, treatment = Treatment, site = Site, plot = PlotID, leaves_tot = LeavesTot, leaves_infec = LeavesInfec, sp = FocalSp) %>%
  mutate(
    treatment = dplyr::recode(treatment, Fungicide = "fungicide", Water = "water"),
    ID = str_remove(ID, "A"),
    age = ifelse(plot < 8, "seedling", "adult")
  )

et2 <- foc %>%
  filter(FocalSp == "Ev") %>%
  select(Month, Treatment, Site, PlotID, ID, LeavesTot, LeavesInfec, FocalSp) %>%
  rename(month = Month, treatment = Treatment, site = Site, plot = PlotID, leaves_tot = LeavesTot, leaves_infec = LeavesInfec, sp = FocalSp) %>%
  mutate(
    treatment = dplyr::recode(treatment, Fungicide = "fungicide", Water = "water"),
    age = ifelse(ID == "A", "adult", "seedling")
  )

et <- full_join(et1, et2)

# Ev
de <- eleaf %>%
  filter(remove == 0) %>%
  mutate(
    sp = "Ev",
    plot = as.numeric(plot)
  ) %>%
  left_join(et)

# Mv
dm <- foc %>%
  filter(FocalSp == "Mv" & (RemoveScans == 0 | is.na(RemoveScans)) & !is.na(LeafArea.pix)) %>%
  select(Month, Treatment, PlantID, Site, PlotID, ID, LeafArea.pix, LesionArea.pix, FocalSp, LeavesTot, LeavesInfec) %>%
  rename(month = Month, treatment = Treatment, plant = PlantID, site = Site, plot = PlotID, leaf_area.pix = LeafArea.pix, lesion_area.pix = LesionArea.pix, sp = FocalSp, leaves_tot = LeavesTot, leaves_infec = LeavesInfec) %>%
  mutate(
    treatment = dplyr::recode(treatment, Fungicide = "fungicide", Water = "water")
  )

# merge
d <- full_join(de, dm) %>% full_join(plots) %>%
  mutate(
    leaves_infec2 = ifelse(leaves_infec == 0 & !is.na(leaf_area.pix), 1, leaves_infec),
    month = factor(month, levels = c("July", "August", "September")),
    week = dplyr::recode(month, July = 0, August = 8, September = 12),
    damage = lesion_area.pix * leaves_infec2 / (leaf_area.pix * leaves_tot),
    density_level = factor(density_level, levels = c("none", "low", "medium", "high")),
    treatment = factor(treatment, levels = c("water", "fungicide"))
  ) %>%
  filter(!is.na(damage))

# rep 0 background across treatments
d2 <- d %>%
  filter(background == "none") %>%
  ungroup() %>%
  select(-background) %>%
  merge(tibble(background = c("Mv seedling", "Ev seedling", "Ev adult")), all = T) %>%
  full_join(filter(d, background != "none"))


#### figures ####

pdf("./output/ev-damage-by-month.pdf",width=16,height=11)
d2 %>%
  filter(sp == "Ev") %>%
  ggplot(aes(x = month, y = damage, colour = treatment)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 2, position = position_dodge(0.2)) +
  facet_grid(background~density_level) +
  scale_colour_manual(values=c("#A5A5A5","#D55E00"), name = "Treatment") +
  theme(axis.text.y=element_text(size=24,colour="black"),
        axis.text.x=element_text(size=24,colour="black", angle = 45, hjust = 1),
        axis.title=element_text(size=26,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(linetype="solid",color="black",fill=NA,size=0.9),
        legend.text=element_text(size=24),
        legend.title=element_text(size=26),
        legend.key=element_blank(),
        plot.title=element_text(size=26,hjust=0.5),
        strip.background=element_blank(),
        strip.text=element_text(size=26)) +
  xlab("Month") + 
  ylab("Proportion leaf area damaged") +
  ggtitle("Elymus virginicus")
dev.off()

pdf("./output/mv-damage-by-month.pdf",width=16,height=11)
d2 %>%
  filter(sp == "Mv") %>%
  ggplot(aes(x = month, y = damage, colour = treatment)) + 
  stat_summary(fun.data = "mean_cl_boot", size = 2, position = position_dodge(0.2)) +
  facet_grid(background~density_level) +
  scale_colour_manual(values=c("#A5A5A5","#D55E00"), name = "Treatment") +
  theme(axis.text.y=element_text(size=24,colour="black"),
        axis.text.x=element_text(size=24,colour="black", angle = 45, hjust = 1),
        axis.title=element_text(size=26,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(linetype="solid",color="black",fill=NA,size=0.9),
        legend.text=element_text(size=24),
        legend.title=element_text(size=26),
        legend.key=element_blank(),
        plot.title=element_text(size=26,hjust=0.5),
        strip.background=element_blank(),
        strip.text=element_text(size=26)) +
  xlab("Month") + 
  ylab("Proportion leaf area damaged") +
  ggtitle("Microstegium vimineum")
dev.off()


#### calculate transmission ####

# find minimum infection
d2 %>%
  filter(damage > 0) %>%
  group_by(sp) %>%
  summarise(min_dam = min(damage))

d2 %>%
  filter(damage > 0) %>%
  ggplot(aes(x = damage)) +
  geom_histogram() +
  facet_wrap(~sp)

# replace NA or 0 values with minimum possible
d2 <- d2 %>%
  mutate(
    damage2 = case_when(damage == 0 & sp == "Mv" ~ 0.0000243,
                        damage == 0 & sp == "Ev" ~ 0.00000380,
                        TRUE ~ damage)
  )

# Calculate an exponential growth rate for each plot
fit_fun <- function(dam, time){
  m <- lm(log(dam) ~ time)
  return(m)
}

dlm <- d2 %>%
  filter((sp == "Ev" & month != "September") | sp == "Mv") %>%
  group_by(sp, site, treatment, plot, background, background_density, density_level) %>%
  summarise(
    intercept = coef(fit_fun(damage2, week))[1],
    r = coef(fit_fun(damage2, week))[2]
  )

d3 <- d2 %>% full_join(dlm) %>%
  mutate(
    pred = exp(intercept + r * week)
  )
  

#### visualize ####

d3 %>%
  filter(sp == "Mv" & treatment == "water") %>%
  ggplot(aes(x = week, y = damage, colour = site)) + 
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_grid(background~density_level)

d3 %>%
  filter(sp == "Mv" & treatment == "fungicide") %>%
  ggplot(aes(x = week, y = damage, colour = site)) + 
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_grid(background~density_level)

d3 %>%
  filter(sp == "Ev" & treatment == "water") %>%
  ggplot(aes(x = week, y = damage, colour = site)) + 
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_grid(background~density_level)

d3 %>%
  filter(sp == "Ev" & treatment == "fungicide") %>%
  ggplot(aes(x = week, y = damage, colour = site)) + 
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_grid(background~density_level)


#### transmission analysis ####

pred_fun <- function(x, ...){
  pred <- predict(x, re.form = NA, ...)
  as.data.frame(pred)
}

dlm2 <- dlm %>%
  filter(!is.na(r)) %>%
  group_by(sp, background) %>%
  nest() %>%
  mutate(m1 = map(.x = data, .f = ~ lmer(r ~ treatment * background_density + (1|site), data = .))) %>%
  mutate(pred = map(.x = m1, ~ pred_fun(.))) %>%
  select(data, sp, background, pred) %>%
  unnest


modesm = lmer(r ~ treatment * background_density + (1|site), data = filter(dlm, background == "Ev seedling" & sp == "Mv"))
summary(modesm)
drop1(modesm, test = "Chisq")
modesm.1 = update(modesm, .~. -treatment:background_density)
drop1(modesm.1, test = "Chisq") # treatment effect
summary(modesm.1) # fungicide reduces transmission

modeam = lmer(r ~ treatment * background_density + (1|site), data = filter(dlm, background == "Ev adult" & sp == "Mv"))
summary(modeam)
drop1(modeam, test = "Chisq")
modeam.1 = update(modeam, .~. -treatment:background_density)
drop1(modeam.1, test = "Chisq") # treatment effect
summary(modeam.1) # fungicide reduces transmission

modmm = lmer(r ~ treatment * background_density + (1|site), data = filter(dlm, background == "Mv seedling" & sp == "Mv"))
summary(modmm)
drop1(modmm, test = "Chisq")
modmm.1 = update(modmm, .~. -treatment:background_density)
drop1(modmm.1, test = "Chisq") # treatment effect
summary(modmm.1) # fungicide reduces transmission

modese = lmer(r ~ treatment * background_density + (1|site), data = filter(dlm, background == "Ev seedling" & sp == "Ev"))
summary(modese)
drop1(modese, test = "Chisq")
modese.1 = update(modese, .~. -treatment:background_density)
drop1(modese.1, test = "Chisq") # no sig effects

modeae = lmer(r ~ treatment * background_density + (1|site), data = filter(dlm, background == "Ev adult" & sp == "Ev"))
summary(modeae)
drop1(modeae, test = "Chisq")
modeae.1 = update(modeae, .~. -treatment:background_density)
drop1(modeae.1, test = "Chisq") # no sig effects

modme = lmer(r ~ treatment * background_density + (1|site), data = filter(dlm, background == "Mv seedling" & sp == "Ev"))
summary(modme)
drop1(modme, test = "Chisq")
modme.1 = update(modme, .~. -treatment:background_density)
drop1(modme.1, test = "Chisq") # no sig effects

# re-format dataset for background type analysis
dlm3 <- dlm %>%
  ungroup() %>%
  mutate(
    background = case_when(background_density == 0 ~ "none",
                           TRUE ~ background)
  ) %>%
  filter(background == "none") %>%
  unique() %>%
  full_join(filter(dlm, background_density > 0))
dim(dlm3)
dim(dlm)

# check that it worked
dlm3 %>%
  group_by(sp, treatment, background) %>%
  summarise(
    mean_r = mean(r, na.rm = T),
    se_r = sd(r, na.rm = T) / sqrt(sum(!is.na(r))),
    n_r = sum(!is.na(r))
  ) 

# stats
mod_m = lmer(r ~ treatment * background + (1|site), data = filter(dlm3, sp == "Mv"))
plot(mod_m) # looks good
library(car)
Anova(lmer(r ~ treatment * background + (1|site), data = filter(dlm3, sp == "Mv")), contrasts = list(treatment = contr.sum, background = contr.sum), type = 3) # treatment effect

mod_e = lmer(r ~ treatment * background + (1|site), data = filter(dlm3, sp == "Ev"))
plot(mod_e) # looks good
Anova(lmer(r ~ treatment * background + (1|site), data = filter(dlm3, sp == "Ev")), contrasts = list(treatment = contr.sum, background = contr.sum), type = 3) # treatment effect
summary(mod_e) # fungicide increases infection

#### visualize ####

# summarise
dlms <- dlm2 %>%
  group_by(sp, treatment, background, background_density) %>%
  summarise(
    mean_r = mean(r, na.rm = T),
    se_r = sd(r, na.rm = T) / sqrt(sum(!is.na(r)))
  ) %>%
  ungroup() 

# rename levels
dlms <- dlms %>%
  mutate(
    background = dplyr::recode(background, "Ev adult" = "adult natives", "Ev seedling" = "1st year natives", "Mv seedling" = "invaders"),
    sp = dplyr::recode(sp, "Ev" = "native", "Mv" = "invader") %>% factor(levels = c("native", "invader"))
    
  )

dlm2 <- dlm2 %>%
  mutate(
    background = dplyr::recode(background, "Ev adult" = "adult natives", "Ev seedling" = "1st year natives", "Mv seedling" = "invaders"),
    sp = dplyr::recode(sp, "Ev" = "native", "Mv" = "invader") %>% factor(levels = c("native", "invader"))
  )

dlm2_sub <- dlm2 %>%
  filter(sp == "invader")

pdf("./output/transmission-by-treatment.pdf",width=8.5,height=6.6)
dlms %>%
  ggplot(aes(x = background_density, y = mean_r)) +
  geom_line(data = dlm2_sub, aes(y = pred, colour = treatment), size = 2) +
  geom_errorbar(aes(ymin = mean_r - se_r, ymax = mean_r + se_r, group = treatment), width = 0, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), size = 5, shape = 21, aes(fill = treatment)) +
  facet_grid(sp~background, scales = "free") +
  scale_colour_manual(values=c("#C7D1AE","#4E3629"), name = "Treatment") +
  scale_fill_manual(values=c("#C7D1AE","#4E3629"), name = "Treatment") +
  theme(axis.text=element_text(size=20,colour="black"),
        axis.title=element_text(size=23,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(linetype="solid",color="black",fill=NA,size=0.9),
        legend.text=element_text(size=20),
        legend.title=element_text(size=23),
        legend.key=element_blank(),
        plot.title=element_text(size=23,hjust=0.5),
        strip.background=element_blank(),
        strip.text=element_text(size=23),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  xlab("Competitor density") + 
  ylab("Disease transmission") + 
  ggtitle("Competitor")
dev.off()




#### old stats and visualization ####

mm1 <- lmer(r ~ background_density * background * treatment + (1|site), data = filter(dlm, sp == "Mv"))
summary(mm1)

drop1(mm1, test = "Chisq")
mm2 <- update(mm1, .~. -background_density:background:treatment)
drop1(mm2, test = "Chisq")
mm3 <- lmer(r ~ background_density + background + treatment + (1|site), data = filter(dlm, sp == "Mv"))
drop1(mm3, test = "Chisq") # treatment

em1 <- lmer(r ~ background_density * background * treatment + (1|site), data = filter(dlm, sp == "Ev"))
summary(em1)

drop1(em1, test = "Chisq")
em2 <- update(em1, .~. -background_density:background:treatment)
drop1(em2, test = "Chisq")
em3 <- update(em2, .~. -background_density:background)
drop1(em3, test = "Chisq")
em4 <- update(em3, .~. -background_density:treatment)
drop1(em4, test = "Chisq")
summary(em4)

# summarise
lms <- dlm %>%
  group_by(background, density_level, sp, treatment) %>%
  summarise(
    meanr = mean(r, na.rm = T),
    ser = sd(r, na.rm = T) / sqrt(sum(!is.na(r)))
  ) %>%
  ungroup() %>%
  filter(density_level != "none" | (density_level == "none" & background == "Ev adult")) %>%
  mutate(
    density_level = as.character(density_level),
    background = case_when(
      density_level == "none" & background == "Ev adult" ~ "none",
      TRUE ~ background),
    background = factor(background, levels = c("none", "Ev adult", "Ev seedling", "Mv seedling")),
    density_level = dplyr::recode(density_level, medium = "med"),
    density_level = factor(density_level, levels = c("none", "low", "med", "high")),
    treatment = dplyr::recode(treatment, water = "Water", fungicide = "Fungicide")
  )

# figure
pdf("./output/damage-by-treatment.pdf",width=11,height=8)
lms %>%
  ggplot(aes(x = density_level, y = meanr, colour = background)) + 
  geom_point(size = 5, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = meanr - ser, ymax = meanr + ser), width = 0, position = position_dodge(0.4)) +
  facet_grid(sp ~ treatment) +
  scale_colour_manual(values=c("black","#002060","#44A5B6","#92D050"), name = "Competitor") +
  theme(axis.text=element_text(size=24,colour="black"),
        axis.title=element_text(size=26,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(linetype="solid",color="black",fill=NA,size=0.9),
        legend.text=element_text(size=24),
        legend.title=element_text(size=26),
        legend.key=element_blank(),
        plot.title=element_text(size=26,hjust=0.5),
        strip.background=element_blank(),
        strip.text=element_text(size=26),
        legend.position = "right") +
  ylab(expression(paste("Damage rate (", week^-1, ")", sep = ""))) + 
  xlab("Competitor density")
dev.off()
