##### info ####

# file: seeds-by-treatment
# author: Amy Kendig
# date last edited: 3/25/19
# goal: see how Mv and Ev seed production is affected by treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lme4)

# run seed data files
source("./code/ev-seeds-data-processing-2018.R")
source("./code/mv-seeds-data-processing-2018.R")

# clear everything except seed data
rm(list = setdiff(ls(), c("eseeds", "mseeds")))

# import other data files
plots <- read_csv("./data/plot-treatments-2018-density-exp.csv")
till <- read_csv("./data/focal-size-infec-jul-2018-density-exp.csv")
foc <- read.csv("~/Google Drive/Microstegium Bipolaris/Data/Field Experiment Data/DensExp_FocalData_Processed_012419.csv")
bgev <- read_csv("~/Google Drive/Microstegium Bipolaris/Data/Field Experiment Data/DensExp_BgTaggedEvData_Processed_012419.csv")


#### edit data ####

# Mv tiller number
dt <- till %>%
  mutate(
    treatment = dplyr::recode(treatment, "Water" = "water", "Fungicide" = "fungicide")
  ) %>%
  filter(sp == "Mv") %>%
  group_by(site, plot, treatment) %>%
  summarise(mean_till = mean(tillers, na.rm = T))

# Ev to remove
rem_ev <- eseeds %>%
  filter(remove == 1) %>%
  mutate(
    plant = paste(site, plot, treatment, age, ID, sep = ".")
  )

# combine Ev across dates
eseeds2 <- eseeds %>%
  group_by(site, plot, treatment, sp, age, ID, focal) %>%
  summarise(
    seeds = sum(seeds, na.rm = T)
  ) %>%
  mutate(
    plant = paste(site, plot, treatment, age, ID, sep = ".")
  ) %>%
  filter(!(plant %in% rem_ev$plant))

# extrac survival data
foc_surv <- foc %>%
  filter(Month == "September" & Focal != "Mv") %>%
  select(Site, PlotID, Treatment, Focal, ID, SurvCor) %>%
  rename(site = Site, plot = PlotID, treatment = Treatment, age = Focal, surv = SurvCor) %>%
  mutate(
    treatment = recode(treatment, Fungicide = "fungicide", Water = "water"),
    age = recode(age, "Ev adult" = "adult", "Ev seedling" = "seedling"),
    sp = "Ev",
    focal = 1
  )

bg_surv <- bgev %>%
  filter(Month == "September") %>%
  select(Site, PlotID, Treatment, ID, SurvCor) %>%
  rename(site = Site, plot = PlotID, treatment = Treatment, surv = SurvCor) %>%
  mutate(
    treatment = recode(treatment, Fungicide = "fungicide", Water = "water"),
    age = case_when(plot < 8 ~ "seedling",
                    plot >7 ~ "adult"),
    ID = str_remove(ID, "A"), 
    sp = "Ev",
    focal = 0
  )

ev_surv <- full_join(foc_surv, bg_surv) %>%
  mutate(
    plant = paste(site, plot, treatment, age, ID, sep = ".")
  ) %>%
  filter(!(plant %in% rem_ev$plant))

# merge plot data with seed data
eseeds2
plots
mseeds
dt
de <- filter(eseeds2, site %in% c("D1", "D2", "D3", "D4")) %>% full_join(ev_surv) %>% full_join(plots)
dm <- full_join(mseeds, plots) %>% left_join(dt)

# see if all have tiller counts
sum(is.na(dm$mean_till))

# estimate per capita seeds for Mv
dm <- dm %>%
  mutate(
    seeds_percap = seeds / 3 * mean_till
  )

# check for missing survival info
sum(is.na(de$surv))
filter(de, is.na(surv)) # no seed numbers anyway
filter(bgev, Site == "D1" & PlotID == 7 & Treatment == "Fungicide" & ID == "R2") %>% data.frame() # couldn't find plant in September, no green earlier
filter(bgev, Site == "D4" & PlotID == 7 & Treatment == "Fungicide" & ID == "R1") %>% data.frame() # couldn't find plant in August or September

# give Ev's 0's if no seeds were collected and they're still alive
filter(de, surv == 1 & is.na(seeds)) # 111 plants

de <- de %>%
  mutate(
    seeds = case_when(surv == 1 & is.na(seeds) ~ 0,
                      TRUE ~ seeds)
  )

# summarize by density
ms <- dm %>%
  group_by(sp, treatment, background, background_density) %>%
  summarise(mean_seeds = mean(seeds_percap, na.rm = T),
            se_seeds = sd(seeds_percap, na.rm = T)/sqrt(sum(!is.na(seeds_percap))))

mp <- ms %>%
  filter(background == "none") %>%
  ungroup() %>%
  select(-background) %>%
  merge(tibble(background = c("Mv seedling", "Ev seedling", "Ev adult")), all = T) %>%
  full_join(filter(ms, background != "none"))

es <- de %>%
  group_by(sp, age, treatment, background, background_density) %>%
  summarise(mean_seeds = mean(seeds, na.rm = T),
            se_seeds = sd(seeds, na.rm = T)/sqrt(sum(!is.na(seeds))))

ep <- es %>%
  filter(background == "none") %>%
  ungroup() %>%
  select(-background) %>%
  merge(tibble(background = c("Mv seedling", "Ev seedling", "Ev adult")), all = T) %>%
  full_join(filter(es, background != "none"))

dat <- full_join(mp, ep)  %>%
  mutate(
    focal = case_when(sp == "Mv" ~ "Mv seedling",
                      TRUE ~ paste(sp, age, sep = " ")),
    treatment = factor(treatment, levels = c("water", "fungicide"))
  )

d <- dm %>%
  select("site", "plot", "treatment", "sp", "background", "background_sp", "background_density", "density_level", "seeds_percap") %>%
  rename(seeds = seeds_percap) %>%
  full_join(de)%>%
  mutate(
    focal = case_when(sp == "Mv" ~ "Mv seedling",
                      TRUE ~ paste(sp, age, sep = " ")),
    treatment = factor(treatment, levels = c("water", "fungicide")),
    background = factor(background, levels = c("none", "Ev seedling", "Ev adult", "Mv seedling"))
  )

d2 <- d %>%
  filter(background == "none") %>%
  select(-background) %>%
  merge(tibble(background = c("Mv seedling", "Ev seedling", "Ev adult")), all = T) %>%
  full_join(filter(d, background != "none")) %>%
  mutate(
    site_plot = paste(site, plot, sep = "_")
  )

# summarize by background type
mb <- dm %>%
  group_by(sp, treatment, background) %>%
  summarise(mean_seeds = mean(seeds_percap, na.rm = T),
            se_seeds = sd(seeds_percap, na.rm = T)/sqrt(sum(!is.na(seeds_percap))))

eb <- de %>%
  group_by(sp, age, treatment, background) %>%
  summarise(mean_seeds = mean(seeds, na.rm = T),
            se_seeds = sd(seeds, na.rm = T)/sqrt(sum(!is.na(seeds))))

bsum <- full_join(mb,eb) %>%
  ungroup() %>%
  mutate(
    focal = case_when(sp == "Mv" ~ "Mv seedling",
                      TRUE ~ paste(sp, age, sep = " ")),
    treatment = factor(treatment, levels = c("water", "fungicide"))
    ) %>% mutate(
  focal = dplyr::recode(focal, "Mv seedling" = "invader", "Ev adult" = "adult native", "Ev seedling" = "1st year native"),
  background = dplyr::recode(background, "Mv seedling" = "invaders", "Ev adult" = "adult natives", "Ev seedling" = "1st year natives") %>% factor(levels = c("none", "1st year natives", "adult natives", "invaders"))
)


#### stats ####

# linear regressions for each focal-background combination
pred_fun <- function(x, ...){
  pred <- predict(x, type = "response", re.form = NA, ...)
  as.data.frame(pred)
}

d3 <- d2 %>%
  filter(!is.na(seeds)) %>%
  group_by(focal, background) %>%
  nest() %>%
  mutate(m1 = map(.x = data, .f = ~ lmer(seeds ~ treatment * background_density + (1|site), data = .))) %>%
  mutate(pred = map(.x = m1, ~ pred_fun(.))) %>%
  select(data, focal, background, pred) %>%
  unnest

# Ev adult regressions
mod_eaea = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Ev adult" & background == "Ev adult")) # can't converge
summary(mod_eaea)
mod_eaea = lmer(seeds ~ treatment * background_density + (1|site_plot), data = filter(d2, focal == "Ev adult" & background == "Ev adult")) 
drop1(mod_eaea, test = "Chisq")
mod_eaea_1 = update(mod_eaea, .~. -treatment:background_density)
drop1(mod_eaea_1, test = "Chisq") # no significant effects
summary(mod_eaea_1)

mod_eaes = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Ev adult" & background == "Ev seedling")) # singular fit
summary(mod_eaes) # no plot variation
mod_eaes = lmer(seeds ~ treatment * background_density + (1|site), data = filter(d2, focal == "Ev adult" & background == "Ev seedling"))
drop1(mod_eaes, test = "Chisq")
mod_eaes_1 = update(mod_eaes, .~. -treatment:background_density)
drop1(mod_eaes_1, test = "Chisq") # treatment effect
summary(mod_eaes_1) # fungicide reduces seeds

mod_eams = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Ev adult" & background == "Mv seedling"))
drop1(mod_eams, test = "Chisq")
mod_eams_1 = update(mod_eams, .~. -treatment:background_density)
drop1(mod_eams_1, test = "Chisq") # no significant effects

# Ev adult mean seeds before replicating 0 plots
mean(filter(de, age == "adult" & treatment == "water")$seeds, na.rm = T) # 72

# Ev seedling regressions
mod_esea = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Ev seedling" & background == "Ev adult")) # singular fit
summary(mod_esea) # no random effects
mod_esea = lm(seeds ~ treatment * background_density , data = filter(d2, focal == "Ev seedling" & background == "Ev adult")) 
drop1(mod_esea, test = "Chisq")
mod_esea_1 = update(mod_esea, .~. -treatment:background_density)
drop1(mod_esea_1, test = "Chisq") # no significant effects

mod_eses = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Ev seedling" & background == "Ev seedling"))
drop1(mod_eses, test = "Chisq")
mod_eses_1 = update(mod_eses, .~. -treatment:background_density)
drop1(mod_eses_1, test = "Chisq") # no significant effects

mod_esms = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Ev seedling" & background == "Mv seedling"))
summary(mod_esms)
drop1(mod_esms, test = "Chisq") # singular fit, probably from low site variation
mod_esms = lmer(seeds ~ treatment * background_density + (1|site_plot), data = filter(d2, focal == "Ev seedling" & background == "Mv seedling"))
drop1(mod_esms, test = "Chisq")
mod_esms_1 = update(mod_esms, .~. -treatment:background_density)
drop1(mod_esms_1, test = "Chisq") # no significant effects, singular fits
mod_esms = lm(seeds ~ treatment * background_density , data = filter(d2, focal == "Ev seedling" & background == "Mv seedling"))
drop1(mod_esms, test = "Chisq")
mod_esms_1 = update(mod_esms, .~. -treatment:background_density)
drop1(mod_esms_1, test = "Chisq") # no significant effects

# mean Ev seedling before replicating 0 plots
mean(filter(d2, age == "seedling")$seeds, na.rm = T) # 9

# Mv regressions
mod_msea = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Mv seedling" & background == "Ev adult"))
drop1(mod_msea, test = "Chisq")
mod_msea_1 = update(mod_msea, .~. -treatment:background_density)
drop1(mod_msea_1, test = "Chisq") # no significant effect
summary(mod_msea)

mod_mses = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Mv seedling" & background == "Ev seedling"))
drop1(mod_mses, test = "Chisq") # significant interaction, but model convergence error
summary(mod_mses) # density increases seed production with water, but not fungicide
mod_mses_1 = update(mod_mses, .~. - treatment:background_density) # this model can't converge
#library(car) # caution: this messes up recode
#Anova(mod_mses, type = 3) #interaction is significant

mod_msms = lmer(seeds ~ treatment * background_density + (1|site/plot), data = filter(d2, focal == "Mv seedling" & background == "Mv seedling")) # singular fit
summary(mod_msms) # site is zero variance
mod_msms = lmer(seeds ~ treatment * background_density + (1|site_plot), data = filter(d2, focal == "Mv seedling" & background == "Mv seedling"))
drop1(mod_msms, test = "Chisq") # singificant interaction
summary(mod_msms) # density decreases seed production with water, but not fungicide
# with water, 6 seeds lost per Mv seedling added

# mean Mv seeds without competitors or disease
mean(filter(d, focal == "Mv seedling" & treatment == "fungicide" & background_density == 0)$seeds, na.rm = T) # 512

# background type instead of density
mod_ms = lmer(seeds ~ treatment * background + (1|site/plot), data = filter(d, focal == "Mv seedling"))
summary(mod_ms)
plot(mod_ms) # variance increases with value
mod_ms_1 = lmer(log(seeds) ~ treatment * background + (1|site/plot), data = filter(d, focal == "Mv seedling"))
plot(mod_ms_1) # better
summary(mod_ms_1)
library(car)
Anova(lmer(log(seeds) ~ treatment * background + (1|site/plot), data = filter(d, focal == "Mv seedling")), contrasts = list(treatment = contr.sum, background = contr.sum), type = 3) # significant background effect

mod_es = lmer(seeds ~ treatment * background + (1|site/plot), data = filter(d, focal == "Ev seedling")) # singular fit
summary(mod_es) # no plot variance
mod_es = lmer(seeds ~ treatment * background + (1|site), data = filter(d, focal == "Ev seedling"))
plot(mod_es) # variance increases with value
mod_es_1 = lmer(log(seeds+1) ~ treatment * background + (1|site), data = filter(d, focal == "Ev seedling"))
plot(mod_es_1) # a little better
Anova(lmer(log(seeds+1) ~ treatment * background + (1|site), data = filter(d, focal == "Ev seedling")), contrasts = list(treatment = contr.sum, background = contr.sum), type = 3)  # significant interaction

mod_ea = lmer(seeds ~ treatment * background + (1|site/plot), data = filter(d, focal == "Ev adult"))
summary(mod_ea)
plot(mod_ea) # variance increases with value
mod_ea_1 = lmer(log(seeds+1) ~ treatment * background + (1|site/plot), data = filter(d, focal == "Ev adult"))
plot(mod_ea_1)
Anova(lmer(log(seeds+1) ~ treatment * background + (1|site/plot), data = filter(d, focal == "Ev adult")), contrasts = list(treatment = contr.sum, background = contr.sum), type = 3) # barely significant background effect


#### visualize ####

# subset d3 for significant effects and rename levels
d3_sub <- d3 %>%
  filter((focal == "Mv seedling" & background != "Ev adult") | (focal == "Ev adult" & background == "Ev seedling")) %>%
  mutate(
    focal = dplyr::recode(focal, "Mv seedling" = "invader", "Ev adult" = "adult native", "Ev seedling" = "1st year natives"),
    background = dplyr::recode(background, "Mv seedling" = "invaders", "Ev adult" = "adult natives", "Ev seedling" = "1st year natives")
  )

# rename for poster
dat <- dat %>%
  mutate(
    focal = dplyr::recode(focal, "Mv seedling" = "invader", "Ev adult" = "adult native", "Ev seedling" = "1st year native"),
    background = dplyr::recode(background, "Mv seedling" = "invaders", "Ev adult" = "adult natives", "Ev seedling" = "1st year natives")
  )
  

pdf("./output/seeds-by-treatment.pdf",width=8.5,height=8.9)
dat %>%
  ggplot(aes(x = background_density, y = mean_seeds)) +
  #geom_line(data = d3_sub, aes(y = pred, colour = treatment), size = 2) +
  geom_errorbar(aes(ymin = mean_seeds - se_seeds, ymax = mean_seeds + se_seeds, group = treatment), width = 0, position = position_dodge(0.3)) +
  geom_point(position = position_dodge(0.3), size = 5, shape = 21, aes(fill = treatment)) +
  facet_grid(focal~background, scales = "free") +
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
        legend.position="bottom",
        legend.direction="horizontal") +
  xlab("Competitor density") + 
  ylab("Per capita seed production") + 
  ggtitle("Competitor")
dev.off()

pdf("./output/seeds-by-treatment-2.pdf",width=8.5,height=6.6)
bsum %>%
  ggplot(aes(x = background, y = mean_seeds)) +
  #geom_line(data = d3_sub, aes(y = pred, colour = treatment), size = 2) +
  geom_errorbar(aes(ymin = mean_seeds - se_seeds, ymax = mean_seeds + se_seeds, group = treatment), width = 0, position = position_dodge(0.3)) +
  geom_point(position = position_dodge(0.3), size = 5, shape = 21, aes(fill = treatment)) +
  facet_wrap(~focal, scales = "free") +
  scale_colour_manual(values=c("#C7D1AE","#4E3629"), name = "Treatment") +
  scale_fill_manual(values=c("#C7D1AE","#4E3629"), name = "Treatment") +
  theme(axis.text.y=element_text(size=20,colour="black"),
        axis.text.x=element_text(size=18,colour="black", angle = 45, hjust = 1),
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
        legend.position="bottom",
        legend.direction="horizontal") +
  xlab("Competitor type") + 
  ylab("Per capita seed production") 
dev.off()

#### estimating competition coefficients ####

# adult Elymus intraspecific competition
dat %>%
  filter(background =="adult natives" & background_density %in% c(4, 8) & focal != "invader" & treatment == "water")

# focal adults
(90.06 - 68.79)/4 # 5 seeds lost per adult added

# focal seedlings
(10.52 - 5.92)/4 # 1 seed lost per adult added


#### simulate density-dependent models ####

a <- -0.7
init <- 750
b <- 2

sim <- tibble(
  density = 0:64
) %>%
  mutate(
    m1 = init / (1 + a * density),
    m2 = init / 1 + (a * density)^b,
    m3 = init * exp(-a * density),
    m4 = init * exp(-a * log(density + 1)),
    m5 = init / (1 + density^a),
    m6 = init / (1 + a * density)^b,
    m7 = 1 + init * (1 - a * density),
    m8 = init * (density + 1)^-a
  ) %>%
  gather(key = "model", value = "seeds", -density)

sim %>% 
  filter(!(model %in% c("m3", "m6", "m7"))) %>%
  ggplot(aes(x = density, y = seeds, colour = model)) +
  geom_line()
