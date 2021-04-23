#### info ####

# file: plot_scale_responses_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/21/21
# goal: plot-level biomass, seeds, and severity


# t-test for each species in each year with plots paired (compares seedling to seedling and adult to adult for Ev)
# use all plots in which that species is planted as a background species
# logit-transform severity
# hedges d to display effect sizes


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(car)
library(tidyverse)
library(effsize)
library(cowplot)

# import data
d1dat <- read_csv("intermediate-data/plot_biomass_seeds_severity_2018_density_exp.csv")
d2dat <- read_csv("intermediate-data/plot_biomass_seeds_severity_2019_density_exp.csv")

# logit adjustment
log_adj = 0.001

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"))

col_pal = c("#FFC000", "#8FD744FF")



#### functions ####

dat_format_fun <- function(res_var, sp_abb, year){

  # select dataset
  if(year == "2018"){
    dat <- d1dat2
  }else{
    dat <- d2dat2
  }

  # filter by species
  dat2 <- dat %>%
    filter(sp == sp_abb)

  # make data wide
  datw <- dat2 %>%
    select(site, plot, treatment, {{res_var}}) %>%
    pivot_wider(names_from = treatment,
                values_from = {{res_var}}) %>%
    drop_na()

  return(datw)

}

t_test_fun <- function(res_var, sp_abb, year){
  
  # format data
  datw <- dat_format_fun(res_var, sp_abb, year)

  # t-test
  mod <- t.test(datw$fungicide, datw$water, paired = T)

  # output
  dato <- tibble(t_val = as.numeric(mod$statistic),
                 t_df = as.numeric(mod$parameter),
                 t_p = mod$p.value)
  return(dato)
  
}

hedges_d_fun <- function(res_var, sp_abb, year){

  # format data
  datw <- dat_format_fun(res_var, sp_abb, year)
  
  # hedge's d
  d <- cohen.d(datw$fungicide, datw$water, hedges.correction = T, paired = T)
  
  # output
  dato <- tibble(hedges_d = d$estimate,
                 d_lower = as.numeric(d$conf.int[1]),
                 d_upper = as.numeric(d$conf.int[2]))
  return(dato)

}


#### edit data ####

# logit-transform proportions
d1dat2 <- d1dat %>%
  mutate(logit_jul_severity = logit(jul_severity, adjust = log_adj),
         logit_late_aug_severity = logit(late_aug_severity, adjust = log_adj),
         logit_sep_severity = logit(sep_severity, adjust = log_adj))

d2dat2 <- d2dat %>%
  mutate(logit_may_severity = logit(may_severity, adjust = log_adj),
         logit_jun_severity = logit(jun_severity, adjust = log_adj),
         logit_jul_severity = logit(jul_severity, adjust = log_adj),
         logit_early_aug_severity = logit(early_aug_severity, adjust = log_adj),
         logit_late_aug_severity = logit(late_aug_severity, adjust = log_adj))


#### visualizations ####

# 2018
ggplot(d1dat2, aes(x = logit_jul_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d1dat2, aes(x = logit_late_aug_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d1dat2, aes(x = logit_sep_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d1dat2, aes(x = biomass.g_m2)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d1dat2, aes(x = seeds)) +
  geom_density() + facet_grid(sp ~ treatment, scales = "free")

# 2019
ggplot(d2dat2, aes(x = logit_may_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d2dat2, aes(x = logit_jun_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d2dat2, aes(x = logit_jul_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d2dat2, aes(x = logit_early_aug_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d2dat2, aes(x = logit_late_aug_severity)) +
  geom_density() + facet_grid(sp ~ treatment)

ggplot(d2dat2, aes(x = biomass.g_m2)) +
  geom_density() + facet_grid(sp ~ treatment, scales = "free")

ggplot(d2dat2, aes(x = seeds)) +
  geom_density() + facet_grid(sp ~ treatment, scales = "free")


#### stats ####

# 2018
d1Sum <- d1dat2 %>%
  pivot_longer(cols = c(jul_severity:seeds, logit_jul_severity:logit_sep_severity),
               names_to = "response",
               values_to = "values") %>%
  filter(!(sp == "Ev" & response == "biomass.g_m2")) %>%
  select(sp, response) %>%
  unique() %>%
  mutate(hedges_d = map2(response, as.character(sp), hedges_d_fun, year = "2018"),
         t_test = map2(response, as.character(sp), t_test_fun, year = "2018")) %>%
  unnest_wider(hedges_d) %>%
  unnest_wider(t_test)
  
# 2019
d2Sum <- d2dat2 %>%
  pivot_longer(cols = c(seeds, biomass.g_m2:logit_late_aug_severity),
               names_to = "response",
               values_to = "values") %>%
  filter(!(sp == "Mv" & response %in% c("logit_may_severity", "may_severity"))) %>%
  select(sp, response) %>%
  unique() %>%
  mutate(hedges_d = map2(response, as.character(sp), hedges_d_fun, year = "2019"),
         t_test = map2(response, as.character(sp), t_test_fun, year = "2019")) %>%
  unnest_wider(hedges_d) %>%
  unnest_wider(t_test)


#### data for figure ####

# combine data
# remove plots without pairs using dat_format_fun
# rename response
# remove unnecessary responses
figdat <- dat_format_fun(sep_severity, "Mv", 2018) %>%
  pivot_longer(cols = c(fungicide, water),
               names_to = "treatment", 
               values_to = "severity") %>%
  full_join(dat_format_fun(biomass.g_m2, "Mv", 2018) %>%
              pivot_longer(cols = c(fungicide, water),
                           names_to = "treatment", 
                           values_to = "biomass.g_m2")) %>%
  full_join(dat_format_fun(seeds, "Mv", 2018) %>%
              pivot_longer(cols = c(fungicide, water),
                           names_to = "treatment", 
                           values_to = "seeds")) %>%
  mutate(sp = "Mv") %>%
  full_join(dat_format_fun(jul_severity, "Ev", 2018) %>%
              pivot_longer(cols = c(fungicide, water),
                           names_to = "treatment", 
                           values_to = "severity") %>%
              full_join(dat_format_fun(seeds, "Ev", 2018) %>%
                          pivot_longer(cols = c(fungicide, water),
                                       names_to = "treatment", 
                                       values_to = "seeds")) %>%
              mutate(sp = "Ev")) %>%
  mutate(year = 2018) %>%
  full_join(dat_format_fun(early_aug_severity, "Mv", 2019) %>%
              pivot_longer(cols = c(fungicide, water),
                           names_to = "treatment", 
                           values_to = "severity") %>%
              full_join(dat_format_fun(biomass.g_m2, "Mv", 2019) %>%
                          pivot_longer(cols = c(fungicide, water),
                                       names_to = "treatment", 
                                       values_to = "biomass.g_m2")) %>%
              full_join(dat_format_fun(seeds, "Mv", 2019) %>%
                          pivot_longer(cols = c(fungicide, water),
                                       names_to = "treatment", 
                                       values_to = "seeds")) %>%
              mutate(sp = "Mv") %>%
              full_join(dat_format_fun(jun_severity, "Ev", 2019) %>%
                          pivot_longer(cols = c(fungicide, water),
                                       names_to = "treatment", 
                                       values_to = "severity") %>%
                          full_join(dat_format_fun(biomass.g_m2, "Ev", 2019) %>%
                                      pivot_longer(cols = c(fungicide, water),
                                                   names_to = "treatment", 
                                                   values_to = "biomass.g_m2")) %>%
                          full_join(dat_format_fun(seeds, "Ev", 2019) %>%
                                      pivot_longer(cols = c(fungicide, water),
                                                   names_to = "treatment", 
                                                   values_to = "seeds")) %>%
                          mutate(sp = "Ev")) %>%
              mutate(year = 2019)) %>%
  mutate(treatment = recode(treatment, "water" = "control (water)"),
         sp = fct_relevel(sp, "Mv"))

sevSum <- d1Sum %>%
  filter((sp == "Mv" & response == "sep_severity") | (sp == "Ev" & response == "jul_severity")) %>%
  mutate(year = 2018) %>%
  full_join(d2Sum %>%
              filter((sp == "Mv" & response == "early_aug_severity") | (sp == "Ev" & response == "jun_severity")) %>%
              mutate(year = 2019)) %>%
  mutate(sig = case_when(t_p < 0.001 ~ "***",
                         t_p < 0.01 ~ "**",
                         t_p < 0.05 ~ "*",
                         t_p < 0.1 ~ "+",
                         TRUE ~ ""),
         txt_size = case_when(sig == "+" ~ "small",
                               TRUE ~ "large"),
         treatment = "control (water)") %>%
  left_join(figdat %>%
              group_by(sp, year) %>%
              summarise(max_val = max(severity, na.rm = T)))

bioSum <- d1Sum %>%
  filter(response == "biomass.g_m2") %>%
  mutate(year = 2018) %>%
  full_join(d2Sum %>%
              filter(response == "biomass.g_m2") %>%
              mutate(year = 2019)) %>%
  mutate(sig = case_when(t_p < 0.001 ~ "***",
                         t_p < 0.01 ~ "**",
                         t_p < 0.05 ~ "*",
                         t_p < 0.1 ~ "+",
                         TRUE ~ ""),
         txt_size = case_when(sig == "+" ~ "small",
                              TRUE ~ "large"),
         treatment = "control (water)") %>%
  left_join(figdat %>%
              group_by(sp, year) %>%
              summarise(max_val = max(biomass.g_m2, na.rm = T)))

seedSum <- d1Sum %>%
  filter(response == "seeds") %>%
  mutate(year = 2018) %>%
  full_join(d2Sum %>%
              filter(response == "seeds") %>%
              mutate(year = 2019)) %>%
  mutate(sig = case_when(t_p < 0.001 ~ "***",
                         t_p < 0.01 ~ "**",
                         t_p < 0.05 ~ "*",
                         t_p < 0.1 ~ "+",
                         TRUE ~ ""),
         txt_size = case_when(sig == "+" ~ "small",
                              TRUE ~ "large"),
         treatment = "control (water)") %>%
  left_join(figdat %>%
              group_by(sp, year) %>%
              summarise(max_val = max(seeds, na.rm = T)))

#### figure ####

# severity
sevFig <- ggplot(figdat, aes(sp, severity * 100, color = treatment)) +
  geom_point(alpha = 0.5, size = 0.75, position = position_dodge(1)) +
  geom_boxplot(fill = NA, outlier.size = 0.75, size = 0.25, position = position_dodge(1)) +
  facet_grid(~ year, switch = "x") +
  geom_text(data = sevSum, color = "black", aes(x = sp, y = max_val * 100, label = sig, size = txt_size)) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_size_manual(values = c(5, 3)) +
  ylab("Lesions (%)") +
  fig_theme
  
# biomass
bioFig <- figdat %>%
  filter(!(sp == "Ev" & year == 2018)) %>%
  ggplot(aes(sp, log(biomass.g_m2), color = treatment)) +
  geom_point(alpha = 0.5, size = 0.75, position = position_dodge(1)) +
  geom_boxplot(fill = NA, outlier.size = 0.75, size = 0.25, position = position_dodge(1)) +
  facet_grid(~ year, space = "free_x", scales = "free_x", switch = "x") +
  geom_text(data = bioSum, color = "black", aes(x = sp, y = log(max_val), label = sig, size = txt_size)) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_size_manual(values = c(5, 3)) +
  ylab(expression(paste("Biomass (log g ", m^2, ")", sep = ""))) +
  fig_theme

# seeds
seedFig <- ggplot(figdat, aes(sp, log(seeds + 1), color = treatment)) +
  geom_point(alpha = 0.5, size = 0.75, position = position_dodge(1)) +
  geom_boxplot(fill = NA, outlier.size = 0.75, size = 0.25, position = position_dodge(1), show.legend = F) +
  facet_grid(~ year, switch = "x") +
  geom_text(data = seedSum, color = "black", aes(x = sp, y = log(max_val + 1), label = sig, size = txt_size), show.legend = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_size_manual(values = c(5, 3)) +
  ylab("Seeds (log)") +
  fig_theme +
  theme(legend.position = c(0.25, 0.95)) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

# combine
pdf("output/plot_scale_responses_2018_2019_density_exp.pdf", width = 7, height = 3)
plot_grid(sevFig, bioFig, seedFig,
          labels = LETTERS[1:3],
          rel_widths = c(1, 0.75, 1), nrow = 1)
dev.off()