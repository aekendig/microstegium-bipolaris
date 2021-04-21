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


# import data
d1dat <- read_csv("intermediate-data/plot_biomass_seeds_severity_2018_density_exp.csv")


# copied and pasted from plot_data_processing_2018_density_exp

# format for t-test
mvSevD1DatS <- mvSevD1DatW %>%
  pivot_wider(names_from = treatment,
              values_from = c(jul_severity, late_aug_severity, sep_severity),
              names_glue = "{.value}_{treatment}")

# t-tests
t.test(mvSevD1DatS$jul_severity_fungicide, mvSevD1DatS$jul_severity_water, paired = T)
t.test(mvSevD1DatS$late_aug_severity_fungicide, mvSevD1DatS$late_aug_severity_water, paired = T)
t.test(mvSevD1DatS$sep_severity_fungicide, mvSevD1DatS$sep_severity_water, paired = T)
# all three sig, biggest difference in Sep


mvBioD1DatW <- mvBioD1Dat2 %>%
  pivot_wider(names_from = treatment,
              values_from = biomass.g_m2)

# t-test
t.test(mvBioD1DatW$fungicide, mvBioD1DatW$water, paired = T)
# not significantly different

mvSeedsD1DatW <- mvSeedsD1Dat2 %>%
  pivot_wider(names_from = treatment,
              values_from = seeds)


# t-test
t.test(mvSeedsD1DatW$fungicide, mvSeedsD1DatW$water, paired = T)
# marginally sig

# format for t-test
evSevD1DatS <- evSevD1DatW %>%
  pivot_wider(names_from = treatment,
              values_from = c(jul_severity, late_aug_severity, sep_severity),
              names_glue = "{.value}_{treatment}") %>%
  drop_na()

# t-tests
t.test(evSevD1DatS$jul_severity_fungicide, evSevD1DatS$jul_severity_water, paired = T)
t.test(evSevD1DatS$late_aug_severity_fungicide, evSevD1DatS$late_aug_severity_water, paired = T)
t.test(evSevD1DatS$sep_severity_fungicide, evSevD1DatS$sep_severity_water, paired = T)
# sig lower in July


# make wide
evSeedsD1DatW <- evSeedsD1Dat2 %>%
  pivot_wider(names_from = treatment,
              values_from = seeds) %>%
  drop_na()

# t-tests
t.test(evSeedsD1DatW$fungicide, evSeedsD1DatW$water, paired = T)
# not sig different