##### info ####

# file: nonlinear-functions
# author: Amy Kendig
# date last edited: 6/5/19
# goal: demonstrations of nonlinear functions that can be used to describe species interactions


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)


#### define functions ####

# References
# Law and Watkinson 1987
# Levine and HilleRisLambers 2009
# Hart et al. 2018

# exponential
ex_fun <- function(y0, x, a) {
  y = y0 * exp(a * x)
  return(y)
}

# Beverton-Holt
bh_fun <- function(y0, x, a) {
  y = y0 / (1 + a * x)
  return(y)
}

# Beverton-Holt with exponent
bhe_fun <- function(y0, x, a, b) {
  y = y0 / (1 + a * x)^b
  return(y)
}

# non-linear isoclines
nli_fun <- function(y0, x, a) {
  y = y0 / (1 + x^a)
  return(y)
}

# linear
ln_fun <- function(y0, x, a) {
  y = y0 - a * x
  return(y)
}

# exponential/log
el_fun <- function(y0, x, a) {
  y = y0 * exp(a * log(x + 1))
  return(y)
}

# another Beverton-Holt with exponent
abhe_fun <- function(y0, x, a, b) {
  y = y0 / (1 + (a * x)^b)
  return(y)
}

# another linear
aln_fun <- function(y0, x, a) {
  y = 1 + y0 * (1 - a * x)
  return(y)
}


#### simulate data ####

# initial y
y0 = 10

# functions
funs = tibble(fun = c("abhe", "aln", "bh", "bhe", "el", "ex", "ln", "nli"),
              fun_name = c("Another Beverton-Holt with exponent",
                           "Another linear",
                           "Beverton-Holt",
                           "Beverton-Holt with exponent",
                           "Exponential/log",
                           "Exponential",
                           "Linear",
                           "Non-linear isoclines"),
              fun_form = c("y = y0 / (1 + (a * x)^b)",
                           "y = 1 + y0 * (1 - a * x)",
                           "y = y0 / (1 + a * x)",
                           "y = y0 / (1 + a * x)^b",
                           "y = y0 * exp(a * log(x + 1))",
                           "y = y0 * exp(a * x)",
                           "y = y0 - a * x",
                           "y = y0 / (1 + x^a)"))

# combine with combinations of x and a
sim_dat = expand.grid(0:64, c(-2:2, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1)) %>%
  merge(funs, all = T) %>%
  as_tibble() %>%
  rename(x = Var1, a = Var2)

# simulate y values
# make plot titles
sim_dat2 <- sim_dat %>%
  mutate(y = case_when(fun == "abhe" ~ abhe_fun(y0, x, a, 0.2),
                       fun == "aln" ~ aln_fun(y0, x, a),
                       fun == "bh" ~ bh_fun(y0, x, a),
                       fun == "bhe" ~ bhe_fun(y0, x, a, 0.2),
                       fun == "el" ~ el_fun(y0, x, a),
                       fun == "ex" ~ ex_fun(y0, x, a),
                       fun == "ln" ~ ln_fun(y0, x, a),
                       fun == "nli" ~ nli_fun(y0, x, a)),
         ptitle = paste(fun_name, fun_form, sep = ":\n"),
         rate = as.factor(a))

# make another parameter combination to look more closely at positive values
sim_dat_pos = expand.grid(0:64, c(0, -0.001, -0.005, -0.01, -0.05, -0.1, -0.2, -0.5, -0.7, -1)) %>%
  merge(funs, all = T) %>%
  as_tibble() %>%
  rename(x = Var1, a = Var2)

# simulate y values
# make plot titles
# switch sign for some functions
sim_dat_pos2 <- sim_dat_pos %>%
  mutate(a = case_when(fun %in% c("el", "ex") ~ abs(a),
                       TRUE ~ a),
         y = case_when(fun == "abhe" ~ abhe_fun(y0, x, a, 0.2),
                       fun == "aln" ~ aln_fun(y0, x, a),
                       fun == "bh" ~ bh_fun(y0, x, a),
                       fun == "bhe" ~ bhe_fun(y0, x, a, 0.2),
                       fun == "el" ~ el_fun(y0, x, a),
                       fun == "ex" ~ ex_fun(y0, x, a),
                       fun == "ln" ~ ln_fun(y0, x, a),
                       fun == "nli" ~ nli_fun(y0, x, a)),
         ptitle = paste(fun_name, fun_form, sep = ":\n"),
         rate = as.factor(a))


#### Figures ####

# number of functions
fun_num <- nrow(funs)

# plotting function
plot_fun <- function(plot_name, dat){
  
  pdf(plot_name)
  for(i in 1:fun_num){
    
    temp_dat <- dat %>%
      filter(fun == funs$fun[i])
    
    print(ggplot(temp_dat, aes(x = x, y = y, color = rate)) +
            geom_line() +
            facet_wrap(~ ptitle) +
            theme_bw() + 
            theme(axis.text = element_text(size = 12, color="black"),
                  axis.title = element_text(size = 14),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.text = element_text(size=12),
                  legend.title = element_text(size=12),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 14, hjust = 0.5)))
  }
  dev.off()
}

# make a plot for each one
plot_fun("./output/nonlinear_functions.pdf", sim_dat2)

# make plot for each without extreme values
plot_fun("./output/nonlinear_functions_trimmed.pdf", sim_dat2 %>%
           filter(y < 50 & y > 0))

# make plot for each without extreme values
plot_fun("./output/nonlinear_functions_positive_trimmed.pdf", sim_dat_pos2 %>%
           filter(y < 50 & y > 0))

