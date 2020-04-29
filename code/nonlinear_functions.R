##### info ####

# file: nonlinear_functions
# author: Amy Kendig
# date last edited: 4/22/20
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
# Liermann and Hilborn 2001 (Eq. 26)

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

# sigmoid Beverton-Holt
sbh_fun <- function(x, r, a, d){
  y = x^(d-1) / (1/r + a * x^d)
  return(y)
}
# a is defined as 1/asymptote for total y (not per capita)
# r and d contribute to the growth rate


#### simulate data ####

# y with no competitors
y0 = 10

# intrinsic growth rate/maximum productivity
r = 0.5

# non-linearity of competition
d = 2

# functions
funs = tibble(fun = c("abhe", "aln", "bh", "bhe", "el", "ex", "ln", "nli", "sbh"),
              fun_name = c("Another Beverton-Holt with exponent",
                           "Another linear",
                           "Beverton-Holt",
                           "Beverton-Holt with exponent",
                           "Exponential/log",
                           "Exponential",
                           "Linear",
                           "Non-linear isoclines",
                           "Sigmoid Beverton-Holt"),
              fun_form = c("y = y0 / (1 + (a * x)^b)",
                           "y = 1 + y0 * (1 - a * x)",
                           "y = y0 / (1 + a * x)",
                           "y = y0 / (1 + a * x)^b",
                           "y = y0 * exp(a * log(x + 1))",
                           "y = y0 * exp(a * x)",
                           "y = y0 - a * x",
                           "y = y0 / (1 + x^a)",
                           "y = x^(d-1) / (1/r + a * x^d)"))

# combine with combinations of x and a
sim_dat = expand.grid(0:64, c(-2:2, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1)) %>%
  merge(funs, all = T) %>%
  as_tibble() %>%
  rename(x = Var1, a = Var2) %>%
  mutate(r = r,
         d = d)

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
                       fun == "nli" ~ nli_fun(y0, x, a),
                       fun == "sbh" ~ sbh_fun(x, r, a, d)),
         ptitle = paste(fun_name, fun_form, sep = ":\n"),
         rate = as.factor(a))

# make another parameter combination to look more closely at positive values
sim_dat_pos = expand.grid(0:64, c(0, -0.001, -0.005, -0.01, -0.05, -0.1, -0.2, -0.5, -0.7, -1)) %>%
  merge(funs, all = T) %>%
  as_tibble() %>%
  rename(x = Var1, a = Var2) %>%
  mutate(r = r,
         d = d)

# simulate y values
# make plot titles
# switch sign for some functions
sim_dat_pos2 <- sim_dat_pos %>%
  mutate(a = case_when(fun %in% c("el", "ex", "sbh") ~ abs(a),
                       TRUE ~ a),
         y = case_when(fun == "abhe" ~ abhe_fun(y0, x, a, 0.2),
                       fun == "aln" ~ aln_fun(y0, x, a),
                       fun == "bh" ~ bh_fun(y0, x, a),
                       fun == "bhe" ~ bhe_fun(y0, x, a, 0.2),
                       fun == "el" ~ el_fun(y0, x, a),
                       fun == "ex" ~ ex_fun(y0, x, a),
                       fun == "ln" ~ ln_fun(y0, x, a),
                       fun == "nli" ~ nli_fun(y0, x, a),
                       fun == "sbh" ~ sbh_fun(x, r, a, d)),
         ptitle = paste(fun_name, fun_form, sep = ":\n"),
         rate = as.factor(a))


#### figures ####

# plotting function
plot_fun <- function(plot_name, dat){
  
  # functions
  dat_funs <- unique(dat$fun)
  
  # number of functions
  fun_num <- length(dat_funs)
  
  pdf(plot_name)
  for(i in 1:fun_num){
    
    temp_dat <- dat %>%
      filter(fun == dat_funs[i])
    
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

plot_fun("./output/nonlinear_functions_positive.pdf", sim_dat_pos2)

# make plot for each without extreme values
plot_fun("./output/nonlinear_functions_trimmed.pdf", sim_dat2 %>%
           filter(y < 50 & y > 0))

plot_fun("./output/nonlinear_functions_positive_trimmed.pdf", sim_dat_pos2 %>% filter(y < 50 & y > 0))


#### simulate only the sigmoid Beverton-Holt ####

# combine with combinations of x and a
sim_dat_sbh = expand_grid(x = 0:64, a = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5)) %>%
  expand_grid(r = c(0.1, 0.5, 1, 5)) %>%
  expand_grid(d = 1:4) %>%
  mutate(y = sbh_fun(x, r, a, d),
         rate = as.factor(a),
         growth = as.factor(r),
         nonlin = as.factor(d),
         ptitle = paste("exponent: ", d, sep = ""))

# number of figures
pages_sbh <- sim_dat_sbh %>% select(d) %>% unique()
fig_sbh <- nrow(pages_sbh)

# figure
pdf("./output/nonlinear_functions_sigmoid_bev_holt.pdf")
for(i in 1:fig_sbh){
  
  # select temporary data
  temp_dat <- filter(sim_dat_sbh, d == pages_sbh$d[i])
  
  # figure
  print(ggplot(temp_dat, aes(x = x, y = y, color = rate)) +
          geom_line() +
          facet_grid(growth ~ ptitle) +
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