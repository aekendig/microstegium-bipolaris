##### info ####

# file: annual_perennial_troubleshooting_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: code cut from other scripts used to troubleshoot


#### annual_perennial_relative_abundance_2019_density_exp ####

# perennial abundance cycles -- why?

# first iteration
abund_dat %>%
  filter(iteration == 1) %>%
  ggplot(aes(time, N, color = species)) +
  geom_line()

param_dat %>%
  filter(iteration == 1) %>%
  ggplot(aes(time, value)) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free")
# perennial adult growing season survival crashed to zero

# find times when survival crashed
filter(param_dat, iteration == 1 & parameter == "U.P" & value < 1)

# look at densities just before and after
filter(abund_dat, iteration == 1 & time %in% c(303:306, 586:589)) %>%
  data.frame()

# perennial survival over a range of densities
tibble(P_dens = 26354.84801,
       A_dens = c(868, 867, 866, 865, 864, 863, 862, 861, 860),
       S_dens = 1156.05843) %>%
  rowwise() %>%
  mutate(U.P = U_P_fun(0, A_dens, S_dens, P_dens, 1),
         U.P.lin = U_P_dens$int_fun[1] + U_P_dens$mv_dens_fun[1] * A_dens + U_P_dens$evS_dens_fun[1] * S_dens + U_P_dens$evA_dens_fun[1] * P_dens,
         U.P.raw = exp(U.P.lin)/(1 + exp(U.P.lin))) %>%
  ggplot(aes(A_dens, U.P)) +
  geom_line()
# the decrease in Microstegium abundance reduced Perennial survival

# is this sensitive to the perennial values?
tibble(A_dens = c(868, 867, 866, 865, 864, 863, 862, 861, 860),
       S_dens = 1156.05843) %>%
  expand_grid(tibble(P_dens = c(26232.93325, 26234, 26235, 26236, 26354.84801))) %>%
  rowwise() %>%
  mutate(U.P = U_P_fun(0, A_dens, S_dens, P_dens, 1)) %>%
  ggplot(aes(A_dens, U.P, color = as.factor(P_dens))) +
  geom_line()
# much more sensitive to the Microstegium abundance because the coefficient is so high:
U_P_dens[1, ]

# why are densities so high?
# compare with iteration 21 of the simulations with disease, perennial abundance asymptotes
param_dat %>%
  filter(iteration == 1 & time > 100) %>%
  full_join(param_dat2 %>%
              filter(iteration == 21 & time > 100)) %>%
  ggplot(aes(time, value, color = as.factor(iteration))) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free")
# over-winter survival of perennial adults is slightly lower

param_dat %>%
  filter(time == 1 & parameter == "w.P") %>%
  mutate(disease = "without disease",
         w.P = round(value, digits = 3)) %>%
  select(iteration, disease, w.P) %>%
  full_join(param_dat2 %>%
              filter(time == 1 & parameter == "w.P") %>%
              mutate(disease = "with disease",
                     w.P = round(value, digits = 3)) %>%
              select(iteration, disease, w.P)) %>%
  arrange(iteration, disease) %>%
  data.frame()
# reduce survival

param_dat %>%
  filter(iteration == 20) %>%
  ggplot(aes(time, value)) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free")
# U.P drops

param_dat %>%
  filter(iteration == 53) %>%
  ggplot(aes(time, value)) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free")

param_dat %>%
  filter(iteration == 53 & parameter == "U.P" & value < 0.75) %>%
  filter(time == min(time))

param_dat %>%
  filter(iteration == 53 & time > 342 & time < 346) %>%
  ggplot(aes(time, value)) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free")
# U.P drops

param_dat %>%
  filter(iteration == 54) %>%
  ggplot(aes(time, value)) +
  geom_line() +
  facet_wrap(~ parameter, scales = "free")
# not sure what is causing fluctuations in annual population