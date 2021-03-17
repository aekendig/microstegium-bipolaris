##### info ####

# file: brms_model_fitting_functions
# author: Amy Kendig
# date last edited: 3/17/21
# goal: functions for general model-fitting routines


#### Priors ####
prior_fun <- function(minX = 0, maxX, shape, scale = 1, dist){
  
  x <- seq(minX, maxX, length.out = 100)
  
  if(dist == "gamma"){
    y <- dgamma(x, shape = shape, scale = scale) # note that this scale is 1/(stan scale)
  } else if(dist == "exponential"){
    y <- dexp(x, shape)
  } else{
    print("unknown distribution")
  }
  
  plot(x, y, type = "l")
  
}


#### Model check ####

mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}


#### Model fit ####

mod_fit_fun <- function(dat, mod, treatCol, ranEff = site, xCol, minX, maxX, yCol){

  outDat <- dat %>%
    select({{treatCol}}) %>%
    unique() %>%
    mutate("{{ranEff}}" := "D5") %>%
    expand_grid(tibble("{{xCol}}" := seq(minX, maxX, length.out = 100))) %>%
    mutate(pred = fitted(mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(mod, newdata = ., allow_new_levels = T)[, "Q97.5"])
  
  outPlot <- ggplot(dat, aes(x = {{xCol}}, y = {{yCol}}, color = {{treatCol}}, fill = {{treatCol}})) +
    geom_ribbon(data = outDat, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_line(data = outDat, aes(y = pred)) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
    stat_summary(geom = "point", fun = "mean", size = 2) +
    theme_bw()
  
  return(list(outDat, outPlot))
  
}
