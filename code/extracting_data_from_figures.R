#### Goal: extract data from figures

# reference: Emmanuel Jjunju, https://www.r-bloggers.com/digitizing-jpeg-graphs-in-r/


#### Set-up ####

# clear all existing data
rm(list=ls())

# libraries
library(jpeg) # to use jpeg images in R
library(zoo)
library(tidyverse)

# digitize functions
ReadAndCal = function(fname){
  ReadImg(fname)
  calpoints <- locator(n=4,type='p',pch=4,col='blue',lwd=2)
  return(calpoints)
}

ReadImg = function(fname){
  img <- readJPEG(fname)
  op <- par(mar=c(0,0,0,0))
  on.exit(par(op))
  plot.new()
  rasterImage(img,0,0,1,1)
}

DigitData = function(col='red',type='p',...){
  type <- ifelse(type=='b','o',type)
  type <- ifelse(type%in%c('l','o','p'),type,'p')
  locator(type=type,col=col,...)
}

Calibrate = function(data,calpoints,x1,x2,y1,y2){
  x 		<- calpoints$x[c(1,2)]
  y 		<- calpoints$y[c(3,4)]
  
  cx <- lm(formula = c(x1,x2) ~ c(x))$coeff
  cy <- lm(formula = c(y1,y2) ~ c(y))$coeff
  
  data$x <- data$x*cx[2]+cx[1]
  data$y <- data$y*cy[2]+cy[1]
  
  return(as.data.frame(data))
}


#### Steps ####

# ReadAndCal opens the jpeg in a plotting window and lets you define points on the x and y axes. You must start by clicking on the left-most x-axis point, then the right-most axis point, followed by the lower y-axis point and finally the upper y-axis point. You don’t need to choose the end points of the axis, only two points on the axis that you know the x or y value for. As you click on each of the 4 points, the coordinates are saved in the object cal.

# DigitData returns you to the figure window, and now you can click on each of the data points you’re interested in retrieving values for. The function will place a dot (colored red in this case) over each point you click on, and the raw x,y coordinates of that point will be saved to the data.points list. When you’re finished clicking points, you need to hit stop/Finish or right-click to stop the data point collection.

# Calibrate converts those raw x,y coordinates into the same scale as the original graph. Feed the function your data.point list, the ReadAndCal list that contains your 4 control points from the first step, and then 4 numeric values that represent the 4 original points you clicked on the x and y axes. These values should be in the original scale of the figure (i.e. read the values off the graph’s tick marks).

#### Digitize Flory et al. 2011 figures ####

# Fig. 2A
(cal_fl11_2a = ReadAndCal("data/lit-figures/Flory_2011_2A.jpg"))
(data_fl11_2a = DigitData(col = 'red'))
df_fl11_2a = Calibrate(data_fl11_2a, cal_fl11_2a, 0, 1, 0, 30)
df_fl11_2a$site = rep(c("1", "2", "3"), each = 4)
df_fl11_2a$measure = rep(c("mean", "mean_se"), 6)
df_fl11_2a$treatment = rep(rep(c("control", "fungicide"), each = 2), 3)

# Fig. 2C
(cal_fl11_2c = ReadAndCal("data/lit-figures/Flory_2011_2C.jpg"))
(data_fl11_2c = DigitData(col = 'red'))
df_fl11_2c = Calibrate(data_fl11_2c, cal_fl11_2c, 0, 1, 0, 30)
df_fl11_2c$site = rep(c("1", "2", "3"), each = 4)
df_fl11_2c$measure = rep(c("mean", "mean_se"), 6)
df_fl11_2c$treatment = rep(rep(c("control", "fungicide"), each = 2), 3)

# Fig. 2D
(cal_fl11_2d = ReadAndCal("data/lit-figures/Flory_2011_2D.jpg"))
(data_fl11_2d = DigitData(col = 'red'))
df_fl11_2d = Calibrate(data_fl11_2d, cal_fl11_2d, 0, 1, 0, 70)
df_fl11_2d$site = rep(c("1", "2", "3"), each = 4)
df_fl11_2d$measure = rep(c("mean", "mean_se"), 6)
df_fl11_2d$treatment = rep(rep(c("control", "fungicide"), each = 2), 3)

# combine
df_fl11_2a2 <- df_fl11_2a %>%
  select(site, treatment, measure, y) %>%
  pivot_wider(names_from = "measure",
              values_from = "y") %>%
  mutate(severity_se = mean_se - mean) %>%
  rename(severity = mean) %>%
  select(-mean_se)

df_fl11_2c2 <- df_fl11_2c %>%
  select(site, treatment, measure, y) %>%
  pivot_wider(names_from = "measure",
              values_from = "y") %>%
  mutate(biomass_g_m2_se = mean_se - mean) %>%
  rename(biomass_g_m2 = mean) %>%
  select(-mean_se)

df_fl11_2d2 <- df_fl11_2d %>%
  select(site, treatment, measure, y) %>%
  pivot_wider(names_from = "measure",
              values_from = "y") %>%
  mutate(seed_heads_se = mean_se - mean) %>%
  rename(seed_heads = mean) %>%
  select(-mean_se)

df_fl11_2 <- df_fl11_2a2 %>%
  full_join(df_fl11_2c2) %>%
  full_join(df_fl11_2d2)


#### Digitize Stricker et al. 2016 figures ####

# Fig. 4AB
(cal_st16_4a = ReadAndCal("data/lit-figures/Stricker_2016_4AB.jpg"))
(data_st16_4a = DigitData(col = 'red'))
df_st16_4a = Calibrate(data_st16_4a, cal_st16_4a, 0, 1, 0, 10)
df_st16_4a$year = rep(c(2012, 2013), each = 16)
df_st16_4a$site = rep(rep(c("1", "2", "3", "4"), each = 4), 2)
df_st16_4a$treatment = rep(rep(c("control", "fungicide"), each = 2), 8)
df_st16_4a$measure = rep(c("mean", "mean_se"), 8)

# Fig. 4CD
(cal_st16_4c = ReadAndCal("data/lit-figures/Stricker_2016_4CD.jpg"))
(data_st16_4c = DigitData(col = 'red'))
df_st16_4c = Calibrate(data_st16_4c, cal_st16_4c, 0, 1, 0, 250)
df_st16_4c$year = rep(c(2012, 2013), each = 16)
df_st16_4c$site = rep(rep(c("1", "2", "3", "4"), each = 4), 2)
df_st16_4c$treatment = rep(rep(c("control", "fungicide"), each = 2), 8)
df_st16_4c$measure = rep(c("mean", "mean_se"), 8)

# Fig. 4EF
(cal_st16_4e = ReadAndCal("data/lit-figures/Stricker_2016_4EF.jpg"))
(data_st16_4e = DigitData(col = 'red'))
df_st16_4e = Calibrate(data_st16_4e, cal_st16_4e, 0, 1, 0, 2000)
df_st16_4e$year = rep(c(2012, 2013), each = 16)
df_st16_4e$site = rep(rep(c("1", "2", "3", "4"), each = 4), 2)
df_st16_4e$treatment = rep(rep(c("control", "fungicide"), each = 2), 8)
df_st16_4e$measure = rep(c("mean", "mean_se"), 8)

# combine
df_st16_4a2 <- df_st16_4a %>%
  select(year, site, treatment, measure, y) %>%
  pivot_wider(names_from = "measure",
              values_from = "y") %>%
  mutate(severity_se = mean_se - mean) %>%
  rename(severity = mean) %>%
  select(-mean_se)

df_st16_4c2 <- df_st16_4c %>%
  select(year, site, treatment, measure, y) %>%
  pivot_wider(names_from = "measure",
              values_from = "y") %>%
  mutate(biomass_g_m2_se = mean_se - mean) %>%
  rename(biomass_g_m2 = mean) %>%
  select(-mean_se)

df_st16_4e2 <- df_st16_4e %>%
  select(year, site, treatment, measure, y) %>%
  pivot_wider(names_from = "measure",
              values_from = "y") %>%
  mutate(seed_heads_se = mean_se - mean) %>%
  rename(seed_heads = mean) %>%
  select(-mean_se)

df_st16_4 <- df_st16_4a2 %>%
  full_join(df_st16_4c2) %>%
  full_join(df_st16_4e2)


#### Save values ####

write_csv(df_fl11_2, "intermediate-data/flory_2011_extracted_figure.csv")

write_csv(df_st16_4, "intermediate-data/stricker_2016_extracted_figure.csv")
