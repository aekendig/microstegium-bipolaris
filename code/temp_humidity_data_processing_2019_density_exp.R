##### info ####

# file: temp_humidity_data_processing_2019_density_exp
# author: Amy Kendig
# date last edited: 1/19/20
# goal: process data from temperature and humdity loggers


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(caTools)


#### import data ####

# make a list of the file names
file_names = list.files(path = "./data/temp-humidity-litter-exp-20191025", 
                        pattern = "[.]csv")

# Import files, all in one list
dat = lapply(paste("./data/temp-humidity-litter-exp-20191025/", file_names, sep = ""), read.csv, header=T)


#### edit data ####

# assign plot to files
for(i in 1:length(file_names)){
  dat[[i]]$id <- sub(".csv", "", file_names[i])
}

# make list into data table
# add and remove columns
# remove data before placement and after collection
# sensors were placed by 3:00 pm on 7/2/19 and collected after 10:00 am on 10/22/19
dat2 <- do.call(rbind.data.frame, dat) %>%
  as_tibble() %>%
  rename(temp = "Temp...C..c.1", rel_hum = "RH.....c.1.2") %>%
  mutate(time = as.POSIXct(Date.Time, format = "%m/%d/%y %H:%M:%S"),
         day = as.Date(Date.Time, format = "%m/%d/%y %H:%M:%S"),
         plot = substr(id, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.factor()) %>%
  select(-c(Date.Time)) %>%
  filter(time >= "2019-07-03 00:00:00" & time <= "2019-10-21 23:59:59")

# daily summaries
day_dat <- dat2 %>%
  group_by(id, plot, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(rel_hum),
            hum_min = min(rel_hum),
            hum_max = max(rel_hum))


#### visualizations ####

ggplot(day_dat, aes(x = day, y = temp_avg, color = plot)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0)
  
ggplot(day_dat, aes(x = day, y = hum_avg, color = plot)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0)

ggplot(filter(dat2, time > "2019-08-01 01:00:00" & time < "2019-08-05 01:00:00"), aes(x = time, y = temp, color = plot)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0)


#### save data #### 
write_csv(dat2, "./intermediate-data/temp_humidity_hourly_2019_density_exp.csv")
write_csv(day_dat, "./intermediate-data/temp_humidity_daily_2019_density_exp.csv")
