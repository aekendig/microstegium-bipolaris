##### info ####

# file: temp_humidity_data_processing_2019_density_exp
# author: Amy Kendig
# date last edited: 8/20/20
# goal: process data from temperature and humdity loggers


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(caTools)


#### import data ####

# make a list of the file names
file_names = list.files(path = "./data/temp-humidity-density-exp-20191025", 
                        pattern = "[.]csv")

# Import files, all in one list
dat = lapply(paste("./data/temp-humidity-density-exp-20191025/", file_names, sep = ""), read.csv, header=T)


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
  mutate(site = substr(id, 1, 2),
         plot = substr(id, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.factor(),
         treatment = "water",
         time = case_when(site == "D3" & plot == 6 ~ as.POSIXlt(Date.Time, format = "%m/%d/%y %H:%M"),
                          TRUE ~ as.POSIXlt(Date.Time, format = "%m/%d/%y %H:%M:%S")),
         day = case_when(site == "D3" & plot == 6 ~ as.Date(Date.Time, format = "%m/%d/%y %H:%M"),
                         TRUE ~ as.Date(Date.Time, format = "%m/%d/%y %H:%M:%S")),
         hum_prop = rel_hum / 100) %>%
#  select(-c(Date.Time)) %>%
  filter(time >= "2019-07-03 00:00:00" & time <= "2019-10-21 23:59:59")

# check days and times
dat2 %>%
  group_by(day) %>%
  summarise(min_time = min(time),
            max_time = max(time))

# check for missing data
filter(dat2, is.na(temp) | is.na(hum_prop)) 
filter(dat2, site == "D3" & plot == 6 & day == as.Date("2019-07-23")) %>% data.frame()
filter(dat2, site == "D3" & plot == 6 & day == as.Date("2019-09-12")) %>% data.frame()
# one hour missing from each day

# remove missing data
dat3 <- dat2 %>%
  filter(!is.na(temp) & !is.na(hum_prop))

# daily summaries
day_dat <- dat3 %>%
  group_by(id, site, plot, treatment, day) %>%
  summarise(temp_avg = mean(temp),
            temp_min = min(temp),
            temp_max = max(temp),
            hum_avg = mean(hum_prop),
            hum_min = min(hum_prop),
            hum_max = max(hum_prop))


#### visualizations ####

ggplot(day_dat, aes(x = day, y = temp_avg, color = site)) +
  geom_point() +
  facet_wrap(~ plot)
  
ggplot(day_dat, aes(x = day, y = hum_avg, color = site)) +
  geom_point() +
  facet_wrap(~ plot)


#### save data #### 
write_csv(dat3, "./intermediate-data/temp_humidity_hourly_2019_density_exp.csv")

