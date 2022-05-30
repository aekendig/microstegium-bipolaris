## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)

# Put data and code in a folder together to be grabbed by make_eml
# Generate metadata files by editing current ones or generating them (see Github page for tutorial)
# Edit and run this script

# will need to put code and data in folders within edi folder to run script
# remove them from folders before final step


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)
library(tidyverse)
library(knitr)


#### import templates ####

# list of data files
dlist <- list.files(path = "./data",
                    pattern = ".csv")

# create high level text files
template_core_metadata(path = "./metadata",
                       license = "CCBY")

# create an attribute file for each data table
template_table_attributes(path = "./metadata",
                          data.path = "./data",
                          data.table = dlist)

# create a categorical value file for each data table
template_categorical_variables(path = "./metadata",
                               data.path = "./data")

# look at units
view_unit_dictionary()


#### data file descriptors ####

dlist

# re-arrange files
dlist2 <- c(dlist[4:5], dlist[1:3])

# description list
ddlist <- NA

ddlist[1] <- "biomass of M. vimineum and E. virginicus"
ddlist[2] <- "establishment and foliar lesions of M. vimineum and E. virginicus"
ddlist[3] <- "heights of plants in Appendix S3 experiment"
ddlist[4] <- "foliar lesions on plants in Appendix S3 experiment"
ddlist[5] <- "leaf wetness of plants in Appendix S3 experiment"

# name list
dnlist <- c("Biomass data",
            "Establishment and lesion data",
            "App. S3 height data", 
            "App. S3 lesion data",
            "App. S3 leaf wetness data")

# print table
# dtable <- data.frame(data = dlist2, description = ddlist)
# kable(dtable)


#### code descriptors ####

# list of code files
clist <- c(list.files(path = "./code",
                      pattern = ".R"),
           list.files(path = "./edi",
                      pattern = ".pdf"))

# remove the eml code
clist <- clist[-3]

# re-arrange code list
clist2 <- clist[c(3, 4, 2, 6, 5, 1, 7)]

# code descripions
cdlist <- c(
  "code to process data for establishment and disease incidence analyses",
  "code to analyze establishment and disease incidence",
  "code to analyze biomass data",
  "code to analyze model residuals",
  "code for figures and values in main text",
  "code for analyses and figures in Appendix S3",
  "description of experiment in Appendix S3"
)

# name list
cnlist <- c("Data processing",
            "Establishment and disease incidence analysis",
            "Biomass analysis",
            "Residuals analysis",
            "Figures and values",
            "App. S3 code", 
            "App. S3 experiment description")

# print table
# ctable <- data.frame(code = clist2, desription = cdlist)
# kable(ctable)


#### make eml ####

# copied data and code from the respective folders and put into edi folder

make_eml(path = "./metadata",
         data.path = "./edi",
         dataset.title = "Invasive grass litter suppresses a native grass species and promotes disease",
         data.table = dlist2,
         data.table.name = dnlist,
         data.table.description = ddlist,
         data.table.quote.character = rep("\"", length(dlist)),
         other.entity = clist2,
         other.entity.name = cnlist,
         other.entity.description = cdlist,
         temporal.coverage = c("2018-05-09", "2018-11-04"),
         geographic.description = "Gainesville, FL, USA",
         geographic.coordinates = c(29.64, -82.36, 29.64, -82.36),
         maintenance.description = "completed", 
         user.id = "aekendig",
         user.domain = "EDI",
         package.id = "edi.1019.1")


#### check warnings ####

eml <- EML::read_eml("./metadata/edi.1019.1.xml")
EML::eml_validate(eml)
