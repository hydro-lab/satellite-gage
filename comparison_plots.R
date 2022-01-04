# Code to compare the computed discharges to the standard gages.

# Written in Matlab and R with RStudio by David Kahler and Mackenzie Martin, 
# Duquesne University, from 2018 to 2021.  The development was supported by 
# the United States Agency for International Development, Southern Africa 
# Regional Mission.  Further information is available at: 
# www.duq.edu/limpopo 
# https://github.com/LimpopoLab 

library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

# BUFFALO CREEK
## USGS Gage

gage <- read_csv("buffalo_usgs_gage.csv")

gage$dt <- gage$datetime
for (i in 1:nrow(gage)) {
     if (gage$tz_cd[i] == "EST") {
          gage$dt[i] <- force_tz(gage$datetime[i], tzone = "EST") # converts to EST
     } else {
          gage$dt[i] <- force_tz(gage$datetime[i], tzone = "EST") - 3600 # converts to EST, even with EDT
     }
}

gage <- gage %>%
     mutate(discharge=as.numeric(discharge_cfs)*0.028316847) %>% # from cfs to m^3/s
     mutate(dt_est=dt-(5*3600)) %>% # forces to EST, will now appear as UTC
     select(dt_est,discharge)

daily <- gage %>%
     mutate(day=as_date(dt_est)) %>%
     group_by(day) %>%
     summarize(daily=mean(discharge, na.rm = TRUE))

## Satellite Gage

