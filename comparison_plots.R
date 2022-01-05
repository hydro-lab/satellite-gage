# Code to compare the computed discharges to the standard gages.

# Written in Matlab and R with RStudio by David Kahler and Mackenzie Martin, 
# Duquesne University, from 2018 to 2021.  The development was supported by 
# the United States Agency for International Development, Southern Africa 
# Regional Mission.  Further information is available at: 
# www.duq.edu/limpopo 
# https://github.com/LimpopoLab 

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(latex2exp)

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

gage <- gage %>%
     mutate(dt=as_date(dt_est)) %>%
     mutate(Measurement="Gage") %>%
     group_by(dt) %>%
     summarize(q=mean(discharge, na.rm = TRUE))

## Satellite Gage
satellite <- read_csv("sat_gage.csv")

satellite <- satellite %>%
     filter(is.na(width)==FALSE) %>%
     mutate(Measurement="Satellite") %>%
     select(-width)

for (i in 1:nrow(satellite)) {
     if (is.na(satellite$q[i])==FALSE) {
          if (satellite$q[i]<0) {
               satellite$q[i] <- NA
          }
     }
}

## COMPARE
gage2 <- gage %>%
     filter(dt>="2017-05-01") %>%
     filter(dt<="2017-10-30")
satellite2 <- satellite %>%
     filter(dt>="2017-05-01") %>%
     filter(dt<="2017-10-30")
discharge <- rbind(gage2,satellite2)
discharge <- rename(discharge,Measurement=Source) # to fix old naming convention, comment if unneeded

# Temporal comparison
ggplot(discharge) +
     geom_point(aes(x=dt,y=q,color=Measurement)) +
     xlab("Date (2017)") +
     ylab(TeX(r'(Discharge ($m^3/s$))')) + # error here is incorrect
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

# One-to-one comparison
discharge2 <- discharge %>%
     pivot_wider(names_from = Measurement, values_from = q) %>%
     filter(is.na(Gage)==FALSE) %>%
     filter(is.na(Satellite)==FALSE)
ggplot(discharge2) +
     geom_point(aes(x=Gage,y=Satellite)) +
     coord_cartesian(xlim = c(0,2), ylim = c(0,2)) +
     xlab("Gage Measurement") +
     ylab("Satellite Measurement") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))
cc <- cor(discharge2$Gage,discharge2$Satellite)
