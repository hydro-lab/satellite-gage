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

## Mutale River

# DWS gage
# Source: https://www.dws.gov.za/Hydrology/Verified/HyDataSets.aspx?Station=A9H029&SiteDesc=RIV
# Date downloaded: 2022 Oct 31
gage <- read_csv("mutale_weir_gage.csv")

gage <- gage %>%
     mutate(dt=ymd(Date)) %>%
     rename(gage=`Q (m^3/s)`) %>%
     select(-QC,-Date) %>%
     pivot_longer(gage,names_to = "measurement")

# Satellite Gage
# Passed from width2flow.R
satellite <- read_csv("sat_gage.csv")

satellite <- satellite %>%
     mutate(dt=as_date(dt)) %>%
     filter(is.na(width)==FALSE) %>%
     rename(satellite=q)

dt <- array(NA,dim = nrow(satellite))
width <- dt
planet <- dt
j <- 1
dt[1] <- satellite$dt[1]
width[1] <- satellite$width[1]
planet[1] <- satellite$satellite[1]
for (i in 2:nrow(satellite)) {
     if (satellite$dt[i] != dt[j]) {
          j <- j + 1
          dt[j] <- satellite$dt[i]
          width[j] <- satellite$width[i]
          planet[j] <- satellite$satellite[i]
     }
}
dt <- dt[1:j]
width <- width[1:j]
planet <- planet[1:j]
satellite <- data.frame(dt,width,planet)
satellite <- satellite %>%
     mutate(dt=as_date(dt)) %>%
     pivot_longer(c(width,planet),names_to = "measurement")

## COMPARE

# Temporal comparison
discharge <- rbind(gage,satellite)
ggplot(discharge) +
     geom_point(aes(x=dt,y=value,color=measurement)) +
     xlab("Date") +
     ylab(TeX(r'(Discharge ($m^3/s$))')) + # error here is incorrect
     ylim(c(0,15)) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(legend.key = element_rect(fill = "white")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

# One-to-one comparison
discharge_wide <- discharge %>%
     pivot_wider(names_from = "measurement",values_from = "value") %>%
     filter(width<30, width>=0)
ggplot(discharge_wide) +
     geom_point(aes(x=gage,y=planet)) +
     #coord_cartesian(xlim = c(0,2), ylim = c(0,2)) +
     xlab(TeX(r'(Gage Discharge ($m^3/s$))')) +
     ylab(TeX(r'(Satellite Discharge ($m^3/s$))')) +
     xlim(c(0,15)) +
     ylim(c(0,15)) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))
cc <- cor(discharge_wide$gage,discharge_wide$planet)
ggplot(discharge_wide) +
     geom_point(aes(x=gage,y=width)) +
     xlab(TeX(r'(Gage Discharge ($m^3/s$))')) +
     ylab("Satellite Width (m)") +
     xlim(c(0,15)) +
     ylim(c(0,30)) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))
cw <- cor(discharge_wide$gage,discharge_wide$width)



