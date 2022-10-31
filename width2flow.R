# Code to convert the widths from the Planet Labs images and compute the 
# discharge based on Manning's equation and the bathymetric profile.

# Written in Matlab and R with RStudio by David Kahler and Mackenzie Martin, 
# Duquesne University, from 2018 to 2021.  The development was supported by 
# the United States Agency for International Development, Southern Africa 
# Regional Mission.  Further information is available at: 
# www.duq.edu/limpopo 
# https://github.com/LimpopoLab 

library(readr)
library(lubridate)
library(ggplot2)
library(latex2exp)
library(dplyr)

# Remember to set working directory
setwd("/Volumes/T7/planet/mutale/width2flow/")

# Channel-specific information
calibration_discharge <- 0.28 # discharge, cubic meters per second
calibration_width <- 4.5 # width, meters
S_0 <- 0.002 # slope, dimensionless

profile <- read_csv("profile.csv") # (2) $location_m and $height_m
if (mean(profile$height_m)>0){ # guesses if the average height is positive, that the height vector points down
     profile$height_m = -profile$height_m # code is based on the height vector to point up
}
p <- ggplot(profile) +
     geom_line(aes(x = location_m, y = height_m)) +
     xlab("Cross-stream distance (m)") +
     ylab("Bathymetry (m)") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))
ggsave("profile.eps", p, device = "eps", dpi = 72)

min_depth <- min(profile$height_m)
W <- profile$location_m
d <- profile$height_m

# Search for the minimum depth that will still generate a single channel (avoid braiding), first identify actual minimum
for (i in 2:(length(W)-1)){
     if (min_depth==d[i]){
          min_depth_index <- i
     }
}
# overwrite if there is a peak, if none, it will not overwrite
for (i in 2:(length(W)-1)){
  if ((d[i]>min_depth)&(d[i-1]<d[i])&(d[i]>d[i+1])){
    min_depth=d[i]
    min_depth_index=i
  }
}
# Find the endpoints of the lowest stage possible in analysis.  First, prepopulate if there it is at the bottom
min_width_1 <- W[min_depth_index]
min_width_index_1 <- min_depth_index
min_width_2 <- W[min_depth_index]
min_width_index_2 <- min_depth_index
for (i in 2:(min_depth_index)){ # approach from index = 1
  if ((d[i-1]>=min_depth)&(min_depth>d[i])){
    min_width_1=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(min_depth-d[i-1])
    min_width_index_1=i
  }
}
for (i in (min_depth_index+2):(length(W))){ # continue from the minimum physical depth to the other side
  if ((d[i-1]<=min_depth)&(min_depth<d[i])){
    min_width_2=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(min_depth-d[i-1])
    min_width_index_2=i-1
  }
}
min_width <- (min_width_2-min_width_1) # minimum width as to avoid braiding

# Find the maximum width that can be used in the code, that is, for which we have bathymetry data.  Again, prepopulate:
max_width_1 <- W[1] # "top" of the channel
max_width_2 <- W[length(W)] # top on the other side
max_depth <- d[1]
max_width_index_1 <- 1
max_width_index_2 <- length(W)
# In the event that the endpoints do not reach the same datum (ideally, they will be at the same height):
# Scenario 1: depth at near bank is lower, adjust far bank
if ((d[1]<d[length(W)])){
  for (i in (min_depth_index_2+1):length(W)){
    if ((d[i-1]<d[1])&(d[1]<d[i])){
      max_width_2=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(d[1]-d[i-1]);
      max_width_index_2=i-1;
    }
  }
}
# Scenario 2: depth at the far bank is lower, adjust near bank
if (d[1]>(d[length(W)])){ 
  max_depth=d[length(W)];
  for (i in (2:(min_width_index_1))){ 
    if ((d[i-1]>d[length(W)])&(d[length(W)]>d[i])){
      max_width_1=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(d[length(W)]-d[i-1]);
      max_width_index_1=i;
    }
  }
}
max_width = max_width_2-max_width_1;

# In order to construct a bathymetry, we need points on both sides of the channel to find level widths.  This will construct those points on each side.
new_far <- array(-9999, dim = c(((min_width_index_1)-(max_width_index_1)+1),2)) # uses the number of "near side" (which start at i = 1) points
new_near <- array(-9999, dim = c(((max_width_index_2)-(min_width_index_2)+1),2)) # uses the number of "far side" points
# Search along the "near" or i=1 to min
for (i in (max_width_index_1):(min_width_index_1-1)){
  for (j in ((min_width_index_2+1):(length(W)))){
    if (((d[j-1]<d[i]))&((d[i]<d[j]))){
      new_far[i-max_width_index_1+1,1]=W[j-1]+(d[i]-d[j-1])*(W[j]-W[j-1])/(d[j]-d[j-1]);
      new_far[i-max_width_index_1+1,2]=d[i];
    }
  }
}
new_far[min_width_index_1-max_width_index_1+1,1] <- min_width_2;
new_far[min_width_index_1-max_width_index_1+1,2] <- min_depth;
# Search along the "far" or min to end of index
for (i in ((min_width_index_2+1):max_width_index_2)){
  for (j in 2:min_width_index_1){
    if (((d[j-1]>d[i]))&((d[i]>d[j]))){
      new_near[i-min_width_index_2,1]=W[j-1]+(d[i]-d[j-1])*(W[j]-W[j-1])/(d[j]-d[j-1]);
      new_near[i-min_width_index_2,2]=d[i];
    }
  }
}
new_near[max_width_index_2+1-min_width_index_2,1] <- min_width_1;
new_near[max_width_index_2+1-min_width_index_2,2] <- min_depth;

# Form new list of positions
prof <- cbind(W, d)
new <- rbind(new_near, new_far)
complet <- rbind(prof, new)
temp <- complet[order(complet[,1]),]
rm(prof,new,complet)
rm(new_near,new_far,W,d)

head <- 0
foot <- 0
rpt <- array(0, dim = nrow(temp))
for (i in 1:nrow(temp)){
     if (i >= 2) {
          if (temp[i,1]==temp[(i-1),1]) { # checks for repeats
               rpt[i] <- rpt[i] + 1
          }
     }
     if ((temp[i,1])<max_width_1){ # checks for anything that is positioned before the maximum width position at the near end
          head=head+1;
     }
     if ((temp[i,1])<=max_width_2){ # counts everything that is positioned before or at the maximum width position at the far end
          foot=foot+1;
     }
}

xsec <- temp[(head+1):foot,]
rpt <- rpt[(head+1):foot]
xsec <- data.frame(xsec[,1],xsec[,2],rpt)
xsec <- xsec %>%
     rename(location_m=xsec...1.,height_m=xsec...2.) %>%
     filter(rpt==0) %>%
     select(-rpt)
rm(temp,rpt)

# Create a table of width intervals (column 4)
levels <- array(-9999, dim = c(nrow(xsec),4))
for (i in 1:(nrow(xsec)-1)) {
     levels[i,1]=xsec$height_m[i] # depth
     levels[i,2]=xsec$location_m[i] # Near bank position
     for (j in i:nrow(xsec)) { # starts at i in case it is a singular minimum point
          if ((xsec$height_m[i])==(xsec$height_m[j])) {
               levels[i,3] <- xsec$location_m[j] # Far bank position
               levels[i,4] <- xsec$location_m[j]-xsec$location_m[i] # width, which will provide the intervals
          }
     }
     if (xsec$height_m[i]<=min_depth) { # at, or below, the lowest
          break
     }
}
levels <- levels[1:i,]

if ((calibration_width>min_width)&&(calibration_width<max_width)){ # check valid calibration values
  for (i in 2:(nrow(levels))){
    if ((levels[i-1,4]>calibration_width)&(calibration_width>=levels[i,4])){
      cal_depth=levels[(i-1),1]-(levels[(i-1),4]-calibration_width)*(levels[(i-1),1]-levels[i,1])/(levels[(i-1),4]-levels[i,4]);
      cal_near=levels[(i-1),2]-(levels[(i-1),4]-calibration_width)*(levels[(i-1),2]-levels[i,2])/(levels[(i-1),4]-levels[i,4]);
      cal_far=levels[(i-1),3]-(levels[(i-1),4]-calibration_width)*(levels[(i-1),3]-levels[i,3])/(levels[(i-1),4]-levels[i,4]);
      for (j in 2:nrow(levels)){
        if ((xsec[(j-1),2]>cal_depth)&(cal_depth>=xsec[j,2])){
          near=j;
        }
      }
      for (j in (near+1):nrow(xsec)){
        if ((xsec[(j-1),2]<=cal_depth)&(cal_depth<xsec[j,2])){
          far=j-1;
        }
      }
      area=(xsec[near,1]-cal_near)*(cal_depth-xsec[near,2])/2;#double check area
      wp=(((xsec[near,1]-cal_near)^2)+(cal_depth-xsec[near,2])^2)^(1/2); # wetted perimeter
      for (j in near:(far-1)){
        area=area+(xsec[(j+1),1]-xsec[j,1])*((cal_depth-xsec[(j+1),2])+(cal_depth-xsec[j,2]))/2;
        wp=wp+((xsec[(j+1),1]-xsec[j,1])^2+(xsec[(j+1),2]-xsec[j,2])^2)^(1/2);
      }
      area=area+(cal_far-xsec[far,1])*(cal_depth-xsec[far,2])/2;
      wp=wp+((cal_far-xsec[far,1])^2+((cal_depth-xsec[far,2])^2)^(1/2));
    }
  }
  
  R_H=area/wp;
  n=area*(R_H^(2/3))*(S_0^(1/2))/calibration_discharge # Manning's n
  
}else{
  print("calibration failed")
}

widths <- read_csv("widths.csv") # (7) $dt, $filename, $ndwi_threshold, $left_m, $right_m, $width_m
# now, the data are in a dataframe called widths.
dt <- widths$dt
width <- widths$width_m
q <- array(NA, dim = nrow(widths))
for (k in 1:(nrow(widths))){
     if (is.na(widths$width_m[k])==FALSE) {
          for (i in 2:nrow(levels)){
               if ((levels[(i-1),4]>widths$width_m[k])&(widths$width_m[k]>=levels[i,4])){
                    wid_depth=levels[i-1,1]-(levels[i-1,4]-widths$width_m[k])*(levels[i-1,1]-levels[i,1])/(levels[i-1,4]-levels[i,4]);
                    wid_near=levels[i-1,2]-(levels[i-1,4]-widths$width_m[k])*(levels[i-1,2]-levels[i,2])/(levels[i-1,4]-levels[i,4]);
                    wid_far=levels[i-1,3]-(levels[i-1,4]-widths$width_m[k])*(levels[i-1,3]-levels[i,3])/(levels[i-1,4]-levels[i,4]);
                    for (j in 2:nrow(levels)){
                         if ((xsec[j-1,2]>wid_depth)&(wid_depth>=xsec[j,2])){
                              near=j;
                         }
                    }
                    for (j in (near+1):nrow(xsec)){
                         if ((xsec[(j-1),2]<=wid_depth)&(wid_depth<xsec[j,2])){
                              far=j-1;
                         }
                    }
                    area=(xsec[near,1]-cal_near)*(wid_depth-xsec[near,2])/2;
                    wp=((xsec[near,1]-cal_near)^2+(wid_depth-xsec[near,2])^2)^(1/2);
                    for (j in near:(far-1)){
                         area <- area+(xsec[j+1,1]-xsec[j,1])*((wid_depth-xsec[j+1,2])+(wid_depth-xsec[j,2]))/2;
                         wp <- wp+((xsec[j+1,1]-xsec[j,1])^2+(xsec[j+1,2]-xsec[j,2])^2)^(1/2);
                    }
                    area <- area+((wid_far-xsec[far,1])*(wid_depth-xsec[far,2])/2);
                    wp <- wp+((wid_far-xsec[far,1])^2+(wid_depth-xsec[far,2])^2)^(1/2);
               }
          }
          R_H=area/wp;
          Q=area*(R_H^(2/3))*(S_0^(1/2))/n;
          if (widths$width_m[k]<min_width){
               R_H=-8;
               area=-8;
               Q=-8;
          }
          if (widths$width_m[k]>max_width){
               R_H=-9;
               area=-9;
               Q=-9;
          }
          q[k]=Q
     }
}

satellite.gage <- data.frame(dt,width,q)
write_csv(satellite.gage, "sat_gage.csv")

ggplot(satellite.gage) +
     geom_point(aes(x=dt,y=q)) +
     ylim(c(0,15)) +
     xlab("Date") +
     ylab(TeX(r'(Discharge ($m^3/s$))')) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

