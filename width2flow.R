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

# Remember to set working directory

# Channel-specific information
calibration_discharge <- 4.567 # discharge, cubic meters per second
calibration_width <- 27.65 # width, meters
S_0 <- 0.0167 # slope, dimensionless

profile <- read_csv("profile.csv") # (2) $location_m and $height_m
if (mean(profile$height_m)>0){ # guesses if the average height is positive, that the height vector points down
     profile$height_m = -profile$height_m # code is based on the height vector to point up
}

wid <- read_csv("width.csv") # (7) $dt, $filename, $ndwi_threshold_3, $ndwi_threshold_2, $left_m, $right_m, $width_m

ggplot(profile) +
     geom_line(aes(x = location_m, y = height_m)) +
     xlab("Cross-stream distance (m)") +
     ylab("Bathymetry (m)") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

min_depth <- min(profile$height_m)

for (i in 2:(length(W)-1)){
  if ((d[i]>min_depth)&(d[i-1]<d[i])&(d[i]>d[i+1])){
    min_depth=d[i];
    min_depth_index=i;
  }
}

# Find the endpoints of the lowest stage possible in analysis.
for (i in 2:(min_depth_index)){
  if ((d[i-1]>=min_depth)&(min_depth>d[i])){
    min_width_1=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(min_depth-d[i-1]);
    min_width_index_1=i;
  }
}
for (i in (min_depth_index+2):(length(W))){
  if ((d[i-1]<=min_depth)&(min_depth<d[i])){
    min_width_2=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(min_depth-d[i-1]);
    min_width_index_2=i-1;
  }
}


min_width <- (min_width_2-min_width_1);

max_width_1=W[1];
max_width_2=W[length(W)];
max_depth <- d[1];

max_width_index_1 <- 1;
max_width_index_2 <- length(W);


for (i in (min_depth_index+2):(length(W))){
  if ((d[i-1]<min_depth)&(min_depth<d[i])){
    min_width_2=W[i-1]+((W[i]-W[i-1])/(d[i]-d[i-1]))*(min_depth-d[i-1]);
    min_width_index_2=i-1;
  }
}

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

new_far <- array(-9999, dim = c(((min_width_index_1)-(max_width_index_1)+1),2))
new_near <- array(-9999, dim = c(((max_width_index_2)-(min_width_index_2)+1),2))

for (i in (max_width_index_1):(min_width_index_1-1)){
  for (j in ((min_width_index_2+1):(length(W)))){
    if (((d[j-1]<d[i]))&((d[i]<d[j]))){
      new_far[i-max_width_index_1+1,1]=W[j-1]+(d[i]-d[j-1])*(W[j]-W[j-1])/(d[j]-d[j-1]);
      new_far[i-max_width_index_1+1,2]=d[i];
    }
  }
}

new_far[min_width_index_1-max_width_index_1+1,1] <- min_width_2;
new_far[min_width_index_1-max_width_index_1+1,2]<- min_depth;
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

prof <- cbind(W, d)
new <- rbind(new_near, new_far)
complet <- rbind(prof, new)
temp <- complet[order(complet[,1]),]

head=0; foot=0;
for (i in 1:nrow(temp)){
  if ((temp[i,1])<max_width_1){
    head=head+1;
  }
  if ((temp[i,1])<=max_width_2){
    foot=foot+1;
  }
}

xsec=temp[(head+1):foot,] 
rm(temp)
levels <- array(-9999, dim = c(nrow(xsec),4))
for (i in 1:nrow(xsec)){
  if (xsec[i,2]<min_depth){
    break
  }
  levels[i,1]=xsec[i,2]; # depth
  levels[i,2]=xsec[i,1]; # Near bank position
  for (j in (i+1):nrow(xsec)){
    if ((xsec[i,2])==(xsec[j,2])){
      levels[i,3] <- xsec[j,1]; # Far bank position
      levels[i,4] <- xsec[j,1]-xsec[i,1]; # width
    }
  }
}

if ((calibration_width>min_width)&&(calibration_width<max_width)){
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
      wp=(((xsec[near,1]-cal_near)^2)+(cal_depth-xsec[near,2])^2)^(1/2);
      for (j in near:(far-1)){
        area=area+(xsec[(j+1),1]-xsec[j,1])*((cal_depth-xsec[(j+1),2])+(cal_depth-xsec[j,2]))/2;
        wp=wp+((xsec[(j+1),1]-xsec[j,1])^2+(xsec[(j+1),2]-xsec[j,2])^2)^(1/2);
      }
      area=area+(cal_far-xsec[far,1])*(cal_depth-xsec[far,2])/2;
      wp=wp+((cal_far-xsec[far,1])^2+((cal_depth-xsec[far,2])^2)^(1/2));
    }
  }
  
  R_H=area/wp;
  n=area*(R_H^(2/3))*(S_0^(1/2))/calibration_discharge;
  
}else{
  print("calibration failed")
}

# now, the data is in a dataframe called wid.
dataout <- array(-9, dim=c(nrow(wid), 4))
#[date,wid]=textread(widths,'%s %f','headerlines',widths_headerlines);
q <- array(0, dim = c(nrow(wid),1))
#fo=fopen(data,'a');
#fprintf(fo,'width, area, hydraulic_radius, discharge\n');
#fprintf(fo,'(m), (m^2), (m), (m^3s^-1)\n');
for (k in 1:(nrow(wid))){
  for (i in 2:nrow(levels)){
    if ((levels[(i-1),4]>wid$width_m[k])&(wid$width_m[k]>=levels[i,4])){
      wid_depth=levels[i-1,1]-(levels[i-1,4]-wid$width_m[k])*(levels[i-1,1]-levels[i,1])/(levels[i-1,4]-levels[i,4]);
      wid_near=levels[i-1,2]-(levels[i-1,4]-wid$width_m[k])*(levels[i-1,2]-levels[i,2])/(levels[i-1,4]-levels[i,4]);
      wid_far=levels[i-1,3]-(levels[i-1,4]-wid$width_m[k])*(levels[i-1,3]-levels[i,3])/(levels[i-1,4]-levels[i,4]);
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
      dataout[k,1] <- area
    }
  }
  R_H=area/wp;
  Q=area*(R_H^(2/3))*(S_0^(1/2))/n;
  if (wid$width_m[k]<min_width){
    R_H=-8;
    area=-8;
    Q=-8;
  }
  if (wid$width_m[k]>max_width){
    R_H=-9;
    area=-9;
    Q=-9;
  }
  
  q[k,1]=Q; 
  dataout[k,2] <- (wid$width_m[k])
  dataout[k,3] <- R_H
  dataout[k,4] <- (wid$filename[k])
  # wid(k),area,R_H,Q);
}

Cr <- c(area)
Hr <- c(R_H)
Ff <- c(n)
Miaw <- c(min_width)
Micd <- c(min_depth)
Maw <- c(max_width)
Mcd <- c(max_depth)
specs <- data.frame(Cr, Hr, Ff, Miaw, Micd, Maw, Mcd)
names(specs)[1] <- "Cross-sectional area at calibration"
names(specs)[2] <- "Hydraulic radius at calibration"
names(specs)[3] <- "Friction factor (Manning)"
names(specs)[4] <- "Minimum allowable width"
names(specs)[5] <- "Minimum corresponding depth"
names(specs)[6] <- "Maximum allowable width"
names(specs)[7] <- "Maximum corresponding depth"

flowoutput <- data.frame(dataout, q)
names(flowoutput)[1] <- "Cross-sectional area at calibration"
names(flowoutput)[2] <- "Width"
names(flowoutput)[3] <- "Hydraulic radius"
names(flowoutput)[4] <- "Date/Filename"
names(flowoutput)[5] <- "Flow (Q)"

write.table(flowoutput, file = "flowoutput.csv", append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)

