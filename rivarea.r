# This method depends on the parameterization of Manning's equation.  
# Prepare a text file that contains the following
# variable definitions:
# calibration_discharge= THE MEASURED DISCHARGE
# calibration_width= THE MEASURED WIDTH
# S_0= THE MEASURED STREAMWISE SLOPE
# profile= THE FILENAME FOR THE CROSS-SECTIONAL PROFILE (desc. below)
# number of headerlines
# widths= THE FILENAME FOR THE WIDTH DATA (desc. below)
# number of headerlines
# specs= specifications output file name
# data= flow output file name
# parameter #calls the above- described file
setwd("/Users/littlesunsh9/Documents/Kahler Lab/planet_order_181828/")
profile <- read.table("bcprofile.txt")
names(profile)[1] <- "W"
names(profile)[2] <- "d"

W <- profile$W
d <- profile$d

# This tests to determine if the depths are input as negative heights or
# positive depths (i.e., orientation of the vertical coordinate).  This
# corrects the values to be negative heights for this analysis

if (mean(d)>0){
    d = -d
    profile$d = -profile$d
}

plot(profile, type = "p", main="Bathymetric Profile", 
     xlab= "Cross-stream-distance (m)" ,
     ylab ="Depth (m)") # Check units

# The river bank at the lower cross-stream-distance will be referred to as
# the near bank and the higher cross-stream-distance will be referred to
# as the far bank.
# Find any irregularities in the bathymetry.  This searches the depths to
# find the heighest peak in the bathymetry and sets the lowest depth value.
min_depth <- min(d);
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
#plot(x = (min_width_1 min_depth), y = (min_width_2 min_depth) );
min_width <- (min_width_2-min_width_1);
#plot(profile, type = "p", main="Bathymetric Profile", 
#  xlab= "Cross-stream-distance (m)" ,
#  ylab ="Depth (m)") # Check units
# min_width_index_1 is the first point below floor from the near bank
# min_width_index_1 is the first point below floor from the far bank

# Find the endpoints of the highest stage possible in analysis.
max_width_1=W[1];
max_width_2=W[length(W)];
max_depth <- d[1];
max_width_index_1 <- 1;
max_width_index_2 <- length(W);
# In the event that the endpoints do not reach the same datum:
# Scenario 1: depth at near bank is lower, adjust far bank

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
#lines([max_width_1 max_width_2],[max_depth max_depth]);

# max_width_index_1 is the first point from the near bank that maps to far
# max_width_index_2 is the first point from the far bank that maps to near

new_far <- array(-9999, dim = c(((min_width_index_1)-(max_width_index_1)+1),2))
new_near <- array(-9999, dim = c(((max_width_index_2)-(min_width_index_2)+1),2))

# Generate width-to-depth function by determination of the widths at each
# change of bank-slope.
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

#plot(W,d,'ob');
#plot(new_near(:,1),new_near(:,2),'xg');
#plot(new_far(:,1),new_far(:,2),'xg');

# Assemble the bank positions to be used for widths.

prof <- cbind(W, d)
new <- rbind(new_near, new_far)
complet <- rbind(prof, new)
temp <- complet[order(complet[,1]),]

head=0; foot=0;
for (i in 1:length(temp)){
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

# Determine area and wetted perimeter for calibration stage: cal_ variables
# are computed and the next lower indeces, near and far, are searched for
# and stored for the numerical integration by the trapazoidal rule and
# wetted perimeter summation by the distance formula.
if ((calibration_width>min_width)&&(calibration_width<max_width)){
    for (i in 2:(length(levels))){
        if ((levels[i-1,4]>calibration_width)&&(calibration_width>=levels[i,4])){
            cal_depth=levels[i-1,1]-(levels[i-1,4]-calibration_width)*(levels[i-1,1]-levels[i,1])/(levels[i-1,4]-levels[i,4]);
            cal_near=levels[i-1,2]-(levels[i-1,4]-calibration_width)*(levels[i-1,2]-levels[i,2])/(levels[i-1,4]-levels[i,4]);
            cal_far=levels[i-1,3]-(levels[i-1,4]-calibration_width)*(levels[i-1,3]-levels[i,3])/(levels[i-1,4]-levels[i,4]);
            for (j in 2:length(levels)){
                if ((xsec[j-1,2]>cal_depth)&&(cal_depth>=xsec[j,2])){
                    near=j;
                }
            }
            for (j in (length(xsec)-length(levels)+1):length(xsec)){
                if ((xsec[j-1,2]<=cal_depth)&&(cal_depth<xsec[j,2])){
                    far=j-1;
                }
            }
            area=(xsec[near,1]-cal_near)*(cal_depth-xsec[near,2])/2;
            wp=((xsec[near,1]-cal_near)^2+(cal_depth-xsec[near,2])^2)^(1/2);
            for (j in near:(far-1)){
                area=area+(xsec[j+1,1]-xsec[j,1])*((cal_depth-xsec[j+1,2])+(cal_depth-xsec[j,2]))/2;
                wp=wp+((xsec[j+1,1]-xsec[j,1])^2+(xsec[j+1,2]-xsec[j,2])^2)^(1/2);
            }
            area=area+(cal_far-xsec[far,1])*(cal_depth-xsec[far,2])/2;
            wp=wp+((cal_far-xsec[far,1])^2+((cal_depth-xsec[far,2])^2)^(1/2));
        }
    }
    R_H=area/wp;
    n=area*(R_H^(2/3))*(S_0^(1/2))/calibration_discharge;
    }else{
        disp('calibration failed');
        }

#f=fopen(specs,'a');
#fprintf(f,'Calibration File\n');
#fprintf(f,'Cross-sectional area at calibration, A=%10.3f\n',area);
#fprintf(f,'Hydraulic radius at calibration, R_{H}=%10.3f\n',R_H);
#fprintf(f,'Friction factor (Manning), n=%6.4f\n',n);
# fprintf(f,'Minimum allowable width:      %8.4f\n',min_width);
# fprintf(f,'Minimum corresponding depth:  %8.4f\n',min_depth);
# fprintf(f,'Maximum allowable width:      %8.4f\n',max_width);
# fprintf(f,'Maximum corresponding depth:  %8.4f\n',max_depth);
# fprintf(f,'In results file, flag: -8 and -9 indicate too small or large width, respectively');
# fclose(f);

# WIDTHS
# Prepare a text file with the date, in a format without spaces, in the
# first column and the cross-stream-distance in the second column.  The
# date is not used by this program; the computed cross-sectional area,
# wetted perimeter, hydraulic radius, and discharge will be output in the
# exact same order as the input widths.

[date,wid]=textread(widths,'%s %f','headerlines',widths_headerlines);
q=zeros(length(wid),1);
#fo=fopen(data,'a');
#fprintf(fo,'width, area, hydraulic_radius, discharge\n');
#fprintf(fo,'(m), (m^2), (m), (m^3s^-1)\n');
for (k in 1:(length(wid))){
    for (i in 2:length(levels)){
        if ((levels(i-1,4)>wid(k))&&(wid(k)>=levels(i,4))){
            wid_depth=levels(i-1,1)-(levels(i-1,4)-wid(k))*(levels(i-1,1)-levels(i,1))/(levels(i-1,4)-levels(i,4));
            wid_near=levels(i-1,2)-(levels(i-1,4)-wid(k))*(levels(i-1,2)-levels(i,2))/(levels(i-1,4)-levels(i,4));
            wid_far=levels(i-1,3)-(levels(i-1,4)-wid(k))*(levels(i-1,3)-levels(i,3))/(levels(i-1,4)-levels(i,4));
            for (j in 2:length(levels)){
                if ((xsec(j-1,2)>wid_depth)&&(wid_depth>=xsec(j,2))){
                    near=j;
                }
            }
            for (j in (length(xsec)-length(levels)+1):length(xsec)){
                if ((xsec[j-1,2]<=wid_depth)&&(wid_depth<xsec[j,2])){
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
    if (wid(k)<min_width){
        R_H=-8;
        area=-8;
        Q=-8;
    }
    if (wid(k)>max_width){
        R_H=-9;
        area=-9;
        Q=-9;
    }
    q(k,1)=Q;
    # fprintf(fo,'%10.3f, %10.3f, %10.3f, %10.3f\n',wid(k),area,R_H,Q);
}
#fclose(fo);
