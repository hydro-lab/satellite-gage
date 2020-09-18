#This program takes in the depth v. cross-stream-distance profile and
# outputs a width-to-area and -to-wetted-perimeter function, evaluated at
# the widths provided.

# Written by David M. Kahler, 12 April 2018, Duquesne University,
# Pittsburgh, PA, 15282, USA.  This program is provided as-is without any
# guarantees or warrenties.  The program may be used or modified as needed.
# It was written in Matlab R2017a.  Please cite the published article in 
# any submitted work.  For errors or questions, please contact: 
# david.m.kahler@gmail.com.

# This method depends on the parameterization of Manning's equation.  
# Prepare a text file with the .m extension that contains the following
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
 parameter #calls the above- described file

# PROFILE
# Prepare a text file with the cross-stream-distance in the first column
# and the measured depth at each distance in the second column.  Remember
# to define filename and headerlines, and test to check format.  This list
# MUST be sorted from smallest to largest width (top to bottom).
[w,d]=textread(profile,'%f %f','headerlines',profile_headerlines);
# This tests to determine if the depths are input as negative heights or
# positive depths (i.e., orientation of the vertical coordinate).  This
# corrects the values to be negative heights for this analysis
if mean(d)>0
    d=-d;
end
figure;
plot(w,d,'-k')
ylabel('Depth (m)');                 % Check units
xlabel('Cross-stream-distance (m)'); % Check units
# The river bank at the lower cross-stream-distance will be referred to as
# the near bank and the higher cross-stream-distance will be referred to
# as the far bank.
hold on;

# Find any irregularities in the bathymetry.  This searches the depths to
# find the heighest peak in the bathymetry and sets the lowest depth value.
min_depth=min(d);
for i=2:(length(w)-1)
    if (d(i)>min_depth)&&(d(i-1)<d(i))&&(d(i)>d(i+1))
        min_depth=d(i);
        min_depth_index=i;
    end
end

# Find the endpoints of the lowest stage possible in analysis.
for i=2:min_depth_index
    if d(i-1)>min_depth&&min_depth>d(i)
        min_width_1=w(i-1)+((w(i)-w(i-1))/(d(i)-d(i-1)))*(min_depth-d(i-1));
        min_width_index_1=i;
    end
end
for i=(min_depth_index+2):length(w)
    if d(i-1)<min_depth&&min_depth<d(i)
        min_width_2=w(i-1)+((w(i)-w(i-1))/(d(i)-d(i-1)))*(min_depth-d(i-1));
        min_width_index_2=i-1;
    end
end
plot([min_width_1 min_width_2],[min_depth min_depth],'-r');
min_width=min_width_2-min_width_1;

# min_width_index_1 is the first point below floor from the near bank
# min_width_index_1 is the first point below floor from the far bank

# Find the endpoints of the highest stage possible in analysis.
max_width_1=w(1);
max_width_2=w(length(w));
max_depth=d(1);
max_width_index_1=1;
max_width_index_2=length(w);
# In the event that the endpoints do not reach the same datum:
# Scenario 1: depth at near bank is lower, adjust far bank
if d(1)<d(length(w))
    for i=(min_depth_index_2+1):length(w)
        if d(i-1)<d(1)&&d(1)<d(i)
            max_width_2=w(i-1)+((w(i)-w(i-1))/(d(i)-d(i-1)))*(d(1)-d(i-1));
            max_width_index_2=i-1;
        end
    end
end
# Scenario 2: depth at the far bank is lower, adjust near bank
if d(1)>d(length(w))
    max_depth=d(length(w));
    for i=2:min_width_index_1
        if d(i-1)>d(length(w))&&d(length(w))>d(i)
            max_width_1=w(i-1)+((w(i)-w(i-1))/(d(i)-d(i-1)))*(d(length(w))-d(i-1));
            max_width_index_1=i;
        end
    end
end
max_width=max_width_2-max_width_1;
plot([max_width_1 max_width_2],[max_depth max_depth],'-r');

# max_width_index_1 is the first point from the near bank that maps to far
# max_width_index_2 is the first point from the far bank that maps to near

# Generate width-to-depth function by determination of the widths at each
# change of bank-slope.
for i=(max_width_index_1):(min_width_index_1-1)
    for j=min_width_index_2+1:(length(w))
        if (d(j-1)<d(i))&&(d(i)<d(j))
            new_far(i-max_width_index_1+1,1)=w(j-1)+(d(i)-d(j-1))*(w(j)-w(j-1))/(d(j)-d(j-1));
            new_far(i-max_width_index_1+1,2)=d(i);
        end
    end
end
new_far(min_width_index_1,1)=min_width_2;
new_far(min_width_index_1,2)=min_depth;
for i=min_width_index_2+1:max_width_index_2
    for j=2:min_width_index_1
        if (d(j-1)>d(i))&&(d(i)>d(j))
            new_near(i-min_width_index_2,1)=w(j-1)+(d(i)-d(j-1))*(w(j)-w(j-1))/(d(j)-d(j-1));
            new_near(i-min_width_index_2,2)=d(i);
        end
    end
end
new_near(max_width_index_2+1-min_width_index_2,1)=min_width_1;
new_near(max_width_index_2+1-min_width_index_2,2)=min_depth;
plot(w,d,'ob');
plot(new_near(:,1),new_near(:,2),'xg');
plot(new_far(:,1),new_far(:,2),'xg');

# Assemble the bank positions to be used for widths. Possibly make a datafram with all d, w, new near, new far
#Preious code goes along near bank and places a point horizontally along the far bank, so that at every point there is a corresponding point
#then goes doen the far bank and does the same. 
#Finds a width at many places along the transect bw two matching endpoints (horizontlally) so we have a width for each end
#w,d was the original set of points, new_near is the new set of points on the near bank, new far is the new ones on the ffar bank
#Sorts those by W (their location along the transect in order)
#[w d means make a new array of W, d, then mutate for new_near, then mutate for new_far]
temp=[w d;new_near;new_far];
temp=sortrows(temp);
head=0; foot=0;
for i=1:length(temp)
    if temp(i,1)<max_width_1
        head=head+1;
    end
    if temp(i,1)<=max_width_2
        foot=foot+1;
    end
end
xsec=temp((head+1):foot,:); clear temp
for i=1:length(xsec)
    if xsec(i,2)<min_depth
        break
    end
    levels(i,1)=xsec(i,2); # depth
    levels(i,2)=xsec(i,1); # Near bank position
    for j=(i+1):length(xsec)
        if xsec(i,2)==xsec(j,2)
            levels(i,3)=xsec(j,1); # Far bank position
            levels(i,4)=xsec(j,1)-xsec(i,1); # width
        end
    end
end

# Determine area and wetted perimeter for calibration stage: cal_ variables
# are computed and the next lower indeces, near and far, are searched for
# and stored for the numerical integration by the trapazoidal rule and
# wetted perimeter summation by the distance formula.
if (calibration_width>min_width)&&(calibration_width<max_width)
    for i=2:length(levels)
        if (levels(i-1,4)>calibration_width)&&(calibration_width>=levels(i,4))
            cal_depth=levels(i-1,1)-(levels(i-1,4)-calibration_width)*(levels(i-1,1)-levels(i,1))/(levels(i-1,4)-levels(i,4));
            cal_near=levels(i-1,2)-(levels(i-1,4)-calibration_width)*(levels(i-1,2)-levels(i,2))/(levels(i-1,4)-levels(i,4));
            cal_far=levels(i-1,3)-(levels(i-1,4)-calibration_width)*(levels(i-1,3)-levels(i,3))/(levels(i-1,4)-levels(i,4));
            for j=2:length(levels)
                if (xsec(j-1,2)>cal_depth)&&(cal_depth>=xsec(j,2))
                    near=j;
                end
            end
            for j=(length(xsec)-length(levels)+1):length(xsec)
                if (xsec(j-1,2)<=cal_depth)&&(cal_depth<xsec(j,2))
                    far=j-1;
                end
            end
            area=(xsec(near,1)-cal_near)*(cal_depth-xsec(near,2))/2;
            wp=((xsec(near,1)-cal_near)^2+(cal_depth-xsec(near,2))^2)^(1/2);
            for j=near:(far-1)
                area=area+(xsec(j+1,1)-xsec(j,1))*((cal_depth-xsec(j+1,2))+(cal_depth-xsec(j,2)))/2;
                wp=wp+((xsec(j+1,1)-xsec(j,1))^2+(xsec(j+1,2)-xsec(j,2))^2)^(1/2);
            end
            area=area+(cal_far-xsec(far,1))*(cal_depth-xsec(far,2))/2;
            wp=wp+((cal_far-xsec(far,1))^2+(cal_depth-xsec(far,2))^2)^(1/2);
        end
    end
    R_H=area/wp;
    n=area*(R_H^(2/3))*(S_0^(1/2))/calibration_discharge;
else
    disp('calibration failed');
end
f=fopen(specs,'a');
fprintf(f,'Calibration File\n');
fprintf(f,'Cross-sectional area at calibration, A=%10.3f\n',area);
fprintf(f,'Hydraulic radius at calibration, R_{H}=%10.3f\n',R_H);
fprintf(f,'Friction factor (Manning), n=%6.4f\n',n);
fprintf(f,'Minimum allowable width:      %8.4f\n',min_width);
fprintf(f,'Minimum corresponding depth:  %8.4f\n',min_depth);
fprintf(f,'Maximum allowable width:      %8.4f\n',max_width);
fprintf(f,'Maximum corresponding depth:  %8.4f\n',max_depth);
fprintf(f,'In results file, flag: -8 and -9 indicate too small or large width, respectively');
fclose(f);

# WIDTHS
# Prepare a text file with the date, in a format without spaces, in the
# first column and the cross-stream-distance in the second column.  The
# date is not used by this program; the computed cross-sectional area,
# wetted perimeter, hydraulic radius, and discharge will be output in the
# exact same order as the input widths.
[date,wid]=textread(widths,'%s %f','headerlines',widths_headerlines);
q=zeros(length(wid),1);
fo=fopen(data,'a');
fprintf(fo,'width, area, hydraulic_radius, discharge\n');
fprintf(fo,'(m), (m^2), (m), (m^3s^-1)\n');
for k=1:length(wid)
    for i=2:length(levels)
        if (levels(i-1,4)>wid(k))&&(wid(k)>=levels(i,4))
            wid_depth=levels(i-1,1)-(levels(i-1,4)-wid(k))*(levels(i-1,1)-levels(i,1))/(levels(i-1,4)-levels(i,4));
            wid_near=levels(i-1,2)-(levels(i-1,4)-wid(k))*(levels(i-1,2)-levels(i,2))/(levels(i-1,4)-levels(i,4));
            wid_far=levels(i-1,3)-(levels(i-1,4)-wid(k))*(levels(i-1,3)-levels(i,3))/(levels(i-1,4)-levels(i,4));
            for j=2:length(levels)
                if (xsec(j-1,2)>wid_depth)&&(wid_depth>=xsec(j,2))
                    near=j;
                end
            end
            for j=(length(xsec)-length(levels)+1):length(xsec)
                if (xsec(j-1,2)<=wid_depth)&&(wid_depth<xsec(j,2))
                    far=j-1;
                end
            end
            area=(xsec(near,1)-cal_near)*(wid_depth-xsec(near,2))/2;
            wp=((xsec(near,1)-cal_near)^2+(wid_depth-xsec(near,2))^2)^(1/2);
            for j=near:(far-1)
                area=area+(xsec(j+1,1)-xsec(j,1))*((wid_depth-xsec(j+1,2))+(wid_depth-xsec(j,2)))/2;
                wp=wp+((xsec(j+1,1)-xsec(j,1))^2+(xsec(j+1,2)-xsec(j,2))^2)^(1/2);
            end
            area=area+(wid_far-xsec(far,1))*(wid_depth-xsec(far,2))/2;
            wp=wp+((wid_far-xsec(far,1))^2+(wid_depth-xsec(far,2))^2)^(1/2);
        end
    end
    R_H=area/wp;
    Q=area*(R_H^(2/3))*(S_0^(1/2))/n;
    if wid(k)<min_width
        R_H=-8;
        area=-8;
        Q=-8;
    end
    if wid(k)>max_width
        R_H=-9;
        area=-9;
        Q=-9;
    end
    q(k,1)=Q;
    fprintf(fo,'%10.3f, %10.3f, %10.3f, %10.3f\n',wid(k),area,R_H,Q);
end
fclose(fo);
