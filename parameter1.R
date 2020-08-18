#Parameter file
cd=4.3; # calibration discharge (the measured discharge), cubic meters per second
cw=27.8; # calibration width (the measured width), meters
S_0=0.0006;# the measured streamwise slope
#INPUT FILES:
profile='bcprofile.txt'; #txt readable file of depths where 
# cross-stream-distance is C1 and measured depth is C2. 
# File must be sorted by cross-stream-distance in ascending order (from 0:width, top:bottom)
profile_headerlines=1;
widths='???'; # Prepare a text file with the date, in a format without spaces, in the
# first column and the cross-stream-distance in the second column.  The
# date is not used by this program; the computed cross-sectional area,
# wetted perimeter, hydraulic radius, and discharge will be output in the
# exact same order as the input widths.
widths_headerlines=1;
# OUTPUT FILES:
    specs='buffalo_cr_specs2.txt'; # specification output
parameters ='buffalo_cr_data2.csv'; #data output - will be comma-separated values

write.table(parameters, file = "parameteroutput.csv", append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)


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