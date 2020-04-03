# This code finds the boundary of the water in a normalized difference water index image based on the histogram of the pixel values.
# this code uses the cropped, single-layer, NDWI image

library(raster)
library(rgdal)
library(rgeos)

# Import raster image, or take it from previous code, set working directory, if needed.
setwd("/Volumes/LaCie2big/widthflowproject/r_test")

fn <- "./20180212cNDWI.tif"
x = raster(fn)
y <- x[]
h = hist(y, breaks=seq(-1,1,by=0.01)) # built-in histogram function.
# Positions: h$mids       number
# Values:    h$counts     integer
bins <- h$mids
v <- h$counts

# Allocate arrays used in analysis
avg <- array(0, dim = c(200,10))
peaks <- array(0, dim = c(200,10))
nop <- array(0, dim = c(1,10))
for (w in 1:10){
  # filter values (v=h$counts) with the averaging window size 2*w+1
  for (i in (w+1):(200-w)){
    avg[i,w] <- ((sum(v[(i-w):(i+w)]))/((2*w)+1))
  }
  # identify and number peaks
  cnt <- 0
  for (j in (w+1):(200-w)){
    if ((avg[j-1,w])<(avg[j,w])){
      if ((avg[j+1,w])<(avg[j,w])){
        cnt <- (cnt+1)
        peaks[j,w] <- cnt
        nop[1,w] <- cnt
      }
    }
  }
}

# set error values for the result vectors in case neither two nor three peaks are found:
threepeak <- -9999
twopeak <- -9999

for (w in 1:10){
  # testing in three peaks
  # due to the order of the w variable, only the 'smoothest' result will be kept
  if ((nop[w])==3){
    # finds the second and third peak
    for (j in 1:200){
      if ((peaks[j,w])==2){
        sec <- j # stores the index of the second peak
      }
      if ((peaks[j,w])==3){
        thr <- j # stores the index of the third peak
      }
    }
    # finds minimum between second and third peak
    m <- max(v) # create variable for minimum, initially set higher than any value
    for (j in (sec):(thr)){
      if ((avg[j,w])<m){
        goal <- j
        m <- avg[j,w]
      }
    }
    threepeak <- (bins[(goal)])
  }
  # test if exactly three peaks were not found
  if ((nop[w])==2){
    # find the position of the first and second (the only) peaks
    for (j in 1:200){
      if ((peaks[j,w])==1){
        fst <- j # stores the index of the second peak
      }
      if ((peaks[j,w])==2){
        sec <- j # stores the index of the third peak
      }
    }
    # finds minimum between first and second peak
    m <- max(v) # create variable for minimum, initially set higher than any value
    for (j in (sec):(thr)){
      if ((avg[j,w])<m){
        goal <- j
        m <- avg[j,w]
      }
    }
    twopeak <- (bins[(goal)])
  }
}

# write to file
date <- Sys.Date()
write.table((c(threepeak, twopeak, date)), file = "wateredge.txt", append = TRUE, sep = ", ", dec = ".")
# output will be in the same order as the input files

# run in command line as:
# r -f peaktest.r
# make sure input and output filenames are coded correctly
# https://cran.r-project.org/doc/manuals/R-intro.html#Invoking-R-from-the-command-line
