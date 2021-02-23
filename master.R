#Run in command line as: Rscript master.R

#Current Update is to inclide a check for existance of the crop boundaries in each raster file else skip

# Raster calcultion to crop Geotiff, also pulls XML metadata: reflectance coefficient (ps:reflectanceCoefficient)
library(rgdal)
library(raster)
library(rgeos)
library(XML)
library(methods)

# remember to set working directory if needed:
setwd("/Users/littlesunsh9/Documents/Kahler Lab/planet_order_181828/")
#Lists for necessary files
#Image list
im <- list.files("/Users/littlesunsh9/Documents/Kahler Lab/planet_order_181828/", pattern = "*AnalyticMS.tif$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)
#Metadata List
g <- list.files("/Users/littlesunsh9/Documents/Kahler Lab/planet_order_181828/", pattern = "*AnalyticMS_metadata.xml$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)

#Inputs from 
#cald <- 4.3; # calibration discharge (the measured discharge), cubic meters per second
#calw <- 27.8; # calibration width (the measured width), meters
#S_0 <- 0.0006;# the measured streamwise slope
#INPUT FILES:
#profile <- read.table('bcprofile.txt'); #txt readable file of depths where 
# cross-stream-distance is C1 and measured depth is C2. 
# File must be sorted by cross-stream-distance in ascending order (from 0:width, top:bottom)
#profile_headerlines=1;
#widths <- read.table('widths.txt');

width <- array(-9, dim=c(length(im),6)) #creates the output file as an array that can be easily read as a table
prowidths <- array(-9, dim=c(length(im),2)) #creates the output file for the parameters for Mannings callibration/use

# LOOP STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for (q in 1:(length(im))){
  #Import raw Planet metadata
  fn <- g[q]
  fl <- xmlParse(fn)
  rc <- setNames(xmlToDataFrame(node=getNodeSet(fl, "//ps:EarthObservation/gml:resultOf/ps:EarthObservationResult/ps:bandSpecificMetadata/ps:reflectanceCoefficient")),"reflectanceCoefficient")
  dm <- as.matrix(rc)
  # 1 Red
  # 2 Green
  # 3 Blue
  # 4 Near infrared
  rc2 <- as.numeric(dm[2]) # Green
  rc4 <- as.numeric(dm[4]) # NIR
  
  # Import raster image, crops to chosen extent
  fn <- im[q]
  pic <- stack(fn)
  # set extent from QGIS analysis:
  # extent format (xmin,xmax,ymin,ymax)
  e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons')
  #crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
  crs(e) <- "+proj=utm +zone=17 +datum=WGS84"
  r <- crop(pic, e)
  # rm(pic) # remove rest of image from RAM
  
  #options(stringsAsFactors = FALSE)
  
  rbrick <- brick(r)
  # calculate normalized difference water index (NDWI).  :
  # calculate NDWI using the green (band 2) and nir (band 4) bands
  ndwi <- ((rc2*r[[2]]) - (rc4*r[[4]])) / ((rc2*r[[2]]) + (rc4*r[[4]]))
  # This formulation follows: Gao, B. (1996). NDWI—A normalized difference water index for remote sensing of vegetation liquid water from space. Remote Sensing of Environment, 53(3), p. 257-266. https://www.sciencedirect.com/science/article/abs/pii/S0034425796000673
  
  # To view, during development:
  #plot(ndwi)
  
  # To export cropped NDWI as a new file and create filename root
  p <- strsplit(im[q], "_3B_AnalyticMS.tif")
  r <- strsplit(p[[1]], "/")
  lr <- tolower(r[[1]])
  len <- length(lr)
  root <- lr[[len]]
  
  writeRaster(x = ndwi,
              filename= paste(root, "cndwi.tif", sep="."),
              format = "GTiff", # save as a tif
              # save as a FLOAT if not default, not integer
              overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.
  
  #attempting the following crop procedure https://gis.stackexchange.com/questions/229356/crop-a-raster-file-in-r
  #attempting the following NDVI as NDWI procedure https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/vegetation-indices-NDVI-in-R/
  
  width[q,1] <- root #for output file: root name of image
  
  # This code finds the boundary of the water in a normalized difference water index image based on the histogram of the pixel values.
  # This code uses the cropped, single-layer, NDWI image
  
  # Import raster image, or take it from previous code, set working directory, if needed.
  
  h = hist(ndwi, breaks=seq(-1,1,by=0.01)) # built-in histogram function.
  # Positions: h$mids       number
  # Values:    h$counts     integer
  bins <- h$mids
  v <- h$counts
  
  #png(filename = root, "hist.png")
  #plot(h)
  #dev.copy(png, filename = paste(root, "hist.png", sep = "."))
  #dev.off()
  setEPS()
  postscript(paste(root, "hist.eps", sep = "."))
  plot(h)
  dev.off()  
  
  
  # Allocate arrays used in analysis
  avg <- array(0, dim = c(200,10))
  peaks <- array(0, dim = c(200,10))
  nop <- array(0, dim = c(1,10))
  for (w in 1:10){
    # filter values (v=h$counts) with the averaging window size 2*w+1
    for (k in (w+1):(200-w)){
      avg[k,w] <- ((sum(v[(k-w):(k+w)]))/((2*w)+1))
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
      for (j in (fst):(sec)){
        if ((avg[j,w])<m){
          goal <- j
          m <- avg[j,w]
        }
      }
      twopeak <- (bins[(goal)])
    }
  }
  # Water's Edge LOOP ENDS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  width[q,2] <-threepeak #for output file: value for the edge of water (3 peak)
  width[q,3] <-twopeak #for output file: value for the edge of water (2 peak)
  
  # write individual water's edge values to the compiled water's edge file
  date <- as.character(Sys.Date())
  ex <- data.frame(threepeak, twopeak, date)
  write.table(ex, file = "wateredge.csv", append = TRUE, sep = ",", dec = ".", col.names = FALSE)
  #write.txt((c(threepeak, twopeak, date)), file = "wateredge.txt", append = TRUE, sep = ", ", dec = ".")
  #write.csv(ex,file = "wateredge.csv")
  # output will be in the same order as the input files
  
  #Finding values along the width of the transect
  #cndwi <- list.files("/Users/littlesunsh9/Documents/planet_order_181828/", pattern = "*.cndwi.tif$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)
  #for (q in 1:(length(cndwi))){
  # fn <- cndwi[q]
  #ndwia <- stack(fn)
  #e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons')
  #crs(e) <- "+proj=utm +zone=17 +datum=WGS84"
  #RDB=(609589.376, 4507801.407)
  #LDB=(609634.607, 4507831.586)
  x1 <- (609589.376)
  x2 <- (609634.607)
  y1 <- (4507801.407)
  y2 <- (4507831.586)
  ma <- (y2-y1)/(x2-x1)
  #this next part will rely on UTM (the coordinates are in meters)
  ra <- (0.1)
  #r is the step size of each point along our width
  t <- sqrt(((x2-x1)^2)+ ((y2-y1)^2))
  f <- ceiling(t/ra)
  pointers <- array(999.999, dim = c(f,2))
  pointers[1,1] <- x1
  pointers[1,2] <- y1
  for (i in 2:f){
    a <- 1
    b <- (-2)*pointers[i-1,1]
    c <- (pointers[i-1,1]^2) - (ra^2)/((ma^2)+1)
    pointers[i,1] <- ((-b)+(sqrt((b^2)-4*a*c)))/(2*a)
    pointers[i,2] <- ((pointers[i-1,2])+(ma*((pointers[i,1])-(pointers[i-1,1]))))
  }
  
  #Check:
  #x <- c(x1,x2)
  #y <- c(y1,y2)
  #plot(x,y)
  
  spat <- SpatialPoints(pointers)
  
  alng <- extract(ndwi, spat, method='simple')
  #plot(alng)
  
  #p <- strsplit(ndwi[q], "cndwi.tif")
  #r <- strsplit(p[[1]], "/")
  #lr <- tolower(r[[1]])
  #len <- length(lr)
  #root <- lr[[len]]
  
  # To export table as a new file
  write.table(alng, file = paste(root, "distwidth.csv", sep="."), append = TRUE, sep = ",", dec = ".", col.names = FALSE)
  
  alng_per <- array(-9, dim=c(f,2)) #allocation for the midpoints
  #when you reach -9 in that array you've reached the end of the midpoints/values found
  restart <- 2 #initial start for i search
  for (j in 1:f){
    cnt=1
    for (i in restart:f){
      if (alng[i]==alng[i-1]){
        cnt=cnt+1
      } 
      else {
        restart <- i+1 #to keep moving forward from the last section without causing a loop at the end of it
        break
      }
    }
    mp <- ((cnt*ra)/2) #ra is the spacing, and mp gives the midpoint of the current distance section
    if (i<(f)){
      alng_per[j,1] <- (((i-1)*ra)-mp) #this makes the first column the location along the index where the midpoint is located, as i will be the end of the section, take away the mp value, and you get the correct location
      alng_per[j,2] <- alng[i-1] 
    }
    else{
      alng_per[j,1] <- ((f*ra)-mp) #this makes the first column the location along the index where the midpoint is located, as i will be the end of the section, take away the mp value, and you get the correct location
      alng_per[j,2] <- alng[i-1] 
    }
    if (i>=(f)){
      break
    } 
  }
  
  for (i in (2:f)){
    if (alng_per[i,2]>threepeak){
      if (alng_per[i-1,2]<threepeak){
        i1 <- alng_per[i-1,1]
        i2 <- alng_per[i,1]
        j1 <- alng_per[i-1,2]
        j2 <- alng_per[i,2]
        n <- threepeak    
        RDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
        break
      }
    }
  }
  for (i in 1:(f-1)){
    if (alng_per[f-i,2]>threepeak){        #expressing the index such that when i = 1, f, and when i = 2, f-1.
      if (alng_per[f-i+1,2]<threepeak){
        i1 <- alng_per[f-i+1,1]
        i2 <- alng_per[f-1,1]
        j1 <- alng_per[f-i+1,2]
        j2 <- alng_per[f-1,2]
        n <- threepeak    
        LDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
        break
      }
    }
  }
  width[q,4] <- LDB #location in meters of bank 1
  width[q,5] <- RDB #location in meters of bank 2
  width[q,6] <- LDB-RDB #gives width in meters
  prowidths[q,1] <- root
  prowidths[q,2] <- LDB-RDB
}
survey <- data.frame(width)
names(survey)[1] <- "filename"
names(survey)[2] <- "waters_edge_value_3peak"
names(survey)[3] <- "waters_edge_value_2peak"
names(survey)[4] <- "LDB_location_m"
names(survey)[5] <- "RDB_location_m"
names(survey)[6] <- "width_m"
write.table(survey, file = "width.csv", append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)

write.table(prowidths, file = "widths.txt", append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)

# run in command line as:
# r -f peaktest.r
# make sure input and output filenames are coded correctly
# https://cran.r-project.org/doc/manuals/R-intro.html#Invoking-R-from-the-command-line