# Code to import image data from downloaded Planet Labs images and record the 
# widths of the river at a pre-programmed location.

# Written in Matlab and R with RStudio by David Kahler and Mackenzie Martin, 
# Duquesne University, from 2018 to 2021.  The development was supported by 
# the United States Agency for International Development, Southern Africa 
# Regional Mission.  Further information is available at: 
# www.duq.edu/limpopo 
# https://github.com/LimpopoLab 

# Run in command line as: Rscript image2width.R

# Raster calculation to crop Geotiff, also pulls XML metadata: reflectance coefficient (ps:reflectanceCoefficient)
library(rgdal)
library(raster)
library(rgeos)
library(XML)
library(methods)
library(sp)
library(parallel)
library(MASS)
library(doParallel)
library(stringr)
library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)

# remember to set working directory if needed

# Lists for necessary files
# Image list
im <- list.files("./", 
                 pattern = "*AnalyticMS.tif$", 
                 full.names = TRUE, 
                 recursive = TRUE, 
                 ignore.case=TRUE, 
                 include.dirs = TRUE)
di <- array(NA, dim = length(im))
for (i in 1:length(im)) {
     a <- str_split(im[i],"/")
     b <- str_split(a[[1]][length(a[[1]])],"_")
     c <- as.character(b[[1]][1])
     d <- as.character(b[[1]][2])
     f <- paste0(c,"T",d)
     di[i] <- ymd_hms(f)
}

# Metadata List
md <- list.files("./", 
                pattern = "*AnalyticMS_metadata.xml$", 
                full.names = TRUE, 
                recursive = TRUE, 
                ignore.case=TRUE, 
                include.dirs = TRUE)
dm <- array(NA, dim = length(md))
for (i in 1:length(md)) {
     a <- str_split(md[i],"/")
     b <- str_split(a[[1]][length(a[[1]])],"_")
     c <- as.character(b[[1]][1])
     d <- as.character(b[[1]][2])
     f <- paste0(c,"T",d)
     dm[i] <- ymd_hms(f)
}

rm(a,b,c,d,f)
id <- array(NA, dim = length(im))  # will match metadata filenames to image filenames and dates
for (i in 1:length(im)) {
     for (j in 1:length(md)) {
          if (di[i]==dm[j]) { # if image date matches metadata date,
               id[i] <- md[j] # store metadata filename matched to image filename and date
          }
     }
}
imagebank <- data.frame(di,im,id)
rm(di,dm,im,md,id,i,j)
imagebank <- imagebank %>% 
     rename(dt=di,md=id) %>% 
     filter(is.na(md)==FALSE) # will contain imagebank data frame with date (dt), image (im), and metadata (md)

# LOOP STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# SINGLE
#widths <- array(NA, dim = c((nrow(imagebank)),7))
#for (q in 1:(2)) { # original loop
# PARALLEL
# Replacing loop with a foreach for parallelization
registerDoParallel(detectCores())
#widths <- foreach (q = 1:2, .combine = 'rbind') %dopar% { # testing loop,
widths <- foreach (q = 1:(nrow(imagebank)), .combine = 'rbind') %dopar% { # parallel computing loop: this changes how data are transferred back from each operation.
     
     #Import raw Planet metadata to get the reflectance coefficients
     fn <- imagebank$md[q]
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
     fn <- imagebank$im[q]
     pic <- stack(fn)
     # set extent from QGIS analysis:
     # extent format (xmin,xmax,ymin,ymax)
     # Buffalo Creek:
     e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons') # Extent needed
     crs(e) <- "+proj=utm +zone=17 +datum=WGS84"
     # Mutale River downstream
     #e <- as(extent( , , , ), 'SpatialPolygons')
     #crs(e) <- "+proj=longlat +datum=WGS84"
     # Set extent from the Planet file !! This is the area from the picture
     test <- as(extent(pic), 'SpatialPolygons') # Extent of image
     crs(test) <- "+proj=utm +zone=17 +datum=WGS84"
     if (gCovers(test,e)) { # returns TRUE if no point in spgeom2 (e, needed) is outside spgeom1 (test, image extent) # used to be (gWithin(e, test, byid = FALSE))
          r <- crop(pic, e)
          rm(pic) # remove rest of image from RAM
          rbrick <- brick(r)
          # calculate NDWI using the green (band 2) and nir (band 4) bands
          ndwi <- ((rc2*r[[2]]) - (rc4*r[[4]])) / ((rc2*r[[2]]) + (rc4*r[[4]]))
          # plot(ndwi) # for viewing during development
          
          # To export cropped NDWI as a new file and create filename root
          p <- strsplit(imagebank$im[q], "_3B_AnalyticMS.tif")
          r <- strsplit(p[[1]], "/")
          lr <- tolower(r[[1]])
          len <- length(lr)
          root <- lr[[len]]
          rm(p,r,lr,len)
          # writeRaster(x = ndwi, ## this does not need to be done, just a nice record.
          #             filename= paste(root, "cndwi.tif", sep="."),
          #             format = "GTiff", # save as a tif, save as a FLOAT if not default, not integer
          #             overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.
          
          output <- array(NA, dim = 7)
          output[1] <- date(as_datetime(imagebank$dt[q]))
          output[2] <- root #for output file: root name of image
          
          # This code finds the boundary of the water in a normalized difference water index 
          # This code uses the cropped, single-layer, NDWI image.  Image based on the histogram of the pixel values.
          
          # Import raster image, or take it from previous code, set working directory, if needed.
          
          h = hist(ndwi, # built-in histogram function.  To find values only.  Plotting is at the end of this loop.
                   breaks=seq(-1,1,by=0.01),
                   plot=FALSE) 
          bins <- h$mids # positions    number
          v <- h$counts # counts        integer
          
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
          threepeak <- -1 # revised error values so the histogram visualization is acceptable; however, after debugging, should go back to -9999
          twopeak <- -1
          
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
               # test in case exactly three peaks were not found
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
          
          # Used in issue #1: Recheck histogram values.  Will comment out after diagnostics
          ndwi_values <- data.frame(ndwi@data@values)
          ndwi_values <- rename(ndwi_values, data=ndwi.data.values)
          h <- ggplot(ndwi_values, aes(x=data)) +
               geom_histogram(breaks = (c(0:200)/100-1), color = "black", fill = "gray", na.rm = TRUE) +
               #geom_vline(aes(xintercept = twopeak), color = "green") +
               geom_vline(aes(xintercept = threepeak), color = "blue") +
               xlab("NDWI") +
               ylab("Count") +
               theme(panel.background = element_rect(fill = "white", colour = "black")) +
               theme(aspect.ratio = 1) +
               theme(axis.text = element_text(face = "plain", size = 12))
          ggsave(paste0(root,"hist.eps"), h, device = "eps", dpi = 72)
          
          # Water's Edge LOOP ENDS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          output[3] <-threepeak #for output file: value for the edge of water (3 peak)
          output[4] <-twopeak #for output file: value for the edge of water (2 peak)
          
          # Buffalo Creek
          # RDB=(609589.376, 4507801.407)
          # LDB=(609634.607, 4507831.586)
          x1 <- (609589.376)
          x2 <- (609634.607)
          y1 <- (4507801.407)
          y2 <- (4507831.586)
          # Mutale River 
          
          # Slopes:
          ma <- (y2-y1)/(x2-x1)
          #this next part will rely on UTM (the coordinates are in meters)
          ra <- (0.1) # r is the step size of each point along our width
          t <- sqrt(((x2-x1)^2)+ ((y2-y1)^2)) # length along search transect
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
          
          spat <- SpatialPoints(pointers)
          
          alng <- extract(ndwi, spat, method='simple')
          # plot(alng, xlab="Position along transect", ylab="NDWI")
          # To export table or NDWI v. position as a new file
          # write.table(alng, file = paste(root, "distwidth.csv", sep="."), append = TRUE, sep = ",", dec = ".", col.names = FALSE)
          
          RDB <- -9999 # preallocate in case of failed search algorithm
          LDB <- -9999
          alng_per <- array(-9, dim=c(f,2)) #allocation for the midpoints
          # when you reach -9 in that array, you've reached the end of the midpoints/values found
          restart <- 2 #initial start for i search
          for (j in 1:f){
               cnt <- 1
               for (i in restart:f) {
                    if (is.na(alng[i])==FALSE) {
                         if (alng[i]==alng[i-1]) { # determines if the next value is equal
                              cnt <- cnt + 1 # counts how many values there are
                         } else {
                              restart <- i + 1 #to keep moving forward from the last section without causing a loop at the end of it
                              break # breaks from current for loop.
                         }
                    }
               }
               mp <- ((cnt*ra)/2) #ra is the spacing, and mp gives the midpoint of the current distance section
               if (is.na(alng[i-1])==FALSE) {
                    if (i<(f)) {
                         alng_per[j,1] <- (((i-1)*ra)-mp) 
                         alng_per[j,2] <- alng[i-1] 
                    } else {
                         alng_per[j,1] <- ((f*ra)-mp) 
                         alng_per[j,2] <- alng[i-1] 
                    }
               }
               if (i>=(f)) {
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
          output[5] <- LDB #location in meters of bank 1
          output[6] <- RDB #location in meters of bank 2
          output[7] <- LDB-RDB #gives width in meters
     }
     #rm(alng_per,avg,dm,e,h,ndwi,nop,peaks,pointers,rbrick,rc,spat,test,a,alng,b,bins,c,cnt,f,fl,fn,goal,i,i1,i2,j,j1,j2,k,LDB,m,ma,mp,n,ra,rc2,rc4,RDB,restart,root,sec,t,thr,threepeak,twopeak,v,w,x1,x2,y1,y2) # this tried to remove vars that didnt exist... oops
     # for single string processing
     # for (i in 1:7) {
     #      widths[q,i] <- output[i]
     # }
     print(output)
}

dt <- as_date(as.numeric(widths[,1]))
filename <- widths[,2]
ndwi_threshold_3 <- as.numeric(widths[,3])
ndwi_threshold_2 <- as.numeric(widths[,4])
left_m <- as.numeric(widths[,5])
right_m <- as.numeric(widths[,6])
width_m <- as.numeric(widths[,7])
widths <- data.frame(dt,filename,ndwi_threshold_3,ndwi_threshold_2,left_m,right_m,width_m)
write_csv(widths, "widths.csv")

