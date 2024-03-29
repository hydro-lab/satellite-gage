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
library(rgdal) # must change to GDAL and PROJ: sf/stars/terra, by 2023
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
setwd("/Volumes/dmk/planet/mutale/data/")

# Lists for necessary files
# Image list
im <- list.files(".", 
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
md <- list.files(".", 
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

# Store geographic data based on crop extent from QGIS analysis:
# extent format (xmin,xmax,ymin,ymax)
## Buffalo Creek:
#e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons') # Extent needed
#crs(e) <- "+proj=utm +zone=17 +datum=WGS84"
## Mutale River downstream
e <- as(extent(245850, 246350, 7478700, 7479200), 'SpatialPolygons')
crs(e) <- "+proj=utm +zone=36 +datum=WGS84" # may need negative y values

# LOOP STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# SINGLE
#widths <- array(NA, dim = c((nrow(imagebank)),7))
#for (q in 1:(2)) { # original loop
# PARALLEL
registerDoParallel(detectCores())
widths <- foreach (q = 1:(nrow(imagebank)), .combine = 'rbind') %dopar% { # parallel computing loop: this changes how data are transferred back from each operation.
     output <- array(NA, dim = 6) # output array - will be filled in if data are valid
     output[1] <- date(as_datetime(imagebank$dt[q])) # to check, use: as_date(output[1])
     
     # To create filename root and establish root for exports
     p <- strsplit(imagebank$im[q], "_3B_AnalyticMS.tif")
     r <- strsplit(p[[1]], "/")
     lr <- tolower(r[[1]])
     len <- length(lr)
     root <- lr[[len]]
     rm(p,r,lr,len)
     output[2] <- root #for output file: root name of image
     # to output files to the exact same image directory, use root <- p[[1]]
     
     #Import raw Planet metadata to get the reflectance coefficients
     fn <- imagebank$md[q]
     fl <- xmlParse(fn)
     rc <- setNames(xmlToDataFrame(node=getNodeSet(fl, "//ps:EarthObservation/gml:resultOf/ps:EarthObservationResult/ps:bandSpecificMetadata/ps:reflectanceCoefficient")),"reflectanceCoefficient")
     dm <- as.matrix(rc)
     rc2 <- as.numeric(dm[2]) # Green
     rc4 <- as.numeric(dm[4]) # NIR
     rm(fl, rc, dm)
     
     # Import raster image, crops to chosen extent
     fn <- imagebank$im[q]
     pic <- stack(fn)
     
     # Set extent from the Planet file !! This is the area from the picture
     test <- as(extent(pic), 'SpatialPolygons') # Extent of image
     #crs(test) <- "+proj=utm +zone=17 +datum=WGS84"
     crs(test) <- "+proj=utm +zone=36 +datum=WGS84"
     if (gCovers(test,e)) { # returns TRUE if no point in spgeom2 (e, needed) is outside spgeom1 (test, image extent) # used to be (gWithin(e, test, byid = FALSE))
          r <- crop(pic, e)
          rm(pic) # remove rest of image from RAM
          rbrick <- brick(r)
          # calculate NDWI using the green (band 2) and nir (band 4) bands
          ndwi <- ((rc2*r[[2]]) - (rc4*r[[4]])) / ((rc2*r[[2]]) + (rc4*r[[4]]))
          # plot(ndwi) # for viewing during development
          
          # writeRaster(x = ndwi, ## this does not need to be done, just a nice record.
          #             filename= paste(root, "cndwi.tif", sep="."),
          #             format = "GTiff", # save as a tif, save as a FLOAT if not default, not integer
          #             overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.
          
          # This code finds the boundary of the water in a normalized difference water index 
          # This code uses the cropped, single-layer, NDWI image.  Image based on the histogram of the pixel values.
          # Import raster image, or take it from previous code, set working directory, if needed.
          
          h = hist(ndwi, # built-in histogram function.  To find values only.  Plotting is at the end of this loop.
                   breaks=seq(-1,1,by=0.01),
                   plot=FALSE) 
          bins <- h$mids # positions    number
          v <- h$counts # counts        integer
          rm(h)
          
          if ((mean(v)) > 1) { # determine if data exist.
               
               # Allocate arrays used in analysis
               maxWindow <- 10 # This is the control on the maximum averaging window AND the size of the following arrays.
               binNumber <- length(bins)
               avg <- array(0, dim = c(binNumber,maxWindow))
               peaks <- array(0, dim = c(binNumber,maxWindow))
               nop <- array(0, dim = c(1,maxWindow))
               for (w in 1:maxWindow){
                    # filter values (v=h$counts) with the averaging window size 2*w+1
                    for (k in (w+1):(binNumber-w)){
                         avg[k,w] <- ((sum(v[(k-w):(k+w)]))/((2*w)+1))
                    }
                    # identify and number peaks
                    cnt <- 0
                    for (j in (w+1):(binNumber-w)){
                         if ((avg[j-1,w])<(avg[j,w])){
                              if ((avg[j+1,w])<(avg[j,w])){
                                   cnt <- (cnt+1)
                                   peaks[j,w] <- cnt
                                   nop[1,w] <- cnt
                              }
                         }
                    }
               }
               
               # AVERAGING VISUALIZATION
               # win <- 6 # the location of the averaging window in avg variable; lower number is small window, larger is smoother
               # troubleshoot <- data.frame(bins,avg[,win])
               # troubleshoot <- rename(troubleshoot, avg = `avg...win.`) # may need to update original variable
               # ggplot(troubleshoot) +
               #      geom_line(aes(x=bins,y=avg)) +
               #      #geom_vline(aes(xintercept=-0.3), color="Blue") +
               #      xlab("NDWI") +
               #      ylab("Count") +
               #      theme(panel.background = element_rect(fill = "white", colour = "black")) +
               #      theme(aspect.ratio = 1) +
               #      theme(axis.text = element_text(face = "plain", size = 12))
               
               # Find first smoothed single-peak avg data
               singleWindow <- -9 # signal value
               for (w in 1:maxWindow) {
                    if (nop[1,w] == 1) {
                         singleWindow <- w
                         peakIndex <- which(peaks[,singleWindow]==1) # which index is the peak
                         peakValue <- bins[peakIndex]
                         break
                    }
               }
               if (singleWindow == -9) { # Only if there is no singleWindow found.
                    #print(root)
                    #print("Averaging failed to find single peak")
                    peakIndex <- 75
                    peakValue <- -0.5
               }
               
               # Find the smooth tail (on the right/positive side of the distribution)
               ndwiSlope <- array(0, dim = c(binNumber)) # derivative of smoothed NDWI histogram/distribution
               for (w in 2:(binNumber-1)) {
                    ndwiSlope[w] <- (avg[w+1,singleWindow] - avg[w-1,singleWindow]) / (bins[w+1] - bins[w-1])
               }
               slopeLimit <- 0.05 * max(abs(ndwiSlope)) # threshold: 5% of max slope
               flatIndex <- 110
               flatValue <- 0.1 # default value to produce a threshold of -0.2
               for (w in (peakIndex+1):binNumber) {
                    if (abs(ndwiSlope[w]) < slopeLimit) {
                         flatIndex <- w
                         flatValue <- bins[flatIndex]
                         break
                    }
               }
               
               # Average the peak and flat values
               ndwiThreshold <- (peakValue+flatValue)/2
               
               # HISTOGRAM VISUALIZATION
               ndwi_values <- data.frame(ndwi@data@values)
               smooth <- data.frame(bins,avg[,singleWindow])
               # h <- ggplot(ndwi_values, aes(x=ndwi.data.values)) +
               #      geom_histogram(breaks = (c(0:200)/100-1), color = "black", fill = "gray", na.rm = TRUE) +
               #      geom_line(data = smooth, aes(x=bins, y= `avg...singleWindow.`), color = "blue") +
               #      geom_vline(aes(xintercept = ndwiThreshold), color = "green") +
               #      xlim(c(-1,1)) +
               #      ylim(c(0,2000)) +
               #      xlab("NDWI") +
               #      ylab("Count") +
               #      theme(panel.background = element_rect(fill = "white", colour = "black")) +
               #      theme(aspect.ratio = 1) +
               #      theme(axis.text = element_text(face = "plain", size = 12))
               # ggsave(paste0(root,".hist.eps"), h, device = "eps", dpi = 72)
               
               # Water's Edge LOOP ENDS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               
               output[3] <- ndwiThreshold
               
               # Buffalo Creek (UTM: 17N)
               #crs <- sp::CRS("+proj=utm +zone=17 +datum=WGS84")
               #x1 <- (609589.376)
               #x2 <- (609634.607)
               #y1 <- (4507801.407)
               #y2 <- (4507831.586)
               # Mutale River downstream (EPSG: 32736, UTM: 36S)
               crs <- sp::CRS("+proj=utm +zone=36 +datum=WGS84")
               x1 <- (246130)
               x2 <- (246066)
               y1 <- (7478894)
               y2 <- (7479006)
               
               # Slopes:
               m <- (y2-y1)/(x2-x1)
               # this is in meters (UTM coordinates)
               dt = 0.1 # step size along transect (m)
               dx = sqrt((dt^2)/((m^2)+1)) # step size in the x-direction
               t <- sqrt(((x2-x1)^2) + ((y2-y1)^2)) # length along search transect
               f <- ceiling(t/dt) # number of steps needed to cross transect
               if (x1 > x2) { # determine which direction the transect should head.
                    dx <- (-1)*dx
               }
               
               transect <- array(999.999, dim = c(f,2)) # [,1] is x, [,2] is y
               transect[1,1] <- x1
               transect[1,2] <- y1
               for (i in 2:f){
                    x <- transect[(i-1),1] # last x value
                    y <- transect[(i-1),2] # last y value
                    transect[i,1] <- x + dx # next x value
                    transect[i,2] <- y + (m * (transect[i,1] - x)) # next y value
               }
               
               endPoints <- SpatialPoints(rbind(c(x1,y1),c(x2,y2)), proj4string = crs)
               # plot(ndwi)
               # points(endPoints)
               spat <- SpatialPoints(transect, proj4string = crs)
               # points(spat)
               
               alng <- extract(ndwi, spat, method='simple')
               
               if (is.na(mean(alng))==FALSE) { # determine if data exist along alng
                    
                    # plot(alng, xlab="Position along transect", ylab="NDWI")
                    # To export table or NDWI v. position as a new file
                    # write.table(alng, file = paste(root, "distwidth.csv", sep="."), append = TRUE, sep = ",", dec = ".", col.names = FALSE)
                    
                    RDB <- -9999 # preallocate in case of failed search algorithm
                    LDB <- -9999
                    alng_per <- array(-9, dim=c(f,2)) # allocation for the midpoints
                    # when you reach -9 in that array, you've reached the end of the midpoints/values found
                    restart <- 2 # initial start for i search
                    for (j in 1:f){
                         cnt <- 1
                         for (i in restart:f) {
                              if ((is.na(alng[i])==FALSE) & (is.na(alng[i-1])==FALSE)) {
                                   if (alng[i]==alng[i-1]) { # determines if the next value is equal
                                        cnt <- cnt + 1 # counts how many values there are
                                   } else {
                                        restart <- i + 1 #to keep moving forward from the last section without causing a loop at the end of it
                                        break # breaks from current for loop.
                                   }
                              }
                         }
                         mp <- ((cnt*dt)/2) #ra is the spacing, and mp gives the midpoint of the current distance section
                         if (is.na(alng[i-1])==FALSE) {
                              if (i<(f)) {
                                   alng_per[j,1] <- (((i-1)*dt)-mp) 
                                   alng_per[j,2] <- alng[i-1] 
                              } else {
                                   alng_per[j,1] <- ((f*dt)-mp) 
                                   alng_per[j,2] <- alng[i-1] 
                              }
                         }
                         if (i>=(f)) {
                              break
                         } 
                    }
                    
                    # Find starting point in the alng_per array:
                    for (i in 1:f) {
                         if (alng_per[i,1] < 0) {
                              fin <- i - 1
                              break
                         }
                    }
                    for (i in (2:fin)){
                         if (alng_per[i,2]>ndwiThreshold){
                              if (alng_per[i-1,2]<ndwiThreshold){
                                   i1 <- alng_per[i-1,1]
                                   i2 <- alng_per[i,1]
                                   j1 <- alng_per[i-1,2]
                                   j2 <- alng_per[i,2]
                                   n <- ndwiThreshold    
                                   RDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
                                   break
                              }
                         }
                    }
                    for (i in 1:(fin-1)){
                         if (alng_per[fin-i,2]>ndwiThreshold){        #expressing the index such that when i = 1, f, and when i = 2, f-1.
                              if (alng_per[fin-i+1,2]<ndwiThreshold){
                                   i1 <- alng_per[fin-i+1,1]
                                   i2 <- alng_per[fin-i,1]
                                   j1 <- alng_per[fin-i+1,2]
                                   j2 <- alng_per[fin-i,2]
                                   n <- ndwiThreshold    
                                   LDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
                                   break
                              }
                         }
                    }
                    output[4] <- LDB #location in meters of bank 1
                    output[5] <- RDB #location in meters of bank 2
                    if ((LDB > 0) & (RDB > 0)) {
                         riverWidth <- abs(LDB - RDB)
                    } else {
                         riverWidth <- -9999
                    }
                    output[6] <- riverWidth #gives width in meters
               } # end of if is.na(mean(alng))==FALSE to determine if data exist along alng
          } # end of if hist is zero/NA - no data exist
     } # end of if gCovers
     
     #rm(alng_per,avg,dm,e,h,ndwi,nop,peaks,pointers,rbrick,rc,spat,test,a,alng,b,bins,c,cnt,f,fl,fn,goal,i,i1,i2,j,j1,j2,k,LDB,m,ma,mp,n,ra,rc2,rc4,RDB,restart,root,sec,t,thr,threepeak,twopeak,v,w,x1,x2,y1,y2) # this tried to remove vars that didnt exist... oops
     # for single string processing
     # for (i in 1:7) {
     #      widths[q,i] <- output[i]
     # }
     print(output)
}

dt <- as_date(as.numeric(widths[,1]))
filename <- widths[,2]
ndwi_threshold <- as.numeric(widths[,3])
left_m <- as.numeric(widths[,4])
right_m <- as.numeric(widths[,5])
width_m <- as.numeric(widths[,6])
widths <- data.frame(dt,filename,ndwi_threshold,left_m,right_m,width_m)
write_csv(widths, "/Users/davidkahler/Documents/R/satellite-gage/widths.csv")
print("Program complete")

