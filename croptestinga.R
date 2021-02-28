#Run in command line as: Rscript master.R

#Current Update is to inclide a check for existance of the crop boundaries in each raster file else skip

# Raster calcultion to crop Geotiff, also pulls XML metadata: reflectance coefficient (ps:reflectanceCoefficient)
library(rgdal)
library(raster)
library(rgeos)
library(XML)
library(methods)
library(sp)

# remember to set working directory if needed:
setwd("/Users/littlesunsh9/Documents/planettestingfolder/")
setwd("/Volumes/LaCie/planet/buffalo/")
#Lists for necessary files
#Image list
im <- list.files("/Users/littlesunsh9/Documents/planettestingfolder/", pattern = "*AnalyticMS.tif$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)
#Metadata List
g <- list.files("/Users/littlesunsh9/Documents/planettestingfolder/", pattern = "*AnalyticMS_metadata.xml$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)

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
    # for testing existance in the workflow
    fn <- "/Volumes/LaCie/planet/buffalo/exist_testing/20180422_153434_0e20/20180422_153434_0e20_3B_AnalyticMS_metadata.xml"
    #fn <- "/Volumes/LaCie/planet/buffalo/exist_testing/20180424_074548_1003/20190514_074548_1003_3B_AnalyticMS_metadata.xml"
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
    fn <- "/Volumes/LaCie/planet/buffalo/exist_testing/20180422_153434_0e20/20180422_153434_0e20_3B_AnalyticMS.tif"
    fn <- "/Volumes/LaCie/planet/buffalo/exist_testing/20180424_074548_1003/20190514_074548_1003_3B_AnalyticMS.tif"
    pic <- stack(fn)
    # set extent from QGIS analysis: !! We need this area for analysis
    # extent format (xmin,xmax,ymin,ymax)
    e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons')
    #crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
    crs(e) <- "+proj=utm +zone=17 +datum=WGS84"
    # Set extent from the Planet file !! This is the area from the picture
    test <- as(extent(pic), 'SpatialPolygons')
    crs(test) <- "+proj=utm +zone=17 +datum=WGS84"
    if (gWithin(e, test, byid = FALSE)) {
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
    }
}


#setwd("/Users/littlesunsh9/Documents/planettestingfolder/")

#test <- "20180411_154228_0f47_3B_AnalyticMS.tif")
#meta <- 
    
    
#    gWithin(spgeom1, spgeom2 = NULL, byid = FALSE, returnDense=TRUE, checkValidity=FALSE)
