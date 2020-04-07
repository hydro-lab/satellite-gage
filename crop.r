# Raster calcultion to crop Geotiff, also pulls XML metadata: reflectance coefficient (ps:reflectanceCoefficient)
library(rgdal)
library(raster)
library(rgeos)
library(XML)
library(methods)

# remember to set working directory if needed:
setwd("/Users/littlesunsh9/Documents/planet_order_181828/")
#Lists for necessary files
#Image list
i <- list.files("/Users/littlesunsh9/Documents/planet_order_181828/", pattern = "*AnalyticMS.tif$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)
#Metadata List
c <- list.files("/Users/littlesunsh9/Documents/planet_order_181828/", pattern = "*AnalyticMS_metadata.xml$", full.names = TRUE, recursive = TRUE, ignore.case=TRUE, include.dirs = TRUE)

# LOOP STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for (q in 1:(length(i))){

fn <- c[q]
fl <- xmlParse(fn)
rc <- setNames(xmlToDataFrame(node=getNodeSet(fl, "//ps:EarthObservation/gml:resultOf/ps:EarthObservationResult/ps:bandSpecificMetadata/ps:reflectanceCoefficient")),"reflectanceCoefficient")
dm <- as.matrix(rc)
# 1 Red
# 2 Green
# 3 Blue
# 4 Near infrared
rc2 <- as.numeric(dm[2]) # Green
rc4 <- as.numeric(dm[4]) # NIR

# Import raster image
fn <- i[q]
pic <- stack(fn)
# set extent from QGIS analysis:
# extent format (xmin,xmax,ymin,ymax)
e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
r <- crop(pic, e)

#options(stringsAsFactors = FALSE)

rbrick <- brick(r)
# calculate normalized difference water index (NDWI).  :
# calculate NDWI using the green (band 2) and nir (band 4) bands
ndwi <- ((rc2*r[[2]]) - (rc4*r[[4]])) / ((rc2*r[[2]]) + (rc4*r[[4]]))
# This formulation follows: Gao, B. (1996). NDWI—A normalized difference water index for remote sensing of vegetation liquid water from space. Remote Sensing of Environment, 53(3), p. 257-266. https://www.sciencedirect.com/science/article/abs/pii/S0034425796000673

# To view, during development
#plot(ndwi)

p <- strsplit(i[q], "_3B_AnalyticMS.tif")
r <- strsplit(p[[1]], "/")
lr <- tolower(r[[1]])
len <- length(lr)
root <- lr[[len]]

# To export NDWI as a new file
writeRaster(x = ndwi,
            filename= paste(root, "ndwi.tif", sep="."),
            format = "GTiff", # save as a tif
            # save as a FLOAT if not default, not integer
            overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.
#attempting the following crop procedure https://gis.stackexchange.com/questions/229356/crop-a-raster-file-in-r
#attempting the following NDVI as NDWI procedure https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/vegetation-indices-NDVI-in-R/

# LOOP ENDS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}
