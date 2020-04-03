# Raster calcultion to crop Geotiff, also pulls XML metadata: reflectance coefficient (ps:reflectanceCoefficient)

library(raster)
library(rgeos)
library(XML)
library(methods)

# remember to set working directory if needed:
setwd("/Volumes/LaCie2big/widthflowproject/r_test")

fn <- "./20180212_155035_0f3c_3B_AnalyticMS_metadata.xml"
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
fn <- "./20180212_155035_0f3c_3B_AnalyticMS.tif"
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

# To view, during development
#plot(ndwi)

# To export NDWI as a new file
writeRaster(x = ndwi,
            filename="test.tif",
            format = "GTiff", # save as a tif
            datatype='INT2S', # save as a INTEGER rather than a float
            overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.
#attempting the following crop procedure https://gis.stackexchange.com/questions/229356/crop-a-raster-file-in-r
#attempting the following NDVI as NDWI procedure https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/vegetation-indices-NDVI-in-R/
