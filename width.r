library(rgdal)
library(raster)
library(rgeos)
library(XML)
library(methods)

# remember to set working directory if needed:
setwd("/Users/littlesunsh9/Documents/planet_order_181828/")

ndwia <- stack("/Users/littlesunsh9/Documents/planet_order_181828/20180411_154228_0f47.cndwi.tif")
e <- as(extent(609555.5999,609709.1999,4507753.099,4507867.5999 ), 'SpatialPolygons')
crs(e) <- "+proj=utm +zone=17 +datum=WGS84"
#RDB=(609589.376, 4507801.407)
#LDB=(609634.607, 4507831.586)
x1 <- (609589.376)
x2 <- (609634.607)
y1 <- (4507801.407)
y2 <- (4507831.586)
m <- (y2-y1)/(x2-x1)
#this next part will rely on UTM (the coordinates are in meters)
r <- (0.1)
#r is the step size of each point along our width
t <- sqrt(((x2-x1)^2)+ ((y2-y1)^2))
f <- ceiling(t/r)
pointers <- array(999.999, dim = c(f,2))
pointers[1,1] <- x1
pointers[1,2] <- y1
for (i in 2:f){
  a <- 1
  b <- (-2)*pointers[i-1,1]
  c <- (pointers[i-1,1]^2) - (r^2)/((m^2)+1)
  pointers[i,1] <- ((-b)+(sqrt((b^2)-4*a*c)))/(2*a)
  pointers[i,2] <- ((pointers[i-1,2])+(m*((pointers[i,1])-(pointers[i-1,1]))))
}

#Check:
#x <- c(x1,x2)
#y <- c(y1,y2)
#plot(x,y)

P <- SpatialPoints(pointers)

alng <- extract(ndwia, P, method='simple')
plot(alng)
