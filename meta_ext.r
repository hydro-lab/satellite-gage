# XML extraction of the reflectance coefficient, ps:reflectanceCoefficient
# Tutorial at: https://www.youtube.com/watch?v=_ZQCdYdhHJI was very helpful for XML parsing.

library(XML)
library(methods)

# remember to set working directory if needed:
setwd("/Volumes/LaCie2big/widthflowproject/r_test")

fn <- "./20180212mXML.xml"
fl <- xmlParse(fn)
rc <- setNames(xmlToDataFrame(node=getNodeSet(fl, "//ps:EarthObservation/gml:resultOf/ps:EarthObservationResult/ps:bandSpecificMetadata/ps:reflectanceCoefficient")),"reflectanceCoefficient")
dm <- as.matrix(rc)
# 1 Red
# 2 Green
# 3 Blue
# 4 Near infrared
rc2 <- as.numeric(dm[2,2]) # Green
rc4 <- as.numeric(dm[4,2]) # NIR
