# This code forms a list of images, writes table to txt
# Useful to check files and determine how many points will be used.

# Written in Matlab and R with RStudio by David Kahler and Mackenzie Martin, 
# Duquesne University, from 2018 to 2021.  The development was supported by 
# the United States Agency for International Development, Southern Africa 
# Regional Mission.  Further information is available at: 
# www.duq.edu/limpopo 
# https://github.com/LimpopoLab 

library(stringr)
library(lubridate)
library(readr)
library(parallel)

im <- list.files("/Volumes/T7/planet/mutale/data/", 
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

registerDoParallel(detectCores())
widths <- foreach (q = 1:(nrow(imagebank)), .combine = 'rbind') %dopar% {
     fn <- imagebank$im[q]
     pic <- stack(fn)
     e <- as(extent(245850, 246350, 7478700, 7479200), 'SpatialPolygons')
     crs(e) <- "+proj=utm +zone=36 +datum=WGS84" # may need negative y values
     test <- as(extent(pic), 'SpatialPolygons') # Extent of image
     crs(test) <- "+proj=utm +zone=36 +datum=WGS84"
     print(gCovers(test,e))
}

write_csv(imagebank, "imagelist.csv")
