# This code forms a list of images, writes table to txt
library(lubridate)
library(readr)

im <- list.files("/Volumes/LaCie2big/BCAnalytic/2017/", 
                 pattern = "*AnalyticMS.tif$", 
                 full.names = TRUE, 
                 recursive = TRUE, 
                 ignore.case=TRUE, 
                 include.dirs = TRUE)

dt <- array("", dim = length(im))
for (i in 1:length(im)) {
      d <- strsplit(im[i], "_")
      dt[i] <- as_date(ymd(d[1]))
}

df <- data.frame(dt, im)
df <- df[order(dt),]

write_csv(df, "imagelist.csv")