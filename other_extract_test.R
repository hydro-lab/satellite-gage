# new function - extract along transect from raster
# extract, in image2width, was not working properly.  Trying this.

# copied from: https://rdrr.io/cran/inlmisc/man/ExtractAlongTransect.html, then annotated
library(inlmisc)

coords <- rbind(c(-100, -90), c(80, 90), c(80, 0), c(40, -40))
plot(coords)
crs <- sp::CRS("+init=epsg:4326")
transect <- sp::SpatialPoints(coords, proj4string = crs)

r <- raster::raster(nrows = 10, ncols = 10,
                    ymn = -80, ymx = 80, crs = crs)
names(r) <- "value"
set.seed(0)
r[] <- runif(raster::ncell(r)) # this appears to generate random values.
r[4, 6] <- NA  # count from top left, replaces value with NA
PlotMap(r) # couldn't get to work
plot(r) # worked

l <- sp::Lines(list(sp::Line(coords)), ID = "Transect")
lines(sp::SpatialLines(list(l), proj4string = crs)) # plots lines
points(transect, pch = 19) # adds points
segs <- ExtractAlongTransect(transect, r)
for (i in seq_along(segs)) points(segs[[i]])

dev.new()
xlab <- "Distance along transect"
ylab <- "Raster value"
xlim <- range(vapply(segs, function(i) {
     range(i@data[, "dist"])
}, c(0, 0)))
ylim <- range(vapply(segs, function(i) {
     range(i@data[, "value"], na.rm = TRUE)
}, c(0, 0)))
PlotGraph(NA, xlab = xlab, ylab = ylab,
          xlim = xlim, ylim = ylim, type = "n")
cols <- GetColors(length(segs), scheme = "bright")
for (i in seq_along(segs))
     lines(segs[[i]]@data[, c("dist", "value")],
           col = cols[i], lwd = 2)
coords <- sp::coordinates(transect)
n <- length(transect)
d <- cumsum(c(0, as.matrix(dist((coords)))[cbind(1:(n - 1), 2:n)]))
abline(v = d, lty = 2)
mtext(sprintf("(%d, %d)", coords[1, 1], coords[1, 2]),
      line = -1, adj = 0, cex = 0.7)
mtext(sprintf("(%d, %d)", coords[n, 1], coords[n, 2]),
      line = -1, adj = 1, cex = 0.7)

graphics.off()
