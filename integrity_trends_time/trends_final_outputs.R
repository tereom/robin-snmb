
# packages 
library("raster")
library("sp")
library("rgdal")

# working directory
setwd("D:/Julian/141_mapas_IE/trends")

# load base raster (to have uniform metadata)
rast <- raster("base.tif")

# load data
data <- read.table("trends.csv",sep=",",header=TRUE)

### linear trends

# linear trend variable
linear <- data$linear

# remove extremely small trends
linear[linear>-0.001&linear<0.001]<-NA

# remove trends with an unisignificant Mann-Kendall test
linear[data$sig>0.05]<-NA

# write trends raster to disk
output <- data.frame(x=data$x,y=data$y,linear)
coordinates(output)=~x+y
gridded(output)=TRUE
output <- raster(output)
projection(output)<- projection(rast)
rf <- writeRaster(output, filename="linear.tif", format="GTiff", overwrite=TRUE)

### reduce noise by removing significant trend-pixels which are isolated 
slope <- output

# minimum size of trend-clumps
size <- 4

# clumped pixels
clumpz <- clump(slope, directions=8)

# frecuency table of clumps
f <- data.frame(freq(clumpz))

# which rows of the data.frame are only represented by clumps under "size" pixels?
str(which(f$count <= size))

# which values do these correspond to?
str(f$value[which(f$count <= size)])

# put these into a vector of clump ID's to be removed
excludeID <- f$value[which(f$count <= size)]

# exclude these
slope[clumpz %in% excludeID] <- NA

# write clumped-significant trends to disk
rf <- writeRaster(slope, filename="linear_clumped.tif", format="GTiff", overwrite=TRUE)

### make shp files of significant-trend polygons

# negative polygons
slope_neg <- slope
slope_neg[slope_neg>=0]<-NA

Rclus <- clump(slope_neg,directions=8)
SPclus <- rasterToPolygons(Rclus, dissolve=TRUE)

# write out a new shapefile (including .prj component)
writeOGR(SPclus, ".", "negative_ie_trend", driver="ESRI Shapefile")

# positive polygons
slope_pos <- slope
slope_pos[slope_pos<=0]<-NA

Rclus <- clump(slope_pos,directions=8)
SPclus <- rasterToPolygons(Rclus, dissolve=TRUE)

# write out a new shapefile (including .prj component)
writeOGR(SPclus, ".", "positive_ie_trend", driver="ESRI Shapefile")
