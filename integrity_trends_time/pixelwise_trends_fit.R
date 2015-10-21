# packages
library("raster")
library("rgdal")
library("zyp")

rm(list = ls())

setwd("E:/work/20150822_IEmaps/maps_20151020/maps/")

# ie files
path = "E:/work/20150822_IEmaps/maps_20151020/maps/"
rasters = list.files(path,pattern="\\.tif$",full.names=TRUE)


# generate raster brick
iebrik <- brick()

# need to re-construct image paths as r reads filenames ina funny order
for (i in 2004:2013)
{
  for (j in 1:12)
  {
   name = paste("out_bn_ie_tabfinal_",toString(i),"_",toString(j),".tif",sep="")
   ierast <- raster(name)
   iebrik <- addLayer(iebrik,ierast)
  }
}


# brik to data table
ie_datatable = data.frame(rasterToPoints(iebrik))

# linear trends
trends_yuep <- zyp.trend.dataframe(indat=ie_datatable,metadata.cols=2, method="yuepilon", conf.intervals=TRUE)

# save to disk
write.table(trends_yuep,"trends.csv",sep=",",row.names=FALSE)
