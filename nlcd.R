#### Code for prepping nlcd data
###  Mitchell O'Neill
##   2/2/2018

library(rgdal)
library(raster)

nlcd <- raster("C:/Users/Localadmin/Documents/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")
#readGDAL("C:/Users/Localadmin/Documents/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

proj4string(nlcd)

nlcd_reclass <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/nlcd_reclass.txt",header = F, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
nlcd_reclass

nlcd2 <- read.table("C:/Users/Localadmin/Desktop/arcStuff/outputs/nlcd/nlcd/w001000.adf",header = F, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)


head(nlcd)
str(nlcd)

#NLCD.2011.Land.Cover.Class
#head(nlcd$NLCD.2011.Land.Cover.Class)
#unique(nlcd$nlcd_2011_landcover_2011_edition_2014_10_10)

#nlcD <- as.data.frame(nlcd$nlcd_2011_landcover_2011_edition_2014_10_10[100:200,100:200])
#head(nlcD)
#summary(nlcD$`nlcd$nlcd_2011_landcover_2011_edition_2014_10_10[100:200, 100:200]`)

#head(nlcd@data)
#str(nlcd@data)
#nlcdata <- nlcd@data
#head(nlcdata@values)
#head(nlcdata@attributes)

nlcdata<- as.data.frame(nlcd@data@attributes)
nlcdata
#head(nlcd@data@attributes[1]$NLCD.2011.Land.Cover.Class)

