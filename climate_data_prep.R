## Prepping climate data
## Mitch O'Neill
## created: 8/8/2017
## last modified: 8/8/2017

library(rgeos)
library(rgdal)
library(raster)

setwd("C:/Users/mwone/Google Drive/WorkingFolder/Predictors/")
states = readOGR(dsn = "states", layer = "US_states")
proj4string(states)

#states.buffer <- gBuffer(states, byid=F, id=NULL, width = .2, quadsegs = 8, capStyle = "ROUND", joinStyle = "ROUND")

EPSG <- make_EPSG() ## creates library of EPSG codes to help assign proj4string
EPSG[grepl("WGS 84$", EPSG$note) == TRUE, ] ## search for WGS 84, #4326
EPSG[grepl("Albers", EPSG$note) == TRUE, ] ## try Conus Albers, NAD83; 5070
EPSG[grepl("5070", EPSG$code) == TRUE, ] ## Conus ALbers

states.t <- spTransform(states, "+init=epsg:5070")

states.buffer3 <- gBuffer(states.t, byid=F, id=NULL, width = 3000, quadsegs = 8, capStyle = "ROUND", joinStyle = "ROUND")


######## clip bioclim data by states buffers 10 and 15 ########
bio_1 <- raster("C:/Users/mwone/Google Drive/WorkingFolder/Predictors/bio/bio_1")


states.buffer3t <- spTransform(states.buffer3, proj4string(bio_1))

fill.in.data <- data.frame("fill in data")
row.names(fill.in.data) <- sapply(slot(states.buffer3t, "polygons"), function(x) slot(x, "ID")) 

buffer3 <- SpatialPolygonsDataFrame(states.buffer3t, data=fill.in.data) #, proj4string= states.buffer3t@proj4string)

writeOGR(buffer3, dsn = "clipped/buffer3km", layer = "buffer3", driver = "ESRI Shapefile")


bio1_3 <- mask(crop(bio_1, states.buffer3t), states.buffer3t)


## writeout bio1_3
writeRaster(bio1_3, filename = "C:/Users/mwone/Google Drive/WorkingFolder/Predictors/clipped/bio1_3.asc",format= "ascii", overwrite=T)





######################
### reading in fishnet

fishnet = readOGR(dsn = "Fishnet", layer = "fishnet")
#plot(fishnet, add=T)
#plot(buffer3, col="red", add=T)
#proj4string(fishnet)
#proj4string(bio1_3)
#fishnet.t <- spTransform(fishnet, proj4string(buffer3))
#
#fishnet_sub <- fishnet.t[buffer3,]
#plot(fishnet_sub, add=T, col="blue")
#
#
#
##So to get all ID's
#state48<- states[48, ]
#plot(state48)
#fishnet.t <- spTransform(fishnet, proj4string(state48))
#
#fishnet_sub <- fishnet.t[state48,]
#plot(fishnet_sub, add=T, col="green")
#fishnet_sub <- mask(fishnet.t, state48)

#################### loop through all bio layers
setwd("C:/Users/mwone/Google Drive/WorkingFolder/Predictors/")

Clim.list <- list.files("climate_data/")

for (j in 1:length(clim.list)) {

  bio.list<- list.files("bio/")
  bio.list <- bio.list[1:19]
  
  for (i in 1:length(bio.list)){
    
    file.in <- paste(clim.list[j], bio.list[i], sep="/")
    dir.create(paste("climate_data_clipped", clim.list[j], sep="/"))
    file.out <- paste("climate_data_clipped", paste(clim.list[j],(paste(bio.list[i], "asc", sep=".")), sep="/"), sep = "/")
  
    bio.i <- raster(file.in)
    
    bio.i <- mask(crop(bio.i, states.buffer3t), states.buffer3t)
    
    writeRaster(bio.i, file.out, format= "ascii")
   
    #clean out garbage
    rm(bio.i)  
    
    print(i)
    
  }
  
  print(clim.list[j])
}


