## Prepping climate data
## Mitch O'Neill
## created: 8/8/2017
## last modified: 8/8/2017

library(rgeos)
library(rgdal)
library(raster)

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files")
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


############################################### 
######### loop through all bio layers #########
###############################################

clim.list <- list.files("Raw/climate_data/") ## list of all files in Raw climate data folder
                                             ## one for each climate model, each with the 19 variables within

for (j in 1:length(clim.list)) { ## iterate through each climate model

  model.folder.path <- paste("Raw/climate_data", clim.list[j], sep="/") 
  ## look within the climate model of the iteration, and   
  bio.list <- list.files(model.folder.path) ## list all of the variables
  bio.list <- bio.list[1:19] ## exlude the metadata files after the 19 bioclimate variables
  
  clipped.model.folder.path <-paste("Products/clipped_climate_data", clim.list[j], sep="/")
  dir.create(clipped.model.folder.path) ## creates a new folder to place the clipped data for the model of the iteration

  for (i in 1:length(bio.list)){ ## loop through each climate variable within the model of the iteration
    
    file.in <- paste(model.folder.path, bio.list[i], sep="/") 
    ## paste together the full file name to bring it into R

    file.out <- paste(clipped.model.folder.path, (paste(bio.list[i], "asc", sep=".")), sep="/")
    ## create filepath for clipped data within the new folder, to export from R
  
    bio.i <- raster(file.in) ## read in file as raster
    
    bio.i <- mask(crop(bio.i, states.buffer3t), states.buffer3t) 
    ## crop raster to extent of united states, and then clip(mask) to the buffer
    
    writeRaster(bio.i, file.out, format= "ascii") ## write out raster as an ESRI raster
   
  print(i) ## keep track of place within model
    
  }
  
  print(clim.list[j]) ## keep track of model within entire loop
}




##########################
### reading in fishnet ###
##########################
fishnet = readOGR(dsn = "Raw/Fishnet", layer = "fishnet")
plot(fishnet)