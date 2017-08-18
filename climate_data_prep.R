#### Prepping climate data
###  Mitch O'Neill
##   created: 8/8/2017
#    last modified: 8/8/2017

library(rgeos)
library(rgdal)
library(raster)

setwd("C:/Users/mwone/Documents/geodata")
states = readOGR(dsn = "states", layer = "US_states") ## read in US shapefile
proj4string(states)


## rgeos expects projected coordinates for buffers, so I need to transform it
EPSG <- make_EPSG() ## creates library of EPSG codes to help assign proj4string
EPSG[grepl("WGS 84$", EPSG$note) == TRUE, ] ## search for WGS 84, #4326
EPSG[grepl("Albers", EPSG$note) == TRUE, ] ## try Conus Albers, NAD83; 5070
EPSG[grepl("5070", EPSG$code) == TRUE, ] ## Conus Albers

states.t <- spTransform(states, "+init=epsg:5070") ## transform to Conus Albers

states.buffer3 <- gBuffer(states.t, byid=F, id=NULL, width = 3000, quadsegs = 8, capStyle = "ROUND", joinStyle = "ROUND")
## US with 3km buffer

## transform buffer to proj4string of the bioclimatic data
states.buffer3t <- spTransform(states.buffer3, proj4string(raster("climate_data/current/bio_1")))

## convert to SpatialPolygonsDataframe so that buffer can be written out
fill.in.data <- data.frame("fill in data") ## filler data to make up dataframe component

row.names(fill.in.data) <- sapply(slot(states.buffer3t, "polygons"), function(x) slot(x, "ID")) 
## row names have to match the ID of the polygon

buffer3 <- SpatialPolygonsDataFrame(states.buffer3t, data=fill.in.data) 
## combine polygon and dataframe into spatial dataframe

writeOGR(buffer3, dsn = "clipped/buffer3km", layer = "buffer3", driver = "ESRI Shapefile")
## write out as ESRI shapefile



############################################### 
######### loop through all bio layers #########
###############################################

clim.list <- list.files("climate_data/") ## list of all files in Raw climate data folder
                                             ## one for each climate model, each with the 19 variables within

for (j in 1:length(clim.list)) { ## iterate through each climate model

  model.folder.path <- paste("climate_data", clim.list[j], sep="/") 
  ## look within the climate model of the iteration, and   
  bio.list <- list.files(model.folder.path) ## list all of the variables
  bio.list <- bio.list[1:19] ## exlude the metadata files after the 19 bioclimate variables
  
  clipped.model.folder.path <-paste("clipped_climate_data", clim.list[j], sep="/")
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
fishnet = readOGR(dsn = "Fishnet", layer = "fishnet")
plot(fishnet)


#################################################
### Repeat clipping with all of North America ###
#################################################

north.america <- as(extent(-170, -52, 7, 70), 'SpatialPolygons')
#model.folder.path <- paste("climate_data", clim.list[j], sep="/") 
## look within the climate model of the iteration, and   

model.folder.path <- "climate_data/current/"
bio.list <- list.files(model.folder.path) ## list all of the variables
bio.list <- bio.list[1:19] ## exlude the metadata files after the 19 bioclimate variables

#clipped2.model.folder.path <-paste("clipped2_climate_data", clim.list[j], sep="/")
clipped2.model.folder.path <-paste("clipped2_climate_data", "current", sep="/")

dir.create(clipped2.model.folder.path) ## creates a new folder to place the clipped data for the model of the iteration

for (i in 1:length(bio.list)){ ## loop through each climate variable within the model of the iteration
  
  file.in <- paste(model.folder.path, bio.list[i], sep="/") 
  ## paste together the full file name to bring it into R
  
  file.out <- paste(clipped2.model.folder.path, (paste(bio.list[i], "asc", sep=".")), sep="/")
  ## create filepath for clipped data within the new folder, to export from R
  
  bio.i <- raster(file.in) ## read in file as raster
  
  bio.i <- crop(bio.i, north.america) 
  ## crop raster to extent of united states, and then clip(mask) to the buffer
  
  writeRaster(bio.i, file.out, format= "ascii") ## write out raster as an ESRI grid
  
  print(i) ## keep track of place within model
  
}





##########################################
######### assessing correlations #########
##########################################

bio.list<- list.files("clipped_climate_data/current", full.names =T) 
## lists all 19 climate rasters in the climate data to be used for fitting

bio.stack <- stack(bio.list) ## creates an object with each of the nine rasters stacked
proj4string(bio.stack) <- proj4string(raster("climate_data/current/bio_1")) 
## assigns reference info based the files where these files came from 

bio.stackM <- as.matrix(bio.stack) ## converts the stacked object into a matrix statistical analysis

corr.matrix <- corr.test(bio.stack, y = NULL, use = "pairwise", method="spearman", adjust="bonferroni", alpha=.05, ci=F)
## matrix of spearman's correlation of the 19 variables, p-value adjusted using Bonferroni correction

write.csv(data.frame(corr.matrix$p), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/p.vals.csv") 
## export p-values in a table format
write.csv(data.frame(corr.matrix$r), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/corr.vals.csv") 
## export correlation coefficients in a table format

##############################################################
######### repeat correlations with north america data ########
##############################################################

bio.list<- list.files("clipped2_climate_data/current", full.names =T) 
## lists all 19 climate rasters in the climate data to be used for fitting

bio.stack <- stack(bio.list) ## creates an object with each of the nine rasters stacked
proj4string(bio.stack) <- proj4string(raster("climate_data/current/bio_1")) 
## assigns reference info based the files where these files came from 

bio.stackM <- as.matrix(bio.stack) ## converts the stacked object into a matrix statistical analysis

corr.matrix <- corr.test(bio.stackM, y = NULL, use = "pairwise", method="spearman", adjust="bonferroni", alpha=.05, ci=F)
## matrix of spearman's correlation of the 19 variables, p-value adjusted using Bonferroni correction

write.csv(data.frame(corr.matrix$p), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/p.valsNA.csv") 
## export p-values in a table format
write.csv(data.frame(corr.matrix$r), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/corr.valsNA.csv") 
## export correlation coefficients in a table format



###################################################################################
#### MESS ANALYSIS ################################################################
###################################################################################

library(dismo)




