#### Prepping climate data
###  Mitch O'Neill
##   created: 8/8/2017
#    last modified: 9/7/2017

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



################################################################ 
######### loop through all bio layers  and clip to L48 #########
################################################################

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

#################################################
### Repeat clipping with all of North America ###
#################################################

ecoregions <- readOGR(dsn="ecoregions", layer = "NA_CEC_Eco_Level1")
## read in CEC ecoregions shapefile

states <- spTransform(states, proj4string(ecoregions)) ## transform states to ecoregions CRS

L48ecoregions <- ecoregions[states,] ## subset polygons which coincide with US

eco <- unique(as.character(L48ecoregions@data$NA_L1NAME)) ## list of ecoregions which coincide w/ L48
eco <- eco[2:11] ## exclude water (great lakes)
L48ecoregions <- ecoregions[as.character(ecoregions$NA_L1NAME) %in% eco, ]
## subset ecoregions on the list of those that coincide with L48
## this convoluted method ensures that discontinuous "islands" of ecoregions of interest are not excluded

L48ecoregions <- gUnaryUnion(L48ecoregions) ## dissolve borders between ecoregions
L48ecoregions <- spTransform(L48ecoregions, "+init=epsg:5070") ## transform to albers = area CONUS

## construct dataframe to append to polygons, because only spatialpolygonsdataframes can be written out
## row names must match polygon id
sapply(slot(L48ecoregions, "polygons"), function(x) slot(x, "ID"))
filler <- "Filler data"
L48d <- data.frame(filler)
L48ecoregions <- SpatialPolygonsDataFrame(L48ecoregions, L48d)
dir.create("L48ecoregions")
writeOGR(L48ecoregions, dsn="L48ecoregions", layer="L48ecoregions", driver = "ESRI Shapefile")

eco <- readOGR(dsn="L48ecoregions", layer="L48ecoregions") ## read in file created above
eco <- spTransform(eco, proj4string(raster("climate_data/current/bio_1"))) ## transofrm to crs of climate data

eco.01    <- gSimplify(eco, tol=0.01, topologyPreserve=T)  ## simplify polygon to save processing time
sapply(slot(eco.01, "polygons"), function(x) slot(x, "ID")) ## find id of polygon
row.names(L48d) <- 0 ## match row.names of filler data to polygon IDs
eco.01    <- SpatialPolygonsDataFrame(eco.01, L48d) ## append filler data to polygons
writeOGR(eco.01, dsn="eco.01", layer="eco.01", driver= "ESRI Shapefile") ## write out file
  
## eco.01 was brought into ArcGIS and a 4km buffer was used to be read in below
L48ecoregions <- readOGR(dsn="C:/Users/mwone/Documents/geodata/eco01buf4k", layer="eco01buf4k")
dir.create("climate_data_clipped_to_ecoregion")
clim.list <- list.files("climate_data/") ## list of all files in Raw climate data folder
## one for each climate model, each with the 19 variables within

for (j in 5:5){
#for (j in 1:length(clim.list)) { ## iterate through each climate model
  
  model.folder.path <- paste("climate_data", clim.list[j], sep="/") 
  ## look within the climate model of the iteration, and   
  bio.list <- list.files(model.folder.path) ## list all of the variables
  bio.list <- bio.list[1:19] ## exlude the metadata files after the 19 bioclimate variables
  
  clipped.model.folder.path <-paste("climate_data_clipped_to_ecoregion", clim.list[j], sep="/")
  dir.create(clipped.model.folder.path) ## creates a new folder to place the clipped data for the model of the iteration
  
  for (i in 1:length(bio.list)){ ## loop through each climate variable within the model of the iteration
    
    file.in <- paste(model.folder.path, bio.list[i], sep="/") 
    ## paste together the full file name to bring it into R
    
    file.out <- paste(clipped.model.folder.path, (paste(bio.list[i], "asc", sep=".")), sep="/")
    ## create filepath for clipped data within the new folder, to export from R
    
    bio.i <- raster(file.in) ## read in file as raster
    
    bio.i <- mask(crop(bio.i, L48ecoregions), L48ecoregions) 
    ## crop raster to extent of united states, and then clip(mask) to the buffer
    
    writeRaster(bio.i, file.out, format= "ascii") ## write out raster as an ESRI grid
    
    print(i) ## keep track of place within model
    
  }
  
  print(clim.list[j])
}
##########################################
######### assessing correlations #########
##########################################

bio.list<- list.files("climate_data_clipped_to_ecoregion/current", full.names =T) 
## lists all 19 climate rasters in the climate data to be used for fitting

bio.stack <- stack(bio.list) ## creates an object with each of the nine rasters stacked
proj4string(bio.stack) <- proj4string(raster("climate_data/current/bio_1")) 
## assigns reference info based the files where these files came from 

bio.stackM <- as.matrix(bio.stack) ## converts the stacked object into a matrix statistical analysis

library(psych)
corr.matrix <- corr.test(bio.stackM, y = NULL, use = "pairwise", method="spearman", adjust="bonferroni", alpha=.05, ci=F)
## matrix of spearman's correlation of the 19 variables, p-value adjusted using Bonferroni correction

write.csv(data.frame(corr.matrix$p), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/p.valsECO.csv") 
## export p-values in a table format
write.csv(data.frame(corr.matrix$r), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/corr.valsECO.csv") 
## export correlation coefficients in a table format
plot(bio.i)
######################################################################################
######################################################################################
######################### JUST GRID CELLS WITH POINTS ################################
######################################################################################
######################################################################################

cellsWpts2 <- readOGR(dsn="cells_with_points2", layer = "cellsWpts2")
## read in shapefile of all grid cells containing points (see points2grids)

bio.list <- list.files("climate_data/current" ) ## list all of the variables
bio.list <- bio.list[1:19] ## exlude the metadata files after the 19 bioclimate variables

clipped.model.folder.path <- ("clipped4_climate_data" )
dir.create(clipped.model.folder.path) ## creates a new folder to place the clipped data for the model of the iteration

for (i in 1:length(bio.list)){ ## loop through each climate variable within the model of the iteration
  
  file.in <- paste("climate_data/current", bio.list[i], sep="/") 
  ## paste together the full file name to bring it into R
  
  file.out <- paste(clipped.model.folder.path, (paste(bio.list[i], "asc", sep=".")), sep="/")
  ## create filepath for clipped data within the new folder, to export from R
  
  bio.i <- raster(file.in) ## read in file as raster
  
  bio.i <- mask(crop(bio.i, cellsWpts2), cellsWpts2) 
  ## crop raster to extent of united states, and then clip(mask) to the buffer
  
  writeRaster(bio.i, file.out, format= "ascii") ## write out raster as an ESRI raster
  
  print(i) ## keep track of place within model
  
}

##########################################
######### assessing correlations #########
##########################################

bio.list <- list.files("clipped4_climate_data/", full.names =T) 
## lists all 19 climate rasters in the climate data to be used for fitting

bio.stack <- stack(bio.list) ## creates an object with each of the nine rasters stacked
proj4string(bio.stack) <- proj4string(raster("climate_data/current/bio_1")) 
## assigns reference info based the files where these files came from 

bio.stackM <- as.matrix(bio.stack) ## converts the stacked object into a matrix statistical analysis

corr.matrix <- corr.test(bio.stackM, y = NULL, use = "pairwise", method="spearman", adjust="bonferroni", alpha=.05, ci=F)
## matrix of spearman's correlation of the 19 variables, p-value adjusted using Bonferroni correction

write.csv(data.frame(corr.matrix$p), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/p.vals4.csv") 
## export p-values in a table format
write.csv(data.frame(corr.matrix$r), "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/corr.vals4.csv") 
## export correlation coefficients in a table format


