#### How to relate points to climate grid cells
###  Mitch O'Neill
##   created 8/18/2017
#    last modified: 9/12/2017

setwd("C:/Users/mwone/Documents/geodata/")
library(rgdal)
library(raster)
library(rgeos)

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_08_22_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
## read in EDDMapS database (most recently prepped version)

edd <- edd[edd$negative == 0 & edd$validInfestedArea == 1 & edd$ReporterFULLName != "US Army Corps of Engineers Ombil Database", ] 
## only include infested area data, exclude absences, exclude values over 10 acres

edd <- edd[edd$grossAreaInAcres <= 1000, ]

####################
sp.list <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList10_11_2017.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)
#head(sp.list)
sp.list <- sp.list$usda.code[sp.list$potential.use2 == 1]

edd <- edd[edd$USDAcode %in% sp.list,]

## extract only the fields that are needed, to save processing time
edd <- data.frame(cbind(1:length(edd$ScientificName), edd$ScientificName, edd$USDAcode, edd$infestedAreaInAcres, edd$negative, edd$Latitude_N, edd$Longitude_N, stringsAsFactors=F))
colnames(edd) <- c("id", "species", "usda", "abundance", "absence", "latitude", "longitude") ## give fields better names

## de-factor the numeric columns that will be used
edd$latitude <- as.numeric(as.character(edd$latitude))
edd$longitude <- as.numeric(as.character(edd$longitude))
edd$id <- as.numeric(as.character(edd$id))
edd$usda <- as.character(edd$usda) ## de-factor the usda code
edd$abundance <- as.numeric(as.character(edd$abundance)) ## de-factor the abundances
edd$abundance[edd$absence == 1] <- 0 ## assign 0 abundance to all absences

coordinates(edd) <- c(7,6) ## specifies long, lat

## all of EDDMapS is supposed to be nad83/wgs84
proj4string(edd) <- CRS("+init=epsg:4326")
#plot(edd)



## read in fishnet (square polygons corresponding to grid cells)
fishnet <- readOGR(dsn = "fishnetPolygon", layer = "fishnetPolygon")
## plot(fishnet) very time consuming

edd <- spTransform(edd, proj4string(fishnet)) ## transform edd points to same proj4string of fishnet

fishnetD <- data.frame(fishnet) ## make dataframe version to append ID column to
cellID <- 1:length(fishnetD$Id) ## make list of unique IDs for each row
fishnetD <- cbind(fishnetD, cellID) ## assign unique ID # to each row
fishnet <- SpatialPolygonsDataFrame(Sr = fishnet, data = fishnetD) 
## convert dataframe back into spatial object, using source spatial object

tab <- over(edd, fishnet) ## overlays points and fishnet squares, yielding rows labeled 
## by the point id (row.name), with information from the cell
## that the point falls in
tab$ptID <- as.numeric(row.names(tab)) ## convert rownames to its own field for easier use
tab$Id <- NULL ## remove meaningless field (all zeroes)
head(tab)
write.csv(tab,"overlayTable1000.csv", row.names=F)


#    ############ subset fishnet to all cells containing points ############
#    ptsPerCell <- data.frame(table(tab$cellID)) ## table with number of points per grid cell polygon
#    ptsPerCell$cellID <- as.numeric(as.character(ptsPerCell$Var1)) # de-factor the var1 and rename it cellID
#    
#    ids <- as.list(ptsPerCell$cellID[ptsPerCell$Freq > 0]) ## make list of cellIDs where there are any points in the cell
#    cellsWpts <- fishnet[fishnet$cellID %in% ids,] ## subset all cells with points in them (using IDs)
#    #plot(cellsWpts)
#    
#    writeOGR(cellsWpts, dsn= "cells_with_points", layer = "cellsWpts", driver = "ESRI Shapefile")
#    ## write out file for later use

#### specify the cell that each point falls in
edd <- data.frame(merge(edd, tab, by.x = "id", by.y = "ptID", all=F))
## for each point (row) in the eddmaps database, merge the cellID corresponding to that
## point from table produced from overlaying the eddmaps points with the fishnet cells
# write.csv(edd, "C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_10_09_2017.csv", row.names=F)
## write out eddmaps database with cell IDs for future use
head(edd)

############### to skip the above code, run the following code: ###############
# edd <- read.table("C:/Users/mwone/Documents/EddmapS data/eddmaps_prepped_08_28_2017.csv", header = T, sep = ",", quote= "\"", 
#                  comment.char= "", stringsAsFactors = F, strip.white = T)
# cellsWpts <- readOGR(dsn="cells_with_points", layer = "cellsWpts")
# tab <- read.table("overlayTable.csv", header=T, sep=",", stringsAsFactors = F)
edd$optional <- NULL ## remove meaningless column
edd$NA. <- NULL
###############################################################################



edd <- edd[!is.na(edd$cellID) & !is.na(edd$usda) & edd$usda != "",]
## exclude rows without a usda code or where the point did not fall in a grid cell
#edd <- edd[edd$absence != 1,] ## remove absence data


rm(ptsPerCell,fishnet,fishnetD,tab,cellID) ## garbage cleaning


#### THIN DATA TO GRID CELL (WITHIN SPECIES)

head(edd) ## create template dataframe matching eddmapps data
id        <- -99
species   <- "TEMPLATE"
usda      <- "TEMPLATE"
abundance <- -99
absence   <- -99
cellID    <- -99
longitude <- -99
latitude  <- -99
edd2 <- data.frame(id,species,usda,abundance,absence,cellID,longitude,latitude)

sp.list <- as.list(as.character(unique(edd$usda)))
## list of all usda codes

for (i in 1:length(sp.list)){ ## loop through all species
  spp.i <- edd[edd$usda == sp.list[i],] ## subset data to the species of the iteration
  cell.ids <- as.list(unique(spp.i$cellID)) ## list all the cells that the species has points in
  print(paste("species number", i, sep=" ")) ## keep track of species of the iteration
  for (j in 1:length(cell.ids)){ ## loop through all cells for that species
    cell.j <- spp.i[spp.i$cellID == cell.ids[j],] ## subset to all points in the cell of the iteration
    keep <- cell.j[cell.j$abundance == max(cell.j$abundance),] ## keep the max value
    keep <- keep[1,] ## only want one record!
    edd2 <- rbind(edd2, keep) ## append the kept row to the master object
  }
#print(i)  
}

edd2 <- edd2[edd2$usda != "TEMPLATE",] ## remove template row
write.csv(edd2, "C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_1000.csv", row.names=F)
## write out for later use

head(edd2) ## where id and row.id differ
