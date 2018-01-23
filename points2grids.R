#### How to relate points to climate grid cells
###  Mitch O'Neill
##   created 8/18/2017
#    last modified: 9/12/2017

setwd("C:/Users/mwone/Documents/geodata/")
library(rgdal)
library(raster)
library(rgeos)

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_12_21_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
head(edd)
## read in EDDMapS database (most recently prepped version)


#sp.list <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/species_summary_12_12_2017.csv", header = T, sep = ",", quote= "\"", 
#                       comment.char= "", stringsAsFactors = F, strip.white = T)
#sp.list <- sp.list$usda.code[sp.list$useable > 0]
#edd <- edd[edd$USDAcode %in% sp.list,]


## extract only the fields that are needed, to save processing time
edd <- data.frame(cbind(1:length(edd$ScientificName), edd$ObjectID, edd$ScientificName, edd$USDAcode, edd$ORD, edd$Latitude_Decimal, edd$Longitude_Decimal, edd$ORDsourceNUM), stringsAsFactors=F)
colnames(edd) <- c("id", "objectid","species", "usda", "abundance", "latitude", "longitude", "source") ## give fields better names
head(edd)

## coords need to be numeric
edd$latitude <- as.numeric(edd$latitude)
edd$longitude <- as.numeric(edd$longitude)
coordinates(edd) <- c(7,6) ## specifies long, lat

## all of EDDMapS is supposed to be nad83/wgs84
proj4string(edd) <- CRS("+init=epsg:4326")
#plot(edd)

## read in fishnet (square polygons corresponding to grid cells)
fishnet <- readOGR(dsn = "fishnetPolygon", layer = "fishnetPolygon")

edd <- spTransform(edd, proj4string(fishnet)) ## transform edd points to same proj4string of fishnet

fishnetD <- data.frame(fishnet) ## make dataframe version to append ID column to
cellID <- 1:length(fishnetD$Id) ## make list of unique IDs for each row
fishnetD <- cbind(fishnetD, cellID) ## assign unique ID # to each row
fishnet <- SpatialPolygonsDataFrame(Sr = fishnet, data = fishnetD) 
## convert dataframe ba/ck into spatial object, using source spatial object

tab <- over(edd, fishnet) ## overlays points and fishnet squares, yielding rows labeled 
## by the point id (row.name), with information from the cell
## that the point falls in
tab$ptID <- as.numeric(row.names(tab)) ## convert rownames to its own field for easier use
tab$Id <- NULL ## remove meaningless field (all zeroes)
head(tab)

write.csv(tab,"overlayTableORD.csv", row.names=F)
#tab <- read.table("overlayTableORD.csv", header = T, sep = ",", quote= "\"", 
#           comment.char= "", stringsAsFactors = F, strip.white = T)

edd <- data.frame(merge(edd, tab, by.x = "id", by.y = "ptID", all=F))
## for each point (row) in the eddmaps database, merge the cellID corresponding to that
## point from table produced from overlaying the eddmaps points with the fishnet cells

edd$optional <- NULL ## remove meaningless column
edd <- edd[!is.na(edd$cellID) & !is.na(edd$usda) & edd$usda != "",]
## exclude rows without a usda code or where the point did not fall in a grid cell
write.csv(edd, "C:/Users/mwone/Documents/EDDMapS data/eddmapsWgrid2.csv", row.names=F)

################## combining 1 and 2 ###################
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmapsWgrid2.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

head(edd)
str(edd)

edd$abundance <- as.numeric(edd$abundance)
edd$source <- as.numeric(edd$source)

summary(as.factor(edd$abundance))
edd$abundance[edd$abundance == 1] <- 2
summary(as.factor(edd$abundance))

#rm(fishnet,fishnetD,tab,cellID) ## garbage cleaning

#### THIN DATA TO GRID CELL (WITHIN SPECIES)
species   <- "TEMPLATE"
max       <- -99
med       <- -99
best      <- -99
no.pts    <- -99
no.best   <- -99
source    <- -99
latitude  <- -99
longitude <- -99
cellID    <- -99
edd2 <- data.frame(species,max,med,best,no.pts,no.best,source,latitude,longitude,cellID, stringsAsFactors = F)

sp.list <- as.list(as.character(unique(edd$usda)))
## list of all usda codes

for (i in 1:length(sp.list)){ ## loop through all species
  spp.i <- edd[edd$usda == sp.list[i],] ## subset data to the species of the iteration
  cell.ids <- as.list(unique(spp.i$cellID)) ## list all the cells that the species has points in
  print(paste("species number", i, sep=" ")) ## keep track of species of the iteration
  for (j in 1:length(cell.ids)){ ## loop through all cells for that species
    cell.j <- spp.i[spp.i$cellID == cell.ids[j],] ## subset to all points in the cell of the iteration
    cell.j <- cell.j[!is.na(cell.j$species),] 
     
    species <- cell.j$usda[1]
    max <- max(cell.j$abundance)
    med <- median(cell.j$abundance)
    best <- median(cell.j$abundance[cell.j$source == min(cell.j$source)])
    no.pts <- length(cell.j$abundance)
    no.best <- length(cell.j$abundance[cell.j$source == min(cell.j$source)])
    source <- min(cell.j$source)
    latitude <- mean(cell.j$latitude)
    longitude <- mean(cell.j$longitude)
    cellID <- unique(cell.j$cellID)
    
  
    keep <- data.frame(species,max,med,best,no.pts,no.best,source,latitude,longitude,cellID, stringsAsFactors = F) 
    edd2 <- rbind(edd2, keep) ## append the kept row to the master object
    #print(j)
  }
print(i)  
}

head(edd2)
edd2 <- edd2[edd2$species != "TEMPLATE",] ## remove template row
sum(edd2$no.pts) == length(edd$id[edd$usda %in% sp.list])

#hist(edd2$med[edd2$no.best>1]/edd2$best[edd2$no.best>1])
#summary(edd2$no.best)
#hist(edd2$no.best)
#edd2[edd2$no.best==2931,]


write.csv(edd2, "C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_12_21_2017.csv", row.names=F)


