#### Bias data sets
###  Mitch O'Neill
##   created 1/31/2018
#    last modified: 2/01/2018

library(rgdal)
library(raster)
setwd("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/")

extentShape = readOGR(dsn = "ArcFiles_2_2_2018/us_shape", layer = "us_shape")
#extentShape2 <- spTransform(extentShape, "+init=epsg:4326")

codes <- read.table("species.codes.1.31.18.csv", header=T, sep=",")

hotspots <- read.table("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/IAS_occurences_final_analysis.csv", 
                         header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
head(hotspots)
hotspots <- hotspots[hotspots$PLANT_CODE %in% codes$x,]

### append new occurencces to the hotspots dataset
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/EDDMaps_plants_10_18_2016.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
## All plant points from EDDMapS, delivered from Rebekah Wallace to Bethany Bradley on October 18, 2016

edd <- edd[edd$negative == 0,]

## makes a column with ScientificName (species name, with subspecies/variety/forma if provided)
edd$ScientificName[edd$SUB_variety == "N/A"] <- paste(edd$SUB_genus[edd$SUB_variety == "N/A"], edd$SUB_species[edd$SUB_variety == "N/A"], sep = " ")
edd$ScientificName[edd$SUB_variety != "N/A"] <- paste(edd$SUB_genus[edd$SUB_variety != "N/A"], edd$SUB_species[edd$SUB_variety != "N/A"], edd$SUB_variety[edd$SUB_variety != "N/A"], sep = " ")

edd <- data.frame(cbind(edd$ScientificName, edd$dateentered,edd$Latitude_Decimal,edd$Longitude_Decimal), stringsAsFactors = F)
colnames(edd) <- c("ScientificName", "date", "lat", "lon")

edd.list <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/speciesList_11_20_2017.csv",
                       header = T, sep = ",", stringsAsFactors = F, strip.white = T)
## download scientific name to USDA code file, from files made manually by Ceara, completed June 1, 2017
## files compiled into one doc by Mitch O'Neill 11/20/2017
## list of species in USDA-sp name file

spp <- as.list(edd.list$eddmap_species)#[!is.na(edd.list$eddmap_species)])
for (sp in spp){
  edd$USDAcode[edd$ScientificName == sp] <- edd.list$USDA_code[edd.list$eddmap_species == sp]
  print(sp)  
}

edd <- edd[!is.na(edd$USDAcode),] 


ssp <- read.table("subspecies2.csv",  header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
## lists subspecies with the species they are to be merged into
head(ssp)
for (i in 1:length(ssp$usda.code)){ ## loop through subspecific codes in my data
  edd$USDAcode[edd$USDAcodeUNMERGED == ssp$usda.code[i]] <- ssp$usda.lump[i]
  ## reassign to lumped species-level code
  print(i)
}

rm(edd.list,ssp,sp,i)
edd <- edd[edd$USDAcode %in% codes$x,]
edd$date <- as.Date(edd$date)
edd <- edd[edd$date > "2012-12-31",]

edd <- data.frame(cbind(edd$USDAcode, edd$lon, edd$lat, rep("EDDMaps2016", NROW(edd))))
colnames(edd) <- colnames(hotspots)

edd$PLANT_CODE <- as.character(edd$PLANT_CODE)
edd$LONGITUDE_DECIMAL <- as.numeric(as.character(edd$LONGITUDE_DECIMAL))
edd$LATITUDE_DECIMAL <- as.numeric(as.character(edd$LATITUDE_DECIMAL))
edd$DATA_SOURCE <- as.character(edd$DATA_SOURCE)

edd$LONGITUDE_DECIMAL[edd$LONGITUDE_DECIMAL > 65 & edd$LONGITUDE_DECIMAL < 130 & edd$Latitude_Decimal >20 & edd$Latitude_Decimal <55] <-
  -1*edd$LONGITUDE_DECIMAL[edd$LONGITUDE_DECIMAL > 65 & edd$LONGITUDE_DECIMAL < 130 & edd$Latitude_Decimal >20 & edd$Latitude_Decimal <55] 

edd <- edd[edd$LONGITUDE_DECIMAL < -65 & edd$LONGITUDE_DECIMAL > -130 & edd$LATITUDE_DECIMAL > 20 & edd$LATITUDE_DECIMAL <55, ]
## removes records I know that I won't use (no coords, coords way outside US, not absence data or valid infested area data)

edd <- edd[!is.na(edd$LATITUDE_DECIMAL) & !is.na(edd$LONGITUDE_DECIMAL),]
coordinates(edd) <- c(2,3) ## specifies long, lat

#all of EDDMapS is supposed to be nad83/wgs84
#EPSG <- make_EPSG() ## creates library of EPSG codes to help assign proj4string
#EPSG[grepl("WGS 84$", EPSG$note) == TRUE, ] ## search for WGS 84 --> 4326
proj4string(edd) <- CRS("+init=epsg:4326")

## clip pts to L48
## load contiguous US shape file
states = readOGR(dsn = "states", layer = "US_states")
## transform states to same proj4string as eddmaps
states <- spTransform(states, "+init=epsg:4326")
edd  <- edd[states, ] ## clip edd_US to L48
edd <- as.data.frame(edd)

edd$PLANT_CODE <- as.character(edd$PLANT_CODE)
edd$DATA_SOURCE <- as.character(edd$DATA_SOURCE)

full_occurences <- rbind(edd,hotspots)
full_occurences <- full_occurences[order(full_occurences$PLANT_CODE),]

coordinates(full_occurences) <- c(2,3)
proj4string(full_occurences) <- CRS("+init=epsg:4326")
full_occurences <- spTransform(full_occurences, proj4string(extentShape))
full_occurences <- full_occurences[extentShape, ] ## 1454711 -- 1452468
full_occurences <- as.data.frame(full_occurences)

### FULL ABUNDANCE BIAS DATASET 
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_1_25_2018.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd <- edd[edd$USDAcode %in% codes$x,]
abun_occurences <- data.frame(cbind(edd$USDAcode,edd$Longitude_Decimal,edd$Latitude_Decimal,rep("EDDMaps2016", NROW(edd))), stringsAsFactors=F)
colnames(abun_occurences) <- colnames(hotspots)
abun_occurences <- abun_occurences[order(abun_occurences$PLANT_CODE),]
abun_occurences$LONGITUDE_DECIMAL <- as.numeric(abun_occurences$LONGITUDE_DECIMAL)
abun_occurences$LATITUDE_DECIMAL <- as.numeric(abun_occurences$LATITUDE_DECIMAL)

coordinates(abun_occurences) <- c(2,3)
proj4string(abun_occurences) <- CRS("+init=epsg:4326")
abun_occurences <- spTransform(abun_occurences, proj4string(extentShape))
abun_occurences <- abun_occurences[extentShape, ]  
abun_occurences <- as.data.frame(abun_occurences)


### ABUNDANCE DATASET FOR MODELLING
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_1_25_2018.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd <- edd[edd$species %in% codes$x,]
abun_modeling <- data.frame(cbind(edd$species,edd$longitude,edd$latitude,rep("EDDMaps2016", NROW(edd))), stringsAsFactors=F)
colnames(abun_modeling) <- colnames(hotspots)
abun_modeling <- abun_modeling[order(abun_modeling$PLANT_CODE),]
abun_modeling$LONGITUDE_DECIMAL <- as.numeric(abun_modeling$LONGITUDE_DECIMAL)
abun_modeling$LATITUDE_DECIMAL <- as.numeric(abun_modeling$LATITUDE_DECIMAL)

coordinates(abun_modeling) <- c(2,3)
proj4string(abun_modeling) <- proj4string(extentShape)
abun_modeling <- abun_modeling[extentShape, ]  
abun_modeling <- as.data.frame(abun_modeling)

write.csv(full_occurences, "MaxEntFiles/full_bias_pts.csv", row.names=F)
write.csv(abun_occurences, "MaxEntFiles/Abun_bias_pts.csv", row.names=F)
write.csv(abun_modeling,  "MaxEntFiles/Abun_model_pts.csv", row.names=F)

rm(full_occurences, abun_occurences, abun_modeling, biasD, codes, full_occurence, hotspots, bias, extentShape, pop, roads, states, spp)


###### OCCURENCE MODELING DATASET ######################
edd <- read.table("MaxEntFiles/full_bias_pts.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

## read in fishnet (square polygons corresponding to grid cells)
fishnet <- readOGR(dsn = "ArcFiles_2_2_2018/fishnet", layer = "fishnet2018")

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

edd$id <- row.names(edd)
edd <- data.frame(merge(edd, tab, by.x = "id", by.y = "ptID", all=F))
## for each point (row) in the eddmaps database, merge the cellID corresponding to that
## point from table produced from overlaying the eddmaps points with the fishnet cells
#write.csv(edd, "full_occurence_wgrids.csv",row.names=F)
#################################################################################
edd <- read.table("full_occurence_wgrids.csv", sep=",", header=T, stringsAsFactors = F)
#colnames(edd) <- c(as.character(edd[1,])

head(edd)
coordinates(edd) <- c(3,4) ## specifies long, lat
proj4string(edd) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

extentShape <- readOGR(dsn = "ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
proj4string(extentShape) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

edd <- edd[extentShape, ]

edd$optional <- NULL ## remove meaningless column
edd <- edd[!is.na(edd$cellID) & !is.na(edd$PLANT_CODE),] ##1452460 to 1452460
## exclude rows without a usda code or where the point did not fall in a grid cell
edd$PLANT_CODE <- as.character(edd$PLANT_CODE)
rm(fishnet,fishnetD,tab,cellID) ## garbage cleaning

edd <- as.data.frame(edd)
head(edd)

#### THIN DATA TO GRID CELL (WITHIN SPECIES)
PLANT_CODE <- "TEMPLATE"
LONGITUDE_DECIMAL <- -99
LATITUDE_DECIMAL <- -99
DATA_SOURCE <- "TEMPLATE"

edd2 <- data.frame(PLANT_CODE, LONGITUDE_DECIMAL, LATITUDE_DECIMAL, DATA_SOURCE, stringsAsFactors = F)

sp.list <- as.list(as.character(unique(edd$PLANT_CODE)))
## list of all usda codes

for (i in 1:length(sp.list)){ ## loop through all species
  spp.i <- edd[edd$PLANT_CODE == sp.list[i],] ## subset data to the species of the iteration
  cell.ids <- as.list(unique(spp.i$cellID)) ## list all the cells that the species has points in
  print(paste("species number", i, sep=" ")) ## keep track of species of the iteration
  for (j in 1:length(cell.ids)){ ## loop through all cells for that species
    cell.j <- spp.i[spp.i$cellID == cell.ids[j],] ## subset to all points in the cell of the iteration
    cell.j <- cell.j[!is.na(cell.j$PLANT_CODE),] 
    
    PLANT_CODE <- cell.j$PLANT_CODE[1]
    LONGITUDE_DECIMAL <- mean(cell.j$LONGITUDE_DECIMAL)
    LATITUDE_DECIMAL <- mean(cell.j$LATITUDE_DECIMAL)
    DATA_SOURCE <- cell.j$DATA_SOURCE[1]
    
    #keep <- data.frame(species,med,no.pts,no.best,source,latitude,longitude,cellID, stringsAsFactors = F) 
    keep <- data.frame(PLANT_CODE, LONGITUDE_DECIMAL, LATITUDE_DECIMAL, DATA_SOURCE, stringsAsFactors = F)
    edd2 <- rbind(edd2, keep) ## append the kept row to the master object
    #print(j)
  }
  print(i)  
}

head(edd2)
full_modeling <- edd2[edd2$PLANT_CODE != "TEMPLATE",] ## remove template row





write.csv(full_modeling,  "MaxEntFiles/full_model_pts.csv", row.names=F)

###############################################################
library(raster)
bias <- stack("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_pop/us_pop/w001001.adf",
              "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_roads/us_roads/w001001.adf")
plot(bias)
proj4string(bias)
plot(bias$w001001.1) ## pop
plot(bias$w001001.2) ## roads

biasD <- cbind(as.data.frame(bias$w001001.1), as.data.frame(bias$w001001.2))
colnames(biasD) <- c("us_pop", "us_roads")
summary(biasD$us_roads)
summary(biasD$us_pop)
biasD$us_pop[is.na(biasD$us_pop) & !is.na(biasD$us_roads)] <- 0
pop <- raster(nrows=nrow(bias), ncols=ncol(bias), ext=extent(bias), crs=crs(bias), vals=biasD$us_pop )

roads <- raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_roads/us_roads/w001001.adf")

nrow(pop)==nrow(roads)
ncol(pop)==ncol(roads)  
extent(pop)==extent(roads) 
proj4string(pop)==proj4string(roads) 

writeRaster(pop,"C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/MaxEntFiles/us_pop_2_5_2018.asc", format="ascii", prj=T, overwrite=T)
writeRaster(roads,"C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/MaxEntFiles/us_roads_2_5_2018.asc", format="ascii", prj=T, overwrite=T)

#######################################################################
#######################################################################
#######################################################################

tabulate <- read.table("Tab_area/tabulate_area.txt",  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
head(tabulate)

tabulate <- tabulate[,3:11]
head(tabulate)
tabulate$total <- rowSums(tabulate, na.rm = FALSE, dims = 1)
head(tabulate)


roads <- raster("MaxEntFiles/us_roads_2_5_2018.asc")
proj4string(roads) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
#"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
#proj4string(roads) <-
roaDs <- as.data.frame(roads)
head(roaDs)

cellid <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=crs(roads), vals=seq(1:ncell(roads)))
plot(cellid)
writeRaster(cellid, "cellIDraster.asc", format="ascii", prj=T)

fishnet <- readOGR(dsn = "ArcFiles_2_2_2018/fishnet", layer = "fishnet2018")
fidraster <- raster("ArcFiles_2_2_2018/FIDraster_2_8/rasterFID.tif")
plot(fidraster)
fids <- as.data.frame(fidraster)
head(fids)
summary(fids)


length(roaDs$us_roads_2_5_2018[!is.na(roaDs$us_roads_2_5_2018)])
roaDs$total[!is.na(roaDs$us_roads_2_5_2018)] <- tabulate$total

totals <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=crs(roads), vals=roaDs$total)
plot(totals)
dev.off()
sqrt(min(tabulate$total))/1000
sqrt(mean(tabulate$total))/1000
sqrt(max(tabulate$total))/1000
hist(tabulate$total) ## 3.7 - 4.4 (4.1)

head(tabulate)
tabulate$prop1 <- tabulate$BIN_1/tabulate$total
tabulate$prop2 <- tabulate$BIN_2/tabulate$total
tabulate$prop3 <- tabulate$BIN_3/tabulate$total
tabulate$prop4 <- tabulate$BIN_4/tabulate$total
tabulate$prop5 <- tabulate$BIN_5/tabulate$total
tabulate$prop6 <- tabulate$BIN_6/tabulate$total
tabulate$prop7 <- tabulate$BIN_7/tabulate$total
tabulate$prop8 <- tabulate$BIN_8/tabulate$total
tabulate$prop9 <- tabulate$BIN_9/tabulate$total

head(fids)
fids$cellID<- as.numeric(row.names(fids))
fids <- fids[!is.na(fids$rasterFID),]
colnames(fids) <- c("FID","cellID")

head(tabulate)
tabulate <- tabulate[,10:19]
tabulate$FID <- as.numeric(row.names(tabulate)) -1 #tabulate$row - 1
head(tabulate)

summary(fids$FID)
summary(tabulate$FID)
nlcd2fid <- merge(fids,tabulate,by="FID")

head(nlcd2fid)
nlcd2fid <- nlcd2fid[order(nlcd2fid$cellID),]
head(nlcd2fid)
plot(nlcd2fid$total~nlcd2fid$cellID)
#newdata <- mtcars[order(mpg),]

roaDs$nlcd1[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop1
roaDs$nlcd2[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop2
roaDs$nlcd3[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop3
roaDs$nlcd4[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop4
roaDs$nlcd5[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop5
roaDs$nlcd6[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop6
roaDs$nlcd7[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop7
roaDs$nlcd8[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop8
roaDs$nlcd9[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$prop9
roaDs$total[!is.na(roaDs$us_roads_2_5_2018)] <- nlcd2fid$total


nlcd1 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd1)
nlcd2 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd2)
nlcd3 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd3)
nlcd4 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd4)
nlcd5 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd5)
nlcd6 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd6)
nlcd7 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd7)
nlcd8 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd8)
nlcd9 <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$nlcd9)
total <- raster(nrows=nrow(roads), ncols=ncol(roads), ext=extent(roads), crs=proj4string(roads), vals=roaDs$total)

plot(total)
#nlcd <- stack(nlcd1,nlcd2,nlcd3,nlcd4,nlcd5,nlcd6,nlcd7,nlcd8,nlcd9)
plot(nlcd1)
plot(nlcd2)
plot(nlcd3)
plot(nlcd4)
plot(nlcd5)
plot(nlcd6)
plot(nlcd7)
plot(nlcd8)
plot(nlcd9)

#dir.create("nlcd")
writeRaster(nlcd1, "nlcd/nlcd1.asc", format="ascii", prj=T, overwrite=T) ## unsuitable
writeRaster(nlcd2, "nlcd/nlcd2.asc", format="ascii", prj=T, overwrite=T) ## develop
writeRaster(nlcd3, "nlcd/nlcd3.asc", format="ascii", prj=T, overwrite=T) ## forest
writeRaster(nlcd4, "nlcd/nlcd4.asc", format="ascii", prj=T, overwrite=T) ## evergreen
writeRaster(nlcd5, "nlcd/nlcd5.asc", format="ascii", prj=T, overwrite=T) ## shrub
writeRaster(nlcd6, "nlcd/nlcd6.asc", format="ascii", prj=T, overwrite=T) ## herbaceous
writeRaster(nlcd7, "nlcd/nlcd7.asc", format="ascii", prj=T, overwrite=T) ## pasture/hay
writeRaster(nlcd8, "nlcd/nlcd8.asc", format="ascii", prj=T, overwrite=T) ## cropland
writeRaster(nlcd9, "nlcd/nlcd9.asc", format="ascii", prj=T, overwrite=T) ## wetland

## read in Worldclim climate rasters clipped to L48; 6 variables manually selected by assessing correlations
bio <- stack("C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_2.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_5.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_6.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_8.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_12.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_15.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_1.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_3.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_4.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_7.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_9.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_10.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_11.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_13.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_14.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_16.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_17.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_18.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_19.asc")
proj4string(bio) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

nlcd <- stack("MaxEntFiles/us_pop_2_5_2018.asc",
              "MaxEntFiles/us_roads_2_5_2018.asc")#,#)#,
#              "nlcd/nlcd1.asc",
#              "nlcd/nlcd2.asc",
#              "nlcd/nlcd3.asc",
#              "nlcd/nlcd4.asc",
#              "nlcd/nlcd5.asc",
#              "nlcd/nlcd6.asc",
#              "nlcd/nlcd7.asc",
#              "nlcd/nlcd8.asc",
#              "nlcd/nlcd9.asc")

#test <- raster("nlcd/nlcd9.asc")

proj4string(nlcd) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

dim(bio)
dim(nlcd)


#extentShape = readOGR(dsn = "ArcFiles_2_2_2018/us_shape", layer = "us_shape")
#plot(extentShape)
#proj4string(bio) <- proj4string(extentShape)
#bio <- mask(crop(bio, extentShape), extentShape)

bio <- mask(crop(bio, nlcd$us_pop_2_5_2018),nlcd$us_pop_2_5_2018)
nlcd <- mask(crop(nlcd, bio$bio_1),bio$bio_1)
#bio <- mask(crop(bio, nlcd),nlcd)
                                                                                                                                                                                                                            dim(bio)
dim(nlcd)

nlcD <- as.data.frame(nlcd)
bioD <- as.data.frame(bio)

nlcor<-cor(nlcD, use="pairwise", method="spearman")
write.csv(nlcor,"nlcor.csv",row.names=F)

ncell(bio)==ncell(nlcd)
nrow(bio)==nrow(nlcd)
ncol(bio)==ncol(nlcd)  
extent(bio)==extent(nlcd) 
proj4string(bio)==proj4string(nlcd) 

ymin(bio) == ymin(nlcd) ## FALSE
ymax(bio) == ymax(nlcd) ## TRUE
xmin(bio) == xmin(nlcd) ## TRUE
xmax(bio) == xmax(nlcd) ## FALSE

#bioD <- as.data.frame(bio)
#stackD <- as.data.frame(stack)
#head(stackD)
#head(bioD)

full_matrix <- cbind(nlcD,bioD) 
head(full_matrix)

length(full_matrix$nlcd1[is.na(full_matrix$nlcd1) & !is.na(full_matrix$bio_1)])
length(full_matrix$nlcd1[!is.na(full_matrix$nlcd1) & is.na(full_matrix$bio_1)])
#writeRaster(stack$layer.1,"checkNLCD.asc", format="ascii",prj=T)
#writeRaster(bio$bio_2 ,"checkBIO.asc", format="ascii",prj=T)

cor.matrix <- cor(full_matrix, use="pairwise", method="spearman")
write.csv(cor.matrix,"nlcd_bio.csv",row.names=F)

plot(nlcd$nlcd1, col="blue")
plot(bio$bio_2, col="orange", add=T)



dir.create("C:/Users/mwone/Documents/Environmental_Data_2_8_2018")
writeRaster(nlcd$nlcd3, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_3.asc", format="ascii")
writeRaster(nlcd$nlcd4, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_4.asc", format="ascii")
writeRaster(nlcd$nlcd5, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_5.asc", format="ascii")
writeRaster(nlcd$nlcd6, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_6.asc", format="ascii")
writeRaster(nlcd$nlcd7, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_7.asc", format="ascii")
writeRaster(nlcd$nlcd8, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_8.asc", format="ascii")
writeRaster(nlcd$nlcd9, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_9.asc", format="ascii")
writeRaster(bio$bio_2, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_2.asc", format="ascii")
writeRaster(bio$bio_5, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_5.asc", format="ascii")
writeRaster(bio$bio_6, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_6.asc", format="ascii")
writeRaster(bio$bio_8, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_8.asc", format="ascii")
writeRaster(bio$bio_12, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_12.asc", format="ascii")
writeRaster(bio$bio_15, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_15.asc", format="ascii")
writeRaster(nlcd$us_roads_2_5_2018, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/us_roads.asc", format="ascii")
writeRaster(nlcd$us_pop_2_5_2018, "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/us_pop.asc", format="ascii")


envi <- stack( "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_8.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/nlcd_9.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_2.asc", 
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_5.asc", 
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_6.asc", 
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_8.asc", 
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_12.asc",
               "C:/Users/mwone/Documents/Environmental_Data_2_8_2018/bio_15.asc" )


abun <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/MaxEntFiles/Abun_model_pts.csv",
                   sep=",", stringsAsFactors = F, header=T)
head(abun)
coordinates(abun)<-c(2,3)
proj4string(abun)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

ext <- extract(envi$nlcd_3, abun)
summary(ext)## no NAs!


######## FINALIZE MAXENT DATASETS
full_modeling   <- read.table("MaxEntFiles/full_model_pts.csv", header=T, sep=",", stringsAsFactors =F)
full_occurences <- read.table("MaxEntFiles/full_bias_pts.csv", header=T, sep=",", stringsAsFactors =F)
abun_occurences <- read.table("MaxEntFiles/Abun_bias_pts.csv", header=T, sep=",", stringsAsFactors =F)
abun_modeling   <- read.table("MaxEntFiles/Abun_model_pts.csv", header=T, sep=",", stringsAsFactors =F)

summary(full_modeling)
summary(full_occurences)
summary(abun_occurences)
summary(abun_modeling)

full_modeling$DATA_SOURCE <- NULL
full_occurences$DATA_SOURCE <- NULL
abun_occurences$DATA_SOURCE <- NULL
abun_modeling$DATA_SOURCE <- NULL


full_occurences$PLANT_CODE <- "BIAS"
abun_occurences$PLANT_CODE <- "BIAS"


coordinates(full_modeling) <- c(2,3)
coordinates(full_occurences) <- c(2,3)
coordinates(abun_occurences) <- c(2,3)
coordinates(abun_modeling) <- c(2,3)

proj4string(full_modeling) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(full_occurences) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(abun_occurences) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(abun_modeling) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

extentShape = readOGR(dsn = "ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
plot(extentShape)
proj4string(extentShape) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

full_modeling <- full_modeling[extentShape, ]     ##306206  to 305949
full_occurences <- full_occurences[extentShape, ] ##1452468 to 1452182
abun_occurences <- abun_occurences[extentShape, ] ##521954  to 521945
abun_modeling <- abun_modeling[extentShape, ]     ##99586   to 99586 


write.csv(full_modeling, "MaxEntFiles/thinned_pts_2_9_2018.csv")
write.csv(full_occurences, "MaxEntFiles/full_pts_2_9_2018.csv", row.names = F)
write.csv(abun_occurences, "MaxEntFiles/full_abundance_pts_2_9_2018.csv", row.names = F)
write.csv(abun_modeling,"MaxEntFiles/thinned_abundance_pts_2_9_2018.csv", row.names = F)

##Remove "optional" column
full_occurences <- read.table("MaxEntFiles/full_pts_2_9_2018.csv", 
                       header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
abun_occurences <- read.table("MaxEntFiles/full_abundance_pts_2_9_2018.csv", 
                              header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
full_occurences$optional <- NULL
abun_occurences$optional <- NULL
write.csv(full_occurences, "MaxEntFiles/full_pts_2_15_2018.csv", row.names = F)
write.csv(abun_occurences, "MaxEntFiles/full_abundance_pts_2_15_2018.csv", row.names = F)




########################################
### separating into individual files
########################################

full_modeling <- read.table("MaxEntFiles/thinned_pts_2_9_2018.csv", header = T, sep = ",", 
                           quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

abun_modeling <- read.table("MaxEntFiles/thinned_abundance_pts_2_9_2018.csv", header = T, sep = ",", 
                           quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

full_modeling$X <- NULL
abun_modeling$X <- NULL
full_modeling$optional <- NULL
abun_modeling$optional <- NULL

sp.list <- unique(full_modeling$PLANT_CODE)
#dir.create("MaxEntFiles/abun")
#dir.create("MaxEntFiles/full")

for(i in 1:length(sp.list)){
  spec.abun <- abun_modeling[abun_modeling$PLANT_CODE == sp.list[i],]
  filename.abun <- paste("MaxEntFiles/abun", paste(sp.list[i], "csv", sep="."), sep="/")
  write.csv(spec.abun, filename.abun, row.names=F)
  

  spec.full <- full_modeling[full_modeling$PLANT_CODE == sp.list[i],]
  filename.full <- paste("MaxEntFiles/full", paste(sp.list[i], "csv", sep="."), sep="/")
  write.csv(spec.full, filename.full, row.names=F)
  
  print(i)
}


#sp.list <- sp.list[order(sp.list)]
#head(sp.list)
#write.csv(sp.list, "MaxEntFiles/speciesList.csv", row.names=F)

#### make list of strings to add to each species' model specification (variables to ignore)
library(glue)
sp.list <- read.table("ordsums2_13_2018.csv", header = T, sep = ",", 
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

sp.list <- data.frame(sp.list$species.code, sp.list$formu2, stringsAsFactors = F)
colnames(sp.list) <- c("code","model")
str(sp.list)
head(sp.list)

var.list <- c("bio_2", "bio_5", "bio_6", "bio_8", "bio_12","bio_15",
              "nlcd_3","nlcd_4","nlcd_5","nlcd_6","nlcd_7","nlcd_8")

              
for (i in 1:NROW(ordsums)){
  ignore <- c("")
  for(j in 1:length(var.list)){
   if(!grepl(var.list[j], sp.list$model[i], fixed=T)){
     ignore <- cbind(ignore,var.list[j])
     }
  }
  sp.list$ignore[i] <- collapse(ignore, sep= " -N ")
  print(i)
}

head(sp.list)
write.csv(sp.list, "MaxEntFiles/speciesList.csv", row.names=F)


###########################################
library(raster)

i <- 0
miny <- -9
maxy <- -9
minx <- -9
maxx <- -9
rows <- -9
cols <- -9
cels <- -9
rasterSum <- data.frame(i,miny,maxy,minx,maxx,rows,cols,cels)

rasterList <- c(list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/", full.names = T),
                "C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc",
                "C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/BIAS_avgABUN.asc")

for (i in 8:15){  ##1-7, 8-15
  lyr <- raster(rasterList[i])
  miny <- ymin(lyr)
  maxy <- ymax(lyr)
  minx <- xmin(lyr)
  maxx <- xmax(lyr)
  rows <- nrow(lyr)
  cols <- ncol(lyr)
  cels <- ncell(lyr)
  rasterSum.i <- data.frame(i,miny,maxy,minx,maxx,rows,cols,cels)
  rasterSum <- rbind(rasterSum,rasterSum.i)
  print(i)
}

rasterSum <- rasterSum[rasterSum$i > 0,]
rasterSum
summary(rasterSum)
max(rasterSum$miny) - min(rasterSum$miny) == 0
max(rasterSum$maxy) - min(rasterSum$maxy) == 0
max(rasterSum$minx) - min(rasterSum$minx) == 0
max(rasterSum$maxx) - min(rasterSum$maxx) == 0 #### difference???
max(rasterSum$cols) - min(rasterSum$cols) == 0
max(rasterSum$rows) - min(rasterSum$rows) == 0
max(rasterSum$cels) - min(rasterSum$cels) == 0

rasterSum$maxxDIFF <- (rasterSum$maxx-mean(rasterSum$maxx))*1000000000
unique(rasterSum$maxxDIFF)
#rasterList
rasterSum$maxx[1]-rasterSum$maxx[12]
rasterSum$maxx[11]-rasterSum$maxx[13]
## Diff of 0.000000001 decimal degrees

# wrong.ext <- stack("C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_3.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_4.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_5.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_6.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_7.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_8.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc", ## renamed in files then changed back
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/BIAS_avgABUN.asc")
# 
# right.ext <- stack("C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/bio_2.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/bio_5.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/bio_6.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/bio_8.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/bio_12.asc",
#                    "C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/bio_15.asc")
# 
# wrong.ext2right <- resample(wrong.ext, right.ext, method="ngb") ##nearest neighbor
# 
writeRaster(w2r$nlcd_3,"C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_3.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$nlcd_4,"C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_4.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$nlcd_5,"C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_5.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$nlcd_6,"C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_6.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$nlcd_7,"C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_7.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$nlcd_8,"C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs/nlcd_8.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$BIAS_avgFULL,"C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc", format="ascii", prj=F, overwrite=T)
writeRaster(w2r$BIAS_avgABUN,"C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/BIAS_avgABUN.asc", format="ascii", prj=F, overwrite=T)
