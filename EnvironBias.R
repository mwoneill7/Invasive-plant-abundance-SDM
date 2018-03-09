#### Environmental bias; MESS; etc.
###  Mitch O'Neill
##   created 9/21/2017
#    last modified 9/25/2017

setwd("C:/Users/mwone/Documents/geodata/")
library(raster)

## read in thinned dataset
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_10acre.csv",  
                  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

## read in species list
spp <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList10_11_2017.csv",  
                  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
sp.list <- spp$usda.code[spp$potential.use_10ac == 1] ## list of useable species

edd <- edd[edd$usda %in% sp.list,] ## subset dataset to useable species
coordinates(edd) <- c(7,8) ## set lat/lon
proj4string(edd) <- CRS("+init=epsg:4326") ## assume wgs84 for all eddmaps data

## read in selected climate variables for L48
bio.stack <- stack("clipped_climate_data/current/bio_2.asc",
                   "clipped_climate_data/current/bio_5.asc",
                   "clipped_climate_data/current/bio_6.asc",
                   "clipped_climate_data/current/bio_8.asc",
                   "clipped_climate_data/current/bio_12.asc",
                   "clipped_climate_data/current/bio_15.asc")

proj4string(bio.stack) <- proj4string(raster("climate_data/current/bio_1"))
## assign proj4string of clipped data to the global dataset they came from

edd <- spTransform(edd, proj4string(bio.stack))
## transorm point data to climate proj4string

edd <- cbind(edd@data, data.frame(extract(bio.stack, edd)))
## extract climate variable values to the points, and bind these to the
## point data records
head(edd)

#     cell.ids <- unique(edd$cellID) ## list of unique cell IDs in the point data

#     edd2 <- edd[1,] ## first row of point dataset, to initiate loop
#     
#     for (j in 2:length(cell.ids)){ ## loop through all cells for that species
#       cell.j <- edd[edd$cellID == cell.ids[j],][1,] ## keep 1 record from cell
#       edd2 <- rbind(edd2, cell.j) ## append the kept row to the master object
#       print(j)
#     }
#     
#     length(edd2$id) == length(unique(edd2$cellID)) ## should be true

## full names of bioclimatic variables, for later titles
letters <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
full.names <- c("mean diurnal range (°C)",
               "max temp of warmest mo. (°C)",
               "min temp of coldest mo. (°C)",
               "mean temp of wettest quar. (°C)",
               "annual precipitation (mm)",
               "precip. seasonality (mm)")

bio.stackD <- as.data.frame(bio.stack)
bio.stackD[,1:4] <- bio.stackD[,1:4]/10
edd[,7:10] <- edd[,7:10]/10


## make dataframe version to make density plots
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.5, cex.axis=1.2, lwd=2, mar=c(5, 4, 4, 2) + 0.1) 
for (i in 1:length(colnames(bio.stackD))){ ## loop through each bioclimatic variable
  d.i <- density(bio.stackD[,i], na.rm=T)  ## density plot for L48
  d.i2 <- density(edd[,(i+6)], na.rm=T)    ## density plot for point data
  #d.i3 <- density(edd2[,(i+6)], na.rm=T)   ## density plot for all grid cells w point data
  
  if(i==1 | i==2 | i==3 | i==6){
    plot(d.i2, col="red", main=letters[i], xlab=full.names[i], sub="")
    #lines(d.i3, col="brown")
    lines(d.i)
  } else {
    plot(d.i, col="black", main=letters[i], xlab=full.names[i], sub="")
    #lines(d.i3, col="brown")
    lines(d.i2, col="red")
  }
  print(i)

}


###################################################################################
###################################################################################
############################### MESS ANALYSIS #####################################
###################################################################################
###################################################################################

library(dismo) ## needed for mess() function

## bio.ref <- stack("clipped4_climate_data/bio_2.asc",
##                  "clipped4_climate_data/bio_5.asc",
##                  "clipped4_climate_data/bio_6.asc",
##                  "clipped4_climate_data/bio_8.asc",
##                  "clipped4_climate_data/bio_12.asc",
##                  "clipped4_climate_data/bio_15.asc") ### stack fitting region
## 
## reference <- as.matrix(bio.ref)  ### make matrix version of stacked reference rasters 

bio.ref <- stack("clipped_climate_data/current/bio_2.asc",
                 "clipped_climate_data/current/bio_5.asc",
                 "clipped_climate_data/current/bio_6.asc",
                 "clipped_climate_data/current/bio_8.asc",
                 "clipped_climate_data/current/bio_12.asc",
                 "clipped_climate_data/current/bio_15.asc") ### stack fitting region

reference <- as.matrix(bio.ref)  ### make matrix version of stacked reference rasters 
prediction <- list.files("clipped_climate_data/")
#prediction <- prediction[prediction != "current"]


dir.create("mess2/")
#setwd("C:/Users/mwone/Documents/geodata/")
for (i in 1:length(prediction)){

  setwd(paste("C:/Users/mwone/Documents/geodata/clipped_climate_data", prediction[i], sep="/"))
  variables <-list.files(getwd())
  bio_2 <- raster(variables[12])
  bio_5 <- raster(variables[15])
  bio_6 <- raster(variables[16])
  bio_8 <- raster(variables[18])  
  bio_12 <- raster(variables[4])
  bio_15 <- raster(variables[7])
  bio.pred <- stack(bio_2, bio_5, bio_6, bio_8, bio_12, bio_15) ### stack prediction region rasters

  mess <- mess(bio.pred, reference, full=FALSE) 
  ### MESS analysis of similarity of prediction region to the reference values

  mess.file <- paste("C:/Users/mwone/Documents/geodata/mess2", paste(prediction[i], "asc", sep=".") , sep="/") 

  writeRaster(mess, mess.file, format= "ascii") ### write out raster for analysis  
  
  print(i)
}

mess <- stack(list.files("mess2/", full.names=T)[1:11])
messD <- as.data.frame(mess)
head(messD)


messD <- messD[messD$bc45 != Inf,]

bc <- length(messD$bc85[messD$bc85 < 0])/length(messD$bc85)
cc <- length(messD$bc85[messD$cc85 < 0])/length(messD$bc85)
gs <- length(messD$bc85[messD$gs85 < 0])/length(messD$bc85)
he <- length(messD$bc85[messD$he85 < 0])/length(messD$bc85)
ni <- length(messD$bc85[messD$in85 < 0])/length(messD$bc85)
cu <- length(messD$current[messD$current < 0])/length(messD$bc85)


mean(c(bc, cc, gs, he, ni))
max(c(bc, cc, gs, he, ni))
min((c(bc, cc, gs, he, ni)))

summary(messD)
head(mess$bc85)



################# visualizing MESS ##################
setwd("C:/Users/mwone/Documents/geodata/")
mess <- stack(list.files("mess2/", full.names=T)[1:11])
proj4string(mess) <- proj4string(raster("climate_data/current/bio_1"))


library(rgdal)
library(rgeos)
states <- readOGR(dsn="states", layer="US_states")
states <- gUnaryUnion(states)
plot(states)
states <- spTransform(states, proj4string(mess))


mess.copy <- mess
mess$bc45[mess$bc45 > 0] <- NA
mess$bc85[mess$bc85 > 0] <- NA
mess$cc45[mess$cc45 > 0] <- NA
mess$cc85[mess$cc85 > 0] <- NA
mess$gs45[mess$gs45 > 0] <- NA
mess$gs85[mess$gs85 > 0] <- NA
mess$he45[mess$he45 > 0] <- NA
mess$he85[mess$he85 > 0] <- NA
mess$in45[mess$in45 > 0] <- NA
mess$in85[mess$in85 > 0] <- NA
mess$current[mess$current > 0] <- NA

# plot(states, col="black") 
# plot(mess$bc45, add=T)
# 
# plot(states, col="black") 
# plot(mess$cc45, add=T)
# 
# plot(states, col="black") 
# plot(mess$gs45, add=T)
# 
# plot(states, col="black") 
# plot(mess$he45, add=T)
# 
# plot(states, col="black") 
# plot(mess$in45, add=T)

par(mfrow = c(2,3))

plot(states, col="black", main = "bc85") 
plot(mess$bc85, add=T) 

plot(states, col="black", main = "cc85") 
plot(mess$cc85, add=T)

plot(states, col="black", main = "gs85") 
plot(mess$gs85, add=T)

plot(states, col="black", main = "he85")
plot(mess$he85, add=T)

plot(states, col="black", main = "in85") 
plot(mess$in85, add=T)



bio.pts <- stack("clipped4_climate_data/bio_2.asc",
                 "clipped4_climate_data/bio_5.asc",
                 "clipped4_climate_data/bio_6.asc",
                 "clipped4_climate_data/bio_8.asc",
                 "clipped4_climate_data/bio_12.asc",
                 "clipped4_climate_data/bio_15.asc") ### grid cells w/pts

plot(bio.ref$bio_2) 
ref2 <- as.matrix(bio.pts)  ### make matrix version of stacked reference rasters
mess2 <- mess(bio.ref, ref2, full=F)
proj4string(mess2) <- proj4string(mess)
mess2[mess2$mess > 0] <- NA

plot(states, col="black", main= "current") 
plot(mess2, add=T)


####################################################################
setwd("C:/Users/mwone/Documents/geodata/")
states <- readOGR(dsn="states", layer="US_states")
ecoregions <- readOGR(dsn="ecoregions", layer = "NA_CEC_Eco_Level1")
## read in CEC ecoregions shapefile

states <- spTransform(states, proj4string(ecoregions)) ## transform states to ecoregions CRS

L48ecoregions <- ecoregions[states,] ## subset polygons which coincide with US
dev.off()
plot(L48ecoregions)

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_09_12_2017.csv",  
                  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

## read in species list
spp <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList08_28_2017.csv",  
                  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
sp.list <- spp$usda.code[spp$potential.use2 == 1] ## list of useable species

edd <- edd[edd$usda %in% sp.list,] ## subset dataset to useable species
coordinates(edd) <- c(7,8) ## set lat/lon
proj4string(edd) <- CRS("+init=epsg:4326") ## assume wgs84 for all eddmaps data

edd <- spTransform(edd, proj4string(L48ecoregions))

table <- over(edd, L48ecoregions)
table(table$NA_L1NAME)
head(table)

bio1 <- raster("climate_data/current/bio_1")
states <- spTransform(states, proj4string(bio1))
bio1 <- mask(crop(bio1,states),states)

L48ecoregion$NA_L1NAME <- as.character(L48ecoregion$NA_L1NAME)
eco <- unique(as.character(L48ecoregions@data$NA_L1NAME))
L48ecoregions <- spTransform(L48ecoregions, proj4string(bio1))

for (i in 2:length(eco)){
  region.i <- L48ecoregions[L48ecoregions$NA_L1NAME == eco[i],] 
  bio.r <- mask(crop(bio1,region.i),region.i)
  bioD <- as.data.frame(bio.r$bio_1)
  print(eco[i])
  print(length(bioD$COUNT[!is.na(bioD$COUNT)]))#gridcells)
  print(i)
  
}





##################################
#### Mess analysis from abun pts to full pts
library(rgdal)
library(raster)
install.packages('dismo')
library(dismo)


ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_3_5.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## best fit models for each species
ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]

formulae <- data.frame(cbind(ordsums$species.code,ordsums$formu), stringsAsFactors = F)
str(formulae)
colnames(formulae) <- c("species","formula")
formulae$bio_2 <- sapply(formulae$formula, grepl, pattern="bio_2", fixed=T)
formulae$bio_5 <- sapply(formulae$formula, grepl, pattern="bio_5", fixed=T)
formulae$bio_6 <- sapply(formulae$formula, grepl, pattern="bio_6", fixed=T)
formulae$bio_8 <- sapply(formulae$formula, grepl, pattern="bio_8", fixed=T)
formulae$bio_12 <- sapply(formulae$formula, grepl, pattern="bio_12", fixed=T)
formulae$bio_15 <- sapply(formulae$formula, grepl, pattern="bio_15", fixed=T)
formulae$nlcd_3 <- sapply(formulae$formula, grepl, pattern="nlcd_3", fixed=T)
formulae$nlcd_4 <- sapply(formulae$formula, grepl, pattern="nlcd_4", fixed=T)
formulae$nlcd_5 <- sapply(formulae$formula, grepl, pattern="nlcd_5", fixed=T)
formulae$nlcd_6 <- sapply(formulae$formula, grepl, pattern="nlcd_6", fixed=T)
formulae$nlcd_7 <- sapply(formulae$formula, grepl, pattern="nlcd_7", fixed=T)
formulae$nlcd_8 <- sapply(formulae$formula, grepl, pattern="nlcd_8", fixed=T)
head(formulae)


spp <- formulae$species

bio <- stack("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc" )

cellIDs <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
#plot(cellIDs)
#plot(bio$bio_2, add=T)

bio <- stack(bio,cellIDs) # "layer" = cellID
# bio$
# abun_all <-
# <-  
#full <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T)
#abun <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/thinned_abundance_pts_2_9_2018.csv", sep=",", header=T)

full_all <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T,
                       stringsAsFactors = F)
full_all$DATA_SOURCE <- NULL



abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- colnames(full_all)
abun_all <- abun_all[abun_all$PLANT_CODE %in% spp,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)


full_all <- rbind(full_all, abun_all)


coordinates(abun_all) <- c(2,3)
ext <- extract(bio, abun_all) ## extract climate values to points
abun_all <- data.frame(cbind(abun_all,ext))
abun_all$optional <- NULL

head(abun_all)
abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)] <- abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)]/10 

coordinates(full_all) <- c(2,3)
ext <- extract(bio$layer, full_all)
full_all <- data.frame(cbind(full_all,ext))
head(full_all)
full_all$optional <-NULL
colnames(full_all) <- c("PLANT_CODE", "CELL", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
#length(unique(full_all$c.714442..733912..710327..708940..699188..699176..757497..715883..))





for (i in 17:length(spp)){
  #formulae[formulae]
  
  abun <- abun_all[abun_all$PLANT_CODE==spp[i],]
  full <- full_all[full_all$PLANT_CODE==spp[i],]
 
  cells <- unique(full$CELL)
  
  bio_id <- data.frame(bio[bio$layer %in% cells,])
  #head(bio_i)
  bio_id[,1:4] <- bio_id[,1:4]/10 
  #head(bio_id)
  
  ids <- data.frame(seq(1:ncell(bio)))
  colnames(ids) <- c("layer")
  
  bio_id <- merge(bio_id, ids, by="layer",all=T)
  bio_i <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
  
  if(formulae$bio_2[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_2)
    names(bio_2) <- "bio_2"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$bio_2 <- NULL
  }
  
  if(formulae$bio_5[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_5)
    names(bio_2) <- "bio_5"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$bio_5 <- NULL
  }
  
  if(formulae$bio_6[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_6)
    names(bio_2) <- "bio_6"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$bio_6 <- NULL
  }
  
  if(formulae$bio_8[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_8)
    names(bio_2) <- "bio_8"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$bio_8 <- NULL
  }
  
  if(formulae$bio_12[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_12)
    names(bio_2) <- "bio_12"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$bio_12 <- NULL
  }
  
  if(formulae$bio_15[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_15)
    names(bio_2) <- "bio_15"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$bio_15 <- NULL
  }
  
  if(formulae$nlcd_3[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_3)
    names(bio_2) <- "nlcd3"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_3 <- NULL
  }
  
  if(formulae$nlcd_4[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_4)
    names(bio_2) <- "nlcd4"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_4 <- NULL
  }
  
  if(formulae$nlcd_5[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_5)
    names(bio_2) <- "nlcd5"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_5 <- NULL
  }
  
  if(formulae$nlcd_6[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_6)
    names(bio_2) <- "nlcd6"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_6 <- NULL
  }
  
  if(formulae$nlcd_7[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_7)
    names(bio_2) <- "nlcd7"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_7 <- NULL
  }
  
  if(formulae$nlcd_8[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_8)
    names(bio_2) <- "nlcd8"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_8 <- NULL
  }
  
  bio_i <- dropLayer(bio_i, 1)
  
  reference <- abun[,2:(NCOL(abun)-3)]
  
  mess <- mess(bio_i, reference, full=FALSE) 
  #plot(mess)
  
  messD <- as.data.frame(mess)
  messD <- data.frame(messD[messD$mess < Inf,])
  colnames(messD) <- "mess"
  
  
  ordsums$negative[ordsums$species.code==spp[i]] <- length(messD$mess[messD$mess < 0])/length(messD$mess)
  
  write.csv(i, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/BIAS_progress.csv", row.names=F)
  print(i)
    
}

#head(full)
#head(abun)
write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_BIAS.csv", row.names=F)

summary(ordsums$negative)
