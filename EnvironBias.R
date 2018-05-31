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
#     #### Mess analysis from abun pts to full pts
#     library(rgdal)
#     library(raster)
#     #install.packages('dismo')
#     library(dismo)
#     
#     
#     ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_FINAL.csv", header = T, sep = ",",  
#                           quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
#     ## best fit models for each species
#     #ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]
#     
#     formulae <- data.frame(cbind(ordsums$species.code,ordsums$formu), stringsAsFactors = F)
#     str(formulae)
#     colnames(formulae) <- c("species","formula")
#     formulae$bio_2 <- sapply(formulae$formula, grepl, pattern="bio_2", fixed=T)
#     formulae$bio_5 <- sapply(formulae$formula, grepl, pattern="bio_5", fixed=T)
#     formulae$bio_6 <- sapply(formulae$formula, grepl, pattern="bio_6", fixed=T)
#     formulae$bio_8 <- sapply(formulae$formula, grepl, pattern="bio_8", fixed=T)
#     formulae$bio_12 <- sapply(formulae$formula, grepl, pattern="bio_12", fixed=T)
#     formulae$bio_15 <- sapply(formulae$formula, grepl, pattern="bio_15", fixed=T)
#     formulae$nlcd_3 <- sapply(formulae$formula, grepl, pattern="nlcd_3", fixed=T)
#     formulae$nlcd_4 <- sapply(formulae$formula, grepl, pattern="nlcd_4", fixed=T)
#     formulae$nlcd_5 <- sapply(formulae$formula, grepl, pattern="nlcd_5", fixed=T)
#     formulae$nlcd_6 <- sapply(formulae$formula, grepl, pattern="nlcd_6", fixed=T)
#     formulae$nlcd_7 <- sapply(formulae$formula, grepl, pattern="nlcd_7", fixed=T)
#     formulae$nlcd_8 <- sapply(formulae$formula, grepl, pattern="nlcd_8", fixed=T)
#     head(formulae)
#     
#     
#     spp <- formulae$species
#     
#     bio <- stack("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
#                  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc" )
#     
#     cellIDs <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
#     #plot(cellIDs)
#     #plot(bio$bio_2, add=T)
#     
#     bio <- stack(bio,cellIDs) # "layer" = cellID
#     # bio$
#     # abun_all <-
#     # <-  
#     #full <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T)
#     #abun <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/thinned_abundance_pts_2_9_2018.csv", sep=",", header=T)
#     
#     full_all <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T,
#                            stringsAsFactors = F)
#     full_all$DATA_SOURCE <- NULL
#     
#     
#     
#     abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
#                       comment.char= "", stringsAsFactors = F, strip.white = T)
#     abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
#     colnames(abun_all) <- colnames(full_all)
#     abun_all <- abun_all[abun_all$PLANT_CODE %in% spp,]
#     abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
#     abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)
#     
#     
#     full_all <- rbind(full_all, abun_all)
#     
#     
#     coordinates(abun_all) <- c(2,3)
#     ext <- extract(bio, abun_all) ## extract climate values to points
#     abun_all <- data.frame(cbind(abun_all,ext))
#     abun_all$optional <- NULL
#     
#     head(abun_all)
#     abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)] <- abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)]/10 
#     
#     coordinates(full_all) <- c(2,3)
#     extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
#     full_all <- full_all[extentShape,] ## 405534, make sure that points not within extent are not used
#     ext <- extract(bio$layer, full_all)
#     full_all <- data.frame(cbind(full_all,ext))
#     head(full_all)
#     full_all$optional <-NULL
#     colnames(full_all) <- c("PLANT_CODE", "CELL", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
#     
#     head(full_all)
#     #length(unique(full_all$c.714442..733912..710327..708940..699188..699176..757497..715883..))
#     
#     
#     
#     
#     
#     for (i in 3:length(spp)){
#       #formulae[formulae]
#       
#       abun <- abun_all[abun_all$PLANT_CODE==spp[i],]
#       full <- full_all[full_all$PLANT_CODE==spp[i],]
#      
#       cells <- unique(full$CELL)
#       
#       bio_id <- data.frame(bio[bio$layer %in% cells,])
#       #head(bio_i)
#       bio_id[,1:4] <- bio_id[,1:4]/10 
#       #head(bio_id)
#       
#       ids <- data.frame(seq(1:ncell(bio)))
#       colnames(ids) <- c("layer")
#       
#       bio_id <- merge(bio_id, ids, by="layer",all=T)
#       bio_i <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
#       
#       if(formulae$bio_2[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_2)
#         names(bio_2) <- "bio_2"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$bio_2 <- NULL
#       }
#       
#       if(formulae$bio_5[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_5)
#         names(bio_2) <- "bio_5"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$bio_5 <- NULL
#       }
#       
#       if(formulae$bio_6[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_6)
#         names(bio_2) <- "bio_6"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$bio_6 <- NULL
#       }
#       
#       if(formulae$bio_8[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_8)
#         names(bio_2) <- "bio_8"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$bio_8 <- NULL
#       }
#       
#       if(formulae$bio_12[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_12)
#         names(bio_2) <- "bio_12"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$bio_12 <- NULL
#       }
#       
#       if(formulae$bio_15[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$bio_15)
#         names(bio_2) <- "bio_15"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$bio_15 <- NULL
#       }
#       
#       if(formulae$nlcd_3[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_3)
#         names(bio_2) <- "nlcd3"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$nlcd_3 <- NULL
#       }
#       
#       if(formulae$nlcd_4[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_4)
#         names(bio_2) <- "nlcd4"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$nlcd_4 <- NULL
#       }
#       
#       if(formulae$nlcd_5[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_5)
#         names(bio_2) <- "nlcd5"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$nlcd_5 <- NULL
#       }
#       
#       if(formulae$nlcd_6[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_6)
#         names(bio_2) <- "nlcd6"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$nlcd_6 <- NULL
#       }
#       
#       if(formulae$nlcd_7[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_7)
#         names(bio_2) <- "nlcd7"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$nlcd_7 <- NULL
#       }
#       
#       if(formulae$nlcd_8[i]){
#         bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_8)
#         names(bio_2) <- "nlcd8"
#         bio_i<- stack(bio_i, bio_2)
#       } else {
#         abun$nlcd_8 <- NULL
#       }
#       
#       bio_i <- dropLayer(bio_i, 1)
#       
#       reference <- abun[,2:(NCOL(abun)-3)]
#       
#       mess <- mess(bio_i, reference, full=FALSE) 
#       #plot(mess)
#       
#       messD <- as.data.frame(mess)
#       colnames(messD) <- "mess"
#       #messD <- data.frame(messD[messD$mess < Inf,])
#       #messD <- data.frame(messD[!is.na(messD$mess),])
#       #colnames(messD) <- "mess"
#       messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
#       #length(messD[messD < Inf])
#       
#       ordsums$mess_min[i] <- min(messD)
#       ordsums$mess_max[i] <- max(messD)
#       ordsums$q5[i]   <- quantile(messD, 0.05)
#       ordsums$q10[i]  <- quantile(messD, 0.10)
#       ordsums$g0[i]   <- length(messD[messD > 0])/length(messD)
#       ordsums$gn5[i]  <- length(messD[messD > -5])/length(messD)
#       ordsums$gn10[i] <- length(messD[messD > -10])/length(messD)
#       ordsums$gn20[i] <- length(messD[messD > -20])/length(messD)
#       ordsums$check[i]<- length(messD)/length(cells)
#       
#       
#       if(length(messD) != length(cells)) {print("DANGER WILL ROBINSON")}
#       
#       #### write out new "full" datasets with abundance points tacked on
#       full$CELL <- NULL
#       file.out <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/species/",spp[i],".csv")
#       write.csv(full, file.out, row.names=F)
#       
#       
#       #ordsums$negative[ordsums$species.code==spp[i]] <- length(messD$mess[messD$mess < 0])/length(messD$mess)
#       #write.csv(i, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/BIAS_progress.csv", row.names=F)
#       print(i)
#         
#     }
#     
#     #head(full)
#     #head(abun)
#     write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_MESS.csv", row.names=F)
#     
#     summary(ordsums$gn5)
#     summary(ordsums$gn10)
#     hist(ordsums$q5, breaks=50)
#     hist(ordsums$q10, breaks=50)
#     
#     length(ordsums$species.code[ordsums$gn10 > .9])
#     117/155
#     min(ordsums$q10[ordsums$gn10 > .9])
#     ## worst one included: 5% of cells ara worse than -50
#     ## 10% of cells are worse than 10 (-9.5)
#     
#     summary(ordsums$mess_max) ## up to 100
#     summary(ordsums$mess_min) ## down to -1758
#     



##################################
####### MESS of full range #######
##################################

library(rgdal)
library(raster)
#install.packages('dismo')
library(dismo)


ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_4_29_2018_updated.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## best fit models for each species
#ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]

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

#full_all <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T,
#                       stringsAsFactors = F)
#full_all$DATA_SOURCE <- NULL



abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental_NEW.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)
abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- c("PLANT_CODE","LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
abun_all <- abun_all[abun_all$PLANT_CODE %in% spp,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)


#full_all <- rbind(full_all, abun_all)
coordinates(abun_all) <- c(2,3)
ext <- extract(bio, abun_all) ## extract climate values to points
abun_all <- data.frame(cbind(abun_all,ext))
abun_all$optional <- NULL

head(abun_all)
abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)] <- abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)]/10 

#coordinates(full_all) <- c(2,3)
#extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
#full_all <- full_all[extentShape,] ## 405534, make sure that points not within extent are not used
#ext <- extract(bio$layer, full_all)
#full_all <- data.frame(cbind(full_all,ext))
#head(full_all)
#full_all$optional <-NULL
#colnames(full_all) <- c("PLANT_CODE", "CELL", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")

#head(full_all)
#length(unique(full_all$c.714442..733912..710327..708940..699188..699176..757497..715883..))


ids <- data.frame(seq(1:ncell(bio)))
colnames(ids) <- c("layer")


for (i in 1:length(spp)){
  #formulae[formulae]
  
  abun <- abun_all[abun_all$PLANT_CODE==spp[i],]

  
  full <- raster(paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",spp[i],".asc",sep=""))
  
  names(full) <- "suitability"
  
  full <- stack(full,bio)
  
  #full$bio_2[full$suitability==0] <- NA
  #full$bio_5[full$suitability==0] <- NA
  #full$bio_6[full$suitability==0] <- NA
  #full$bio_8[full$suitability==0] <- NA
  #full$bio_12[full$suitability==0] <- NA
  #full$bio_15[full$suitability==0] <- NA
  #full$nlcd_3[full$suitability==0] <- NA
  #full$nlcd_4[full$suitability==0] <- NA
  #full$nlcd_5[full$suitability==0] <- NA
  #full$nlcd_6[full$suitability==0] <- NA
  #full$nlcd_7[full$suitability==0] <- NA
  #full$nlcd_8[full$suitability==0] <- NA
  
  
  ##########################
  cells <- unique(full$layer[full$suitability==1])
  
  bio_id <- data.frame(bio[bio$layer %in% cells,])
  #head(bio_i)
  bio_id[,1:4] <- bio_id[,1:4]/10 
  #head(bio_id)
  
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
    names(bio_2) <- "nlcd_3"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_3 <- NULL
  }
  
  if(formulae$nlcd_4[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_4)
    names(bio_2) <- "nlcd_4"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_4 <- NULL
  }
  
  if(formulae$nlcd_5[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_5)
    names(bio_2) <- "nlcd_5"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_5 <- NULL
  }
  
  if(formulae$nlcd_6[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_6)
    names(bio_2) <- "nlcd_6"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_6 <- NULL
  }
  
  if(formulae$nlcd_7[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_7)
    names(bio_2) <- "nlcd_7"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_7 <- NULL
  }
  
  if(formulae$nlcd_8[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$nlcd_8)
    names(bio_2) <- "nlcd_8"
    bio_i<- stack(bio_i, bio_2)
  } else {
    abun$nlcd_8 <- NULL
  }
  
  bio_i <- dropLayer(bio_i, 1)
  ####################################
  #full<-dropLayer(full,c(1,14))

  reference <- abun[,2:(NCOL(abun)-3)]
  
  mess <- mess(bio_i, reference, full=FALSE) 
  #plot(mess)
  
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  #messD <- data.frame(messD[messD$mess < Inf,])
  #messD <- data.frame(messD[!is.na(messD$mess),])
  #colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  #length(messD[messD < Inf])
  
  #ordsums$ALLVARSmess_min[i] <- min(messD)
  #ordsums$ALLVARSmess_max[i] <- max(messD)
  #ordsums$ALLVARSq5[i]   <- quantile(messD, 0.05)
  #ordsums$ALLVARSq10[i]  <- quantile(messD, 0.10)
  #ordsums$ALLVARSg0[i]   <- length(messD[messD > 0])/length(messD)
  #ordsums$ALLVARSgn5[i]  <- length(messD[messD > -5])/length(messD)
  ordsums$MESSgn10[i] <- length(messD[messD > -10])/length(messD)
  ordsums$MESSg.e.n10[i] <- length(messD[messD >= -10])/length(messD)
  #ordsums$MESSgn20[i] <- length(messD[messD > -20])/length(messD)
  #ordsums$ALLVARScheck[i]<- length(messD)/length(cells)
  
  
  #if(length(messD) != length(cells)) {print("DANGER WILL ROBINSON")}
  
  #### write out new "full" datasets with abundance points tacked on
  #full$CELL <- NULL
  #file.out <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/species/",spp[i],".csv")
  #write.csv(full, file.out, row.names=F)
  
  
  #ordsums$negative[ordsums$species.code==spp[i]] <- length(messD$mess[messD$mess < 0])/length(messD$mess)
  #write.csv(i, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/BIAS_progress.csv", row.names=F)
  print(i)
  
}



length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.9 & ordsums$kappa >0 & ordsums$kappaP < 0.05])  ## 64
length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.85 & ordsums$kappa >0 & ordsums$kappaP < 0.05]) ## 73
length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.80 & ordsums$kappa >0 & ordsums$kappaP < 0.05]) ## 91
length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.75 & ordsums$kappa >0 & ordsums$kappaP < 0.05]) ## 95

length(ordsums$species.code[ordsums$MESSgn10 > 0.95 & ordsums$kappa >0 & ordsums$kappaP < 0.05])  ## 64


hist(ordsums$kappa[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.9 & ordsums$kappa >0 & ordsums$kappaP < 0.05], main="All variables", xlab="kappa")
median(ordsums$AUC)


length(ordsums$species.code[ordsums$MESSgn10 >.9])
length(ordsums$species.code[ordsums$kappaP <0.05 & ordsums$kappa > 0])

hist(ordsums$MESSgn10, xlim=c(0.1,1), main="All variables", xlab="proportion not extrapolated")

hist(ordsums$AUC, xlim=c(0.3,1), main="All variables", xlab="AUC", breaks=100)

write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_29_2018.csv", row.names=F)
hist(ordsums$MESSgn10, breaks=30)

#########################################################
############ ALL VARS
#########################################################
library(rgdal)
library(raster)
#install.packages('dismo')
library(dismo)


ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_29_2018.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## best fit models for each species
#ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]
spp <- ordsums$species

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
bio$bio_2 <- bio$bio_2/10
bio$bio_5 <- bio$bio_5/10
bio$bio_6 <- bio$bio_6/10
bio$bio_8 <- bio$bio_8/10

#cellIDs <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
#plot(cellIDs)
#plot(bio$bio_2, add=T)

#bio <- stack(bio,cellIDs) # "layer" = cellID
# bio$
# abun_all <-
# <-  
#full <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T)
#abun <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/thinned_abundance_pts_2_9_2018.csv", sep=",", header=T)

#full_all <- read.table("C:/Users/Localadmin/Documents/Google_Drive_Overflow_3_9_2018/MaxEntFiles/full_model_pts.csv", sep=",", header=T,
#                       stringsAsFactors = F)
#full_all$DATA_SOURCE <- NULL



abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental_NEW.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)
abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- c("PLANT_CODE","LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
abun_all <- abun_all[abun_all$PLANT_CODE %in% spp,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)


#full_all <- rbind(full_all, abun_all)
coordinates(abun_all) <- c(2,3)
ext <- extract(bio, abun_all) ## extract climate values to points
abun_all <- data.frame(cbind(abun_all,ext))
abun_all$optional <- NULL

head(abun_all)
abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)] <- abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)]/10 


#ids <- data.frame(seq(1:ncell(bio)))
colnames(ids) <- c("layer")


for (i in 11:length(spp)){
  #formulae[formulae]
  if(ordsums$kappaP[i] < 0.05 & ordsums$kappa[i]> 0 & ordsums$AUC > 0.7){
  abun <- abun_all[abun_all$PLANT_CODE==spp[i],]
  
  
  full <- raster(paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",spp[i],".asc",sep=""))
  
  names(full) <- "suitability"
  
  full <- stack(full,bio)
  

  
  ##########################
  #cells <- unique(full$layer[full$suitability==1])
  
  full$bio_2[full$suitability == 0] <- NA
  full$bio_5[full$suitability == 0] <- NA
  full$bio_6[full$suitability == 0] <- NA
  full$bio_8[full$suitability == 0] <- NA
  full$bio_12[full$suitability == 0] <- NA
  full$bio_15[full$suitability == 0] <- NA
  full$nlcd_3[full$suitability == 0] <- NA
  full$nlcd_4[full$suitability == 0] <- NA
  full$nlcd_5[full$suitability == 0] <- NA
  full$nlcd_6[full$suitability == 0] <- NA
  full$nlcd_7[full$suitability == 0] <- NA
  full$nlcd_8[full$suitability == 0] <- NA
  
  #bio_id <- data.frame(bio[bio$layer %in% cells,])
  ##head(bio_i)
  #bio_id[,1:4] <- bio_id[,1:4]/10 
  ##head(bio_id)
  #
  #bio_id <- merge(bio_id, ids, by="layer",all=T)
  ##bio_i <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_id$)
  #
  #bio_i <- dropLayer(bio_i,1)
  ####################################
 full<-dropLayer(full,1)
  
  reference <- abun[,2:(NCOL(abun)-3)]
  
  mess <- mess(full, reference, full=FALSE) 
  #plot(mess)
  
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  #messD <- data.frame(messD[messD$mess < Inf,])
  #messD <- data.frame(messD[!is.na(messD$mess),])
  #colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  #length(messD[messD < Inf])
  
  #ordsums$ALLVARSmess_min[i] <- min(messD)
  #ordsums$ALLVARSmess_max[i] <- max(messD)
  #ordsums$ALLVARSq5[i]   <- quantile(messD, 0.05)
  #ordsums$ALLVARSq10[i]  <- quantile(messD, 0.10)
  #ordsums$ALLVARSg0[i]   <- length(messD[messD > 0])/length(messD)
  #ordsums$ALLVARSgn5[i]  <- length(messD[messD > -5])/length(messD)
  ordsums$MESSgn10ALL[i] <- length(messD[messD > -10])/length(messD)
  ordsums$MESSg.e.n10ALL[i] <- length(messD[messD >= -10])/length(messD)
  #ordsums$MESSgn20[i] <- length(messD[messD > -20])/length(messD)
  #ordsums$ALLVARScheck[i]<- length(messD)/length(cells)
  
  
  #if(length(messD) != length(cells)) {print("DANGER WILL ROBINSON")}
  
  #### write out new "full" datasets with abundance points tacked on
  #full$CELL <- NULL
  #file.out <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/species/",spp[i],".csv")
  #write.csv(full, file.out, row.names=F)
  
  
  #ordsums$negative[ordsums$species.code==spp[i]] <- length(messD$mess[messD$mess < 0])/length(messD$mess)
  #write.csv(i, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/BIAS_progress.csv", row.names=F)
  print(i)
}
}


length(ordsums$species[ordsums$MESSg.e.n10ALL > .9& ordsums$kappaP < 0.05 & ordsums$kappa> 0 & ordsums$AUC > 0.7])
hist(ordsums$MESSg.e.n10ALL)
write.csv(ordsums,"bloopALLMESS.csv")




######################################
######## COMPARATIVE MESS ANALYSES
#### FULL VARS VS SELECTED VARS
#### ABUN TO OCC
#### OCC TO RANGE
#### ABUN TO RANGE

library(rgdal)
library(raster)
#install.packages('dismo')
library(dismo)


ordsums <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/old_ordsums/ordsums_4_19_2018_updated.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## best fit models for each species
#ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]

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

bio$bio_2 <- bio$bio_2/10
bio$bio_5 <- bio$bio_5/10
bio$bio_6 <- bio$bio_6/10
bio$bio_8 <- bio$bio_8/10

cellIDs <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
#plot(cellIDs)
#plot(bio$bio_2, add=T)

bio <- stack(bio,cellIDs) # "layer" = cellID



occ.list <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/species", full=T)


occ_full<- read.table(occ.list[1], sep=",",header=T,stringsAsFactors = F) 
for (i in 2:length(occ.list)) {
 
 occ_i <- read.table(occ.list[i], sep=",",header=T,stringsAsFactors = F) 
 occ_full <- rbind(occ_full,occ_i)
 print(i)
}


head(occ_full)
coordinates(occ_full)<-c(2,3)
ext <- extract(bio$layer, occ_full)
occ_full <- cbind(occ_full,ext)
occ_full <- data.frame(occ_full)
occ_full$optional<-NULL
colnames(occ_full) <- c("PLANT_CODE","CELLID", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")

  
  occ_i <- occ_full[occ_full$PLANT_CODE == spp[1],]
  occ_thinned <- occ_i[1,]
  cellids <- unique(occ_i$CELLID)
  for (j in 2:length(cellids)){
    occ_i_j <- occ_i[occ_i$CELLID == cellids[j],][1,]
    occ_thinned <- rbind(occ_thinned, occ_i_j)
  }
  
  
for(i in 2:length(unique(occ_full$PLANT_CODE))){
  occ_i <- occ_full[occ_full$PLANT_CODE == spp[i],]
  cellids <- unique(occ_i$CELLID)
  for (j in 1:length(cellids)){
    occ_i_j <- occ_i[occ_i$CELLID == cellids[j],][1,]
    occ_thinned <- rbind(occ_thinned, occ_i_j)
  }

 print(i)  
}
occBACKUP <- occ_thinned
#NROW(occ_thinned)/NROW(occ_full)
NROW(abun_all)/NROW(occ_thinned) ## 32%
  
coordinates(occ_thinned)  <- c(3,4) 
ext <- extract(bio, occ_thinned)  
occ_thinned <-cbind(occ_thinned,ext)
head(occ_thinned)
#colnames(occ_thinned)


abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)
abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- c("PLANT_CODE","LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
abun_all <- abun_all[abun_all$PLANT_CODE %in% spp,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)
#full_all <- rbind(full_all, abun_all)
coordinates(abun_all) <- c(2,3)
ext <- extract(bio, abun_all) ## extract climate values to points
abun_all <- data.frame(cbind(abun_all,ext))
abun_all$optional <- NULL

head(abun_all)
#abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)] <- abun_all[,(ncol(abun_all)-14):(ncol(abun_all)-11)]/10 
summary(abun_all)
summary(occ_thinned)


#coordinates(full_all) <- c(2,3)
#extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
#full_all <- full_all[extentShape,] ## 405534, make sure that points not within extent are not used
#ext <- extract(bio$layer, full_all)
#full_all <- data.frame(cbind(full_all,ext))
#head(full_all)
#full_all$optional <-NULL
#colnames(full_all) <- c("PLANT_CODE", "CELL", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")

#head(full_all)
#length(unique(full_all$c.714442..733912..710327..708940..699188..699176..757497..715883..))


#ids <- data.frame(seq(1:ncell(bio)))
#colnames(ids) <- c("layer")


bioD <- as.data.frame(bio)
summary(bioD)


ordsums$MESS_abun2occ_SEL <- -99
ordsums$MESS_abun2range_SEL <- -99
ordsums$MESS_occ2range_SEL <- -99
ordsums$MESS_occ2range_ALL <- -99
ordsums$MESS_abun2range_ALL <- -99
ordsums$MESS_abun2occ_ALL <- -99


#redo 1,5; 1-var ones
for (i in 2:length(spp)){
  #formulae[formulae]
  if(is.na(ordsums$MESS_abun2occ_SEL[i] )){ 
  abun <- abun_all[abun_all$PLANT_CODE==spp[i],]
    abun_vars <- abun[,2:(NCOL(abun)-3)]
    
  occ <- data.frame(occ_thinned[occ_thinned$PLANT_CODE==spp[i],])
    occ_vars <- occ[,3:(NCOL(occ)-4)]
    
  range <- raster(paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",spp[i],".asc",sep=""))
  names(range) <- "suitability"

  ###########################
  #### USING ALL 12 VARIABLES
  ###########################

 
  ### abun to occ
  bio_occ <- bio
  bio_occ$bio_2[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$bio_5[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$bio_6[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$bio_8[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$bio_12[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$bio_15[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$nlcd_3[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$nlcd_4[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$nlcd_5[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$nlcd_6[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$nlcd_7[!(bio_occ$layer %in% occ$layer)] <- NA
  bio_occ$nlcd_8[!(bio_occ$layer %in% occ$layer)] <- NA
  
  #dropLayer(bio_occ$layer,13)
  bio_occ<- dropLayer(bio_occ, 13)
    
  mess <- mess(bio_occ, abun_vars, full=FALSE) 
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  ordsums$MESS_abun2occ_ALL[i] <- length(messD[messD > -10])/length(messD)

  #summary(as.data.frame(bio_occ))
  #summary(abun_vars)
  
    
  ## abun to range
  range <- stack(range,bio)
  range$bio_2[range$suitability==0] <- NA
  range$bio_5[range$suitability==0] <- NA
  range$bio_6[range$suitability==0] <- NA
  range$bio_8[range$suitability==0] <- NA
  range$bio_12[range$suitability==0] <- NA
  range$bio_15[range$suitability==0] <- NA
  range$nlcd_3[range$suitability==0] <- NA
  range$nlcd_4[range$suitability==0] <- NA
  range$nlcd_5[range$suitability==0] <- NA
  range$nlcd_6[range$suitability==0] <- NA
  range$nlcd_7[range$suitability==0] <- NA
  range$nlcd_8[range$suitability==0] <- NA
  range <- dropLayer(range,1)
  range <- dropLayer(range,13)
  
  mess <- mess(range, abun_vars, full=F)
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  ordsums$MESS_abun2range_ALL[i] <- length(messD[messD > -10])/length(messD)
  
  #  plot(mess)
  ## occ to range
  mess <- mess(range, occ_vars, full=F)
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  ordsums$MESS_occ2range_ALL[i] <- length(messD[messD > -10])/length(messD)
  
  #nlayers(range)
  #NCOL(occ_vars)

  
  
  #############################
  #### USING SELECTED VARIABLES
  #############################
  
  
  ##########################
  #cells <- unique(range$layer[full$suitability==1])
  
  # bio_id <- data.frame(bio[bio$layer %in% cells,])
  # bio_id[,1:4] <- bio_id[,1:4]/10 

  #bio_id <- merge(bio_id, ids, by="layer",all=T)
  #bio_i <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
  ###############################
  
  
  ##range raster
  range <- as.data.frame(range)
  bio_occ <- as.data.frame(bio_occ)
  range2 <- cellIDs
  occ2 <- cellIDs
  
  if(formulae$bio_2[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_2)
    names(bio_2) <- "bio_2"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_2)
    names(bio_2) <- "bio_2"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$bio_2 <- NULL
    occ_vars$bio_2 <- NULL
  }
  
  if(formulae$bio_5[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_5)
    names(bio_2) <- "bio_5"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_5)
    names(bio_2) <- "bio_5"
    occ2<- stack(occ2, bio_2)
    
  } else {
    abun_vars$bio_5 <- NULL
    occ_vars$bio_5 <- NULL
  }
  
  if(formulae$bio_6[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_6)
    names(bio_2) <- "bio_6"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_6)
    names(bio_2) <- "bio_6"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$bio_6 <- NULL
    occ_vars$bio_6 <- NULL
  }
  
  if(formulae$bio_8[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_8)
    names(bio_2) <- "bio_8"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_8)
    names(bio_2) <- "bio_8"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$bio_8 <- NULL
    occ_vars$bio_8 <- NULL
  }
  
  if(formulae$bio_12[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_12)
    names(bio_2) <- "bio_12"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_12)
    names(bio_2) <- "bio_12"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$bio_12 <- NULL
    occ_vars$bio_12 <- NULL
  }
  
  if(formulae$bio_15[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_15)
    names(bio_2) <- "bio_15"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_15)
    names(bio_2) <- "bio_15"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$bio_15 <- NULL
    occ_vars$bio_15 <- NULL
  }
  
  if(formulae$nlcd_3[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_3)
    names(bio_2) <- "nlcd_3"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_3)
    names(bio_2) <- "nlcd_3"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$nlcd_3 <- NULL
    occ_vars$nlcd_3 <- NULL
  }
  
  if(formulae$nlcd_4[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_4)
    names(bio_2) <- "nlcd_4"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_4)
    names(bio_2) <- "nlcd_4"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$nlcd_4 <- NULL
    occ_vars$nlcd_4 <- NULL
  }
  
  if(formulae$nlcd_5[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_5)
    names(bio_2) <- "nlcd_5"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_5)
    names(bio_2) <- "nlcd_5"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$nlcd_5 <- NULL
    occ_vars$nlcd_5 <- NULL
  }
  
  if(formulae$nlcd_6[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_6)
    names(bio_2) <- "nlcd_6"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_6)
    names(bio_2) <- "nlcd_6"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$nlcd_6 <- NULL
    occ_vars$nlcd_6 <- NULL
  }
  
  if(formulae$nlcd_7[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_7)
    names(bio_2) <- "nlcd_7"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_7)
    names(bio_2) <- "nlcd_7"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$nlcd_7 <- NULL
    occ_vars$nlcd_7 <- NULL
  }
  
  if(formulae$nlcd_8[i]){
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_8)
    names(bio_2) <- "nlcd_8"
    range2<- stack(range2, bio_2)
    
    bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_8)
    names(bio_2) <- "nlcd_8"
    occ2<- stack(occ2, bio_2)
  } else {
    abun_vars$nlcd_8 <- NULL
    occ_vars$nlcd_8 <- NULL
  }
  
  range2 <- dropLayer(range2, 1)
  occ2 <- dropLayer(occ2, 1)
  
  abun_vars <- data.frame(abun_vars)
  occ_vars <- data.frame(occ_vars)
  
  if(nlayers(occ2)== 1){
    
    range2b <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=as.vector(range2))
    occ2b <-   raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=as.vector(occ2))
    #colnames(abun_vars) <- "layer"
    #colnames(occ_vars) <- "layer"
    range2 <- stack(range2,range2b)
    occ2 <- stack(occ2,occ2b)
    abun_vars$layer <- abun_vars[,1]
    occ_vars$layer <- occ_vars[,1]
    
  }
  ####################################MESS_abun2range_ALL

  ### abun to occ
  mess <- mess(occ2, abun_vars, full=FALSE) 
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  ordsums$MESS_abun2occ_SEL[i] <- length(messD[messD > -10])/length(messD)
  
  ## abun to range
  mess <- mess(range2, abun_vars, full=FALSE) 
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  ordsums$MESS_abun2range_SEL[i] <- length(messD[messD > -10])/length(messD)
  
  ## occ to range  
  mess <- mess(range2, occ_vars, full=FALSE) 
  messD <- as.data.frame(mess)
  colnames(messD) <- "mess"
  messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
  ordsums$MESS_occ2range_SEL[i] <- length(messD[messD > -10])/length(messD)
  } 
  
  print(i)
  
}


summary(ordsums)
ordsums$MESS_abun2occ_SEL[ordsums$MESS_abun2occ_SEL== -99] <- NA
ordsums$MESS_abun2range_SEL[ordsums$MESS_abun2range_SEL== -99] <- NA
ordsums$MESS_occ2range_SEL[ordsums$MESS_occ2range_SEL== -99] <- NA
summary(ordsums)

length(ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05 & ordsums$AUC > 0.7 &
                              ordsums$MESS_abun2occ_ALL > 0.9]) ## 59
length(ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05 & ordsums$AUC > 0.7 &
                              ordsums$MESS_abun2range_ALL > 0.9]) ## 36
length(ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05 & ordsums$AUC > 0.7 &
                              ordsums$MESS_abun2occ_SEL > 0.9]) ## 93
length(ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05 & ordsums$AUC > 0.7 &
                              ordsums$MESS_abun2range_SEL > 0.9]) ## 70

spplist59 <- ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05 & ordsums$AUC > 0.7 &
                                    ordsums$MESS_abun2occ_ALL > 0.9]
spplist70 <- ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05 & ordsums$AUC > 0.7 &
                                    ordsums$MESS_abun2range_SEL > 0.9]
write.csv(spplist59, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/spplist59.csv", row.names=F)
write.csv(spplist70, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/spplist70.csv", row.names=F)

hist(ordsums$MESS_occ2range_ALL[ordsums$AUC >0.95])
hist(ordsums$MESS_occ2range_SEL[ordsums$AUC >0.95])

dev.off()
plot(ordsums$MESS_abun2occ_ALL~ordsums$MESSgn10)
abline(h=.9)
abline(v=.9)

writr.table

plot(ordsums$MESS_abun2occ_ALL~ordsums$MESSgn10)
abline(h=.9)
abline(v=.9)
plot(ordsums$MESS_abun2range_SEL~ordsums$MESSgn10)
abline(h=.9)
abline(v=.9)

write.csv(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_4_2018.csv", row.names=F)

#length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.9 & ordsums$kappa >0 & ordsums$kappaP < 0.05])  ## 64
#length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.85 & ordsums$kappa >0 & ordsums$kappaP < 0.05]) ## 73
#length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.80 & ordsums$kappa >0 & ordsums$kappaP < 0.05]) ## 91
#length(ordsums$species.code[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.75 & ordsums$kappa >0 & ordsums$kappaP < 0.05]) ## 95
#
#length(ordsums$species.code[ordsums$MESSgn10 > 0.95 & ordsums$kappa >0 & ordsums$kappaP < 0.05])  ## 64
#
#
#hist(ordsums$kappa[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.9 & ordsums$kappa >0 & ordsums$kappaP < 0.05], main="All variables", xlab="kappa")
#median(ordsums$AUC)
#
#
#length(ordsums$species.code[ordsums$MESSgn10 >.9])
#length(ordsums$species.code[ordsums$kappaP <0.05 & ordsums$kappa > 0])
#
#hist(ordsums$MESSgn10, xlim=c(0.1,1), main="All variables", xlab="proportion not extrapolated")
#
#hist(ordsums$AUC, xlim=c(0.3,1), main="All variables", xlab="AUC", breaks=100)

write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_29_2018.csv", row.names=F)
hist(ordsums$MESSgn10, breaks=30)



######################
### bias diagrams
######################

occ_all <- as.data.frame(occ_thinned)
head(occ_all)
summary(bioD$bio_2)

full.names <- c("mean diurnal range (°C)",
                "max temp of warmest mo. (°C)",
                "min temp of coldest mo. (°C)",
                "mean temp of wettest quar. (°C)",
                "annual precipitation (mm)",
                "precip. seasonality (mm)",
                "deciduous/mixed",
                "pine forest",
                "shrubland",
                "grassland",
                "pasture/hay",
                "cultivated crops")



png("C:/Users/Localadmin/Google Drive/GRC_figs/envi_bias.png",width=1500,height=800)
par(mfrow=c(3,4), cex.main=1.2, cex.lab=1.5, cex.axis=1.2, lwd=2, mar=c(5, 4, 4, 2) + 0.1 ) 
for (i in 1:12){ ## loop through each bioclimatic variable
  d.i <- density(bioD[,i], na.rm=T)  ## density plot for L48
  d.i2 <- density(occ_all[,(i+2)], na.rm=T)    ## density plot for point data
  #d.i3 <- density(edd2[,(i+6)], na.rm=T)   ## density plot for all grid cells w point data
  
  
    plot(d.i, col="black", main=letters[i], xlab=full.names[i], sub="", xlim = c(min(d.i$x,d.i2$x),max(d.i$x,d.i2$x)), ylim= c(min(d.i$y,d.i2$y),max(d.i$y,d.i2$y)))
    #lines(d.i3, col="brown")
    lines(d.i2, col="red")

    
  print(i)
  
}
dev.off()
# par(mfrow=c(3,4), cex.main=1.2, cex.lab=1.5, cex.axis=1.2, lwd=2, mar=c(5, 4, 4, 2) + 0.1) 
# for (i in 1:12){ ## loop through each bioclimatic variable
#   d.i <- density(bioD[,i], na.rm=T)  ## density plot for L48
#   d.i2 <- density(occ_full[,(i+2)], na.rm=T)    ## density plot for point data
#   #d.i3 <- density(edd2[,(i+6)], na.rm=T)   ## density plot for all grid cells w point data
#   
#   if(i==1 | i==2 | i==3 | i==6){
#     plot(d.i2, col="red", main=letters[i], xlab=full.names[i], sub="")
#     #lines(d.i3, col="brown")
#     lines(d.i)
#   } else {
#     plot(d.i, col="black", main=letters[i], xlab=full.names[i], sub="")
#     #lines(d.i3, col="brown")
#     lines(d.i2, col="red")
#   }
#   print(i)
#   
# }

#########################################################
########### non-extrap impact ranges ####################
#########################################################
library(raster)
library(dismo)

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

bio$bio_2 <- bio$bio_2/10
bio$bio_5 <- bio$bio_5/10
bio$bio_6 <- bio$bio_6/10
bio$bio_8 <- bio$bio_8/10

#cellIDs <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))

#bio <- stack(bio,cellIDs) # "layer" = cellID



ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_4_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  
#ordsums2 <- ordsums[ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > 0.9 & ordsums$kappa > 0 & ordsums$kappaP < 0.05,]  
ordsums2 <- ordsums

filenames <- data.frame(ordsums2$species.code)
colnames(filenames) <- "x"

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

filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN/",filenames$x,".asc")

abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)


abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- c("PLANT_CODE","LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
abun_all <- abun_all[abun_all$PLANT_CODE %in% filenames$x,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)
#full_all <- rbind(full_all, abun_all)
coordinates(abun_all) <- c(2,3)
ext <- extract(bio, abun_all) ## extract climate values to points
abun_all <- data.frame(cbind(abun_all,ext))
abun_all$optional <- NULL

dir.create( "C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/")

for (i in 108:length(spp)){
  if(ordsums$hiAbunSQKM[i]>0 & !is.na(ordsums$hiAbunSQKM[i])){

    abun <- abun_all[abun_all$PLANT_CODE==spp[i],]
    abun_vars <- abun[,2:(NCOL(abun)-2)]
    
    range <- raster(filenames$files[i])
    names(range) <- "suitability"

    range <- stack(range,bio)
    range$bio_2[range$suitability==0] <- NA
    range$bio_5[range$suitability==0] <- NA
    range$bio_6[range$suitability==0] <- NA
    range$bio_8[range$suitability==0] <- NA
    range$bio_12[range$suitability==0] <- NA
    range$bio_15[range$suitability==0] <- NA
    range$nlcd_3[range$suitability==0] <- NA
    range$nlcd_4[range$suitability==0] <- NA
    range$nlcd_5[range$suitability==0] <- NA
    range$nlcd_6[range$suitability==0] <- NA
    range$nlcd_7[range$suitability==0] <- NA
    range$nlcd_8[range$suitability==0] <- NA
    range <- dropLayer(range,1)
    #range <- dropLayer(range,13)
    
    #############################
    range2 <- range$nlcd_8
    range <- as.data.frame(range)
    

    
    if(formulae$bio_2[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_2)
      names(bio_2) <- "bio_2"
      range2<- stack(range2, bio_2)
      
      
    } else {
      abun_vars$bio_2 <- NULL
      
    }
    
    if(formulae$bio_5[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_5)
      names(bio_2) <- "bio_5"
      range2<- stack(range2, bio_2)
      
      
    } else {
      abun_vars$bio_5 <- NULL
   
    }
    
    if(formulae$bio_6[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_6)
      names(bio_2) <- "bio_6"
      range2<- stack(range2, bio_2)
 
    } else {
      abun_vars$bio_6 <- NULL

    }
    
    if(formulae$bio_8[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_8)
      names(bio_2) <- "bio_8"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$bio_8 <- NULL

    }
    
    if(formulae$bio_12[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_12)
      names(bio_2) <- "bio_12"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$bio_12 <- NULL

    }
    
    if(formulae$bio_15[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_15)
      names(bio_2) <- "bio_15"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$bio_15 <- NULL

    }
    
    if(formulae$nlcd_3[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_3)
      names(bio_2) <- "nlcd_3"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$nlcd_3 <- NULL

    }
    
    if(formulae$nlcd_4[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_4)
      names(bio_2) <- "nlcd_4"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$nlcd_4 <- NULL

    }
    
    if(formulae$nlcd_5[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_5)
      names(bio_2) <- "nlcd_5"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$nlcd_5 <- NULL

    }
    
    if(formulae$nlcd_6[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_6)
      names(bio_2) <- "nlcd_6"
      range2<- stack(range2, bio_2)
      

    } else {
      abun_vars$nlcd_6 <- NULL

    }
    
    if(formulae$nlcd_7[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_7)
      names(bio_2) <- "nlcd_7"
      range2<- stack(range2, bio_2)

    } else {
      abun_vars$nlcd_7 <- NULL

    }
    
    if(formulae$nlcd_8[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_8)
      names(bio_2) <- "nlcd_8"
      range2<- stack(range2, bio_2)

    } else {
      abun_vars$nlcd_8 <- NULL

    }
    
    range2 <- dropLayer(range2, 1)
    #abun_vars <- data.frame(abun_vars)

    
    if(nlayers(range2)== 1){
      
      range2b <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=as.vector(range2))

      #colnames(abun_vars) <- "layer"
      #colnames(occ_vars) <- "layer"
      range2 <- stack(range2,range2b)

      abun_vars$layer <- abun_vars[,1]

      
    }
    ################################################
    
    ## abun to range
    mess <- mess(range2, abun_vars, full=FALSE) 
    messD <- as.data.frame(mess)
    messD <- cbind(messD, as.data.frame(bio$bio_2))
    #colnames(messD) <- "mess"
    #messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
    #ordsums$MESS_abun2range_SEL[i] <- length(messD[messD > -10])/length(messD)
    messD$mess[messD$mess > -10 & messD$mess != Inf] <- 1
    messD$mess[messD$mess == Inf & !is.na(messD$bio_2)] <- 0
    messD$mess[messD$mess <= -10 & !is.na(messD$bio_2)] <- 0
    messD$mess[is.na(messD$bio_2)] <- NA

    
    mess <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=messD$mess)

    writeRaster(mess,paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/",spp[i],".asc"), overwrite=T)
}
    print(i)
  
}
