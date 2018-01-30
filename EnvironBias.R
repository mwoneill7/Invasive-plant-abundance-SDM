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
