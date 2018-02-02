#### Bias data sets
###  Mitch O'Neill
##   created 1/31/2018
#    last modified: 2/01/2018

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/")

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

library(rgdal) ## for spatial data
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


### FULL ABUNDANCE BIAS DATASET 
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_1_25_2018.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd <- edd[edd$USDAcode %in% codes$x,]
abun_occurences <- data.frame(cbind(edd$USDAcode,edd$Longitude_Decimal,edd$Latitude_Decimal,rep("EDDMaps2016", NROW(edd))), stringsAsFactors=F)
colnames(abun_occurences) <- colnames(hotspots)
abun_occurences <- abun_occurences[order(abun_occurences$PLANT_CODE),]

### ABUNDANCE DATASET FOR MODELLING
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_1_25_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd <- edd[edd$species %in% codes$x,]
abun_modeling <- data.frame(cbind(edd$species,edd$longitude,edd$latitude,rep("EDDMaps2016", NROW(edd))), stringsAsFactors=F)
colnames(abun_modeling) <- colnames(hotspots)
abun_modeling <- abun_modeling[order(abun_modeling$PLANT_CODE),]

dir.create("MaxEntFiles")
write.csv(full_occurences, "MaxEntFiles/full_bias_pts.csv", row.names=F)
write.csv(abun_occurences, "MaxEntFiles/Abun_bias_pts.csv", row.names=F)
write.csv(full_modeling,  "MaxEntFiles/full_model_pts.csv", row.names=F)
write.csv(abun_modeling,  "MaxEntFiles/Abun_model_pts.csv", row.names=F)
