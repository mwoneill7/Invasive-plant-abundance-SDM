#### Environmental bias; MESS; etc.
###  Mitch O'Neill
##   created 9/21/2017
#    last modified 9/25/2017

setwd("C:/Users/mwone/Documents/geodata/")
library(raster)

## read in thinned dataset
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_09_12_2017.csv",  
                  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

## read in species list
spp <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList08_28_1017.csv",  
                  header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
sp.list <- spp$usda.code[spp$potential.use2 == 1] ## list of useable species

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

cell.ids <- unique(edd$cellID) ## list of unique cell IDs in the point data

edd2 <- edd[1,] ## first row of point dataset, to initiate loop

for (j in 2:length(cell.ids)){ ## loop through all cells for that species
  cell.j <- edd[edd$cellID == cell.ids[j],][1,] ## keep 1 record from cell
  edd2 <- rbind(edd2, cell.j) ## append the kept row to the master object
  print(j)
}

length(edd2$id) == length(unique(edd2$cellID)) ## should be true

## full names of bioclimatic variables, for later titles
fullnames <- c("Mean Diurnal Range",
               "Max Temp of Warmest Month",
               "Min Temp of Coldest Month",
               "Mean Temp of Wettest Quarter",
               "Annual Precip",
               "Precip Seasonality")

bio.stackD <- as.data.frame(bio.stack)
## make dataframe version to make density plots
 
for (i in 1:length(colnames(bio.stackD))){ ## loop through each bioclimatic variable
  d.i <- density(bio.stackD[,i], na.rm=T)  ## density plot for L48
  d.i2 <- density(edd[,(i+6)], na.rm=T)    ## density plot for point data
  d.i3 <- density(edd2[,(i+6)], na.rm=T)   ## density plot for all grid cells w point data
  plot(d.i2, col="red", main=fullnames[i], sub="")
  lines(d.i3, col="brown")
  lines(d.i)

  print(i)

}
