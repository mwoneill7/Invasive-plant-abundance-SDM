#### Full workflow for identifying absence data, presence data, dealing with spp names, etc.
###  Mitch O'Neill
##   created 6/5/2017
#    last modified: 8/22/2017

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files")

edd <- read.table("file:///C:/Users/mwone/Documents/EDDMapS data/EDDMaps_plants_10_18_2016.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
## All plant points from EDDMapS, delivered from Rebekah Wallace to Bethany Bradley on October 18, 2016
 
edd <- edd[edd$SUB_species != "spp" & edd$SUB_species != "spp." & edd$SUB_species != "Spp" & edd$SUB_species != "Spp." & 
            edd$SUB_species != "NULL" & edd$SUB_species!= "N/A" & edd$SUB_species != "n/a" & edd$SUB_species != "NA" & !is.na(edd$SUB_species),]

## makes a column with ScientificName (species name, with subspecies/variety/forma if provided)
edd$ScientificName[edd$SUB_variety == "N/A"] <- paste(edd$SUB_genus[edd$SUB_variety == "N/A"], edd$SUB_species[edd$SUB_variety == "N/A"], sep = " ")
edd$ScientificName[edd$SUB_variety != "N/A"] <- paste(edd$SUB_genus[edd$SUB_variety != "N/A"], edd$SUB_species[edd$SUB_variety != "N/A"], edd$SUB_variety[edd$SUB_variety != "N/A"], sep = " ")
edd$ReporterFULLName <- paste(edd$ReporterFName, edd$ReporterLName, sep = " ")
 
### Address the confounding of canopycover and percentcover
### from Rebekah Wallace
## IF: Canopycover IS NOT NULL and percentcover IS NULL; THEN Canopycover is the area covered by the invasive species
## IF: Percentcover IS NOT NULL; THEN percentcover is the area covered by the invasive species and canopycover
## (NULL or NOT NULL) would be the native vegetation's cover.
edd$invasivecover[edd$canopycover != "NULL" & edd$percentcover == "NULL"] <- edd$canopycover[edd$canopycover != "NULL" & edd$percentcover == "NULL"]
edd$invasivecover[edd$percentcover != "NULL"] <- edd$percentcover[edd$percentcover != "NULL"]
edd$invasivecover[edd$canopycover == "NULL" &  edd$percentcover == "NULL"] <- "NULL"

length(edd$invasivecover[edd$canopycover != "NULL" & edd$percentcover == "NULL"] )
#### 264871 canopy
length(edd$invasivecover[edd$percentcover != "NULL"] )
#### 184536 percent

edd <- edd[edd$negative == 0,] ## remove asence data because it is spatially biased and not usefull

#### GET RID OF UNNECESSARY COLUMNS
colnames(edd)

edd <- data.frame(edd$ScientificName,
            edd$infestedAreaInAcres, 
            edd$grossAreaInAcres,
            edd$invasivecover,
            edd$stemcount,
            edd$NumberObserved,
            edd$TXabundance,
            edd$AbundanceText,
            edd$PercentPlotCovered,
            edd$density,
            edd$ReporterFULLName,
            edd$ReporterOrg,
            edd$Latitude_Decimal,
            edd$Longitude_Decimal, stringsAsFactors=F)

colnames(edd) <- c("ScientificName",
                 "infestedAreaInAcres",
                 "grossAreaInAcres",
                 "invasivecover",
                 "stemcount",
                 "NumberObserved",
                 "TXabundance",
                 "AbundanceText",
                 "PercentPlotCovered",
                 "density",
                 "ReporterFULLName",
                 "ReporterOrg",
                 "Latitude_Decimal",
                 "Longitude_Decimal" )

### flag valid abundance data (0= not valid; 1= valid); only local extent data used as abundance data in this study

### local extent data was read in as character data because "NULL" is recognized as such in R
### convert to numeric class; note that "NULL" will be coerced to NA, which will cause a warning message
edd$infestedAreaInAcres <- as.numeric(edd$infestedAreaInAcres)
edd$grossAreaInAcres <- as.numeric(edd$grossAreaInAcres)

edd$validInfestedArea[(edd$grossAreaInAcres >= edd$infestedAreaInAcres) & (!is.na(edd$infestedAreaInAcres)) & edd$infestedAreaInAcres > 0] <- 1
## whenever gross area is larger than or equal infested are (as required by definition), the record is valid
edd$validInfestedArea[is.na(edd$infestedAreaInAcres)] <- 0 ## if infested area is missing, record is not valid
edd$validInfestedArea[is.na(edd$grossAreaInAcres)] <- 0    ## if gross area is missing, record is not valid
edd$validInfestedArea[edd$infestedAreaInAcres <= 0] <- 0   ## if infested area is 0, or less than zero, record is not valid **
edd$validInfestedArea[edd$ReporterFULLName == "US Army Corps of Engineers Ombil Database"] <- 0
edd$validInfestedArea[edd$grossAreaInAcres < edd$infestedAreaInAcres] <- 0
## wherever gross area is less than infested area, the record is not valid

summary(as.factor(edd$validInfestedArea))
## ensure that all records were classified as valid or not valid abundance data

## GET RID OF POINTS WITH NO ABUNDANCE DATA

edd$TXabundance[(edd$TXabundance == "None" | edd$TXabundance == "0" | edd$TXabundance == "NULL"|
                 edd$TXabundance == "" | edd$TXabundance == "Select Abundance" )] <- NA
  # where TXabundance implies absence or is NULL
  
edd$density[(edd$density == "None" | edd$density == "0 PLANTS" | edd$density == "Trace, NONE FOUND" | edd$density == "NULL" | 
             edd$density == "*" | edd$density == "Select Abundance" | edd$density == "" | edd$density == " " )] <- NA
  # where the entry for density implies absence or is NULL  
  
edd$AbundanceText[edd$AbundanceText == "None                          " | edd$AbundanceText == "NULL" ] <- NA
  # where entry for AbundanceText implies absences or is NULL
  
edd$invasivecover[edd$invasivecover == "00%" | edd$invasivecover == "none" | edd$invasivecover == "0" | edd$invasivecover == "0%        " | edd$invasivecover == "0         " | edd$invasivecover == "Select One" |
     edd$invasivecover == "NULL" | edd$invasivecover == " " | edd$invasivecover == "N/A" | edd$invasivecover == "unknown" | edd$invasivecover == "" | edd$invasivecover == "          " ] <- NA
  # where canopycover/percentcover implies absence
  
edd$PercentPlotCovered[edd$PercentPlotCovered == "0" | edd$PercentPlotCovered == "0,0" | edd$PercentPlotCovered == "00%" | edd$PercentPlotCovered == " 0%" |
     edd$PercentPlotCovered == "NULL"] <- NA # where percentplotcovered implies absence or is NULL
  
edd$stemcount[edd$stemcount == "0" | edd$stemcount == "0.000000" | edd$stemcount == "NULL" | edd$stemcount == "unspecified" | 
     edd$stemcount == "N/A" | edd$stemcount == " " | edd$stemcount == "" | edd$stemcount == "No Value"] <- NA
  # where stemcount implies absence or is NULL
  
edd$NumberObserved[edd$NumberObserved == "0" |
     edd$NumberObserved == "NULL" | edd$NumberObserved== "unknown" | edd$NumberObserved == "" | edd$NumberObserved == " " ] <- NA

edd <- edd[!(edd$validInfestedArea == 0 & is.na(edd$NumberObserved) & is.na(edd$invasivecover) & is.na(edd$stemcount) & 
             is.na(edd$TXabundance) & is.na(edd$AbundanceText) & is.na(edd$PercentPlotCovered) & is.na(edd$density)),]

# 1,529,053
# 1,528,622
# 911,422

## SUBSET ALL POINTS IN CONUS
## prep coordinate fields; rgdal requires numeric coordinates
edd$Latitude_Decimal <- as.numeric(edd$Latitude_Decimal)
edd$Longitude_Decimal <- as.numeric(edd$Longitude_Decimal)

summary(edd$Latitude_Decimal)
summary(edd$Longitude_Decimal)
## note: 441 nulls
edd <- edd[!is.na(edd$Latitude_Decimal) & !is.na(edd$Longitude_Decimal),]

## all points should be long: -65 to -130, lat: 20 to 55

#### one has a positive Longitude but is reported by Fairfax County Park Authority, VA
#### Switching sign will place it in Fairfax.
#edd$Longitude_Decimal[edd$Longitude_Decimal > 0 & edd$ReporterOrg == "Fairfax County Park Authority, VA"] <-  
#  -1*edd$Longitude_Decimal[edd$Longitude_Decimal > 0 & edd$ReporterOrg == "Fairfax County Park Authority, VA"]

## WRONG SIGN FOR LONGITUDE
probs <- edd[edd$Longitude_Decimal > 65 & edd$Longitude_Decimal < 130 & edd$Latitude_Decimal >20 & edd$Latitude_Decimal <55, ]#& edd$validInfestedArea == 1 & !is.na(edd$USDAcode), ] 
## NONE WITH WRONG SIGN FOR LATITUDE
probs <- edd[edd$Latitude_Decimal > -55 & edd$Latitude_Decimal < -20 & edd$Longitude_Decimal < -65 & edd$Longitude_Decimal > -130, ]
## NONE HAVE SWITCHING ERROR
probs <- edd[edd$Longitude_Decimal >20 & edd$Longitude_Decimal < 55 & edd$Latitude >-130 & edd$Longitude -65,]
rm(probs)

## fix positive latitudes
edd$Longitude_sign_change <- 0
edd$Longitude_sign_change[edd$Longitude_Decimal > 65 & edd$Longitude_Decimal < 130 & edd$Latitude_Decimal >20 & edd$Latitude_Decimal <55] <- 1
edd$Longitude_Decimal[edd$Longitude_Decimal > 65 & edd$Longitude_Decimal < 130 & edd$Latitude_Decimal >20 & edd$Latitude_Decimal <55] <-
  -1*edd$Longitude_Decimal[edd$Longitude_Decimal > 65 & edd$Longitude_Decimal < 130 & edd$Latitude_Decimal >20 & edd$Latitude_Decimal <55] 

library(rgdal) ## for spatial data
edd <- edd[!is.na(edd$Longitude_Decimal) & !is.na(edd$Latitude_Decimal) & edd$Longitude_Decimal < -65 & edd$Longitude_Decimal > -130 & edd$Latitude_Decimal > 20 & edd$Latitude_Decimal <55, ]
## removes records I know that I won't use (no coords, coords way outside US, not absence data or valid infested area data)
## 911422 to 782480
head(edd)

coordinates(edd) <- c(14,13) ## specifies long, lat

#all of EDDMapS is supposed to be nad83/wgs84
EPSG <- make_EPSG() ## creates library of EPSG codes to help assign proj4string
EPSG[grepl("WGS 84$", EPSG$note) == TRUE, ] ## search for WGS 84 --> 4326
proj4string(edd) <- CRS("+init=epsg:4326")

## clip pts to L48
## load contiguous US shape file
states = readOGR(dsn = "states", layer = "US_states")
proj4string(states) ## wgs84
## transform states to same proj4string as eddmaps
states <- spTransform(states, "+init=epsg:4326")
plot(states) ## plot eddmaps onto US
plot(edd, pch=20, add=T) ## plot points on top
edd  <- edd[states, ] ## clip edd_US to L48
plot(edd, pch=20, add=T, col="purple") 
## replot subsetted points in new color to ensure that it worked
dev.off()

edd <- data.frame(edd, stringsAsFactors = F) #### bring back into dataframe.format
edd$optional <- NULL ## remove meaningless column

## none by US Army Corps
length(edd$stemcount[edd$ReporterFULLName == "US Army Corps of Engineer" ])

## Area is continuous, check unique vals for others
unique(edd$PercentPlotCovered)
unique(edd$TXabundance)
unique(edd$AbundanceText)
unique(edd$stemcount)
unique(edd$NumberObserved)
unique(edd$invasivecover)
summary(as.factor(edd$density))

density <- data.frame(summary(as.factor(edd$density[!is.na(edd$density) & is.na(edd$PercentPlotCovered) &
                   is.na(edd$TXabundance) & is.na(edd$stemcount) & edd$validInfestedArea == 0 &
                   is.na(edd$NumberObserved) & is.na(edd$invasivecover) & is.na(edd$AbundanceText)])))

length(edd$density[!is.na(edd$density) & is.na(edd$PercentPlotCovered) &
                     is.na(edd$TXabundance) & is.na(edd$stemcount) & edd$validInfestedArea == 0 &
                     is.na(edd$NumberObserved) & is.na(edd$invasivecover) & is.na(edd$AbundanceText)])

edd$density <- as.character(edd$density)
edd$ReporterFULLName <- as.character(edd$ReporterFULLName)
edd$invasivecover <- as.character(edd$invasivecover)
edd$NumberObserved <- as.character(edd$NumberObserved)
edd$stemcount <- as.character(edd$stemcount)
edd$TXabundance <- as.character(edd$TXabundance)
edd$AbundanceText <- as.character(edd$AbundanceText)
edd$PercentPlotCovered <- as.character(edd$PercentPlotCovered)


# notfound <- c("Arctotheca prostrata",
# "Thymus praecox ssp. britannicus",
# "Campanula cervicaria",
# "Cecropia palmata",
# "Heteropterys brachiata",
# "Clinopodium menthifolium ssp. ascendens",
# "Dichrostachys cinerea ssp. africana",
# "Viburnum opulus ssp. opulus",
# "Cirsium clavatum var. osterhoutii",
# "Lepidium chalepense",
# "Nitellopsis obtusa",
# "Petasites fragrans",
# "Prunus subhirtella var. pendula")

# nf <- edd[edd$ScientificName %in% notfound,]
# summary(as.factor(nf$ScientificName))
# length(edd$ScientificName[edd$ScientificName == "Arctotheca calendula"])


### update species names (via USDA codes)
edd.list <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/speciesList_11_20_2017.csv",
                       header = T, sep = ",", stringsAsFactors = F, strip.white = T)
## download scientific name to USDA code file, made manually by Ceara, completed June 1, 2017
## compiled by Mitch O'Neill 11/20/2017

head(edd.list)

## list of species in USDA-sp name file
spp <- as.list(edd.list$eddmap_species)#[!is.na(edd.list$eddmap_species)])

for (sp in spp){
  edd$USDAcode[edd$ScientificName == sp] <- edd.list$USDA_code[edd.list$eddmap_species == sp]
  print(sp)  
}

edd <- edd[!is.na(edd$USDAcode),] #738742 to 738433

##### merge subspecies #####
edd$USDAcodeUNMERGED <- edd$USDAcode
## keep a copy of the column before merging, so that I can go back if necessary

ssp <- read.table("subspecies2.csv",  header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
## lists subspecies with the species they are to be merged into
head(ssp)
for (i in 1:length(ssp$usda.code)){ ## loop through subspecific codes in my data
  edd$USDAcode[edd$USDAcodeUNMERGED == ssp$usda.code[i]] <- ssp$usda.lump[i]
  ## reassign to lumped species-level code
  print(i)
}


### REMOVE AQUATICS
aqua <- read.table("aquatics.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
## list of aquatic species that haven't been eliminated yet

edd <- edd[!(edd$USDAcode %in% aqua$usda.code),]


write.csv(edd, "C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_11_25_2017.csv", row.names = F)

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_11_25_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

########### COUNT ############

edd$count[!is.na(edd$stemcount)] <- edd$stemcount[!is.na(edd$stemcount)]
edd$count[!is.na(edd$NumberObserved)] <- edd$NumberObserved[!is.na(edd$NumberObserved)] 
edd$count[!is.na(edd$TXabundance)] <- edd$TXabundance[!is.na(edd$TXabundance)] 

#count <- edd[!is.na(edd$count),] 
edd$countN <- as.numeric(edd$count)
edd$countC[is.na(edd$countN)] <- edd$count[is.na(edd$countN)]

b1 <- c("Single Plant", "Single stem", "One", "one")
b2t10 <- c("1 - 10 Plants", "2-5", "1-10", "0-10", "6-10", "2-10", "2-10                                              ")
b11t100 <- c("10 - 100 Plants", "10-20", "20-50", "50+", "10-100", "> 15", "<100                                              ",
             "20 - 99", "20-99", "20-99 ", "11-15", "11-100", "<100", "11-100                                            ")
b101t1000 <- c("100 - 500 Plants", "500 Plants +", "100 - 1000 Plants", "100 - 999", "100-500",
               "100-999", "100-999 ", "101-1000", "100-1,000", "101-1,000", "100-1000", "1,000")
bG1000 <- c(">1000", ">1000 ", "More than 1000", ">1,000", "> 1000 Plants")
bL20 <- c("<20", "<20 ", "Less than 20")

few      <- c("Few", "few")
rare     <- "Rare"
uncommon <- "Uncommon"
common   <- "Common"
abundant <- "Abundant"
many     <- c("A Lot", "many", "many SELECTED", "many ", "Many")  

edd$countQUANT[edd$countC %in% b1] <- "1 plant" 
edd$countQUANT[edd$countC %in% b2t10] <- "2-10 plants"
edd$countQUANT[edd$countC %in% b11t100] <- "11-100 plants"
edd$countQUANT[edd$countC %in% b101t1000] <- "101-1000 plants"
edd$countQUANT[edd$countC %in% bG1000] <- ">1000 plants"
edd$countQUANT[edd$countC %in% bL20] <- "<20 plants"
edd$countQUANT[edd$countN == 1] <- "1 plant" 
edd$countQUANT[edd$countN >= 2 & edd$countN <= 10] <- "2-10 plants"
edd$countQUANT[edd$countN >= 11 & edd$countN <= 100] <- "11-100 plants"
edd$countQUANT[edd$countN >= 101 & edd$countN <= 1000] <- "101-1000 plants"
edd$countQUANT[edd$countN >= 1001] <- ">1000 plants"

summary(as.factor(edd$countQUANT)) 

edd$countQUAL[edd$countC %in% few]      <- "few"
edd$countQUAL[edd$countC %in% rare]     <- "rare"
edd$countQUAL[edd$countC %in% uncommon] <- "uncommon"
edd$countQUAL[edd$countC %in% common]   <- "common"
edd$countQUAL[edd$countC %in% abundant] <- "abundant"
edd$countQUAL[edd$countC %in% many]     <- "many"

############## COVER ##############

edd$coverN<-as.numeric(edd$invasivecover)
edd$coverC[is.na(edd$coverN)] <- edd$invasivecover[is.na(edd$coverN)]
edd$coverP <- as.numeric(gsub("%", "", edd$coverC))


trace <- c("Trace (Less than 1%)",
           "(T) Trace = less than 1%",
           "<1%       ",
           "<2%                                               ",
           "0 - 2.5%  ",
           "<2.5%     ",
           "COVER C <2%",
           "trace-low ",
           "trace",
           "Trace",
           "Trace     ",
           "trace     " )

low <- c("(L) Low = 1%-5%",
         "Low (1.0 - 5.0%)",
         "1 = <5%, n",
         "1-5%      ",
         "1 - 5%    ",
         "1-5       ", 
         "<5%       ",
         "<5%, few (",
         "<5%, singl",
         "<5%, numer",
         "low-medium",
         "low-med   ",
         "low       ",
         "Low       ",
         "low",
         "Lo        ",
         "Low" )


mod <- c("Moderate (5.1 - 25%)",
         "5-25%     ",
         "2 = 5 - 25",
         "5.1-25%   ",
         "5 - 25%   ",
         "6-25%     ",
         "6 To 25   ",
         "6-25       ",
         "(M) = 5.1%-25%",
         "Moderate: 6-25%",
         "6-25      ",
         "COVER C 2 to 25%",
         "2-25%                                             ",
         "0-25%     ",
         "<25%",
         "2.5 - 20% ",
         "Moderate  ",
         "Moderate ",
         "moderate",
         "moderate  ",
         "Moderate",
         "Medium",
         "medium    ",
         "med-high  " )

high <- c("26-50%                                            ",
          "3 = 25 - 5",
          "26-50%    ",
          "25-50%    ",
          "25-50%",
          "26-50     ",
          "High: 26-50%",
          "COVER C 26 to 50%",
          "COVER C 25 to 50%",
          "24-75%    ",
          "25.1-100% ",
          "High (25.1-100%)",
          "(H) High = 25.1%-100%",
          ">25%      ",
          "High",
          "high      ",
          "High      ",
          "HIgh      ",
          "high",
          "20 - 50%  ")

majority <- c("50-75%",
              "50-75%    ",
              "51-75%    ",
              "4 = 50 - 7",
              "75-95%    ",
              "95-100%   ",
              "5 = 75 - 1",
              "76-100%   ",
              "75-100%",
              "75-100%   ",
              ">75%      ",
              "COVER C >50%",
              "COVER C >50% Geo",
              "51-100    ",
              "Majority: 51-100%",
              "50 - 80%  ",
              "80 - 97.5%" )

edd$cover[edd$coverC %in% trace] <- "trace"
edd$cover[edd$coverC %in% low] <- "low"
edd$cover[edd$coverC %in% mod] <- "medium"
edd$cover[edd$coverC %in% high] <- "high"
edd$cover[edd$coverC %in% majority] <- "majority"

edd$cover[edd$coverN == 1] <- "trace" ##observations under 1 are ignored
edd$cover[edd$coverN > 1 & edd$coverN <= 5] <- "low"
edd$cover[edd$coverN > 5 & edd$coverN <= 25] <- "medium"
edd$cover[edd$coverN > 25 & edd$coverN <= 50] <- "high"
edd$cover[edd$coverN > 50 & edd$coverN <= 100] <- "majority"

edd$cover[edd$coverP <= 1] <- "trace"  
edd$cover[edd$coverP > 1 & edd$coverP <= 5] <- "low"
edd$cover[edd$coverP > 5 & edd$coverP <= 25] <- "medium"
edd$cover[edd$coverP > 25 & edd$coverP <= 50] <- "high"
edd$cover[edd$coverP > 50 & edd$coverP <= 100] <- "majority"

summary(as.factor(edd$cover))


### DENSITY
edd$densityN <- as.numeric(edd$density)
edd$density[!is.na(edd$densityN)] <- NA
edd$densityN <- NULL

for(i in 1:length(edd$density[!is.na(edd$density)])){
  
  edd$dPERC[!is.na(edd$density)][i] <- grepl("%", edd$density[!is.na(edd$density)][i])
 
  if(grepl("PLANT", edd$density[!is.na(edd$density)][i], ignore.case=T) |
     grepl("ROSETTE", edd$density[!is.na(edd$density)][i], ignore.case=T)){
       edd$dPLANT[!is.na(edd$density)][i] <- TRUE
     } else {
       edd$dPLANT[!is.na(edd$density)][i] <- FALSE
     }
 
  if(grepl("low", edd$density[!is.na(edd$density)][i], ignore.case=T) |
     grepl("med", edd$density[!is.na(edd$density)][i], ignore.case=T) |
     grepl("high", edd$density[!is.na(edd$density)][i], ignore.case=T) |
     grepl("maj", edd$density[!is.na(edd$density)][i], ignore.case=T) |
     grepl("moderate", edd$density[!is.na(edd$density)][i], ignore.case=T) |
     grepl("trace", edd$density[!is.na(edd$density)][i], ignore.case=T)){#)
      edd$dTLMH[!is.na(edd$density)][i] <- TRUE
  } else {
      edd$dTLMH[!is.na(edd$density)][i] <- FALSE
    }
   
  print(i)
}

#write.csv(edd, "C://Users/mwone/Google Drive/eddDENSITY.csv",row.names=F)
edd <- read.table("C://Users/mwone/Google Drive/eddDENSITY.csv", header = T, sep = ",", quote= "\"", 
                   comment.char= "", stringsAsFactors = F, strip.white = T)

edd$countQUANT[edd$countQUANT == ">20 plants"] <- "<20 plants" ## this is really <20 from a typo

edd <- edd[!is.na(edd$countQUANT) | !is.na(edd$cover) | edd$dTLMH == TRUE
           | edd$dPLANT == TRUE | edd$dPERC == TRUE | edd$validInfestedArea == 1 | !is.na(edd$AbundanceText),]

edd <-edd[!is.na(edd$ScientificName),]

edd$Dcount[edd$dPLANT == TRUE & is.na(edd$countQUANT) & !is.na(edd$density)] <- edd$density[edd$dPLANT == TRUE & is.na(edd$countQUANT) & !is.na(edd$density)]
edd$countFromDensity<- 0
edd$countFromDensity[edd$dPLANT == TRUE & is.na(edd$countQUANT) & !is.na(edd$density)]  <- 1

edd$Dcover[edd$dPERC == TRUE & is.na(edd$cover) & !is.na(edd$density)] <- edd$density[edd$dPERC == TRUE & is.na(edd$cover)  & !is.na(edd$density)]
edd$Dcover[edd$dTLMH == TRUE & is.na(edd$cover) & !is.na(edd$density)] <- edd$density[edd$dTLMH == TRUE & is.na(edd$cover) & !is.na(edd$density)]
edd$coverFromDensity <- 0
edd$coverFromDensity[edd$dPERC == TRUE & is.na(edd$cover) & !is.na(edd$density)] <-1
edd$coverFromDensity[edd$dTLMH == TRUE & is.na(edd$cover) & !is.na(edd$density)] <-1

#unique(edd$Dcount)
b1 <- c("Trace, 1 PLANT","Single Plant","1 Plant","1 PLANT")
b2t10 <- c("1 - 10 Plants","2 PLANTS","7 PLANTS","5 PLANTS" ,"3 PLANTS","8 PLANTS","4 PLANTS","10 PLANTS","6 PLANTS","9 PLANTS","Trace, 2 PLANTS","TD;  4 PLANTS","TD;  SEVERAL PLANTS","Trace, 5 PLANTS","Trace, 10 PLANTS")
b11t100 <- c("LD;  75-100  PLANTS","20 X 50 FT.  50 PLANTS","TD;  DOZ PLANTS","TD;  1 DOZ PLANTS","20-30 PLANTS","TD;  DOZEN PLANTS","Trace, 12 PLANTS","10 - 100 Plants","23 PLANTS","20 ROSETTES","14 PLANTS","15 PLANTS","20 PLANTS","40 PLANTS","13 PLANTS","12 PLANTS","50 PLANTS","30 PLANTS","25 PLANTS","11 PLANTS","17 PLANTS","18 PLANTS","67 PLANTS","100 PLANTS","90 PLANTS","45 PLANTS")
b101t1000 <- c("100 - 1000 Plants" )
bG1000 <- c("> 1000 Plants")
bL20 <- c("1-20 plants")
#few <- c("LD;  FEW PLANTS OVER .25","FEW PLANTS")

edd$countQUANT[edd$Dcount %in% b1] <- "1 plant" 
edd$countQUANT[edd$Dcount %in% b2t10] <- "2-10 plants"
edd$countQUANT[edd$Dcount %in% b11t100] <- "11-100 plants"
edd$countQUANT[edd$Dcount %in% b101t1000] <- "101-1000 plants"
edd$countQUANT[edd$Dcount %in% bG1000] <- ">1000 plants"
edd$countQUANT[edd$Dcount %in% bL20] <- "<20 plants"
#edd$countQUAL[edd$Dcount %in% few] <- "few"

trace <- c("Trace","Trace, 5 PLANTS","Trace, 50 FT","Trace, NONE","Trace, .10 ACRE","Trace, 10X10","Trace, 1 PLANT","Trace, 2 PLANTS","Trace, 12 PLANTS","Trace, 10 PLANTS","Trace, 40X15","Trace, 100X100","Trace, 150X100","Trace, POINT SOURCE","TRACE, 0.01","TRACE, 0.25","TRACE, 1.75","TRACE, 0.1","TRACE, 1","Trace: Less than 1%","TRACE, 0.5","<2%","TRACE")
low <- c("LOW, 0.01","1 to 5 %","0 - 5%","0% - 5 %","0 % - 5 %","0% - 5%","0%-5%","Low, 100 FT","Low, 100X75","Low, 50X50","Low, N/A","Low, 1.5 ACRE","Low, 30 FT","Low, 15X200","Low, 50 FT","Low, 3 ACRE","Low, 40X20","Low, 30X10","Low, 10X10","Low, 30X30","LOW, 0.1","Low, POINT SOURCE","LOW, 0.25","LOW, 0.5","LOW, 1.5","LOW, 5","Low, .10 - .25 ACRE","Low, .50 - 1 ACRES","LOW, 2.5","Low, .25 - .50 ACRES","LOW, 10","LOW, 0.2","LOW, 1","LOW, 2","Low: 1 to 5%","Low","0-10%","LOW","<5%")
mod <- c("Moderate","6% - 25 %","6%-25%","6% - 25%","6 - 25%","moderate,small","moderate","Medium, 2.5 ACRE","Medium, 100+ ACRES","Medium, 25X150","Medium, 0.5 ACRE","Medium, 20X40","Medium, 50X100","Medium, 20X20","Medium, 4 ACRE","Medium, 60X30","Medium, 1 AC","Medium, 0.25 AC","Medium, 200X40","Medium, 5X10","Medium, 10X10","Medium, 30X200","Medium, 50X20","Medium, 60X20","Medium, 100X100","Medium, 1 ACRE","MEDIUM, 0.01","Medium, .10 AC","MEDIUM, 3.5","MEDIUM, 0.5","MEDIUM, 0.25","MEDIUM, 0.1","MEDIUM, 2.5","Medium, POINT SOURCE","Medium, .10 - .25 ACRE","MEDIUM, 2","MEDIUM, 1","Medium, .50 - 1 ACRES","MEDIUM, 0.33","MEDIUM, 0.75","MEDIUM, 5","MEDIUM, 3","MEDIUM, 4","MEDIUM, 1.5","MEDIUM, 0.2","Medium, .25 - .50 ACRES","Moderate: 6 to 25%","Medium","11-25%","2 - 25%","MEDIUM","5-25%",">25%","10%","20%")
high <- c("High, 1 ACRE","50%","26 % - 50 %","26% - 50%","26% - 50 %","26% -50%","High, 100 FT","High, 60 FT","High, 300 FT","High, 0.5 AC","High, 0.25 AC","High, 5X5","High, 0.15 ACRES","High, 60X50","High, 40X40","High, 100X100","High, 50X40","High, 50X20","High, 5 ACRE +","High, .5 ACRE","High, 20 ACRES","High, 50X50","High, 40X300","HIGH, 1.5","HIGH, 2.25","HIGH, 0.01","HIGH, 10","HIGH, 0.1","HIGH, 0.33","HIGH, 3","High, .10 - .25 ACRE","HIGH, 1.25","High, .50 - 1 ACRES","HIGH, 0.25","High, POINT SOURCE","High, .25 - .50 ACRES","HIGH, 1","HIGH, 0.75","HIGH, 0.5","HIGH, 25","HIGH, 5","HIGH, 2","HIGH, 0.51","HIGH, 40","HIGH, 300","HIGH, 7","HIGH, 0.2","High: 26 to 50%","High: 26-100%","High","26-50%","26 - 50%","HIGH","30%","40%")
majority <- c("Majority: 51% or greater","51-75%","76-100%",">50%","51-100%","70%","90%","80%")

edd$cover[edd$Dcover %in% trace] <- "trace"
edd$cover[edd$Dcover %in% low] <- "low"
edd$cover[edd$Dcover %in% mod] <- "medium"
edd$cover[edd$Dcover %in% high] <- "high"
edd$cover[edd$Dcover %in% majority] <- "majority"

edd <- edd[!is.na(edd$cover) | !is.na(edd$countQUANT) | edd$validInfestedArea == 1 | !is.na(edd$AbundanceText),]

edd$countORD[edd$countQUANT == "1 plant"] <- 1
edd$countORD[edd$countQUANT == "2-10 plants" | edd$countQUANT == "<20 plants"] <- 2
edd$countORD[edd$countQUANT == "11-100 plants"] <- 3
edd$countORD[edd$countQUANT == "101-1000 plants" | edd$countQUANT == ">1000 plants"] <- 4 
hist(edd$countORD)
hist(edd$countORD[edd$countFromDensity == 0])
hist(edd$countORD[edd$countFromDensity == 1])
summary(as.factor(edd$Dcount))
#18127+27160
#unique(edd$density[!is.na(edd$density) & is.na(edd$countQUANT) & edd$dPLANT==T & is.na(edd$Dcount)])

edd$coverORD[edd$cover == "trace"] <- 1
edd$coverORD[edd$cover == "low"] <- 2
edd$coverORD[edd$cover == "medium"] <- 3
edd$coverORD[edd$cover == "high" | edd$cover == "majority"] <- 4
hist(edd$coverORD)
hist(edd$coverORD[edd$coverFromDensity == 1])
hist(edd$coverORD[edd$coverFromDensity == 0])
unique(edd$AbundanceText)

edd$textORD[edd$AbundanceText == "Single Plant                  "] <- 1
edd$textORD[edd$AbundanceText == "Scattered Plants              " | edd$AbundanceText == "Linearly scattered            "] <- 2
edd$textORD[edd$AbundanceText == "Scattered Dense Patches       "] <- 3
edd$textORD[edd$AbundanceText == "Dominant cover                " | edd$AbundanceText == "Dense Monoculture             "] <- 4
hist(edd$textORD)

cor.test(edd$textORD,edd$coverORD, method="spearman") ## .34
cor.test(edd$textORD,edd$countORD, method="spearman") ## .64
cor.test(edd$countORD,edd$coverORD, method="spearman") ## .36
cor.test(edd$coverORD,edd$infestedAreaInAcres, method="spearman")  ## .44
cor.test(edd$textORD,edd$infestedAreaInAcres, method="spearman")   ## .41
cor.test(edd$countORD,edd$infestedAreaInAcres, method="spearman")  ## .63

plot(edd$coverORD,log(edd$infestedAreaInAcres))
plot(edd$textORD, log(edd$infestedAreaInAcres))
plot(edd$countORD,log(edd$infestedAreaInAcres))
#write.csv(edd2, "C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_11_29_2017.csv", row.names=F)

edd$ORD[!is.na(edd$coverORD)] <- edd$coverORD[!is.na(edd$coverORD)]
edd$ORD[!is.na(edd$textORD)] <- edd$textORD[!is.na(edd$textORD)]
edd$ORD[!is.na(edd$countORD)] <- edd$countORD[!is.na(edd$countORD)]

hist(edd$coverORD)
hist(edd$countORD)
hist(edd$textORD)
hist(as.numeric(edd$ORD))

#######################################

length(edd$countORD[!is.na(edd$countORD) & edd$countFromDensity==1])/length(edd$countORD[!is.na(edd$countORD) ])
## 73 percent of counts from density
length(edd$coverORD[!is.na(edd$coverORD) & edd$coverFromDensity==1])/length(edd$coverORD[!is.na(edd$coverORD) ])
## 29 percent of covers from density

edd$ORD[!is.na(edd$countORD)] <- edd$countORD[!is.na(edd$countORD)]

edd$logAREA[edd$validInfestedArea == 1] <- log(edd$infestedAreaInAcres[edd$validInfestedArea == 1])

quantile(edd$logAREA[edd$ORD==1],na.rm=T)
#exp(-3.863046)
plot(edd$logAREA~as.factor(edd$countORD))
plot(edd$logAREA~as.factor(edd$coverORD))
plot(edd$logAREA~as.factor(edd$textORD))
plot(edd$logAREA[!is.na(edd$textORD) | !is.na(edd$countORD)]~as.factor(edd$ORD[!is.na(edd$textORD) | !is.na(edd$countORD)]))

exp(quantile(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==1], .5, na.rm=T))
exp(quantile(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==2], .5, na.rm=T))
exp(quantile(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==3], .5, na.rm=T))
exp(quantile(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==4], .5, na.rm=T))

exp(mean(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==1 & !is.na(edd$logAREA)]))
exp(mean(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==2 & !is.na(edd$logAREA)]))
exp(mean(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==3 & !is.na(edd$logAREA)]))
exp(mean(edd$logAREA[(!is.na(edd$textORD) | !is.na(edd$countORD)) & edd$ORD==4 & !is.na(edd$logAREA)]))

##########################################################
edd$infestedAreaInAcres[edd$validInfestedArea == 0]<-NA
edd$infORD[edd$infestedAreaInAcres <= 20/43560] <- 1
edd$infORD[edd$infestedAreaInAcres <= 500/43560 & edd$infestedAreaInAcres> 20/43560] <- 2
edd$infORD[edd$infestedAreaInAcres <= .25 & edd$infestedAreaInAcres> 500/43560] <- 3
edd$infORD[edd$infestedAreaInAcres > .25] <- 4

cor.test(edd$infORD,edd$countORD,method="pearson")
cor.test(edd$infORD,edd$textORD,method="pearson")
cor.test(edd$infORD,edd$coverORD,method="pearson")
cor.test(edd$infORD,as.numeric(edd$ORD),method="pearson")
#############################################################
hist(edd$infORD)
hist(log(edd$infestedAreaInAcres[edd$infORD==4]))
hist(edd$infORD[is.na(edd$ORD)])
hist(edd$textORD[!is.na(edd$infORD)])
hist(edd$ORD)
hist(edd$textORD)

plot(edd$logAREA[!is.na(edd$textORD) | !is.na(edd$countORD)]~as.factor(edd$ORD[!is.na(edd$textORD) | !is.na(edd$countORD)]))
abline(h=log(10/43560))
abline(h=log(250/43560))
abline(h=log(.25))

plot(edd$logAREA~as.factor(edd$ORD))
abline(h=log(20/43560))
abline(h=log(500/43560))
abline(h=log(.25))

plot(factor(edd$ORD, levels=c(1,2,3,4),ordered=T)~factor(edd$infORD, levels=c(1,2,3,4), ordered=T))

edd$ORD[is.na(edd$ORD) & !is.na(edd$infORD)] <- edd$infORD[is.na(edd$ORD) & !is.na(edd$infORD)]
summary(edd$ORD)
hist(edd$ORD)

###############################

edd$ORD[!is.na(edd$infORD)] <- edd$infORD[!is.na(edd$infORD)]
edd$ORD[!is.na(edd$coverORD)] <- edd$coverORD[!is.na(edd$coverORD)]
edd$ORD[!is.na(edd$textORD)] <- edd$textORD[!is.na(edd$textORD)]
edd$ORD[!is.na(edd$countORD)] <- edd$countORD[!is.na(edd$countORD)]

edd$ORDsource[!is.na(edd$infORD)] <- "area"
edd$ORDsource[!is.na(edd$coverORD)] <- "cover"
edd$ORDsource[!is.na(edd$textORD)] <- "text"
edd$ORDsource[!is.na(edd$countORD)] <- "count"


summary(as.factor(edd$ORDsource))
269141/663262   ## 41% cover
229523/663236   ## 35% text
123118/(663262) ## 19% count
41480/(663262)  ## 6% area


