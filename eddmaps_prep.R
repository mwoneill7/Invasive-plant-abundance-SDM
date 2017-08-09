#### Full workflow for identifying absence data, presence data, dealing with spp names, etc.
###  Mitch O'Neill
##   6/27/2017

###### updating species names, etc. of EDDMapS dataset
#####  code by: Mitch O'Neill
####   created: 6/5/2017
###    last modified: 6/16/2017

setwd("C:/Users/mwone/Google Drive")

edd <- read.table("Fall 2016/EDDMapS_data_expl_2017/EDDMaps_plants_10_18_2016.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
## this is the full invasive plant dataset EDDMaps, obtained by Bethany Bradley on October 18, 2016

head(edd)

## makes a column with ScientificName (species name, with subspecies/variety/forma if provided)
edd$ScientificName[edd$SUB_variety == "N/A"] <- paste(edd$SUB_genus[edd$SUB_variety == "N/A"], edd$SUB_species[edd$SUB_variety == "N/A"], sep = " ")
edd$ScientificName[edd$SUB_variety != "N/A"] <- paste(edd$SUB_genus[edd$SUB_variety != "N/A"], edd$SUB_species[edd$SUB_variety != "N/A"], edd$SUB_variety[edd$SUB_variety != "N/A"], sep = " ")


### Address the confounding of canopycover and percentcover
### from Rebekah Wallace
## IF: Canopycover IS NOT NULL and percentcover IS NULL; THEN Canopycover is the area covered by the invasive species
## IF: Percentcover IS NOT NULL; THEN percentcover is the area covered by the invasive species and canopycover
## (NULL or NOT NULL) would be the native vegetation's cover.
edd$invasivecover[edd$canopycover != "NULL" & edd$percentcover == "NULL"] <- edd$canopycover[edd$canopycover != "NULL" & edd$percentcover == "NULL"]
edd$invasivecover[edd$percentcover != "NULL"] <- edd$percentcover[edd$percentcover != "NULL"]
edd$invasivecover[edd$canopycover == "NULL" &  edd$percentcover == "NULL"] <- "NULL"
summary(edd$invasivecover)
unique(edd$invasivecover)



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
## **absences/zero abundance are addressed separately, below.
edd$validInfestedArea[edd$grossAreaInAcres < edd$infestedAreaInAcres] <- 0
## wherever gross area is less than infested area, the record is not valid

unique(edd$validInfestedArea)
summary(edd$validInfestedArea)
## ensure that all records were classified as valid or not valid abundance data


### update species names (via USDA codes)
sp.list <- read.table("NSF_GSS_shared/Synonomies/speciesList06_01_2017.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T)
## download scientific name to USDA code file, compiled manually by Ceara, completed June 1, 2017

## list of species in USDA-sp name file
spp <- as.list(sp.list$eddmap_species[!is.na(sp.list$eddmap_species)])

for (sp in spp){
  length(edd$ScientificName[is.na(edd$invasivecover)])
  
  edd$USDAcode[edd$ScientificName == sp] <- sp.list$USDA_code[sp.list$eddmap_species == sp & !is.na(sp.list$eddmap_species)]
  print(sp)  
}



## flag absences (0=occurence, 1 = potential absence)
edd$potential.absence <- 0 ## default is that a row is not absence data

## any records that have an entry in the abundance column that resembles absence/zero abundance
## assigned 1 (possible absence)
edd$potential.absence[edd$TXabundance == "None" | edd$TXabundance == "0" |
                        # where TXabundance implies absence
                        edd$grossAreaInAcres == 0 | edd$infestedAreaInAcres == 0 | 
                        # records were infested acres, gross acres, or both are zero
                        edd$density == "None" | edd$density == "0 PLANTS" | edd$density == "Trace, NONE FOUND" |
                        # where the entry for density implies absences
                        edd$AbundanceText == "None                          " |
                        # where entry for AbundanceText implies absences
                        edd$invasivecover == "00%" | edd$invasivecover == "none" | edd$invasivecover == "0" | edd$invasivecover == "0%        " | edd$invasivecover == "0         " |
                        # where canopycover/percentcover implies absence
                        edd$PercentPlotCovered == "0" | edd$PercentPlotCovered == "0,0" | edd$PercentPlotCovered == "00%" | edd$PercentPlotCovered == " 0%" |
                        # where percentplotcovered implies absence 
                        edd$stemcount == "0" | edd$stemcount == "0.000000" |
                        # where stemcount implies absence
                        edd$NumberObserved == "0"]  <- 1
## where 0 plants were observed


edd$potential.absence[edd$negative == 1] <- 2 ## these rows are directly recorded as absences
## and are thus interpreted as such. Some cleaning may be helpful





## out of potential absences, check for contradictions of other columns
edd$contradiction[edd$potential.absence == 1 | edd$potential.absence == 2] <- 1 ## contradicted by default, unless I show that there is no contradiction

edd$contradiction[(edd$potential.absence == 1 | edd$potential.absence == 2 ) & ## where each abundance field must NOT contain nonzero abundance data
                    
                    (edd$TXabundance == "None" | edd$TXabundance == "0" |
                       edd$TXabundance == "NULL"| edd$TXabundance == "" | edd$TXabundance == "Select Abundance" ) &
                    # where TXabundance implies absence or is NULL
                    
                    (edd$grossAreaInAcres == 0 | edd$infestedAreaInAcres == 0 | is.na(edd$infestedAreaInAcres) | is.na(edd$grossAreaInAcres)) &
                    # records were infested acres, gross acres, or both are zero or NULL
                    
                    (edd$density == "None" | edd$density == "0 PLANTS" | edd$density == "Trace, NONE FOUND" |
                       edd$density == "NULL" | edd$density == "*" | edd$density == "Select Abundance" | edd$density == "" | edd$density == " " ) &
                    # where the entry for density implies absence or is NULL  
                    
                    (edd$AbundanceText == "None                          " | edd$AbundanceText == "NULL" ) &
                    # where entry for AbundanceText implies absences or is NULL
                    
                    (edd$invasivecover == "00%" | edd$invasivecover == "none" | edd$invasivecover == "0" | edd$invasivecover == "0%        " | edd$invasivecover == "0         " | edd$invasivecover == "Select One" |
                       edd$invasivecover == "NULL" | edd$invasivecover == " " | edd$invasivecover == "N/A" | edd$invasivecover == "unknown" | edd$invasivecover == "" | edd$invasivecover == "          " ) &
                    # where canopycover/percentcover implies absence
                    
                    (edd$PercentPlotCovered == "0" | edd$PercentPlotCovered == "0,0" | edd$PercentPlotCovered == "00%" | edd$PercentPlotCovered == " 0%" |
                       edd$PercentPlotCovered == "NULL" ) & # where percentplotcovered implies absence or is NULL
                    
                    (edd$stemcount == "0" | edd$stemcount == "0.000000" |
                       edd$stemcount == "NULL" | edd$stemcount == "unspecified" | edd$stemcount == "N/A" | edd$stemcount == " " | edd$stemcount == "" ) &
                    # where stemcount implies absence or is NULL
                    
                    (edd$NumberObserved == "0" |
                       edd$stemcount == "NULL" | edd$stemcount == "unknown" | edd$stemcount == "" | edd$stemcount == " " )] <- 0

edd[edd$negative == 1 & edd$contradiction ==1, ]
## there is one contradiction, but upon further inspecting that record, the entry in the comment field clear indicates that this is an absence



## out of potential absences, check if the reporter (or if NULL, the organization) uses the negative column
edd$ReporterFULLName <- paste(edd$ReporterFName, edd$ReporterLName, sep = " ")
## make a reporter full name field by pasting FName and LName

length(edd$ReporterFULLName[edd$ReporterFULLName == "NULL NULL"])
## NOTE: the 18 records with out Reporter names end up getting the unique name "NULL NULL"
## As none of these records use the "negative field", they end up passing through this quality control test
## (i.e. while records where the reporter has been shown to use the negative column in other records will be 
## excluded, these reporter-less records will not be excluded)
## NOTE: Reporter Names were used because about half of the observations do not use organization
## PERHAPS LOOK INTO



#####################################################################################################
class(edd$dateentered)
head(as.Date(edd$dateentered))
head((edd$dateentered))
edd$dateenteredD <- as.Date(edd$dateentered)
summary(edd$dateenteredD)
min(edd$dateenteredD)

## loop through reporters to assess negative column use
rep.list <- as.list(unique(edd$ReporterFULLName[edd$potential.absence == 1]))


for (i in 1:length(rep.list)){
  
  
  
  ## if the reporter uses the absence colum
  if(length(edd$ReporterFULLName[edd$negative == 1 & edd$ReporterFULLName == rep.list[i]]) > 0){
    
    
    ## then a 1 is entered in the negative use column for all of that reporter's observations, AND
    edd$negUse[edd$ReporterFULLName == rep.list[i]] <- 1 
    
    
    ## for all potential absences from the reporter of the iteration, 
    edd$negUse2[edd$potential.absence == 1 & edd$ReporterFULLName == rep.list[i] & 
                  ## The ones entered earlier than than the first use of the negative column for that reporter are flagged as 0         
                  (edd$dateenteredD < min(edd$dateenteredD[edd$negative == 1 & edd$ReporterFULLName == rep.list[i]]))] <- 0 
    
    
    ## for all potential absences from the reporter of the iteration, 
    edd$negUse2[edd$potential.absence == 1 & edd$ReporterFULLName == rep.list[i] & 
                  ## The ones not entered earlier than than the first use of the negative column for that reporter are flagged as 1         
                  (edd$dateenteredD >= min(edd$dateenteredD[edd$negative == 1 & edd$ReporterFULLName == rep.list[i]]))] <- 1 
    
  } else { 
    
    ## if negative column is not used, then a 0 is entered in the negative use column
    edd$negUse[edd$ReporterFULLName == rep.list[i]] <- 0
    
  }
  
  
  
  
  print(i)
  
}




#### Review of absence data cleaning column definitions
####
####
#### potential.absences : (0) records with no indication of absence
####                      (1) records with an indication of absence in an abundance column
####                      (2) records with a (1) in the negative field
####
#### contradictions     : of the potential absences == 1 OR 2;
####                      (0) records with no contradictions
####                      (1) records with contradictions
####
#### negUse             : (0) the reporter of this record does not use the absence column
####                      (1) the reporter of this record uses the absence column
####
#### negUse2            : (0) this potential absence record was entered before the reporter used the column
####                      (1) this potential absence record was not entered before the reporter used the column
####

#### absences are:    potential.absences = 2   AND   contradictions = 0
####                  potential.absences = 1   AND   contradictions = 0   AND   negUse = 0
####                  potential.absences = 1   AND   contradictions = 0   AND   negUse = 1   AND   negUse2 = 0

#### using entry date does not save any data

edd$AbsenceData[(edd$potential.absence == 2) | 
                  (edd$potential.absence == 1 & edd$contradiction == 0 & edd$negUse == 0) |
                  (edd$potential.absence == 1 & edd$contradiction == 0 & edd$negUse == 1 & edd$negUse2 == 0)] <- 1

write.csv(edd, "eddmaps_prepped.csv", row.names = F)


## prep coordinate fields
edd$Latitude_N <- as.numeric(edd$Latitude_Decimal)
edd$Longitude_N <- as.numeric(edd$Longitude_Decimal)
summary(edd$Latitude_N)
summary(edd$Longitude_N)
edd1 <- edd[edd$validInfestedArea == 1 | edd$AbsenceData == 1, ]
edd1 <- edd1[!is.na(edd1$ScientificName), ]
summary(edd1$Latitude_N)
summary(edd1$Longitude_N)
edd1[edd1$Longitude_N > 0, ] # edd$LonN[edd$LongitudeN == 77.30791] <- -77.30791 ### one has a positive Longitude, but it says fairfax VA and switching sign will place it in Fairfax.
### CHeck out funky lay/lons



     
head(edd)
library(rgdal)
coordinates(edd) <- c(11,10) ## specifies long, lat, CHANGE NUMBERS


EPSG <- make_EPSG() ## creates library of EPSG codes to help assign proj4string
EPSG[grepl("WGS 84$", EPSG$note) == TRUE, ] ## search for WGS 84 --> 4326

#assume that all EDDMapS data is wgs??
proj4string(edd) <- CRS("+init=epsg:4326")

## check to see points on states file
## load contiguous US shape file
setwd("C://Users/mwone/Downloads/")
states = readOGR(dsn = "states", layer = "US_states")
proj4string(states) ## wgs84
## transform states to same proj4string as eddmaps
states2 <- spTransform(states, "+init=epsg:4326")
edd_US  <- edd[states2, ] ## clip edd_US to L48
plot(states2) ## plot eddmaps onto US
plot(edd_US, pch=20, add=T)
