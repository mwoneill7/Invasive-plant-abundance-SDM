#### Combining eddmaps species list with w/ invasive listings
###  Mitch O'Neill
##   created: 6/13/2017; 
#    last modified 8/9/2017

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files")

#### eddmaps species list (all species with infested area downloaded from EDDMAPS October 2016)
edd.list <- read.table("Raw/speciesList06_01_2017.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
edd.list$edd <- 1


#### list of species from the occurence-based hotspots analysis (Allen and Bradley 2016)
occ.list <- read.table("C:/Users/mwone/Google Drive/NSF_GSS_shared/Hotspots_and_Abundance/occ_data_with_2017_codes.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
occ.list$occ <- 1
head(occ.list)

#### USDA federal/State listings (Downloaded from PLANTS 5/30/2017)
usda.list <- read.table("Raw/PLANTS.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
usda.list$usdaCHK <- 1
#### subsets plants with federal and state listing status, not other invasive plant lists available in USDA database
head(usda.list)

plants.i  <- usda.list[grep("L48(I)",   usda.list$Native.Status, fixed=TRUE), ]
plants.n  <- usda.list[grep("L48(N)",   usda.list$Native.Status, fixed=TRUE), ]
plants.in <- usda.list[grep("L48(I,N)", usda.list$Native.Status, fixed=TRUE), ]
plants.i$L48  <- "I"
plants.n$L48  <- "N"
plants.in$L48 <- "I/N"
plants.nativity <- rbind(plants.i, plants.n, plants.in)
head(plants.nativity)

#### IPAUS species list (Downloaded in 2014; database last updated May 2012)
ipa.list <- read.table("Raw/ipa2.csv", header = T, sep = ",", 
                       stringsAsFactors = F, quote= "\"", comment.char= "", strip.white = T)
ipa.list$IPA <- 1
head(ipa.list)



## full eddmaps data download from 
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_08_09_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)


ed.list <- as.list(unique(edd$USDAcode))

usda.code <- "TEMPLATE"
scientific.name <- "TEMPLATE"
comment <- "TEMPLATE"
state   <- -99
federal <- -99
other   <- -99
ipaus   <- -99
hotspot <- -9
nativity<- "TEMPLATE"
usdaCHK <- -99
no.abs  <- -99
no.abun <- -99
no.big <- -99

ed.listing <- data.frame(scientific.name, usda.code, comment, state, federal, other, ipaus, hotspot, no.abun, no.abs, no.big, nativity, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)


for (i in 1:length(ed.list)){
  scientific.name <- edd.list$PLANTS_species[edd.list$USDA_code == ed.list[i]][1]
  comment <- edd.list$Comment[edd.list$USDA_code == ed.list[i]][1]
  usda.code <- ed.list[i]
  usda.i <- usda.list[usda.list$Accepted.Symbol == ed.list[i],]
  usda.i <- usda.i[!is.na(usda.i$Accepted.Symbol), ]
  
  
  # the length of where any are not NULL is > 0, then listing =1    
  if (length(usda.i$Accepted.Symbol[usda.i$State.Noxious.Status != ""]) > 0) {
    state <- 1} else {state <- 0}
  
  if (length(usda.i$Accepted.Symbol[usda.i$Federal.Noxious.Status == "NW"]) > 0) {
    federal <- 1} else  {federal <- 0}
  
  if (length(usda.i$Accepted.Symbol[usda.i$Invasive != ""]) > 0) {
    other <- 1}  else {other <- 0}
  
  if (length(ipa.list$habit[ipa.list$USDA_code_2017 == ed.list[i]]) > 0) {
    ipaus <- 1}  else {ipaus <- 0}

  if (length(occ.list$habit[occ.list$NewCode == ed.list[i]]) > 0){
    hotspot <- 1} else {hotspot <- 0}

  if(length(plants.nativity$L48[plants.nativity$Accepted.Symbol == ed.list[i]]) > 0) {
    L48 <- plants.nativity$L48[plants.nativity$Accepted.Symbol == ed.list[i]] 
  } else {L48 = "N/A"}
  
  if(length(usda.i$Accepted.Symbol[usda.i$usdaCHK == 1]) > 0) {
    usdaCHK <- 1} else {usdaCHK <- 0}
  
    
  edd.abun <- edd[edd$USDAcode == ed.list[i] & edd$validInfestedArea == 1, ]
  no.abun <- length(edd.abun$USDAcode[!is.na(edd.abun$USDAcode)])
  
  edd.abs <- edd[edd$USDAcode == ed.list[i] & edd$negative == 1, ]
  no.abs <- length(edd.abs$USDAcode[!is.na(edd.abs$USDAcode)])
  
  edd.big <-  edd[edd$USDAcode == ed.list[i] & (edd$infestedAreaInAcres > 20 | edd$ReporterFULLName == "US Army Corps of Engineers Ombil Database"), ]
  no.big <- length(edd.big$USDAcode[edd.big$USDAcode == 1])
  
  ed.listing.i <- data.frame(scientific.name, usda.code, comment, state, federal, other, ipaus, hotspot, no.abun, no.abs, no.big, L48, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)
  names(ed.listing.i) <- c("scientific.name", "usda.code", "comment", "state", "federal", "other", "ipaus", "hotspot", "no.abun", "no.abs", "no.big", "nativity", "usdaCHK") ## "no.abun", "no.abs", "usdaCHK")  
  ed.listing <- rbind(ed.listing, ed.listing.i)
  
  print(i)
  
}

ed.listing <- ed.listing[ed.listing$usda.code != "TEMPLATE", ]
ed.listing <- ed.listing[!is.na(ed.listing$usda.code),]
unique(ed.listing$usda.code)
summary(ed.listing)
head(ed.listing)






## ensure that synonyms are checked in ipa.list

no.ipaus <- ed.listing$scientific.name[ed.listing$ipaus == 0]

for (i in 1:length(no.ipaus)){
  list.i <- usda.list$Scientific.Name[usda.list$Accepted.Symbol == ed.listing$usda.code[ed.listing$scientific.name == no.ipaus[i]]]
  a <- 0
  
  for (j in 1:length(list.i)){
    a.j  <- length(ipa.list$scientific_name[ipa.list$scientific_name == list.i[j]])
    a <- a + a.j
    }
  
  if (a > 0) {
    ed.listing$ipaus2[ed.listing$scientific.name == no.ipaus[i]] <- 1
  }
  print(i)
}


## two spp that are on ipaus were not in ipaus file, so i'm entering their ipaus status manually
ed.listing$ipaus[!is.na(ed.listing$ipaus2) & ed.listing$ipaus2 == 1] <- 1
ed.listing[ed.listing$state == 0 & ed.listing$federal == 0 & ed.listing$hotspot == 1 & ed.listing$other == 0 & ed.listing$ipaus == 0, ]
ed.listing$ipaus[ed.listing$scientific.name == "Centaurea x moncktonii"] <- 1


#######################################################################
#######################      FLAG AQUATICS     ########################
#######################################################################

aqua <- read.table("Products/aquatics.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")

for (i in 1:length(aqua$usda.code)){
  ed.listing$Aquatic[ed.listing$usda.code == aqua$usda.code[i]] <- 1
  print(i)
}

summary(edd.list$Aquatic)



