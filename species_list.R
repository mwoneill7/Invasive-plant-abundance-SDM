#### Combining eddmaps species list with w/ invasive listings
###  Mitch O'Neill
##   created: 6/13/2017; 
#    last modified 8/9/2017

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files")

#### eddmaps species list (all species with infested area downloaded from EDDMAPS October 2016)
edd.list <- read.table("Raw/speciesList06_01_2017.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
edd.list$edd <- 1## to mark that this species was on the edd.list (for upcoming merges)


#### list of species from the occurence-based hotspots analysis (Allen and Bradley 2016)
occ.list <- read.table("C:/Users/mwone/Google Drive/NSF_GSS_shared/Hotspots_and_Abundance/occ_data_with_2017_codes.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
occ.list$occ <- 1 ## to mark that this species was on the occurence-based hotspots analysis (for upcoming merges)
head(occ.list)

#### USDA federal/State listings (Downloaded from PLANTS 5/30/2017)
usda.list <- read.table("Raw/PLANTS.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
usda.list$usdaCHK <- 1
head(usda.list) ## to mark that this species was in the PLANTS database (for upcoming merges)

#### nativity is only listed on first species name for a species code, AND
#### nativity for each place (canada, alaska, etc) is lumped into one messy field, SO
#### extract all rows with each of the three possibilities (I, N, I/N)
plants.i  <- usda.list[grep("L48(I)",   usda.list$Native.Status, fixed=TRUE), ]
plants.n  <- usda.list[grep("L48(N)",   usda.list$Native.Status, fixed=TRUE), ]
plants.in <- usda.list[grep("L48(I,N)", usda.list$Native.Status, fixed=TRUE), ]

### make a cleaner nativity column for L48 based on what their L48 status is
plants.i$L48  <- "I"
plants.n$L48  <- "N"
plants.in$L48 <- "I/N"

### bind rows back into a single dataframe
plants.nativity <- rbind(plants.i, plants.n, plants.in)
head(plants.nativity)

#### IPAUS species list (Downloaded in 2014; database last updated May 2012)
ipa.list <- read.table("Raw/ipa2.csv", header = T, sep = ",", 
                       stringsAsFactors = F, quote= "\"", comment.char= "", strip.white = T)
ipa.list$IPA <- 1 ## marks that the species appeared on IPA for (upcoming merge) 
head(ipa.list)


## prepped eddmaps database (see eddmaps_prep.R) 
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_08_09_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
ed.list <- as.list(unique(edd$USDAcode)) ## list of all usda-codes in prepped file



## create template dataframe for a master file with info for each species
usda.code <- "TEMPLATE" ## for usda code
scientific.name <- "TEMPLATE" ## for accepted scientific name (plants)
comment <- "TEMPLATE" ## pulls any comments from the eddmaps-usdacode-conversion file
state   <- -99 ## state status?
federal <- -99 ## federal status
other   <- -99 ## listed by another source? (Weeds of the U.S.; PLANTS)
ipaus   <- -99 ## listed on ipaus?
hotspot <- -99 ## was this plant in the occurence hotspot analysis (Allen and Bradley 2016)
nativity<- "TEMPLATE" ## nativity of plant in L48
usdaCHK <- -99 ## doublt check that usda code is in PLANTS database
no.abs  <- -99 ## number of absence points
no.abun <- -99 ## number of points with valid infested area
no.big <- -99  ## number of those observations with unusually high infested area
min.abun <- -99 ## minimum non-zero infested area value
max.abun <- -99 ## maximum infested area value, excluding unusually high ones

ed.listing <- data.frame(scientific.name, usda.code, comment, state, federal, other, ipaus, hotspot, no.abun, no.abs, no.big, min.abun, max.abun, nativity, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)
## make template object

# loops through list of all species with usda codes in prepped database
for (i in 1:length(ed.list)){
  
  scientific.name <- edd.list$PLANTS_species[edd.list$USDA_code == ed.list[i]][1]
  ## pulls accepted scientific name (PLANTS) for usdacode
  
  comment <- edd.list$Comment[edd.list$USDA_code == ed.list[i]][1]
  ## pulls comments from eddmaps-usda conversion file
  
  usda.code <- ed.list[i] ## pulls usda code of iteration
  
  usda.i <- usda.list[usda.list$Accepted.Symbol == ed.list[i],] ## subsets usda database for all synonyms of usda code
  usda.i <- usda.i[!is.na(usda.i$Accepted.Symbol), ] ## removes nonsence rows
  
  
  ## if any of the synonyms (subset of usda database) have state listing, then state = 1    
  if (length(usda.i$Accepted.Symbol[usda.i$State.Noxious.Status != ""]) > 0) {
    state <- 1} else {state <- 0}
 
  ## if any of the synonyms (subset of usda database) have federal listing, then federal = 1 
  if (length(usda.i$Accepted.Symbol[usda.i$Federal.Noxious.Status == "NW"]) > 0) {
    federal <- 1} else  {federal <- 0}
  
  ## if any of the synonyms (subset of usda database) have any other listing (invasive), then other = 1 
  if (length(usda.i$Accepted.Symbol[usda.i$Invasive != ""]) > 0) {
    other <- 1}  else {other <- 0}
  
  ## if species appears on ipa, then ipaus=1
  if (length(ipa.list$habit[ipa.list$USDA_code_2017 == ed.list[i]]) > 0) {
    ipaus <- 1}  else {ipaus <- 0}

  ## if species appears on hotspots list (2016) then hotspot = 1
  if (length(occ.list$habit[occ.list$NewCode == ed.list[i]]) > 0){
    hotspot <- 1} else {hotspot <- 0}

  ## if the species is in the nativity object (i.e. if usda provided an L48 nativity)
  if(length(plants.nativity$L48[plants.nativity$Accepted.Symbol == ed.list[i]]) > 0) {
    L48 <- plants.nativity$L48[plants.nativity$Accepted.Symbol == ed.list[i]] 
  } else {L48 = "N/A"} ## if species is not in nativity object, then usda did not provide nativity
  
  if(length(usda.i$Accepted.Symbol[usda.i$usdaCHK == 1]) > 0) {
    usdaCHK <- 1} else {usdaCHK <- 0} ## verify that species was in PLANTS
  
  ## subset all records for species of iteration with valid infested area
  edd.abun <- edd[edd$USDAcode == ed.list[i] & edd$validInfestedArea == 1, ]
  no.abun <- length(edd.abun$USDAcode[!is.na(edd.abun$USDAcode)])
  ## number of records, exluding NA rows
  
  ## subset all records for species of iteration with negative=1
  edd.abs <- edd[edd$USDAcode == ed.list[i] & edd$negative == 1, ]
  no.abs <- length(edd.abs$USDAcode[!is.na(edd.abs$USDAcode)])
  ## number of records, exluding NA rows
  
  ## subset all records that have valid infested area,
  no.big <- length(edd.abun$USDAcode[!is.na(edd.abun$USDAcode) & (edd.abun$infestedAreaInAcres > 20 | edd.abun$ReporterFULLName == "US Army Corps of Engineers Ombil Database")]) 
  ## number of records thats are unusually large (>20) or reported by the large-value reporter, exluding NA rows
  
  min.abun <- min(edd.abun$infestedAreaInAcres[!is.na(edd.abun$USDAcode)]) ## minimum non-zero abundance
  max.abun <- max(edd.abun$infestedAreaInAcres[!is.na(edd.abun$USDAcode) & edd.abun$infestedAreaInAcres < 20 & edd.abun$ReporterFULLName != "US Army Corps of Engineers Ombil Database"])
  ## maximum infested area, excluding unusually large values and excluding large-value reporter
  
  ## compiles all of this info into dataframe, matching format of template
  ed.listing.i <- data.frame(scientific.name, usda.code, comment, state, federal, other, ipaus, hotspot, no.abun, no.abs, no.big, min.abun, max.abun, L48, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)
  names(ed.listing.i) <- c("scientific.name", "usda.code", "comment", "state", "federal", "other", "ipaus", "hotspot", "no.abun", "no.abs", "no.big", "min.abun", "max.abun", "nativity", "usdaCHK") ## "no.abun", "no.abs", "usdaCHK")  
  ed.listing <- rbind(ed.listing, ed.listing.i) ## concatenate row from iteration to template (and rows from past iterations)
  
  print(i) ## keeps track of place
  
}


head(usda.list)

ed.listing <- ed.listing[ed.listing$usda.code != "TEMPLATE", ] ## removes template rows
ed.listing <- ed.listing[!is.na(ed.listing$usda.code),] ## removes NA rows
unique(ed.listing$usda.code)
summary(ed.listing)
head(ed.listing)

## ensure that synonyms are checked in ipa.list

## extract all species in eddmaps that were not on ipaus, to check synonyms
no.ipaus <- ed.listing$scientific.name[ed.listing$ipaus == 0]

for (i in 1:length(no.ipaus)){ ## loop through those species
  ## extract list of synonyms
  list.i <- usda.list$Scientific.Name[usda.list$Accepted.Symbol == ed.listing$usda.code[ed.listing$scientific.name == no.ipaus[i]]]
  a <- 0 ## this will tally all of the synonyms on the IPA; reset to zero for each new set of synonyms
  
  for (j in 1:length(list.i)){ ## loop through synonyms
    a.j  <- length(ipa.list$scientific_name[ipa.list$scientific_name == list.i[j]]) 
    ## if the synonym is in ipaus, then a.j will equal 1, otherwise it will equal zero
    a <- a + a.j ## adds whether or not synonym was in IPA (yes=1, no=0)
    }
  
  if (a > 0) { ## if any of the synonyms were on ipaus
    ed.listing$ipaus[ed.listing$scientific.name == no.ipaus[i]] <- 1
  } ## then reassign ipaus to a 1
  
  print(i) ## keep track of place
}


#######################################################################
#######################      FLAG AQUATICS     ########################
#######################################################################

aqua <- read.table("Products/aquatics.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
## list of aquatic species that haven't been eliminated yet

ed.listing$Aquatic <- 0

for (i in 1:length(aqua$usda.code)){ ## loop through aquatics
  ed.listing$Aquatic[ed.listing$usda.code == aqua$usda.code[i]] <- 1 ## flag any matches
  print(i) ## keep place
}

summary(as.factor(ed.listing$Aquatic))

## flag species I know I can't use
ed.listing$potential.use <- 1 ## default is that I can use them, then remove them based on certain criteria

ed.listing$potential.use[ed.listing$Aquatic == 1 | ## can't use aquatics
                        (ed.listing$state + ed.listing$ipaus + ed.listing$federal + ed.listing$other) < 1 |
                          ## can't use species that haven't been recognized as problematic
                        (ed.listing$no.abs + ed.listing$no.abun < 20) | ed.listing$no.abun < 15 |
                          ## data deficient
                        (ed.listing$nativity == "N") ] <- 0  
                          ## can't use native plants

### flag records where there is enough data right now, but if large values are removed, there may not be.
ed.listing$not.enough.reasonable.values <- 0
ed.listing$not.enough.reasonable.values[(ed.listing$no.abs + ed.listing$no.abun - ed.listing$no.big < 20) | 
                                        (ed.listing$no.abun - ed.listing$no.big < 15)] <- 1

write.csv(ed.listing, "Products/SpeciesList08_10_2017.csv", row.names = F)
#### output file
#### manual work to follow; i.e.
#### where nativity = N/A manually assign nativity using MoBot, etc.