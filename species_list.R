#### Combining eddmaps species list with w/ invasive listings
###  Mitch O'Neill
##   created: 6/13/2017; 
#    last modified 8/9/2017

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files")

#### eddmaps species list (all species with infested area downloaded from EDDMAPS October 2016)
edd.list <- read.table("SpeciesList_11_20_2017.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
## practice with ssp

edd.list$edd <- 1## to mark that this species was on the edd.list (for upcoming merges)
head(edd.list)

#### list of species from the occurence-based hotspots analysis (Allen and Bradley 2016)
occ.list <- read.table("C:/Users/mwone/Google Drive/NSF_GSS_shared/Hotspots_and_Abundance/occ_data_with_2017_codes.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
occ.list$occ <- 1 ## to mark that this species was on the occurence-based hotspots analysis (for upcoming merges)
head(occ.list)

#### USDA federal/State listings (Downloaded from PLANTS 5/30/2017)
usda.list <- read.table("PLANTS.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
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
ipa.list <- read.table("ipa2.csv", header = T, sep = ",", 
                       stringsAsFactors = F, quote= "\"", comment.char= "", strip.white = T)
ipa.list$IPA <- 1 ## marks that the species appeared on IPA for (upcoming merge) 
head(ipa.list)


## prepped eddmaps database (see eddmaps_prep.R) 
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_prepped_12_18_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
ed.list <- as.list(unique(edd$USDAcode)) ## list of all usda-codes in prepped file

## create template dataframe for a master file with info for each species
usda.code <- "TEMPLATE" ## for usda code
scientific.name <- "TEMPLATE" ## for accepted scientific name (plants)
#comment <- "TEMPLATE" ## pulls any comments from the eddmaps-usdacode-conversion file
state   <- -99 ## state status?
federal <- -99 ## federal status
other   <- -99 ## listed by another source? (Weeds of the U.S.; PLANTS)
ipaus   <- -99 ## listed on ipaus?
hotspot <- -99 ## was this plant in the occurence hotspot analysis (Allen and Bradley 2016)
nativity<- "TEMPLATE" ## nativity of plant in L48
usdaCHK <- -99 ## doublt check that usda code is in PLANTS database
#no.abs  <- -99 ## number of absence points
no.abun <- -99 ## number of points with valid infested area
#no.big <- -99  ## number of those observations with unusually high infested area
#min.abun <- -99 ## minimum non-zero infested area value
#max.abun <- -99 ## maximum infested area value, excluding unusually high ones
#L48 <- "TEMPLATE"

ed.listing <- data.frame(scientific.name, usda.code, state, federal, other, ipaus, hotspot, no.abun, nativity, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)
## make template object

# loops through list of all species with usda codes in prepped database
for (i in 1:length(ed.list)){
  
  #scientific.name <- edd.list$scientific.name[edd.list$usda.code == ed.list[i]][1]
  scientific.name <- edd.list$PLANTS_species[edd.list$USDA_code == ed.list[i]][1]
  ## pulls accepted scientific name (PLANTS) for usdacode

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

  no.abun <- length(edd$USDAcode[edd$USDAcode == ed.list[i]])  
  #################
  # ## subset all records for species of iteration with valid infested area
  # edd.abun <- edd[edd$USDAcode == ed.list[i] & edd$validInfestedArea == 1, ]
  # no.abun <- length(edd.abun$USDAcode[!is.na(edd.abun$USDAcode)])
  # ## number of records, exluding NA rows
  #
  # ## subset all records for species of iteration with negative=1
  # edd.abs <- edd[edd$USDAcode == ed.list[i] & edd$negative == 1, ]
  # no.abs <- length(edd.abs$USDAcode[!is.na(edd.abs$USDAcode)])
  # ## number of records, exluding NA rows
  #
  # # subset all records that have valid infested area,
  # no.big <- length(edd.abun$USDAcode[!is.na(edd.abun$USDAcode) & (edd.abun$infestedAreaInAcres > 20 | edd.abun$ReporterFULLName == "US Army Corps of Engineers Ombil Database")]) 
  # ## number of records thats are unusually large (>20) or reported by the large-value reporter, exluding NA rows
  #
  # min.abun <- min(edd.abun$infestedAreaInAcres[!is.na(edd.abun$USDAcode)]) ## minimum non-zero abundance
  # max.abun <- max(edd.abun$infestedAreaInAcres[!is.na(edd.abun$USDAcode) & edd.abun$infestedAreaInAcres < 20 & edd.abun$ReporterFULLName != "US Army Corps of Engineers Ombil Database"])
  # ## maximum infested area, excluding unusually large values and excluding large-value reporter
  #################
  
  ## compiles all of this info into dataframe, matching format of template
  ed.listing.i <- data.frame(scientific.name, usda.code, state, federal, other, ipaus, hotspot, no.abun, L48, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)
  names(ed.listing.i) <- c("scientific.name", "usda.code", "state", "federal", "other", "ipaus", "hotspot", "no.abun", "nativity", "usdaCHK") ## "no.abun", "no.abs", "usdaCHK")  
  ed.listing <- rbind(ed.listing, ed.listing.i) ## concatenate row from iteration to template (and rows from past iterations)
  
  print(i) ## keeps track of place
  
}


ed.listing <- ed.listing[ed.listing$usda.code != "TEMPLATE", ] ## removes template rows
#ed.listing <- ed.listing[!is.na(ed.listing$usda.code),] ## removes NA rows
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

summary(ed.listing$no.abun)
ed.listing <- ed.listing[ed.listing$no.abun >= 40,]

# #######################################################################
# #######################      FLAG AQUATICS     ########################
# #######################################################################
# 
#  aqua <- read.table("aquatics.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
#  ## list of aquatic species that haven't been eliminated yet
#  
#  ed.listing$Aquatic <- 0
#  
#  for (i in 1:length(aqua$usda.code)){ ## loop through aquatics
#    ed.listing$Aquatic[ed.listing$usda.code == aqua$usda.code[i]] <- 1 ## flag any matches
#    print(i) ## keep place
#  }
#  
#  summary(as.factor(ed.listing$Aquatic))
#  ed.listing <- ed.listing[ed.listing$Aquatic == 0,]
#  #ed.listing <- ed.listing[ed.listing$no.abun >= 20,] ## 939 to 407
# 
## flag species I know I can't use
#ed.listing$potential.use <- 1 ## default is that I can use them, then remove them based on certain criteria
#
#ed.listing$potential.use[ed.listing$Aquatic == 1 | ## can't use aquatics
#                        (ed.listing$state + ed.listing$ipaus + ed.listing$federal + ed.listing$other) < 1 |
#                          ## can't use species that haven't been recognized as problematic
#                        (ed.listing$no.abs + ed.listing$no.abun < 20) | ed.listing$no.abun < 15 |
#                          ## data deficient
#                        (ed.listing$nativity == "N") ] <- 0  
#                          ## can't use native plants
# 
# ### flag records where there is enough data right now, but if large values are removed, there may not be.
# ed.listing$data <- 0
# ed.listing$not.enough.reasonable.values[(ed.listing$no.abs + ed.listing$no.abun - ed.listing$no.big < 20) | 
#                                         (ed.listing$no.abun - ed.listing$no.big < 15)] <- 1

########################
 
#write.csv(ed.listing, "Species_List_11_21_2017.csv", row.names = F)
#### output file
#### manual work to follow; i.e.
#### where nativity = N/A manually assign nativity using MoBot, etc.

#ed.listing<- read.table("Species_List_11_21_2017.csv",header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
invListing <- read.table("master list.csv",header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "") 
head(invListing)
invListing$Number.of.States[is.na(invListing$Number.of.States)]<-0
invListing$Federal.Noxious[is.na(invListing$Federal.Noxious)]<-0

spp <- ed.listing$usda.code

for (i in 1:length(spp)){
  if(length(invListing$NewCode[invListing$NewCode==spp[i]])>=1){
    ed.listing$stateNEW[ed.listing$usda.code==spp[i]] <- max(invListing$Number.of.States[invListing$NewCode == spp[i]])
    ed.listing$FedNEW[ed.listing$usda.code==spp[i]] <- max(invListing$Federal.Noxious[invListing$NewCode == spp[i]])
    
  }
print(i)
}

nativity2 <- read.table("Nativity2.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
for (i in (1:length(nativity2$usda.code))){
  ed.listing$nativity[ed.listing$usda.code == nativity2$usda.code[i]] <- nativity2$Nativity[nativity2$usda.code == nativity2$usda.code[i]]
}

head(ed.listing)
summary(as.factor(ed.listing$nativity))
ed.listing <- ed.listing[ed.listing$nativity == "I",]
summary(ed.listing$stateNEW)
summary(ed.listing$no.abun)

#summary(ed.listing$FedNEW)
#
#probs <- ed.listing[ed.listing$state == 1 & is.na(ed.listing$stateNEW),]
#probsp <- probs$usda.code
#
#write.csv(probs, "C://Users/mwone/Desktop/probs.csv", row.names=F)
#
#probs <- ed.listing[ed.listing$state == 1 & ed.listing$stateNEW == 0,]
#probs <- probs[!is.na(probs$usda.code),]
#write.csv(probs, "C://Users/mwone/Desktop/probs2.csv", row.names=F)
#
#probs <- ed.listing[ed.listing$federal == 1 & is.na(ed.listing$FedNEW),]
#probs <- probs[!is.na(probs$usda.code),]
#write.csv(probs, "C://Users/mwone/Desktop/probs2.csv", row.names=F)
#
#probs <- ed.listing[ed.listing$federal == 1 & ed.listing$FedNEW==0,]
#probs <- probs[!is.na(probs$usda.code),]
#write.csv(probs, "C://Users/mwone/Desktop/probs2.csv", row.names=F)

#probs <- ed.listing[ed.listing$state == 1 & ed.listing$stateNEW == 0,]

#probs <- ed.listing[ed.listing$state == 1 & is.na(ed.listing$stateNEW),]
#max(c(NA,7)) ## LOOK INTO LATER
#ed.listing2 <- ed.listing[(!is.na(ed.listing$stateNEW) & ed.listing$stateNEW >1) | (ed.listing$FedNEW ==1 & !is.na(ed.listing$FedNEW)) |
#                            ed.listing$other==1 | ed.listing$ipaus==1 ,]
#ed.listing2 <- ed.listing2[ed.listing2$no.abun > 20 & ed.listing2$nativity != "N",]
#head(ed.listing2)

ed.listing$useable <- 0
ed.listing$useable[ed.listing$other==1] <- .75
ed.listing$useable[ed.listing$ipaus==1 | (!is.na(ed.listing$FedNEW) & ed.listing$FedNEW == 1) | 
                     (!is.na(ed.listing$stateNEW) & ed.listing$stateNEW > 0)] <- 1
#ed.listing$useable[ed.listing$other==1 & ed.listing$nativity == "I/N"] <- 0.25
#ed.listing$useable[ed.listing$nativity == "I/N" & (ed.listing$ipaus==1 | ed.listing$state==1 | ed.listing$federal ==1)] <- 0.5
summary(as.factor(ed.listing$useable))
#summary(as.factor(edd$ORD[edd$USDAcode == "URDI" ]))
ed.listing[ed.listing$useable == .75,]
#ed.listing$useable[ed.listing$no.abun < 20] <- 0
#ed.listing[ed.listing$useable == .25,]
#unique(ed.listing$ipaus)
#unique(ed.listing$stateNEW)
#summary(ed.listing2$no.abun)
#hist(ed.listing2$no.abun,breaks=100)

write.csv(ed.listing,"species_summary_12_12_2017.csv",row.names=F)
ed.listing <- read.table("species_summary_12_12_2017.csv", header = T, sep = ",", quote= "\"", 
                   comment.char= "", stringsAsFactors = F, strip.white = T)


## thinned, number in each category
edd2 <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_12_18_2017.csv", header = T, sep = ",", quote= "\"", 
                   comment.char= "", stringsAsFactors = F, strip.white = T)
edd2$abundance <- ceiling(edd2$best)


sp <- ed.listing$usda.code[ed.listing$useable>0]

for (i in 1:length(sp)){
  edd.sp <- edd2[edd2$species == sp[i],]
  ed.listing$no.1[ed.listing$usda.code == sp[i]] <- length(edd.sp$species[edd.sp$abundance == 1])
  ed.listing$no.2[ed.listing$usda.code == sp[i]] <- length(edd.sp$species[edd.sp$abundance == 2])
  ed.listing$no.3[ed.listing$usda.code == sp[i]] <- length(edd.sp$species[edd.sp$abundance == 3])
  ed.listing$no.4[ed.listing$usda.code == sp[i]] <- length(edd.sp$species[edd.sp$abundance == 4])
  #ed.listing$no.01[ed.listing$usda.code == sp[i]]<- length(edd.sp$species[edd.sp$abundance == 1 & edd.sp$cov01 == 0]) 
  
  if(length(unique(edd.sp$abundance))==4){
  ed.listing$rarest.bin[ed.listing$usda.code == sp[i]] <- min(summary(as.factor(edd.sp$abundance)))
  #ed.listing$second.bin[ed.listing$usda.code == sp[i]] <- min(summary(as.factor(edd.sp$abundance)))
  } else {
  ed.listing$rarest.binc[ed.listing$usda.code == sp[i]] <- 0
  }
  
  
  print(i)
}


length(ed.listing$usda.code[ed.listing$rarest.bin >= 10 & !is.na(ed.listing$rarest.bin)])
## 116 species (median of best)

species.list <- ed.listing[ed.listing$rarest.bin >= 10 & !is.na(ed.listing$rarest.bin),]
write.csv(ed.listing, "species_list_12_18_17.csv",row.names=F)

edd3 <- edd2[edd2$species %in% species.list$usda.code,]
write.csv(edd3, "edd_final_12_18_17.csv",row.names=F)




##############################################
lost <- ed.listing[ed.listing$rarest.bin <10,]
edd3 <- edd2[edd2$species %in% species.list$usda.code,]
write.csv(edd3, "edd_final_12_15_17.csv",row.names=F)
write.csv(ed.listing, "species_list_12_15_17.csv",row.names=F)
unique(edd3$species)
head(edd3)


coordinates(edd3) <- c(5,4)
plot(edd3, pch=".")

#sppF <- species.list$usda.code
probs <- probs[probs$species %in% sppF,]
unique(probs$species)

##############################################
cov01 <- ed.listing[ed.listing$no.01 <10 & ed.listing$no.1>=10 & ed.listing$no.2>=10 & ed.listing$no.3>=10 & ed.listing$no.4>=10,]
cov01[!is.na(cov01$scientific.name),]
## cov01 saves one species (VEIN9)
#summary(ed.listing)

nrow(ed.listing[ed.listing$rarest.bin < 10 & !is.na(ed.listing$rarest.bin),])
nrow(ed.listing[ed.listing$rarest.bin >= 5 & !is.na(ed.listing$rarest.bin),])
 species.list <- ed.listing[ed.listing$rarest.bin >= 10 & !is.na(ed.listing$rarest.bin),]
species.list$n <- species.list$no.01 + species.list$no.2 + species.list$no.3 + species.list$no.4
hist(species.list$n, breaks=50)
hist(species.list$rarest.bin, breaks=50)
summary(species.list$n)


length(ed.listing$rarest.bin[ed.listing$rarest.bin >= 5 & ed.listing$useable > 0])  ## 272 (270)
length(ed.listing$rarest.bin[ed.listing$rarest.bin >= 10]) ## 217 (213)




##############################################














ed.listing[ed.listing$rarest.bin >= 5 & ed.listing$useable < 1 & ed.listing$useable > 0,]
## same as below, but also urtica dioica & taraxacum officionale
ed.listing[ed.listing$rarest.bin >= 10 & ed.listing$useable < 1 & ed.listing$useable > 0,]
## Bromus inermis and Oplismenus hirtellus
 
## none rely on the "other" column
#ed.listing[ed.listing$other==1 & ed.listing$ipaus==0 & ed.listing$FedNEW == 0 & ed.listing$stateNEW ==0,]
hist(ed.listing$no.abun[ed.listing$rarest.bin >= 5 & ed.listing$useable >.5],breaks=100) #& ed.listing$no.abun < 120
range(ed.listing$no.abun[ed.listing$rarest.bin >= 5 & ed.listing$useable >.5],na.rm=T) ## 20 42080
quantile(ed.listing$no.abun[ed.listing$rarest.bin >= 5 & ed.listing$useable >.5],na.rm=T)
hist(ed.listing$no.abun[ed.listing$rarest.bin < 5 & ed.listing$useable >.5],breaks=100)
ed.listing[ed.listing$no.abun > 500 & ed.listing$rarest.bin <5,]
length(ed.listing$no.abun[ed.listing$no.abun > 100 & ed.listing$rarest.bin <5])


length(ed.listing$rarest.bin[ed.listing$rarest.bin == ed.listing$no.1])/383 # 64%
length(ed.listing$rarest.bin[ed.listing$rarest.bin == ed.listing$no.2])/383 # 12%
length(ed.listing$rarest.bin[ed.listing$rarest.bin == ed.listing$no.3])/383 # 14%
length(ed.listing$rarest.bin[ed.listing$rarest.bin == ed.listing$no.4])/383 # 30%

#ed.listing[ed.listing$usda.code == "PUMOL",]
#ed.listing$no.1<-NULL
#ed.listing$no.2<-NULL
#ed.listing$no.3<-NULL
#ed.listing$no.4<-NULL


ed <- ed.listing[ed.listing$useable == 1 | ed.listing$usda.code == "OPHI",]

length(ed$usda.code[ed$rarest.bin >=5 & ed$no.1 >0 & ed$no.2 >0 & ed$no.3>0 & ed$no.4 >0])
length(ed$usda.code[ed$rarest.bin >=5 & ed$no.1 >0 & ed$no.2 >0 & ed$no.3>0 & ed$no.4 >0 & ed$no.abun >= 50])
length(ed$usda.code[ed$rarest.bin >=10 & ed$no.1 >0 & ed$no.2 >0 & ed$no.3>0 & ed$no.4 >0])
length(ed$usda.code[ed$rarest.bin >=10 & ed$no.abun >= 50  & ed$no.1 >0 & ed$no.2 >0 & ed$no.3>0 & ed$no.4 >0])

#ed2 <- ed[ed$rarest.bin < 5,]

#ed2 <- data.frame(cbind(ed2$usda.code, ed2$no.1, ed2$no.2, ed2$no.3, ed2$no.4), stringsAsFactors = F)
#colnames(ed2)<- c("sp","one","two","three","four")
#write.csv(ed2, "C://Users/mwone/Desktop/ed2.csv", row.names=F)

ed.listing$one5[ed.listing$no.1 <5] <- 0
ed.listing$two5[ed.listing$no.2 <5] <- 0
ed.listing$three5[ed.listing$no.3 <5] <- 0
ed.listing$four5[ed.listing$no.4 <5] <- 0

ed.listing$one5[ed.listing$no.1 >= 5] <- 1
ed.listing$two5[ed.listing$no.2 >= 5] <- 1
ed.listing$three5[ed.listing$no.3 >= 5] <- 1
ed.listing$four5[ed.listing$no.4 >= 5] <- 1

ed.listing$cats5 <- ed.listing$one5 +ed.listing$two5 +ed.listing$three5 + ed.listing$four5
summary(as.factor(ed.listing$cats5[ed.listing$useable == 1]))
ed.listing[ed.listing$cats5==3 & ed.listing$useable == 1,]

ed.listing$one10[ed.listing$no.1 <10] <- 0
ed.listing$two10[ed.listing$no.2 <10] <- 0
ed.listing$three10[ed.listing$no.3 <10] <- 0
ed.listing$four10[ed.listing$no.4 <10] <- 0

ed.listing$one10[ed.listing$no.1 >= 10] <- 1
ed.listing$two10[ed.listing$no.2 >= 10] <- 1
ed.listing$three10[ed.listing$no.3 >= 10] <- 1
ed.listing$four10[ed.listing$no.4 >= 10] <- 1

ed.listing$cats10 <- ed.listing$one10 +ed.listing$two10 +ed.listing$three10 + ed.listing$four10
summary(as.factor(ed.listing$cats10[ed.listing$useable == 1]))

summary(as.factor(ed.listing$cats5[ed.listing$useable == 1 & ed.listing$no.abun >= 50]))
summary(as.factor(ed.listing$cats10[ed.listing$useable == 1 & ed.listing$no.abun >= 50]))

summary(as.factor(ed.listing$cats5[ed.listing$useable == 1 & ed.listing$no.1 == ed.listing$rarest.bin]))
summary(as.factor(ed.listing$cats10[ed.listing$useable == 1 & ed.listing$no.1 == ed.listing$rarest.bin]))

summary(as.factor(ed.listing$cats5[ed.listing$useable == 1 & ed.listing$no.abun >= 50& ed.listing$no.1 == ed.listing$rarest.bin]))
summary(as.factor(ed.listing$cats10[ed.listing$useable == 1 & ed.listing$no.abun >= 50& ed.listing$no.1 == ed.listing$rarest.bin]))

#ed2$one[ed2$one >= 5] <- 1
#ed2$two[ed2$two >= 5] <- 1
#ed2$three[ed2$three >= 5] <- 1
#ed2$four[ed2$four >= 5] <- 1

ed.listing$no.1b <- ed.listing$no.1 + ed.listing$no.2

length(ed.listing$usda.code[ed.listing$no.1b >=5 & ed.listing$no.3 >= 5 & ed.listing$no.4 >=5 & ed.listing$useable==1 & ed.listing$no.abun >= 50])
length(ed.listing$usda.code[ed.listing$no.1b >=10 & ed.listing$no.3 >= 10 & ed.listing$no.4 >= 10 & ed.listing$useable==1 & ed.listing$no.abun >= 50])

ed.listing[ed.listing$scientific.name=="Arctotheca calendula",]

hist(as.numeric(edd$abundance), main = "After thinning", xlab = "abundance rank")







############################################################################
edd2 <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_08_28_2017.csv", header = T, sep = ",", quote= "\"", 
                   comment.char= "", stringsAsFactors = F, strip.white = T)

ed.listing <- read.table("SpeciesList08_10_2017.csv", header = T, sep = ",", quote= "\"", 
                         comment.char= "", stringsAsFactors = F, strip.white = T)

for (i in 1:length(ed.listing$usda.code)){
  
  ## subset all records for species of iteration with valid infested area
  edd.abun <- edd2[edd2$usda == ed.listing$usda.code[i] & edd2$abundance > 0, ]
  ed.listing$no.abunG[ed.listing$usda.code == ed.listing$usda.code[i]] <- 
    length(edd.abun$usda[!is.na(edd.abun$usda)])
  ## number of records, exluding NA rows
  
  ## subset all records for species of iteration with negative=1
  edd.abs <- edd2[edd2$usda == ed.listing$usda.code[i] & edd2$absence == 1, ]
  ed.listing$no.absG[ed.listing$usda.code == ed.listing$usda.code[i]] <- 
    length(edd.abs$usda[!is.na(edd.abs$usda)])
  ## number of records, exluding NA rows
  print(i)
}

### new potential use column using number of grrid cells rather than number of points
ed.listing$potential.use2 <- 1 ## default is that I can use them, then remove them based on certain criteria

ed.listing$potential.use2[ed.listing$Aquatic == 1 | ## can't use aquatics
                            (ed.listing$state + ed.listing$ipaus + ed.listing$federal + ed.listing$other) < 1 |
                            ## can't use species that haven't been recognized as problematic
                            (ed.listing$no.absG + ed.listing$no.abunG < 20) | ed.listing$no.abunG < 15 |
                            ## data deficient
                            (ed.listing$nativity == "N") ] <- 0  
## can't use native plants


write.csv(ed.listing, "SpeciesList08_28_2017.csv", row.names=F)
## write out for later use


############################################################################
########## add the number of unique grid cells to the master file ##########
##################### USING MAX THRESHOLD OF 10 ACRES ######################
############################################################################

edd2 <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_10acre.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

ed.listing <- read.table("SpeciesList09_26_2017.csv", header = T, sep = ",", quote= "\"", 
                        comment.char= "", stringsAsFactors = F, strip.white = T)

for (i in 1:length(ed.listing$usda.code)){
  
## subset all records for species of iteration with valid infested area
edd.abun <- edd2[edd2$usda == ed.listing$usda.code[i] & edd2$abundance > 0, ]
ed.listing$no.abunG10[ed.listing$usda.code == ed.listing$usda.code[i]] <- 
  length(edd.abun$usda[!is.na(edd.abun$usda)])
## number of records, exluding NA rows

print(i)
}

### new potential use column using number of grrid cells rather than number of points
ed.listing$potential.use_10ac <- 1 ## default is that I can use them, then remove them based on certain criteria

ed.listing$potential.use_10ac[ed.listing$Aquatic == 1 | ## can't use aquatics
                           (ed.listing$state + ed.listing$ipaus + ed.listing$federal + ed.listing$other) < 1 |
                           ## can't use species that haven't been recognized as problematic
                           (ed.listing$no.abunG10 < 20) |
                           ## data deficient
                           (ed.listing$nativity == "N") ] <- 0  
                           ## can't use native plants


write.csv(ed.listing, "SpeciesList10_11_2017.csv", row.names=F)
## write out for later use

############################################################################
########## add the number of unique grid cells to the master file ##########
##################### USING MAX THRESHOLD OF 02 ACRES ######################
############################################################################

edd2 <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_2acre.csv", header = T, sep = ",", quote= "\"", 
                   comment.char= "", stringsAsFactors = F, strip.white = T)

ed.listing <- read.table("SpeciesList10_11_2017.csv", header = T, sep = ",", quote= "\"", 
                         comment.char= "", stringsAsFactors = F, strip.white = T)

for (i in 1:length(ed.listing$usda.code)){
  
  ## subset all records for species of iteration with valid infested area
  edd.abun <- edd2[edd2$usda == ed.listing$usda.code[i] & edd2$abundance > 0, ]
  ed.listing$no.abunG2[ed.listing$usda.code == ed.listing$usda.code[i]] <- 
    length(edd.abun$usda[!is.na(edd.abun$usda)])
  ## number of records, exluding NA rows
  
  print(i)
}
summary(ed.listing$no.abunG2)
### new potential use column using number of grrid cells rather than number of points
ed.listing$potential.use_2ac <- 1 ## default is that I can use them, then remove them based on certain criteria

ed.listing$potential.use_2ac[ed.listing$Aquatic == 1 | ## can't use aquatics
                            (ed.listing$state + ed.listing$ipaus + ed.listing$federal + ed.listing$other) < 1 |
                            ## can't use species that haven't been recognized as problematic
                            (ed.listing$no.abunG2 < 20) |
                            ## data deficient
                            (ed.listing$nativity == "N") ] <- 0  
## can't use native plants


write.csv(ed.listing, "SpeciesList10_11_2017.csv", row.names=F)
## write out for later use

#################### Add Species Trait Data ####################
masterList <- read.table("file:///C:/Users/mwone/Google Drive/NSF_GSS_shared/Hotspots_and_Abundance/master list.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

speciesList <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList08_28_1017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
head(speciesList)
head(masterList)

traits <- data.frame(cbind(masterList$NewCode, masterList$habit, masterList$duration_final), stringsAsFactors=F)
head(traits)

colnames(traits)<-c("usda.code", "habit", "duration")
traits$usda.code[1] <- "FEBR7"
head(traits)

spList <- merge(speciesList, traits, by="usda.code", all.x=T)

#write.csv(spList, "SpeciesList09_26_2017.csv", row.names=F)
speciesList <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList10_11_2017.csv", header = T, sep = ",", quote= "\"", 
                          comment.char= "", stringsAsFactors = F, strip.white = T)

speciesList <- speciesList[speciesList$potential.use_10ac == 1,]


speciesList$habitF <- as.factor(speciesList$habit)
speciesList$durationF <- as.factor(speciesList$duration)
summary(speciesList$durationF)
summary(speciesList$habitF)

speciesList$both <- as.factor(paste(speciesList$duration, speciesList$habit, sep=" "))
summary(speciesList$both)

hist(speciesList$no.abunG10, breaks=25, main=" ", 
     xlab= "no. of grid cells with record", ylab= "no. of species",
     cex.axis=1.2, cex.lab=1.5, mar=c(5, 4, 4, 2))

# summary(spList2$no.abunG)
# str(spList2$no.abs)
# length(spList2$no.absG[spList2$no.absG > 0]) # 48
# summary(spList2$no.absG[spList2$no.absG > 0])# 1-1,700; mean 326
# hist(spList2$no.absG[spList2$no.absG > 0], breaks=10)
# 
# summary(spList2$no.absG[spList2$no.absG > 0]/(spList2$no.abunG[spList2$no.absG > 0]+ spList2$no.absG[spList2$no.absG > 0]))
# hist(spList2$no.absG[spList2$no.absG > 0]/(spList2$no.abunG[spList2$no.absG > 0]+ spList2$no.absG[spList2$no.absG > 0]))
# ##absence heavy is a prob
# 
# edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_09_12_2017.csv",  
#                   header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
# 
# ## read in species list
# spp <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/SpeciesList08_28_2017.csv",  
#                   header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
# sp.list <- spp$usda.code[spp$potential.use2 == 1] ## list of useable species
# 
# edd <- edd[edd$usda %in% sp.list,] ## subset dataset to useable species
# coordinates(edd) <- c(7,8) ## set lat/lon
# proj4string(edd) <- CRS("+init=epsg:4326") ## assume wgs84 for all eddmaps data
#
# states <- readOGR(dsn="C:/Users/mwone/Documents/geodata/states", layer="US_states")
# states <- spTransform(states, proj4string(edd))
# plot(states)
# plot(edd, col="red", pch=19, cex=.3, add=T)
# plot(edd[edd$absence==1,], col="black", pch=19, cex=.3, add=T)
# 
# proj4string(edd)
# plot(states)
#
# splist<- spp[spp$potential.use2 == 1,]
# splist$a <- (splist$no.abs-splist$no.absG)/splist$no.abs
# splist$p <- (splist$no.abun-splist$no.abunG)/splist$no.abun
# 
# summary(splist$a)
# summary(splist$p)
# hist(splist$a)
# hist(splist$p)
# 
# edd <- read.table("C:/Users/mwone/Documents/EddmapS data/eddmaps_prepped_08_22_2017.csv", header = T, sep = ",", quote= "\"", 
#                   comment.char= "", stringsAsFactors = F, strip.white = T)
# 
# negs <- edd[edd$negative ==1,]
# negs <- negs[negs$USDAcode %in% sp.list,]
# summary(as.factor(negs$ReporterFULLName))
# negs$Comments[negs$Comments != "NULL"]
# 
# head(spp)
# hist(spp$no.abunG[spp$potential.use2 == 1], breaks=30, main="Histogram of number of cells with records per species", xlab="no. of cells with records", ylab="no. of species")
