### Combining eddmaps species list with w/ invasive listings
##  Mitch O'Neill
#   6/13/2017; 7/18/2017

setwd("C:/Users/mwone/Google Drive")


#### eddmaps species list (all species with infested area downloaded from EDDMAPS October 2016)
edd.list <- read.table("NSF_GSS_shared/Synonomies/speciesList06_01_2017.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
edd.list$edd <- 1



#### list of species from the occurence-based hotspots analysis (Allen and Bradley 2016)
occ.list <- read.table("NSF_GSS_shared/Hotspots_and_Abundance/occ_data_with_2017_codes.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
occ.list$occ <- 1
head(occ.list)



#### IPAUS species list (Downloaded in 2014; database last updated May 2012)
#~ipa.list <- read.table("NSF_GSS_shared/Hotspots_and_Abundance/IPAUS species list with USDA codes.csv", header = T, sep = ",", 
#~                       stringsAsFactors = F, quote= "\"", comment.char= "", strip.white = T)
#~ipa.list$IPA <- 1
#~head(ipa.list)
#~
#~#### read in file with list of IPAUS species that don't match
#~ipa.list$Scientific.Name <- ipa.list$scientific_name
#~
#~ipa.prob.list <- read.table("NSF_GSS_shared/Synonomies/IPAUS_prob_spp.csv", header = T, sep = ",", 
#~                            stringsAsFactors = F, quote= "\"", comment.char= "", strip.white = T)
#~
#~spp.list1 <- as.list(ipa.prob.list$ipa_species)
#~
#~for (i in 1:length(spp.list1)) {
#~  ipa.list$Scientific.Name[ipa.list$scientific_name == spp.list1[i]] <- ipa.prob.list$name[ipa.prob.list$ipa_species == spp.list1[i]]
#~  print(i)
#~}
#~


#### USDA federal/State listings (Downloaded from PLANTS 5/30/2017)
usda.list <- read.table("NSF_GSS_shared/Synonomies/PLANTS.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
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




## CHECK WHERE IPA SPECIES NAMES OCCUR MORE THAN ONCE ON USDA-list before merge
#    list.check <- as.list(unique(ipa.list$Scientific.Name))
#    
#    sp.name <- "TEMPLATE"
#    matches <- -9
#    check <- data.frame (sp.name, matches)
#    
#    for (i in 1:length(list.check)){
#      matches.i <- usda.list[usda.list$Scientific.Name == list.check[i], ]
#      matches <- length(matches.i$Scientific.Name)
#      sp.name <- as.character(list.check[i])
#      check.i <- data.frame(sp.name, matches)
#      check <- rbind(check, check.i)
#      
#    }
#    
#    check.probs <- check[check$matches != 1 & check$sp.name != "TEMPLATE", ]

## there are some potential problems with using Scientific name to merge ipa.list and usda.list, for some species on
## ipa.list appear more than once on usda.list

## loop through all species on ipa.list, extract USDA codes, flag manuals
## manuals using author names on the IPA list!

#     for (i in 1:length(ipa.list$Scientific.Name)) {
#       
#       n <- length(usda.list$Accepted.Symbol[usda.list$Scientific.Name == ipa.list$Scientific.Name[i]])
#       
#       if (n == 0) { ## if the species is not in the PLANTS database, the new USDA code is "not found"
#         ipa.list$USDA_code_2017[i] <- "not found" } 
#       if (n == 1)  { ## if the species is in the PLANTS database just once
#         ipa.list$USDA_code_2017[i] <- usda.list$Accepted.Symbol[usda.list$Scientific.Name == ipa.list$Scientific.Name[i]] }
#       ## pulls the new USDA code from the PLANTS database for the species in this iteration of the loop
#       if (n > 1)  { ## if the species is not in the PLANTS database more than once
#         ipa.list$USDA_code_2017[i] <- "manual" } ## flag to be dealt with manually
#       
#       print(i)
#       
#     }
#write.csv(ipa.list, file="ipa2---.csv", row.names=F)


#    unique(ipa.list$Accepted_Code[ipa.list$USDA_code_2017 == "not found"])
#    ## if only value is "" the "not founds are all plants that don't occur in USDA PLANTS
#    ipaWnoUSDA <- ipa.list[ipa.list$USDA_code_2017 == "not found",]
#    ## keep for later to cross ref with eddmaps species with no USDA

## after manually editing the USDA codes, read it back in
ipa.list <- read.table("file:///C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/SpList/ipa2.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")


#ipa.list$merge.column <- paste(ipa.list$USDA_code_2017, ipa.list$Scientific.Name, sep = " ")
#usda.list$merge.column <- paste(usda.list$Accepted.Symbol, usda.list$Scientific.Name, sep = " ")
#usda.ipa <- merge(ipa.list, usda.list, by = "merge.column", all = TRUE)
#write.csv(usda.ipa, file="C:/Users/mwone/Desktop/merge1.csv", row.names=F)
#ed.list[295]



edd <- read.table("C:/Users/mwone/Google Drive/NSF_GSS_shared/Watch_Lists/eddmaps_prepped.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
edd <- edd[edd$AbsenceData == 1 | edd$validInfestedArea == 1, ]

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
###### USDA CODES, take out once edd_prepped is updated with these codes ########
#################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

## list of species in USDA-sp name file
spp <- as.list(edd.list$eddmap_species[!is.na(edd.list$eddmap_species)])

for (sp in spp){
  length(edd$ScientificName[is.na(edd$invasivecover)])
  
  edd$USDAcode[edd$ScientificName == sp] <- edd.list$USDA_code[edd.list$eddmap_species == sp & !is.na(edd.list$eddmap_species)]
  print(sp)  
}

edd <- edd[!is.na(edd$USDAcode) & edd$USDAcode != "", ]

check <- edd[edd$Infestedarea == "NULL" &
             edd$Grossarea == "NULL" &
             edd$invasivecover == "NULL" &
             edd$TXabundance == "NULL" &
             edd$density == "NULL" &
             edd$stemcount == "NULL" &
             edd$Aggressiveness == "NULL" &
             edd$NumberObserved == "NULL", ]
check <- check[!is.na(check$USDAcode), ]
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
edd <- edd[!is.na(edd$USDAcode) & edd$USDAcode != "", ]


ed.list <- as.list(unique(edd.list$USDA_code[edd.list$USDA_code != ""]))

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

ed.listing <- data.frame(scientific.name, usda.code, comment, state, federal, other, ipaus, hotspot, no.abun, no.abs, nativity, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)


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

edd.abs <- edd[edd$USDAcode == ed.list[i] & edd$AbsenceData == 1, ]
no.abs <- length(edd.abs$USDAcode[!is.na(edd.abs$USDAcode)])

 
  
  ed.listing.i <- data.frame(scientific.name, usda.code, comment, state, federal, other, ipaus, hotspot, no.abun, no.abs, L48, usdaCHK, stringsAsFactors = F) ## no.abun, no.abs, usdaCHK, stringsAsFactors = F)
  names(ed.listing.i) <- c("scientific.name", "usda.code", "comment", "state", "federal", "other", "ipaus", "hotspot", "no.abun", "no.abs", "nativity", "usdaCHK") ## "no.abun", "no.abs", "usdaCHK")  
  ed.listing <- rbind(ed.listing, ed.listing.i)
  print(i)
}

ed.listing <- ed.listing[ed.listing$usda.code != "TEMPLATE", ]
unique(ed.listing$usda.code)
summary(ed.listing)




###################################################################################################
####### complete nad27 file only once
nad27 <- read.table("file:///C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/Datums/nadprops.csv", sep = ",", stringsAsFactors = F, header = T)
unique(nad27$usda.code2017)


nad.list <- nad27$spp

for (i in 1:length(nad.list)) {

## if the usda code has not been manually enterred, pull it from somewhere else
   if (nad27$usda.code2017[nad27$spp == nad.list[i]] == "") {
  
  ## pull from plants
  n <- length(usda.list$Accepted.Symbol[usda.list$Scientific.Name == ipa.list$Scientific.Name[i]])
          
          if (n == 0) { ## if the species is not in the PLANTS database, the new USDA code is "not found"
            nad27$usda.code2017[nad27$spp == nad.list[i]] <- "not found" 
            } 
          if (n == 1 & length(edd.list$USDA_code[edd.list$PLANTS_species == nad.list[i]]) > 0)  { ## if the species is in the PLANTS database just once
            nad27$usda.code2017[nad27$spp == nad.list[i]]<- usda.list$Accepted.Symbol[usda.list$Scientific.Name == nad.list[i]]
            }
          if (n > 1)  { ## if the species is not in the PLANTS database more than once
            nad27$usda.code2017[nad27$spp == nad.list[i]] <- "manual" 
            } ## flag to be dealt with manually
          }
  
       ## pull from ipa2
   if (nad27$usda.code2017[nad27$spp == nad.list[i]] == "" & length(ipa.list$USDA_code_2017[ipa.list$scientific_name == nad.list[i]]) > 0) {
       nad27$usda.code2017[nad27$spp == nad.list[i]] <- ipa.list$USDA_code_2017[ipa.list$scientific_name == nad.list[i]]
   } 
       ## pull from  edd.list
   if (nad27$usda.code2017[nad27$spp == nad.list[i]] == "" & length(edd.list$USDA_code[edd.list$PLANTS_species == nad.list[i]]) > 0) {
       nad27$usda.code2017[nad27$spp == nad.list[i]] <- edd.list$USDA_code[edd.list$eddmap_species == nad.list[i]]
   } 
   if (nad27$usda.code2017[nad27$spp == nad.list[i]] == "" & length(edd.list$USDA_code[edd.list$PLANTS_species == nad.list[i]]) > 0) {
       nad27$usda.code2017[nad27$spp == nad.list[i]] <- edd.list$USDA_code[edd.list$PLANTS_species == nad.list[i]]
      }
 print(i)
 
}

### write out to finish manually
write.csv(nad27, "C:/Users/mwone/Desktop/GoogleDriveOverflow/nad27b.csv")

#############################################################################

nad27b<- read.table("C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/Datums/nad27b.csv", sep = ",", stringsAsFactors = F, header=T)
nad.codes <- nad27b$usda.code2017

for (i in 1:length(nad.codes)) {
  ed.listing$nad27[ed.listing$usda.code == nad.codes[i]] <- nad27b$nad27[nad27b$usda.code2017 == nad.codes[i]]
  
}

summary(ed.listing)

write.csv(ed.listing, "C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/speciesList.csv")




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


## two spp that are on ipaus were not in ipaus file, so i'm enterring their ipaus status manually
chk <- ed.listing[!is.na(ed.listing$ipaus2) & ed.listing$ipaus2 == 1, ]
ed.listing$ipaus[!is.na(ed.listing$ipaus2) & ed.listing$ipaus2 == 1] <- 1
ed.listing[ed.listing$state == 0 & ed.listing$federal == 0 & ed.listing$hotspot == 1 & ed.listing$other == 0 & ed.listing$ipaus == 0, ]
ed.listing$ipaus[ed.listing$scientific.name == "Centaurea x moncktonii"] <- 1


## rate useability: 0= unuseable, 1= possible, 2= probable, 3= useable

ed.listing$use[ed.listing$no.abun < 15 | (ed.listing.i$no.abs + ed.listing$no.abun) < 20] <- 0
## data deficient
ed.listing$use[ed.listing$state == 0 & ed.listing$other == 0 & ed.listing$ipaus == 0 & ed.listing$federal == 0 ] <- 0
## not listed as invasive
ed.listing$use[is.na(ed.listing$use) & (ed.listing$state == 1 | ed.listing$other == 1 | ed.listing$ipaus == 1 | ed.listing$federal == 1) &
               ed.listing$nativity == "I"] <- 3
## enough data, listing, introduced
ed.listing$use[is.na(ed.listing$use) & (ed.listing$state == 1 | ed.listing$other == 1 | ed.listing$ipaus == 1 | ed.listing$federal == 1) &
                 ed.listing$nativity == "N/A"] <- 2
ed.listing$use[is.na(ed.listing$use) & ed.listing$hotspot == 1] <- 2
ed.listing$use[is.na(ed.listing$use) & ed.listing$hotspot == 0] <- 1

length(ed.listing$use[ed.listing$use == 2])

summary(ed.listing$use)


write.csv(ed.listing, "C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/splist.csv", row.names=F )


####################################################################################################



ed.listing <- read.table("file:///C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/SpList/splistUPDATED.csv", sep = ",", header=T, stringsAsFactors = F)
summary(ed.listing)




ed.listing$use2[ed.listing$no.abun < 15 | (ed.listing$no.abs + ed.listing$no.abun) < 20] <- 0
## data deficient
ed.listing$use2[ed.listing$state == 0 & ed.listing$other == 0 & ed.listing$ipaus == 0 & ed.listing$federal == 0 ] <- 0
## not listed as invasive
ed.listing$use2[is.na(ed.listing$use2) & (ed.listing$state == 1 | ed.listing$other == 1 | ed.listing$ipaus == 1 | ed.listing$federal == 1) &
                 ed.listing$nativity == "I"] <- 3
## enough data, listing, introduced
ed.listing$use2[is.na(ed.listing$use2) & (ed.listing$state == 1 | ed.listing$other == 1 | ed.listing$ipaus == 1 | ed.listing$federal == 1) &
                 ed.listing$nativity == "N/A"] <- 2
ed.listing$use2[is.na(ed.listing$use2) & ed.listing$hotspot == 1] <- 2
ed.listing$use2[is.na(ed.listing$use2) & ed.listing$hotspot == 0] <- 1

ed.listing$useNEW <- ed.listing$use2/3

spp.manual <- ed.listing$usda.code[!is.na(ed.listing$useb)]

for (i in 1:length(spp.manual)){
  ed.listing$useNEW[ed.listing$usda.code == spp.manual[i]] <- ed.listing$useb[ed.listing$usda.code == spp.manual[i]]
  print(i)
}

write.csv(ed.listing, "C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/splistNEW.csv", row.names=F )




## who knows whats going on below this point  


#ed.listing3 <- ed.listing
#head(ed.listing3)
#
#
#abunYES <- ed.listing3[(ed.listing3$no.abun >= 15 & (ed.listing3$no.abun + ed.listing3$no.abs) >= 20) & 
#                (ed.listing3$hotspot == 1 | 
#                (ed.listing3$nativity == "I" & 
#                (ed.listing3$federal == 1 | ed.listing3$state == 1 | ed.listing3$other == 1 | ed.listing3$ipaus == 1))), ] #<- 1
#
#abunMANUAL <- ed.listing3[(ed.listing3$no.abun >= 15 & (ed.listing3$no.abun + ed.listing3$no.abs) >= 20) & 
#                          ed.listing3$hotspot == 0 & ed.listing3$nativity != "I" & 
#                          !(ed.listing3$federal == 0 & ed.listing3$state == 0 & ed.listing3$other == 0 & ed.listing3$ipaus == 0), ] #<- 1
#
#
#
#abunPROB <- ed.listing3[(ed.listing3$no.abun >= 15 & (ed.listing3$no.abun + ed.listing3$no.abs) >= 20) &
#                        (ed.listing3$hotspot == 1 | 
#                        ed.listing3$L48 == "I" ), ] 
#                        #(ed.listing3$federal == 1 | ed.listing3$state == 1 | ed.listing3$other == 1 | ed.listing3$ipaus == 1))), ] #<- 1
#
#
#
#length(ed.listing3$usda.code[(ed.listing3$no.abun >= 15 & (ed.listing3$no.abun + ed.listing3$no.abs) >= 20)])



#write.csv(abunYES, "abunYES2.csv", row.names = F)
#
#write.csv(abunPROB, "probableSPPlist.csv", row.names = F)
#
#head(abunyesOLD)
#
#abunyesOLD <- read.table("C:/Users/mwone/Google Drive/NSF_GSS_shared/Hotspots_and_Abundance/tentative_abundance_spp_list.csv", sep=",", stringsAsFactors = F, header = T)
#
#abun.mixed <- merge(abunYES, abunyesOLD, by = "usda.code", all=TRUE)
#write.csv(abun.mixed, "abunYES3.csv", row.names = F)

aepo <- edd[edd$USDAcode == "AEPO", ]
alpe4 <- edd[edd$USDAcode == "ALPE4", ]
alpe4 <- alpe4[alpe4$AbsenceData == 1,]
alpe4 <- alpe4[!is.na(alpe4$USDAcode), ]


length(abunYES$usda.code[abunYES$no.abs >= 1])
## 107
length(abunYES$usda.code[abunYES$no.abs >= .25*(abunYES$no.abs + abunYES$no.abun)])
length(abunYES$usda.code[abunYES$no.abs >= .33*(abunYES$no.abs + abunYES$no.abun)])
length(abunYES$usda.code[abunYES$no.abs >= .50*(abunYES$no.abs + abunYES$no.abun)])






#######################################################################
setwd("C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/SpList/speciesDecisions/")
edd.list <- read.table("SpeciesList71217.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
aqua <- read.table("aquatics.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")

for (i in 1:length(aqua$usda.code)){
  edd.list$Aquatic[edd.list$usda.code == aqua$usda.code[i]] <- 1
  print(i)
}
summary(edd.list$Aquatic)
edd.list$use[edd.list$Aquatic ==1] <- 0
write.csv(edd.list, "SpeciesList71417.csv", row.names = F)


##### merge subspecies
edd$USDAcodeUNMERGED <- edd$USDAcode

ssp <- read.table("file:///C:/Users/mwone/Desktop/GoogleDriveOverflow/WorkingFolder/SpList/speciesDecisions/subspecies.csv",  header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")

for (i in 1:length(ssp$usda.code)){
  edd$USDAcode[edd$USDAcodeUNMERGED == ssp$usda.code[i]] <- ssp$usda.lump[i]
  print(i)
}


#######################################################################
################ Don't use my absences, just negatives ################
edlist <- read.table("C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/SpList/speciesDecisions/SpeciesList71817.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
#edd <- read.table("C:/Users/mwone/Desktop/GoogleDriveOverFlow/eddmaps_prepped_6_30_2017.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")

ed.list <- as.list(unique(edlist$usda.code[edlist$usda.code != ""]))


for (i in 1:length(ed.list)){
  edd.abs <- edd[edd$USDAcode == ed.list[i] & edd$negative == 1, ]
  edlist$no.neg[edlist$usda.code == ed.list[i]] <- length(edd.abs$USDAcode[!is.na(edd.abs$USDAcode)])
  print(i)
}

edlist$no.pot.abs <- edlist$no.abs - edlist$no.neg
edlist$enoughdata[(edlist$no.neg + edlist$no.abun) >= 20 & edlist$no.abun >= 15] <- 1
edlist$enoughdata[(edlist$no.neg + edlist$no.abun) < 20 | edlist$no.abun < 15] <- 0

summary(edlist$enoughdata)
#write.csv(edlist, "C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/SpList/speciesDecisions/SpeciesList71717.csv", row.names = F)


###########################################################################################


edd$LatN <- as.numeric(edd$Latitude_Decimal)
edd$LonN <- as.numeric(edd$Longitude_Decimal)
### one has a positive Longitude, but it says fairfax VA and switching sign will place it in Fairfax.
edd$LonN[edd$LonN == 77.30791] <- -77.30791 

edd$coord <- "bad"
edd$coord[(edd$LatN > 24 & edd$LatN < 50) & (edd$LonN < -65 & edd$LonN > -125)] <- "good"

edd$EXP <- 0
edd$EXP[edd$validInfestedArea == 1 & edd$coord == "good"& edd$ReporterFULLName != "US Army Corps of Engineers Ombil Database" ] <- 1



for (i in 1:length(ed.list)){
  edd.i <- edd[edd$USDAcode == ed.list[i] & edd$negative == 0 & edd$EXP == 1, ]
  ## exp = 1 means the lat/lon falls in the bounding box, not reported by USACE, valid infested area
  edd.i <- edd.i[!is.na(edd.i$USDAcode),]
  
  edlist$infestedBelow5[edlist$usda.code == ed.list[i]] <-     length(edd.i$infestedAreaInAcres[edd.i$infestedAreaInAcres <= 5] )
  edlist$infestedBelow100[edlist$usda.code == ed.list[i]] <-   length(edd.i$infestedAreaInAcres[edd.i$infestedAreaInAcres <= 100])
  edlist$DensityBelowOne5[edlist$usda.code == ed.list[i]] <-   length(edd.i$infestedAreaInAcres[edd.i$infestedAreaInAcres/edd.i$grossAreaInAcres == 1 & edd.i$infestedAreaInAcres <= 5])
  edlist$DensityBelowOne100[edlist$usda.code == ed.list[i]] <- length(edd.i$infestedAreaInAcres[edd.i$infestedAreaInAcres/edd.i$grossAreaInAcres == 1 & edd.i$infestedAreaInAcres <= 100])
  print(i)
}



edlist$enoughdataI5   <- 0
edlist$enoughdataI100 <- 0
edlist$enoughdataD5   <- 0
edlist$enoughdataD100 <- 0

edlist$enoughdataI5[edlist$infestedBelow5 + edlist$no.neg >= 20 & edlist$infestedBelow5 >= 15]           <- 1
edlist$enoughdataI100[edlist$infestedBelow100 + edlist$no.neg >= 20 & edlist$infestedBelow100 >= 15]     <- 1
edlist$enoughdataD5[edlist$DensityBelowOne5 + edlist$no.neg >= 20 & edlist$DensityBelowOne5 >= 15]       <- 1
edlist$enoughdataD100[edlist$DensityBelowOne100 + edlist$no.neg >= 20 & edlist$DensityBelowOne100 >= 15] <- 1

edlist <- edlist[!is.na(edlist$use),]


write.csv(edlist, "C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/SpList/speciesDecisions/SpeciesList71817.csv", row.names = F)
write.csv(edd, "C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/eddmaps_prepped_07_18_17.csv", row.names = F)
edd <- read.csv("C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/eddmaps_prepped_07_18_17.csv", header = T, sep = ",", stringsAsFactors = F, strip.white = T, quote= "\"", comment.char= "")
head(edd)
write.csv(edd, "C:/Users/mwone/Desktop/GoogleDriveOverFlow/WorkingFolder/eddmaps_prepped_07_18_17.csv", row.names = F)


