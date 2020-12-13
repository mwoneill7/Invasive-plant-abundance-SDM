#### Making clean files to share in publication
###  created: 9/11/2019
##   last modified: 9/11/2019
#    Mitch O'Neill




#### files needed
## 1. environmental layers
## 2. list of species ---- sp.list <- unique(speciesdata$code)
## 3. bias file ---- should I make it in the code
## 4. ordsums ---- should I make it in the code?

edd <- read.table("file:///C:/Users/mwone/Documents/abs_pres_occbias_figs/edd_w_environmental_NEW.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
head(edd)


spp <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/speciesList_1_25_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)


spp <- data.frame(spp$scientific.name,spp$usda.code)
head(spp)
#edd$species.name <- "TEST"

edd2 <- merge(edd,spp, by.x="species", by.y="spp.usda.code", all.x=T, all.y=F)
head(edd2)

#for(i in 1:length(edd$species)){
#  
#  edd$species.name[i] <- spp$spp.scientific.name[spp$spp.usda.code==edd$species[i]]
#  
#  print(i)
#}




edd2 <- data.frame(edd2$species, edd2$spp.scientific.name, edd2$abundance, edd2$bio_2, 
                   edd2$bio_5, edd2$bio_6,  edd2$bio_8, edd2$bio_12, edd2$bio_15, edd2$nlcd_3, 
                   edd2$nlcd_4, edd2$nlcd_5, edd2$nlcd_6, edd2$nlcd_7, edd2$nlcd_8)
head(edd2)

colnames(edd2) <- c("usda_code","scientific_name","abundance", "bio_2", "bio_5", "bio_6", "bio_8", 
                    "bio_12", "bio_15", "nlcd_3", "nlcd_4", "nlcd_5", "nlcd_6", "nlcd_7", "nlcd_8")
head(edd2)
length(unique(edd2$usda_code))
summary(edd2$usda_code, max.sum=200)


list <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/spplist70.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
head(list)

edd2 <- edd2[edd2$usda_code %in% list$x,]

length(unique(edd2$usda_code))
write.csv(edd2,"file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/Cleaned_up/abundance.csv", row.names=F)
