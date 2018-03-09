### Code to run Maxent models
## find 95% MTP threshold values 
## and create binary maps

################################################################################
## Setup
###############################################################################
library(maptools)
library(raster)
library(rgdal)
library(colorRamps)

### memory ################################################################################
#memory.limit()     #16341
#memory.size()      #640.54
#memory.size(max=T) #695.69
#memory.size(max=F) #642.95

#options(java.parameters = "-Xmx14g" )
#memory.limit(size=NA)

#install.packages("colorRamps")
## download MaxEnt: https://biodiversityinformatics.amnh.org/open_source/maxent/
## note that this code was developed with MaxEnt v.3.3.1 and has not been tested with the new release

# set where yout maxent jar file is
#maxent.location='C:/Users/mwone/Google Drive/Ordinal/Occurence/MaxEnt/maxent.jar'
#maxent.location='C:/Users/Localadmin/Google Drive/Ordinal/Occurence/MaxEnt/maxent.jar'
#maxent.location='C:/Users/mwone/Desktop/MaxEnt/maxent.jar'
#########################################################################################

maxent.location='C:/Users/Localadmin/Documents/MaxEnt/maxent.jar'
## put all environmental data in one folder-- all files must be ascii, with EXACT same resolution and extent

#########################################################################################

## define directory of environmental data for fitting
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs",sep="")

sp.list <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/MaxEntFiles/spForms3_4_2018.csv", header = T, sep = ",", 
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
head(sp.list)
#sp.list<- sp.list[order(sp.list$code),]
#head(sp.list)

for(i in 1:length(sp.list$code)){
  
  ## set up model options, linear and quadratic features only, use sampling bias surface
  spp_bias_log=paste0("java -jar ",maxent.location, 
                      " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature ",
                      sp.list$ignore[i], " biastype=3")
  ##spp_bias_log=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature -N pop_2_5_us -N road_2_5_us biastype=3")
  
  ### ALL SPECIES ###
  ## define point location samples (all species in one .csv), must be set up with three columns (species, lon, lat, in that order)
  samples=paste("samplesfile=", paste("C:/Users/Localadmin/Documents/Maxent_modeling/FULL/species", 
                                      paste(sp.list$code[i], "csv", sep="."), sep="/"), sep="")
  ## define output directory
  directory<- paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/FINAL_OUTPUT_FULL", sp.list$code[i], sep="/")
  dir.create(directory)
  output=paste("outputdirectory=", directory, sep="")
  
  bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc",sep="")
  ## average bias ascii
  
  ## call model to run
  system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))
  
  
  
  ### REPEAT FOR ONLY ABUNDANCE POINTS ###
  ## define point location samples (all species in one .csv), must be set up with three columns (species, lon, lat, in that order)
  samples=paste("samplesfile=", paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/species", 
                                      paste(sp.list$code[i], "csv", sep="."), sep="/"), sep="")
  ## define output directory
  directory <- paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/FINAL_OUTPUT_ABUN", sp.list$code[i], sep="/")
  dir.create(directory)
  output=paste("outputdirectory=", directory, sep="")
  
  bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/BIAS_avgABUN.asc",sep="")
  ## average bias ascii
  
  ## call model to run
  system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))
  
  print(i)
  write.csv(i,"progressMAXENT.csv",row.names=F)
}


####################################################################################
## Find 95% threshold for each species AND get training/test sample sizes
###################################################################################
## maxent output represents relative occurence rate (ROR)
## we want a binary map, so we need to know what ROR value to use to define "suitable" vs "unsuitable"
## all maxent surfaces must sum to 1, so the threshold value can vary from species to species
## therefore we need to calculate the threshold by species (below)

## create list of modeled species
#d=read.table("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/IAS_occurences_final_analysis.csv",sep=",",header=T)

#setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/")

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_FINAL.csv", header = T, sep = ",",  
                                 quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

spp=ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP < 0.05]
#colnames(spp)="PLANT_CODE"
#spp$PLANT_CODE=as.character(spp$PLANT_CODE)
#head(spp)
#str(spp)
#dim(spp)


dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis")
dataset_type <- c("FULL","ABUN")

for (type in dataset_type) {
  
  
  #create dataframe to hold 95% MTP logistic threshold values for each run
  spp_val=data.frame("Species",15,0,"CtySppRich_Bias",100000, 100000,stringsAsFactors=F) 
  colnames(spp_val)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n","train_n")
  #spp_val$Species=as.character(spp_val$Species)
  #spp_val$Bias=as.character(spp_val$Bias)
  
  ## get all output files (to deal with any skipped spp and variable run number) C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/FINAL_OUTPUT_ABUN
  setwd(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/",type,"/FINAL_OUTPUT_",type))
  #a=list.files(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE)
  
  for (i in spp){        #loop for each species
    for (j in 0:9){      #loop for each iteration within a species
      
      #spp_iter=paste(i,"_",j,"_samplePredictions.csv",sep="")
      
      #if(spp_iter %in% a == TRUE) {
      file=paste(i,"/",i,"_",j,"_samplePredictions.csv",sep="") #iteration j point-level model output 
      
      d=read.table(file, sep=",",header=T)
      d$run=j
      dt=d[d$Test.or.train == "train",]  #subset iteration data to training data only
      dtest=d[d$Test.or.train == "test",]  #subset iteration data to testing data only
      
      samp=dim(dt)
      samp1=dim(dtest)
      
      train_n=samp[1] #get training sample size for this iteration
      test_n=samp1[1]  #get test sample size for this iteration
      
      tp=round(0.95*samp[1], 0) #number of points for 95% threshold in each iteration
      
      dt=dt[order(dt$Logistic.prediction, decreasing=T),] #order the logistic output largest to smallest
      
      sp_runs=data.frame(i, j, dt$Logistic.prediction[tp], "CtySppRich_Bias", test_n, train_n) #create dataframe with species, run, logistic value of the point that bounds the 95% percentile, bias model, test and train n
      colnames(sp_runs)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n", "train_n")
      
      spp_val=rbind(spp_val,sp_runs) 
      #}

            
      #print(j)
    }
    ### get AUC scores
    aucs=read.table(paste(i,"/maxentResults.csv",sep=""), header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
    if(type=="FULL"){
       ordsums$AUCtest_full[ordsums$species.code == i] <- aucs$Test.AUC[11]
       ordsums$AUCtrain_full[ordsums$species.code == i] <- aucs$Training.AUC[11]
       #print(i)
      } else {
       ordsums$AUCtest_abun[ordsums$species.code == i] <- aucs$Test.AUC[11]
       ordsums$AUCtrain_abun[ordsums$species.code == i] <- aucs$Training.AUC[11]
        
      }
    print(i) 
  }
  
  
  ## trim junk row from spp_val
  #head(spp_val)
  #tail(spp_val)
  #dim(spp_val)
  #spp_val1=spp_val[2:8900,]
  spp_val1=spp_val[spp_val$Species!="Species",]
  
  
  ## calculate mean threshold per species per range from spp_val1
  thresh=as.data.frame(tapply(spp_val1$X95_MTP_thresh_logistic, list(spp_val1$Species, spp_val1$Bias), mean))
  thresh$Species=rownames(thresh)
  #head(thresh)
  #dim(thresh)
  
  test=as.data.frame(tapply(spp_val1$test_n, list(spp_val1$Species, spp_val1$Bias), mean))
  test$Species=rownames(test)
  colnames(test)=c("test_n","Species")
  #head(test)
  #dim(test)
  
  train=as.data.frame(tapply(spp_val1$train_n, list(spp_val1$Species, spp_val1$Bias), mean))
  train$Species=rownames(train)
  colnames(train)=c("train_n","Species")
  #head(train)
  #dim(train)
  
  samples=cbind(test, train)
  #head(samples)
  #colnames(samples)=c("test_n","Species","train_n","Species1")
  #samples$match[samples$Species == samples$Species1]="TRUE"
  #print(table(samples$match))
  
  samples1=subset(samples, select=c("Species","test_n","train_n"))
  samples1$total_n=samples1$test_n + samples1$train_n
  samples1$Bias="CtySppRich"
  samples1$Climate="Current"
  #head(samples1)
  #dim(samples1)
  
  write.table(thresh, file=paste0("C:/Users/Localadmin/Documents/Maxent_modeling/thresholds",type,".csv"),sep=",",row.names=F)
  write.table(samples1,file=paste0("C:/Users/Localadmin/Documents/Maxent_modeling/sample_sizes",type,".csv"),sep=",",row.names=F) 
  

  
  #########################################################################
  ### Create binary maps for each species
  #########################################################################
  
  ## use the thresholds above to go from ROR to binary
  
  ### create dataset to hold suitable cell counts for each species 
  #range_size=data.frame(-100,-100,"NA","NA")
  #colnames(range_size)=c("Var1","Freq","spp","bias")
  #range_size
  
  dir.create(paste0("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",type))
  
  
  ###### create binary rasters based on 95% MTP
  
  for (i in unique(thresh$Species)){  
    file=paste(i,"/",i,"_avg.asc",sep="")
    d1=readAsciiGrid(file) 
    
    d1@data[d1@data >= thresh$CtySppRich_Bias[thresh$Species==i]]=1
    d1@data[d1@data < thresh$CtySppRich_Bias[thresh$Species==i]]=0
    
    out_file=paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",type,"/",i,"_",type,".asc",sep="")
    write.asciigrid(d1, out_file)
    
    tab=as.data.frame(table(d1@data))
    tab$spp=i
    #tab$bias="CtySppRich_Bias"
    
    #range_size=rbind(range_size,tab)
    #image(d1, col=c("grey","blue"))
    
    if(type=="FULL"){
      ordsums$range_size_full[ordsums$species.code==i] <- tab$Freq[tab$Var1==1]
      ordsums$range_check_full[ordsums$species.code==i] <- sum(tab$Freq)
    } else {
      ordsums$range_size_abun[ordsums$species.code==i] <- tab$Freq[tab$Var1==1]
      ordsums$range_check_abun[ordsums$species.code==i] <- sum(tab$Freq)
    }
      
    
    rm(d1)
    gc()
    print(i)
  }
  
  #head(range_size)
  #dim(range_size)
  #range_size1=range_size[2:1793,]
  #range_size1=range_size[range_size$Var1 == -100,]
  
  #range_size2=range_size1[range_size1$Var1 ==1,]
  #head(range_size2)
  #dim(range_size2)
  
  print(type)
}

summary(ordsums$AUCtest_full)
summary(ordsums$AUCtrain_full)
summary(ordsums$AUCtest_abun)
summary(ordsums$AUCtrain_abun)
plot(ordsums$AUCtest_full~ordsums$AUCtest_abun)
abline(0,1)

summary(ordsums$range_check_abun)
summary(ordsums$range_check_full)
summary(ordsums$range_size_abun)
summary(ordsums$range_size_full)

summary(ordsums$range_size_abun/ordsums$range_size_full)
write.csv(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_maxent.csv", row.names=F)

####################################################
##### compare Maxent Full VS Abun w/in and across sp
####################################################

setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis")

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_maxent.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

ordsums <- ordsums[ordsums$kappa > 0 & ordsums$kappaP <0.05,]
summary(ordsums$AUCtest_abun)


library(raster)
library(rgdal)

#########################################
#species     <- "TEMPLATE"
#auc         <- -99

#ncell_abun  <- -99
#lat_abun    <- -99
#lon_abun    <- -99

#ncell_full  <- -99
#lat_full    <- -99
#lon_full    <- -99

#abun1_full1 <- -99
#abun0_full1 <- -99
#abun1_full0 <- -99

#compare <- data.frame(species,ncell_abun,lat_abun,lon_abun,ncell_full,lat_full,lon_full,abun1_full1,abun0_full1,abun1_full0)
##########################################

for(i in 1:length(ordsums$species.code)){
  
  abun <- raster(paste0("ABUN/",ordsums$species.code[i],"_ABUN.asc"))
  full <- raster(paste0("FULL/",ordsums$species.code[i],"_FULL.asc"))

  comp <- as.data.frame(stack(abun,full)) 
  colnames(comp) <- c("abun","full")
  comp <- comp[!is.na(comp$abun),]
  
  ordsums$abun1_full1[i] <- length(comp$abun[comp$abun == 1 & comp$full == 1 ])
  ordsums$abun0_full1[i] <- length(comp$abun[comp$abun == 0 & comp$full == 1 ])
  ordsums$abun1_full0[i] <- length(comp$abun[comp$abun == 1 & comp$full == 0 ])  
  
  print(i)
  #compare.i <- data.frame(species,auc,ncell_abun,lat_abun,lon_abun,ncell_full,lat_full,lon_full,abun1_full1,abun0_full1,abun1_full0)
  #compare <- rbind(compare,compare.i)
}

ordsums$full_in_common <- ordsums$abun1_full1/ordsums$range_size_full
ordsums$full_w_o_abun <- ordsums$abun0_full1/ordsums$range_size_full
ordsums$abun_in_common <- ordsums$abun1_full1/ordsums$range_size_abun
ordsums$abun_w_o_full <- ordsums$abun1_full0/ordsums$range_size_abun


ordsums$full_chk <- ordsums$full_in_common + ordsums$full_w_o_abun
summary(ordsums$full_chk)
ordsums$abun_chk <- ordsums$abun_in_common + ordsums$abun_w_o_full
summary(ordsums$abun_chk)

############ exploration ####################
summary(ordsums$full_in_common) 
summary(ordsums$full_w_o_abun)

summary(ordsums$abun_in_common) 
summary(ordsums$abun_w_o_full)


# for full models: 3.5% - 100%  (38.3%) in common
#                  0%   - 96.5% (61.7%) without 1

# for abun models: 70.8%- 100%  (95.9%) in common
#                  0%   - 29.2% (04.1%) without 1   

ordsums[ordsums$full_in_common < .1,] # RICO3 SEPU7 VIVI
ordsums[ordsums$abun_w_o_full > .20,] # ELRE4 RAFI SCTE TUFA VEFO

extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
proj4string(extentShape) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

abun <- raster(paste0("ABUN/","ELRE4","_ABUN.asc"))
full <- raster(paste0("FULL/","ELRE4","_FULL.asc"))

#par(mfrow = c(3,1))
plot(abun)
plot(full)
plot(extentShape, col="grey")

abunpts <- read.table(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/species/", "ELRE4", ".csv"), header=T, stringsAsFactors = F, sep=",")
fullpts <- read.table(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/species/", "ELRE4", ".csv"), header=T, stringsAsFactors = F, sep=",")
head(abunpts)
coordinates(abunpts) <- c(2,3)
coordinates(fullpts) <- c(2,3)
proj4string(abunpts) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
proj4string(fullpts) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
plot(fullpts, add=T)
plot(abunpts, pch=19, col="red", add=T)


## percent point w abundance data
for (i in 1:length(ordsums$species.code)){
  ordsums$no.pts.full[i] <- NROW(read.table(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/species/", ordsums$species.code[i], ".csv")))
  ordsums$no.pts.abun[i] <- NROW(read.table(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/species/", ordsums$species.code[i], ".csv")))
  
  print(i)
}


ordsums$prop_pts_w_abun <- ordsums$no.pts.abun/ordsums$no.pts.full
plot(ordsums$full_in_common~ordsums$prop_pts_w_abun)
plot(ordsums$abun_in_common~ordsums$prop_pts_w_abun)

cor.test(ordsums$full_in_common,ordsums$prop_pts_w_abun,method="pearson") ##.22
cor.test(ordsums$abun_in_common,ordsums$prop_pts_w_abun,method="spearman") ##.17
cor.test(ordsums$full_in_common,ordsums$AUCtest_full,method="spearman")    ##.27
cor.test(ordsums$abun_in_common,ordsums$AUCtest_full,method="spearman")    ##.05
cor.test(ordsums$full_in_common,ordsums$AUCtest_abun,method="spearman")    ##-.24
cor.test(ordsums$abun_in_common,ordsums$AUCtest_abun,method="spearman")    ##.24
cor.test(ordsums$full_in_common,ordsums$no.pts,method="spearman")          ##.01
cor.test(ordsums$abun_in_common,ordsums$no.pts,method="spearman")          ##.02

hist(ordsums$prop_pts_w_abun)
ordsums[ordsums$prop_pts_w_abun >1,]



########################################################
## Create Species Richness Map 
######################################################


library(raster)
library(rgdal)
library(ordinal)
library(reshape2)
### use abundance models to predict within maxent (full) ranges

edd <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## eddmaps abundance data (thinned) with environmental data

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_FINAL.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## best fit models for each species
ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]

formulae <- data.frame(cbind(ordsums$species.code,ordsums$formu), stringsAsFactors = F)
str(formulae)
colnames(formulae) <- c("species","formula")

spp <- formulae$species[order(formulae$species)]
edd <- edd[edd$species %in% spp,]
length(unique(edd$species))

bio <- stack("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc" )
## assign proj4string to that of one of the raw files before clipping
proj4string(bio) <-  "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
## environmental variables
bioD <- as.data.frame(bio)
str(bioD)

setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/FULL")
## location of maxent ranges
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/")
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN")
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN")

for (i in 1:length(spp)){
  
  species <- edd[edd$species == spp[i],]
  species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3))
  
  M <- clm(formula(formulae$formula[formulae$species==spp[i]]), data = species)
  #AIC(M)
  
  range <- raster(paste0(spp[i],"_FULL.asc"))
  #plot(range)
  rangeD <- as.data.frame(range)
  colnames(rangeD) <- c("presence")
  
  p <- predict(M, bioD, type="class")
  p <- melt(p)
  p <- p[,1]
  
  p <- data.frame(cbind(p, rangeD$presence))
  colnames(p) <- c("abundance","presence")
  p$abundance_presence <- p$abundance*p$presence
  
  p$hi_abundance <- p$abundance_presence/3
  p$hi_abundance[p$hi_abundance < 1] <- 0
  
  
  abun <- raster(nrows=nrow(range),ncols=ncol(range),ext=extent(range),crs=proj4string(bio),vals=p$abundance_presence)
  hi_abun <- raster(nrows=nrow(range),ncols=ncol(range),ext=extent(range),crs=proj4string(bio),vals=p$hi_abundance)

  writeRaster(abun, paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/",spp[i],".asc"),overwrite=T)
  writeRaster(hi_abun,paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/",spp[i],".asc"),overwrite=T)
  
  print(i)
  
}

head(species)

  #summary(as.factor(p$hi_abundance))
  plot(abun)
  plot(hi_abun)
  extent(bio)==extent(range)
  ncell(bio)==ncell(range)
  nrow(bio)==nrow(range)
  ncol(bio)==ncol(range)
  
  head(rangeD)
  #summary(as.factor(rangeD$ACJA_FULL))

ordsums[ordsums$species.code=="ACJA",]

#(adapt for abundance sums once you have prediction maps from ordinal regression)
### Sum binary maps 

filenames = list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN",pattern="*.asc", full.names=TRUE)
ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)

out_file="C:/Users/Localadmin/Documents/Maxent_modeling/ORDINAL/richness_map_HI_ABUN.asc"
writeRaster(overlap, filename=out_file)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/ORDINAL/richness_map_HI_ABUN.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()

##### repeat for full points Maxent
filenames = list.files("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/FULL/",pattern="*.asc", full.names=TRUE)
ldf=stack(filenames)

overlap=calc(ldf, sum)

out_file="C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_FULL.asc"
writeRaster(overlap, filename=out_file)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_FULL.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()

### repeat for abundance points Maxent
filenames = list.files("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/ABUN/",pattern="*.asc", full.names=TRUE)
ldf=stack(filenames)

overlap2=calc(ldf, sum)

out_file="C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.asc"
writeRaster(overlap2, filename=out_file)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap2, col.regions=colors)
dev.off()




###############################################################################################################
### Hotspot maps
################################################################################################################



current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/ORDINAL/richness_map_HI_ABUN.asc")
quantile(current_rich$richness_map_HI_ABUN)

rich <- as.data.frame(current_rich)
rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN < 29] <- 0
rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN >= 29] <- 1

hi_abun_hotspot <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_HI_ABUN)
plot(hi_abun_hotspot)

####################################################
current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_FULL.asc")
quantile(current_rich$richness_map_FULL)

rich <- as.data.frame(current_rich)
head(rich)
rich$richness_map_FULL[rich$richness_map_FULL < 69] <- 0
rich$richness_map_FULL[rich$richness_map_FULL >= 69] <- 1


full_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_FULL)
plot(full_hotspot)

####################################################
current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.asc")
quantile(current_rich$richness_map_ABUN)

rich <- as.data.frame(current_rich)
head(rich)
rich$richness_map_ABUN[rich$richness_map_ABUN < 43] <- 0
rich$richness_map_ABUN[rich$richness_map_ABUN >= 43] <- 1


abun_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_ABUN)
plot(abun_hotspot)



############################################
comp <- as.data.frame(stack(abun_hotspot,full_hotspot))
head(comp)
comp$x[comp$layer.1==0 & comp$layer.2==0 ] <- 0 #neither
comp$x[comp$layer.1==1 & comp$layer.2==0 ] <- 1 #just abun
comp$x[comp$layer.1==1 & comp$layer.2==1 ] <- 3 #just full
comp$x[comp$layer.1==0 & comp$layer.2==1 ] <- 6 #both

length(comp$x[comp$x == 3])/length(comp$layer.1[comp$layer.1 == 1])
length(comp$x[comp$x == 3])/length(comp$layer.1[comp$layer.2 == 1])

comp_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp$x)
## based on richness maps here, adapt to summed abundance maps
plot(comp_hotspot)

##########################################
comp2 <- as.data.frame(stack(hi_abun_hotspot,full_hotspot))
head(comp2)
comp2$x[comp2$layer.1==0 & comp2$layer.2==0 ] <- 0 #neither
comp2$x[comp2$layer.1==1 & comp2$layer.2==0 ] <- 1 #just abun
comp2$x[comp2$layer.1==1 & comp2$layer.2==1 ] <- 3 #just full
comp2$x[comp2$layer.1==0 & comp2$layer.2==1 ] <- 6 #both

length(comp2$x[comp2$x == 3])/length(comp2$layer.1[comp2$layer.1 == 1])
length(comp2$x[comp2$x == 3])/length(comp2$layer.1[comp2$layer.2 == 1])
.5*(length(comp2$x[comp2$x == 3])/length(comp2$layer.1[comp2$layer.1 == 1]) + length(comp2$x[comp2$x == 3])/length(comp2$layer.1[comp2$layer.2 == 1]))
comp_hotspot2<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp2$x)
## based on richness maps here, adapt to summed abundance maps
plot(comp_hotspot2)
colors=rev(heat.colors(100))
spplot(comp_hotspot2, col=colors)
## 98% overlap


################## average proportion of abundant species
################## plot continuous deviation from that
#comp2$prop <- comp2$layer.1/comp2$layer.2
hi_abun_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/ORDINAL/richness_map_HI_ABUN.asc")
full_max_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_FULL.asc")

comp3 <- as.data.frame(stack(hi_abun_rich,full_max_rich))
head(comp3)
comp3$prop<- comp3$richness_map_HI_ABUN/comp3$richness_map_FULL

hist(comp3$prop)
mean(comp3$prop, na.rm=T)
median(comp3$prop, na.rm=T)
sd(comp3$prop, na.rm=T)

comp3$dev <- (comp3$prop - mean(comp3$prop,na.rm=T))/sd(comp3$prop, na.rm=T)
summary(comp3$dev)
hist(comp3$dev)
quantile(comp3$dev,na.rm=T)

comp3$sum[comp3$prop>=1] <- 2
comp3$sum[comp3$prop<1 & comp3$prop>-1] <- 1
comp3$sum[comp3$prop<=-1] <- 0

comp3$sum2[comp3$dev>=2] <- 2
comp3$sum2[comp3$dev<2 & comp3$dev>-2] <- 1
comp3$sum2[comp3$dev<=-2] <- 0


abun_prop <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp3$prop)
plot(abun_prop)

abun_deviation <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp3$dev)
#colors=rev(heat.colors(100))
#spplot(abun_deviation, col=colors)
plot(abun_deviation)

abun_deviation_sum<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp3$sum)
plot(abun_deviation_sum)
#colors=rev(heat.colors(100))
#spplot(abun_deviation_sum, col=colors)


abun_deviation_sum2<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp3$sum2)
plot(abun_deviation_sum2)

#########################################################
##### NOT YET ADAPTED
##### Current Richness
current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_FULL.asc")
#proj4string(current_rich) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs""+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
quantile(current_rich, probs=seq(0,1,0.05)) #get quantiles of richness

## set up the matrix required for the function reclassify
from=matrix(quantile(current_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to=matrix(quantile(current_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes=matrix(seq(0,1,0.05), nrow=21,ncol=1)
rcl=cbind(from,to,becomes)
colnames(rcl)=c("from","to","becomes")
rcl

for (i in 2:21){
  rcl[i,1]= rcl[i-1,2]+0.00000001
}
rcl[1,1]=0
rcl

# reclassify richness raster to quantiles
curr_rich_quant <- reclassify(current_rich, rcl)

## plot quantiles map
pdf(file="C:/Users/Localadmin/Documents/MaxEnt_modeling/quantilesFULL.pdf",width=11,height=8.5)
cuts=quantile(curr_rich_quant, probs=seq(0,1,0.05))
colors=rev(heat.colors(length(cuts)))
colors=matlab.like2(length(cuts))
legend=c(0,0.25,0.5,0.75,1)
legend_place=c(0,0.25,0.5,0.75,1)
print(spplot(curr_rich_quant,  col.regions=colors, at=cuts, colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) ))
dev.off()

### create binary hotspot map with top 25%
## create reclass table
from1=matrix(quantile(curr_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to1=matrix(quantile(curr_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes1=matrix(c(rep(0,16),rep(1,5)), nrow=21,ncol=1)
rcl1=cbind(from1,to1,becomes1)
colnames(rcl1)=c("from","to","becomes")
rcl1


for (i in 2:21){
  rcl1[i,1]= rcl1[i-1,2]+0.001
}
rcl1

#reclassify to hotspots
curr_rich_hot <- reclassify(curr_rich_quant, rcl1)

## plot hotspots map
pdf(file="C:/Users/Localadmin/Documents/MaxEnt_modeling/hotspotFULL.pdf",width=11,height=8.5)
cuts=c(0,0.5,1)
colors=c("grey70","red")
legend=c("","Hotspot")
legend_place=c(0,0.75)
print(spplot(curr_rich_hot,  col.regions=colors, at=cuts,colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) )) #
dev.off()

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/quantile_richness_FULL.asc"
writeRaster(curr_rich_quant, filename=out_file)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/hotspot_richness_FULL.asc"
writeRaster(curr_rich_hot, filename=out_file)

##################################
#### repeat for abun maxent models
##### Current Richness
current_rich=raster("Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.asc")
proj4string(current_rich) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs""+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
quantile(current_rich, probs=seq(0,1,0.05)) #get quantiles of richness

## set up the matrix required for the function reclassify
from=matrix(quantile(current_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to=matrix(quantile(current_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes=matrix(seq(0,1,0.05), nrow=21,ncol=1)
rcl=cbind(from,to,becomes)
colnames(rcl)=c("from","to","becomes")
rcl

for (i in 2:21){
  rcl[i,1]= rcl[i-1,2]+0.00000001
}
rcl[1,1]=0
rcl

# reclassify richness raster to quantiles
curr_rich_quant <- reclassify(current_rich, rcl)

## plot quantiles map
pdf(file="C:/Users/Localadmin/Documents/MaxEnt_modeling/quantilesABUN.pdf",width=11,height=8.5)
cuts=quantile(curr_rich_quant, probs=seq(0,1,0.05))
colors=rev(heat.colors(length(cuts)))
colors=matlab.like2(length(cuts))
legend=c(0,0.25,0.5,0.75,1)
legend_place=c(0,0.25,0.5,0.75,1)
print(spplot(curr_rich_quant,  col.regions=colors, at=cuts, colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) ))
dev.off()

### create binary hotspot map with top 25%
## create reclass table
from1=matrix(quantile(curr_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to1=matrix(quantile(curr_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes1=matrix(c(rep(0,16),rep(1,5)), nrow=21,ncol=1)
rcl1=cbind(from1,to1,becomes1)
colnames(rcl1)=c("from","to","becomes")
rcl1


for (i in 2:21){
  rcl1[i,1]= rcl1[i-1,2]+0.001
}
rcl1

#reclassify to hotspots
curr_rich_hot <- reclassify(curr_rich_quant, rcl1)

## plot hotspots map
pdf(file="C:/Users/Localadmin/Documents/MaxEnt_modeling/hotspotABUN.pdf",width=11,height=8.5)
cuts=c(0,0.5,1)
colors=c("grey70","red")
legend=c("","Hotspot")
legend_place=c(0,0.75)
print(spplot(curr_rich_hot,  col.regions=colors, at=cuts,colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) )) #
dev.off()

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/quantile_richness_ABUN.asc"
writeRaster(curr_rich_quant, filename=out_file)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/hotspot_richness_ABUN.asc"
writeRaster(curr_rich_hot, filename=out_file)



#########################################################################################
## MaxEnt: Sampling bias model, raw output
#########################################################################################
## this is one way to estimate sampling bias
## the model estimates where invasive plants are likely to have been sampled based on "ease of access" variables
## use ALL points (including duplicates within a climate grid cell) for ALL species in your analysis for this step
## we call this target group sampling

## set up model options, only non-climate predictors in the environmental directory, linear and quadratic features only
bias_TG=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=raw noremoveduplicates noaskoverwrite replicates=10 nothreshold nohinge noautofeature")
#bias_TG=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=raw noremoveduplicates noaskoverwrite replicates=10 nothreshold nohinge noautofeature -N ann_precip -N tmax_jul -N tmin_jan")

## define output directory
#output=paste("outputdirectory=","C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/TG_allPts",sep="")
output=paste("outputdirectory=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL",sep="")

## define directory of environmental data for fitting
#environmental=paste("environmentallayers=","C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/US_ASCIIs",sep="")
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/bias_ASCIIs",sep="")
## put in all of the environmental data

## define point location samples (all species in one .csv, give all the same name, e.g., bias)
## otherwise each species will be treated seperately (good for later, but want bias across ALL species here)
#samples=paste("samplesfile=","C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/IAS_bias_points_maxent.csv",sep="")
samples=paste("samplesfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/TG_allPts_FULL/full_pts_2_15_2018.csv",sep="")

## call model to run
system(paste(bias_TG,output,environmental,samples, "autorun"))

## fair warning: maxEnt dumps a ton of files into the output directory
## the output of this model will be your bias layer for the species-level models
## this is NOT what we did in the hotspots paper, so don't get confused

#########################################################################################
## MaxEnt: Current Climate, species richness bias, logistic output
#########################################################################################

## set up model options, linear and quadratic features only, use sampling bias surface
spp_bias_log=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature")

## define output directory
output=paste("outputdirectory=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/Bias_Logistic_FULL",sep="")

## define directory of environmental data for fitting
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/bias_ASCIIs",sep="")

## define point location samples (all species in one .csv)
## must be set up with three columns (species, lon, lat, in that order)
samples=paste("samplesfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/TG_allPts_FULL/full_pts_2_15_2018.csv",sep="")

## define sampling bias surface
## this should be the average ascii file from your bias model above
bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/RENAME.asc",sep="")

## call model to run
system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))






#########################################################################################
## REPEAT FOR JUST ABUN POINTS
#########################################################################################

## RAW SAMPLING BIAS
## this is one way to estimate sampling bias
## the model estimates where invasive plants are likely to have been sampled based on "ease of access" variables
## use ALL points (including duplicates within a climate grid cell) for ALL species in your analysis for this step
## we call this target group sampling

## set up model options, only non-climate predictors in the environmental directory, linear and quadratic features only
bias_TG=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=raw noremoveduplicates noaskoverwrite replicates=10 nothreshold nohinge noautofeature")
#bias_TG=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=raw noremoveduplicates noaskoverwrite replicates=10 nothreshold nohinge noautofeature -N ann_precip -N tmax_jul -N tmin_jan")

## define output directory
#output=paste("outputdirectory=","C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/TG_allPts",sep="")
output=paste("outputdirectory=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN",sep="")

## define directory of environmental data for fitting
#environmental=paste("environmentallayers=","C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/US_ASCIIs",sep="")
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/bias_ASCIIs",sep="")
## put in all of the environmental data

## define point location samples (all species in one .csv, give all the same name, e.g., bias)
## otherwise each species will be treated seperately (good for later, but want bias across ALL species here)
#samples=paste("samplesfile=","C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/IAS_bias_points_maxent.csv",sep="")
samples=paste("samplesfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/TG_allPts_ABUN/full_abundance_pts_2_15_2018.csv",sep="")

## call model to run
system(paste(bias_TG,output,environmental,samples, "autorun"))

## fair warning: maxEnt dumps a ton of files into the output directory
## the output of this model will be your bias layer for the species-level models
## this is NOT what we did in the hotspots paper, so don't get confused

#########################################################################################
## MaxEnt: Current Climate, species richness bias, logistic output
#########################################################################################

## set up model options, linear and quadratic features only, use sampling bias surface
spp_bias_log=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature")
#"java -jar ",maxent.location, " nowarnings noprefixes -E responsecurveS jackknife outputformat=raw      noremoveduplicates noaskoverwrite replicates=10 nothreshold nohinge noautofeature"
## define output directory
output=paste("outputdirectory=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/Bias_Logistic_ABUN",sep="")

## define directory of environmental data for fitting
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/bias_ASCIIs",sep="")

## define point location samples (all species in one .csv)
## must be set up with three columns (species, lon, lat, in that order)
samples=paste("samplesfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/TG_allPts_ABUN/full_abundance_pts_2_15_2018.csv",sep="")

## define sampling bias surface
## this should be the average ascii file from your bias model above
bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/RENAME.asc",sep="")

## call model to run
system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))


#########################################################################################
## MaxEnt: Current Climate, species richness bias, logistic output
#########################################################################################

## set up model options, linear and quadratic features only, use sampling bias surface
spp_bias_log=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature")

## define output directory
output=paste("outputdirectory=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/Bias_Logistic_FULL",sep="")

## define directory of environmental data for fitting
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/bias_ASCIIs",sep="")

## define point location samples (all species in one .csv)
## must be set up with three columns (species, lon, lat, in that order)
samples=paste("samplesfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/TG_allPts_FULL/full_pts_2_15_2018.csv",sep="")

## define sampling bias surface
## this should be the average ascii file from your bias model above
bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/RENAME.asc",sep="")

## call model to run
system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))




















####################################################################################
## Find 95% threshold for each species AND get training/test sample sizes
###################################################################################
## maxent output represents relative occurence rate (ROR)
## we want a binary map, so we need to know what ROR value to use to define "suitable" vs "unsuitable"
## all maxent surfaces must sum to 1, so the threshold value can vary from species to species
## therefore we need to calculate the threshold by species (below)

## create list of modeled species
d=read.table("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/IAS_occurences_final_analysis.csv",sep=",",header=T)
spp=as.data.frame(unique(d$PLANT_CODE))
colnames(spp)="PLANT_CODE"
spp$PLANT_CODE=as.character(spp$PLANT_CODE)
head(spp)
str(spp)
dim(spp)

#create dataframe to hold 95% MTP logistic threshold values for each run
spp_val=data.frame("Species",15,0,"CtySppRich_Bias",100000, 100000,stringsAsFactors=F) 
colnames(spp_val)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n","train_n")
#spp_val$Species=as.character(spp_val$Species)
#spp_val$Bias=as.character(spp_val$Bias)

## get all output files (to deal with any skipped spp and variable run number)
setwd("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/CtyRich_Logistic")
a=list.files(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE)

for (i in unique(spp$PLANT_CODE)){        #loop for each species
  for (j in 0:9){                        #loop for each iteration within a species
    
    spp_iter=paste(i,"_",j,"_samplePredictions.csv",sep="")
    
    if(spp_iter %in% a == TRUE) {
      file=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/CtyRich_Logistic/",i,"_",j,"_samplePredictions.csv",sep="") #iteration j point-level model output 
      
      d=read.table(file, sep=",",header=T)
      d$run=j
      dt=d[d$test.or.train == "train",]  #subset iteration data to training data only
      dtest=d[d$test.or.train == "test",]  #subset iteration data to testing data only
      
      samp=dim(dt)
      samp1=dim(dtest)
      
      train_n=samp[1] #get training sample size for this iteration
      test_n=samp1[1]  #get test sample size for this iteration
      
      tp=round(0.95*samp[1], 0) #number of points for 95% threshold in each iteration
      
      dt=dt[order(dt$Logistic.prediction, decreasing=T),] #order the logistic output largest to smallest
      
      sp_runs=data.frame(i, j, dt$Logistic.prediction[tp], "CtySppRich_Bias", test_n, train_n) #create dataframe with species, run, logistic value of the point that bounds the 95% percentile, bias model, test and train n
      colnames(sp_runs)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n", "train_n")
      
      spp_val=rbind(spp_val,sp_runs) 
    }
    
    print(j)
  }
  print(i)
}

## trim junk row from spp_val
head(spp_val)
tail(spp_val)
dim(spp_val)
spp_val1=spp_val[2:8900,]

## calculate mean threshold per species per range from spp_val1
thresh=as.data.frame(tapply(spp_val1$X95_MTP_thresh_logistic, list(spp_val1$Species, spp_val1$Bias), mean))
thresh$Species=rownames(thresh)
head(thresh)
dim(thresh)

test=as.data.frame(tapply(spp_val1$test_n, list(spp_val1$Species, spp_val1$Bias), mean))
test$Species=rownames(test)
colnames(test)=c("test_n","Species")
head(test)
dim(test)

train=as.data.frame(tapply(spp_val1$train_n, list(spp_val1$Species, spp_val1$Bias), mean))
train$Species=rownames(train)
colnames(train)=c("train_n","Species")
head(train)
dim(train)

samples=cbind(test, train)
head(samples)
colnames(samples)=c("test_n","Species","train_n","Species1")
samples$match[samples$Species == samples$Species1]="TRUE"
table(samples$match)

samples1=subset(samples, select=c("Species","test_n","train_n"))
samples1$total_n=samples1$test_n + samples1$train_n
samples1$Bias="CtySppRich"
samples1$Climate="Current"
head(samples1)
dim(samples1)

write.table(thresh, file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/IAS_95Perc_thresholds_CtySppRichBias_Current_FINAL.csv",sep=",",row.names=F)
write.table(samples1,file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/IAS_Model_SampleSize_Current_FINAL.csv",sep=",",row.names=F) 

#########################################################################
### Create binary maps for each species
#########################################################################

## use the thresholds above to go from ROR to binary

### create dataset to hold suitable cell counts for each species 
range_size=data.frame(-100,-100,"NA","NA")
colnames(range_size)=c("Var1","Freq","spp","bias")
range_size

###### create binary rasters based on 95% MTP

for (i in unique(thresh$Species)){  
  file=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/CtyRich_Logistic/",i,"_avg.asc",sep="")
  d1=readAsciiGrid(file) 
  
  d1@data[d1@data >= thresh$CtySppRich_Bias[thresh$Species==i]]=1
  d1@data[d1@data < thresh$CtySppRich_Bias[thresh$Species==i]]=0
  
  out_file=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/binary/sppRichCtyBias_current",i,"_SppRichBias_FivePercMTP.asc",sep="")
  write.asciigrid(d1, out_file)
  
  tab=as.data.frame(table(d1@data))
  tab$spp=i
  tab$bias="CtySppRich_Bias"
  
  range_size=rbind(range_size,tab)
  
  #image(d1, col=c("grey","blue"))
  
  rm(d1)
  gc()
  print(i)
}

head(range_size)
dim(range_size)
range_size1=range_size[2:1793,]

range_size2=range_size1[range_size1$Var1 ==1,]
head(range_size2)
dim(range_size2)






##################################################################
### Create Projections with Future Climate (example with 1 GCM)
##################################################################
##################################################################
### BCC-CSM1-1 model, RCP4.5

environmental="C:/Users/mwone/Google Drive/Ordinal/Occurence/Climate_Data/Analysis/rcp45/BCC_CSM1_1 "

for (i in unique(spp$PLANT_CODE))  { # loop returns warnings if a species has <10 iterations, but runs fine on however many lambdas files the species *does* have
  for (j in 0:9){
    file1=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/CtyRich_Logistic/",i,"_",j,".lambdas ",sep="") #iteration j lambda file 
    output=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/Future/rcp45/BCC_CSM1_1/",i,"_",j,sep="")
    projection <- paste0("java -cp ",maxent.location, " density.Project ", file1, environmental, output, sep="") # enter any args here
    system(paste(projection, "autorun"))     
    print(j)
    print(i)
    
  }
}

#################### Create Mean Future Projection from 10 iterations

for (i in unique(spp$PLANT_CODE)){        #loop for each species
  
  a1=list.files("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/Future/rcp45/BCC_CSM1_1",pattern=paste(i,"*", sep=""), full.names=TRUE)
  
  proj=stack(a1)
  
  avg=calc(proj, mean)   
  
  out_file=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/Future/rcp45/BCC_CSM1_1/Mean/",i,"_BCC_CSM1_1_rcp45.asc", sep="")
  writeRaster(avg, filename=out_file) 
  
  rm(a1, proj, avg)
  gc()
  
  print(i)
}


####################################################################
### Create binary maps and calculate range size (future, 1 GCM)
####################################################################

thresh=read.table("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/IAS_95Perc_thresholds_CtySppRichBias_Current_FINAL.csv",sep=",",header=T)

### create dataset to hold suitable cell counts for each species 
range_size=data.frame(-100,-100,"NA","NA","NA")
colnames(range_size)=c("Var1","Freq","spp","model","RCP")
range_size

###### create binary rasters based on 95% MTP

for (i in unique(thresh$Species)){  
  file=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/Future/rcp45/BCC_CSM1_1/Mean/",i,"_BCC_CSM1_1_rcp45.asc",sep="")
  d1=readAsciiGrid(file) 
  
  d1@data[d1@data >= thresh$CtySppRich_Bias[thresh$Species==i]]=1
  d1@data[d1@data < thresh$CtySppRich_Bias[thresh$Species==i]]=0
  
  out_file=paste("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/binary/BCC_CSM1_1_rcp45/",i,"_BCC_CSM1_1_rcp45.asc",sep="")
  write.asciigrid(d1, out_file)
  
  tab=as.data.frame(table(d1@data))
  tab$spp=i
  tab$model="BCC_CSM1_1"
  tab$RCP="rcp_4_5"
  
  range_size=rbind(range_size,tab)
  
  #image(d1, col=c("grey","blue"))
  
  rm(d1)
  gc()
  print(i)
}

head(range_size)
dim(range_size)
range_size1=range_size[2:1781,]

range_size2=range_size1[range_size1$Var1 ==1,]
head(range_size2)
dim(range_size2)

########################################################
## Create Species Richness Map (adapt for abundance sums once you have prediction maps from ordinal regression)
######################################################

### Sum binary maps 
filenames = list.files("C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/FINAL_OUTPUT/binary/sppRichCtyBias_current",pattern="*.asc", full.names=TRUE)
ldf=stack(filenames)

overlap=calc(ldf, sum)

out_file="C:/Users/mwone/Google Drive/Ordinal/Occurence/Hotspots_output/Richness_CtyRichBias_FivePercMTP_Current.asc"
writeRaster(overlap, filename=out_file)

pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/SppRichness_CtyRichBias_5PercMTP_Current.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()

##############################################################################################
## Calculate and Plot RCP 4.5 richness mean and sd from future model projections
##############################################################################################

## this calculates the mean richness across 13 future GCM maps of species richness
## you can adapt this code to find the mean abundance predictions across your GCMs
## assumes each GCM abundance map is ascii format

filenames = list.files("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/RCP45_recalc",pattern="*_rcp45.asc", full.names=TRUE)
filenames
length(filenames)

ldf=stack(filenames)

mean_45=calc(ldf, mean)
sd_45=calc(ldf, sd)

out_file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Ensemble Maps/Richness_Mean_RCP45_recalc.asc"
outfile1="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Ensemble Maps/Richness_SD_RCP45_recalc.asc"
writeRaster(mean_45, filename=out_file)
writeRaster(sd_45, filename=outfile1)

ls()
rm(ldf, filenames, mean_45, sd_45, outfile1, out_file)
gc()

####################################################################################
## Plot RCP 4.5 Ensemble Mean and SD Richness 
###################################################################################

## plots the results of the future ensemble calculations above

mean_45=readAsciiGrid("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Ensemble Maps/Richness_Mean_RCP45_recalc.asc") 
proj4string(mean_45) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

sd_45=readAsciiGrid("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Ensemble Maps/Richness_SD_RCP45_recalc.asc") 
proj4string(sd_45) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

srg=readOGR("C:/Users/Jenica/Dropbox/National Invasives SDMs/Intermediate Data/Ecoregions_L48_wgs84.shp", layer="Ecoregions_L48_wgs84")
proj4string(srg)

#setEPS()
#postscript("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/SppRichness_RCP45_Ensemble_Mean_recalc.eps",width=3.5,height=3)
pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/SppRichness_RCP45_Ensemble_Mean_recalc.pdf",width=11,height=8.5)
plot_poly=list("sp.lines", as(srg, "SpatialLines"), col="grey20")
colors=rev(heat.colors(100))
spplot(mean_45, col.regions=colors, sp.layout=list(plot_poly))
dev.off()

setEPS()
postscript("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Figures/PNAS/SppRichness_SD_RCP45_recalc.eps",width=3.5,height=3)
#pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/SppRichness_SD_RCP45_recalc.pdf",width=11,height=8.5)
plot_poly=list("sp.lines", as(srg, "SpatialLines"), col="grey20")
colors=blue2red(20)
spplot(sd_45, col.regions=colors, sp.layout=list(plot_poly))
dev.off()

###############################################################################################################
### Hotspot maps
################################################################################################################

## based on richness maps here, adapt to summed abundance maps

##### Current Richness
current_rich=raster("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Richness_CtyRichBias_FivePercMTP_Current_recalc.asc")
proj4string(current_rich) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

quantile(current_rich, probs=seq(0,1,0.05)) #get quantiles of richness

## set up the matrix required for the function reclassify
from=matrix(quantile(current_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to=matrix(quantile(current_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes=matrix(seq(0,1,0.05), nrow=21,ncol=1)
rcl=cbind(from,to,becomes)
colnames(rcl)=c("from","to","becomes")
rcl

for (i in 2:21){
  rcl[i,1]= rcl[i-1,2]+0.00000001
}
rcl[1,1]=0
rcl

# reclassify richness raster to quantiles
curr_rich_quant <- reclassify(current_rich, rcl)

## plot quantiles map
pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Hotspot_quantiles_Current_recalc.pdf",width=11,height=8.5)
cuts=quantile(curr_rich_quant, probs=seq(0,1,0.05))
colors=rev(heat.colors(length(cuts)))
colors=matlab.like2(length(cuts))
legend=c(0,0.25,0.5,0.75,1)
legend_place=c(0,0.25,0.5,0.75,1)
print(spplot(curr_rich_quant,  col.regions=colors, at=cuts, colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) ))
dev.off()

### create binary hotspot map with top 25%
## create reclass table
from1=matrix(quantile(curr_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to1=matrix(quantile(curr_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes1=matrix(c(rep(0,16),rep(1,5)), nrow=21,ncol=1)
rcl1=cbind(from1,to1,becomes1)
colnames(rcl1)=c("from","to","becomes")
rcl1


for (i in 2:21){
  rcl1[i,1]= rcl1[i-1,2]+0.001
}
rcl1

#reclassify to hotspots
curr_rich_hot <- reclassify(curr_rich_quant, rcl1)

## plot hotspots map
pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Hotspot_binary_Current_recalc.pdf",width=11,height=8.5)
cuts=c(0,0.5,1)
colors=c("grey70","red")
legend=c("","Hotspot")
legend_place=c(0,0.75)
print(spplot(curr_rich_hot,  col.regions=colors, at=cuts,colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) )) #
dev.off()

out_file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Richness_quantiles_current_recalc.asc"
writeRaster(curr_rich_quant, filename=out_file)

out_file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Richness_hotspots_current_recalc.asc"
writeRaster(curr_rich_hot, filename=out_file)

######### Future Ensemble Average 
## same process with future ensemble average richness
fut_rich=raster("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Ensemble Maps/Richness_Mean_RCP45_recalc.asc")
proj4string(fut_rich) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

quantile(fut_rich, probs=seq(0,1,0.05))

from=matrix(quantile(fut_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to=matrix(quantile(fut_rich, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes=matrix(seq(0,1,0.05), nrow=21,ncol=1)
rcl=cbind(from,to,becomes)
colnames(rcl)=c("from","to","becomes")
rcl

for (i in 2:21){
  rcl[i,1]= rcl[i-1,2]+0.00000001
}
rcl[1,1]=0
rcl

fut_rich_quant <- reclassify(fut_rich, rcl)

pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Hotspot_quantiles_future_recalc.pdf",width=11,height=8.5)
cuts=quantile(fut_rich_quant, probs=seq(0,1,0.05))
colors=rev(heat.colors(length(cuts)))
colors=matlab.like2(length(cuts))
legend=c(0,0.25,0.5,0.75,1)
legend_place=c(0,0.25,0.5,0.75,1)
print(spplot(fut_rich_quant,  col.regions=colors, at=cuts, colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) ))
dev.off()

### create binary hotspot map with top 25%
from1=matrix(quantile(fut_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
to1=matrix(quantile(fut_rich_quant, probs=seq(0,1,0.05)), nrow=21,ncol=1)
becomes1=matrix(c(rep(0,16),rep(1,5)), nrow=21,ncol=1)
rcl1=cbind(from1,to1,becomes1)
colnames(rcl1)=c("from","to","becomes")
rcl1


for (i in 2:21){
  rcl1[i,1]= rcl1[i-1,2]+0.001
}
rcl1

fut_rich_hot <- reclassify(fut_rich_quant, rcl1)

pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Hotspot_binary_future_recalc.pdf",width=11,height=8.5)
cuts=c(0,0.5,1)
colors=c("grey70","red")
legend=c("","Hotspot")
legend_place=c(0,0.75)
print(spplot(fut_rich_hot,  col.regions=colors, at=cuts,colorkey=list(width=2,labels=list(at=legend_place,labels=legend)) )) #
dev.off()

out_file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Richness_quantiles_future_recalc.asc"
writeRaster(fut_rich_quant, filename=out_file)

out_file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Richness_hotspots_future_recalc.asc"
writeRaster(fut_rich_hot, filename=out_file)


##### Hotspot change

fut_rich_hot1=fut_rich_hot
fut_rich_hot1[fut_rich_hot1 >0]=2

unique(fut_rich_hot1) #0 and 2

hot_change=fut_rich_hot1 + curr_rich_hot

#setEPS()
#postscript("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Hotspot_change_recalc.eps",width=3.5,height=3)
pdf(file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/Hotspot_change_recalc.pdf",width=11,height=8.5) 
# 0 means pixel NOT hotspot in either time period
# 1 means pixel hotspot in current climate only
# 2 means pixel hotspot in future climate only
# 3 means pixel hotspot in current and future 

plot_poly=list("sp.lines", as(srg, "SpatialLines"), col="grey20")
cuts=c(-.5, 0.5, 1.5, 2.5, 3.5)
#legend=c("Non-hotspot","Current climate hotspot","Future climate hotspot","Current & future hotspot")
legend=c("None","Current","Future","Both")    

print(spplot(hot_change,col.regions=c("grey90","royalblue","red", "orange"), at=cuts, sp.layout=list(plot_poly),
             colorkey=list(width=2,labels=list(at=c(0,1,2,3),labels=legend)))) 

dev.off()

out_file="C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Hotspots_change_recalc.asc"
writeRaster(hot_change, filename=out_file)

hot_change1=raster("C:/Users/Jenica/Dropbox/National Invasives SDMs/Analysis/richness_asciis/Hotspots_change_recalc.asc")

freq(hot_change1)


