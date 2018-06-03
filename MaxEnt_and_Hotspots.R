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

sp.list <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/spForms3_4_2018.csv", header = T, sep = ",", 
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
head(sp.list)
#sp.list<- sp.list[order(sp.list$code),]
#head(sp.list)
spp_bias_log = paste0("java -jar ",maxent.location, 
                      " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 noproduct nothreshold nohinge writeplotdata noautofeature biastype=3")
## linear + quadratic terms, thinning, logistic, 10 reps, bias surface

for(i in 1:length(sp.list$code)){
   ################# all vs abun points code ############################
#  ## set up model options, linear and quadratic features only, use sampling bias surface
#  spp_bias_log=paste0("java -jar ",maxent.location, 
#                      " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature ",
#                      sp.list$ignore[i], " biastype=3")
#  ##spp_bias_log=paste0("java -jar ",maxent.location, " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 nothreshold nohinge writeplotdata noautofeature -N pop_2_5_us -N road_2_5_us biastype=3")
#  
#   ### ALL SPECIES ###
#   ## define point location samples (all species in one .csv), must be set up with three columns (species, lon, lat, in that order)
#   samples=paste("samplesfile=", paste("C:/Users/Localadmin/Documents/Maxent_modeling/FULL/species", 
#                                      paste(sp.list$code[i], "csv", sep="."), sep="/"), sep="")
#   ## define output directory
#   directory<- paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/FINAL_OUTPUT_FULL", sp.list$code[i], sep="/")
#   dir.create(directory)
#   output=paste("outputdirectory=", directory, sep="")
#   
#   bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc",sep="")
#   ## average bias ascii
#   
#   ## call model to run
#   system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))
#  
#  
#  
#   ### REPEAT FOR ONLY ABUNDANCE POINTS ###
#   ## define point location samples (all species in one .csv), must be set up with three columns (species, lon, lat, in that order)
#   samples=paste("samplesfile=", paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/species", 
#                                       paste(sp.list$code[i], "csv", sep="."), sep="/"), sep="")
#   ## define output directory
#   directory <- paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/FINAL_OUTPUT_ABUN", sp.list$code[i], sep="/")
#   dir.create(directory)
#   output=paste("outputdirectory=", directory, sep="")
#   
#   bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/BIAS_OUTPUT_ABUN/BIAS_avgABUN.asc",sep="")
#   ## average bias ascii
#   
#   ## call model to run
#   system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))
 
  
   ### all variables ############ 
    
  ## define point location samples (all species in one .csv), must be set up with three columns (species, lon, lat, in that order)
  samples=paste("samplesfile=", paste("C:/Users/Localadmin/Documents/Maxent_modeling/species", 
                                      paste(sp.list$code[i], "csv", sep="."), sep="/"), sep="")
  ## define output directory
  directory<- paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/FINAL_OUTPUT", sp.list$code[i], sep="/")
  dir.create(directory)
  output=paste("outputdirectory=", directory, sep="")
  bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc",sep="")
  ## average bias ascii
  ## call model to run
  system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))
  
  
  print(i)
  write.csv(i,"progressMAXENT.csv",row.names=F)
}

rm(sp.list,bias,directory,environmental,i,maxent.location,output,samples,spp_bias_log)
gc()

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

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_5_2018.csv", header = T, sep = ",",  
                                 quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

spp=ordsums$species.code
#colnames(spp)="PLANT_CODE"
#spp$PLANT_CODE=as.character(spp$PLANT_CODE)
#head(spp)
#str(spp)
#dim(spp)


#dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis")
#dataset_type <- c("FULL","ABUN","FULL_all_vars")
#for (type in dataset_type){

  
  #create dataframe to hold 95% MTP logistic threshold values for each run
  spp_val=data.frame("Species",15,0,"CtySppRich_Bias",100000, 100000,stringsAsFactors=F) 
  colnames(spp_val)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n","train_n")
  #spp_val$Species=as.character(spp_val$Species)
  #spp_val$Bias=as.character(spp_val$Bias)
  
  ## get all output files (to deal with any skipped spp and variable run number) C:/Users/Localadmin/Documents/MaxEnt_modeling/ABUN/FINAL_OUTPUT_ABUN
  setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/FINAL_OUTPUT/")
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
    ordsums$AUC[ordsums$species.code == i] <- aucs$Test.AUC[11]
    
    #### for past model runs ####

    #if(type=="FULL"){
    #   ordsums$AUCtest_full[ordsums$species.code == i] <- aucs$Test.AUC[11]
    #   ordsums$AUCtrain_full[ordsums$species.code == i] <- aucs$Training.AUC[11]
    #}
    #
    #if(type=="ABUN"){
    #   ordsums$AUCtest_abun[ordsums$species.code == i] <- aucs$Test.AUC[11]
    #   ordsums$AUCtrain_abun[ordsums$species.code == i] <- aucs$Training.AUC[11]
    #} 
    #if(type=="FULL_all_vars"){
    #   ordsums$AUCtest_allvars[ordsums$species.code == i] <- aucs$Test.AUC[11]
    #   ordsums$AUCtrain_allvars[ordsums$species.code == i] <- aucs$Training.AUC[11]
    #}  
    #######################
    
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
  
  write.table(thresh, file= "C:/Users/Localadmin/Documents/Maxent_modeling/summaries/thresholds.csv" ,sep=",",row.names=F)
  write.table(samples1,file= "C:/Users/Localadmin/Documents/Maxent_modeling/summaries/sample_sizes.csv" ,sep=",",row.names=F) 
  

  
  #########################################################################
  ### Create binary maps for each species
  #########################################################################
  
  ## use the thresholds above to go from ROR to binary
  
  ### create dataset to hold suitable cell counts for each species 
  #range_size=data.frame(-100,-100,"NA","NA")
  #colnames(range_size)=c("Var1","Freq","spp","bias")
  #range_size
  
  #dir.create("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/")
  
  
  ###### create binary rasters based on 95% MTP
  
  for (i in unique(thresh$Species)){  
    file=paste(i,"/",i,"_avg.asc",sep="")
    d1=readAsciiGrid(file) 
    
    d1@data[d1@data >= thresh$CtySppRich_Bias[thresh$Species==i]]=1
    d1@data[d1@data < thresh$CtySppRich_Bias[thresh$Species==i]]=0
    
    out_file=paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",i,".asc",sep="")
    write.asciigrid(d1, out_file)
    

    
    
    ##### for past model runs ###############
    #  tab=as.data.frame(table(d1@data))
    #  tab$spp=i
    #  tab$bias="CtySppRich_Bias"
    #  range_size=rbind(range_size,tab)
    #  image(d1, col=c("grey","blue"))
    #  if(type=="FULL"){
    #    ordsums$range_size_full[ordsums$species.code==i] <- tab$Freq[tab$Var1==1]
    #    ordsums$range_check_full[ordsums$species.code==i] <- sum(tab$Freq)
    #  }
    #  if(type=="ABUN") {
    #    ordsums$range_size_abun[ordsums$species.code==i] <- tab$Freq[tab$Var1==1]
    #    ordsums$range_check_abun[ordsums$species.code==i] <- sum(tab$Freq)
    #  }
    #  if(type=="FULL_all_vars") {
    #    ordsums$range_size_allvars[ordsums$species.code==i] <- tab$Freq[tab$Var1==1]
    #    ordsums$range_check_allvars[ordsums$species.code==i] <- sum(tab$Freq)
    #  }
    #########################  
    
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
  
  #print(type)
  #}

summary(ordsums$AUC)
length(ordsums$AUC[ordsums$AUC > .7])


###################################
## EXPLORATION OF MODEL PERFORMANCE
###################################
plot(ordsums$AUCtest_full~ordsums$AUCtest_abun, xlab="AUC abundance pts", ylab="AUC full pts" )
abline(0,1, col="blue", lwd=3)
abline(h=0.7, col="grey")
abline(v=0.7, col="grey")

par(mfrow=c(2,1))
#hist(ordsums$AUCtest_abun, xlim=c(.45,1.0), main="", xlab="AUC abundance pts")
hist(ordsums$AUCtest_full, xlim=c(.45,1.0), main="", xlab="AUC abun vars")
hist(ordsums$AUCtest_allvars, xlim=c(.45,1.0), main="", xlab="AUC all pts")
length(ordsums$AUCtest_allvars[ordsums$AUCtest_allvars >.7 & ordsums$AUCtest_full <.7]) ## using all vars saves 27 species!
length(ordsums$AUCtest_allvars[ordsums$AUCtest_allvars >.7 & ordsums$kappa > 0 & ordsums$kappaP < 0.05 & ordsums$gn10 > .9]) ## 86
length(ordsums$AUCtest_allvars[ordsums$AUCtest_allvars >.7 & ordsums$kappa > 0 & ordsums$kappaP < 0.05 & ordsums$gn5 > .95]) ## 55

summary(ordsums$range_check_abun)
summary(ordsums$range_check_full)
summary(ordsums$range_size_abun)
summary(ordsums$range_size_full)
summary(ordsums$range_size_abun/ordsums$range_size_full)
hist(ordsums$range_size_abun/ordsums$range_size_full)
length(ordsums$species.code[ordsums$range_size_abun/ordsums$range_size_full > 1])
ordsums[ordsums$range_size_abun/ordsums$range_size_full > 1,]


write.csv(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_3_23_2018.csv", row.names=F)



####################################################
##### compare Maxent Full VS Abun w/in and across sp   ~~~ NO LONGER USED
####################################################

setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis")

#ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_maxent.csv", header = T, sep = ",",  
#                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

#ordsums <- ordsums[ordsums$kappa > 0 & ordsums$kappaP <0.05,]
summary(ordsums$AUCtest_abun)


library(raster)
library(rgdal)


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

############ exploration
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
hist(ordsums$prop_pts_w_abun)


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




### exploraton 3/20
length(ordsums$species.code[ordsums$AUCtest_full >.7]) ### 121/155
length(ordsums$species.code[ordsums$AUCtest_abun >.7]) ### 148/155
length(ordsums$species.code[ordsums$AUCtest_full >.75]) ### 100/155

length(ordsums$species.code[ordsums$kappa > 0 & ordsums$kappaP <0.05]) ## 125/155
length(ordsums$species.code[ordsums$kappa > 0.2]) ## 67


length(ordsums$species.code[ordsums$gn10 > .9]) ### 117
length(ordsums$species.code[ordsums$gn5 > .95]) ### 78

plot(ordsums$gn10, ordsums$full_in_common)
cor.test(ordsums$gn10, ordsums$full_in_common, method="spearman") ## .6
length(ordsums$species.code[ordsums$full_in_common > .5])



summary(ordsums$prop_pts_w_abun)
summary(ordsums$full_in_common) 
hist(ordsums$prop_pts_w_abun)
hist(ordsums$full_in_common)

cor.test(ordsums$prop_pts_w_abun, ordsums$full_in_common, method="spearman") # .3
cor.test(ordsums$AUCtest_full, ordsums$full_in_common, method="spearman") # .3

cor.test(ordsums$AUCtest_full, ordsums$no.vars, method="spearman") # .08
cor.test(ordsums$AUCtest_full, ordsums$prop_pts_w_abun, method="spearman") # .5
cor.test(ordsums$gn10, ordsums$AUCtest_full, method="spearman") # .4
125/155
117/155
148/155
121/155
length(ordsums$species.code[ordsums$gn10 > .9 & ordsums$AUCtest_full >.7  & ordsums$kappa > 0 & ordsums$kappaP <0.05]) ### 77
length(ordsums$species.code[ordsums$gn5 > .95 & ordsums$AUCtest_full >.75 & ordsums$kappa > 0.2]) ### 17


####################################################
##### compare Maxent Full VS All w/in and across sp    ~~~ NO LONGER USED
####################################################

for(i in 1:length(ordsums$species.code)){
  
  all <- raster(paste0("FULL_all_vars/",ordsums$species.code[i],"_FULL_all_vars.asc"))
  full <- raster(paste0("FULL/",ordsums$species.code[i],"_FULL.asc"))
  
  comp <- as.data.frame(stack(all,full)) 
  colnames(comp) <- c("all","full")
  comp <- comp[!is.na(comp$a),]
  
  ordsums$all1_full1[i] <- length(comp$all[comp$all == 1 & comp$full == 1 ])
  ordsums$all0_full1[i] <- length(comp$all[comp$all == 0 & comp$full == 1 ])
  ordsums$all1_full0[i] <- length(comp$all[comp$all == 1 & comp$full == 0 ])  
  
  ordsums$range_size_all[i] <- length(comp$all[comp$all==1])
  
  print(i)
  #compare.i <- data.frame(species,auc,ncell_abun,lat_abun,lon_abun,ncell_full,lat_full,lon_full,abun1_full1,abun0_full1,abun1_full0)
  #compare <- rbind(compare,compare.i)
}

ordsums$all_in_common <- ordsums$all1_full1/ordsums$range_size_all
ordsums$all_w_o_full<- ordsums$all1_full0/ordsums$range_size_all
ordsums$full_in_common2 <- ordsums$all1_full1/ordsums$range_size_full
ordsums$full_w_o_all <- ordsums$all0_full1/ordsums$range_size_full


ordsums$full_chk <- ordsums$full_in_common2 + ordsums$full_w_o_all
summary(ordsums$full_chk)
ordsums$abun_chk <- ordsums$all_in_common + ordsums$all_w_o_full
summary(ordsums$all_chk)
#### EXPLORE

summary(ordsums$all_in_common)
summary(ordsums$all_w_o_full)  
summary(ordsums$full_in_common2)
summary(ordsums$full_w_o_all)  ## more often in full but not in all

summary(ordsums$range_size_all)
summary(ordsums$range_size_full)

summary((ordsums$range_size_full-ordsums$range_size_all)/ordsums$range_size_full)
hist((ordsums$range_size_all)/ordsums$range_size_full)
### median 2/3
### median : 50% larger,    35% smaller


write.csv(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_3_23_2018.csv", row.names=F)



####################################################
########## PREDICT ORDINAL ABUNDANCE WITHIN RANGES##
####################################################
library(maptools)
library(raster)
library(rgdal)
library(colorRamps)
library(ordinal)
library(reshape2)
library(rgeos)
length(list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN"))
edd <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## eddmaps abundance data (thinned) with environmental data

#ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_3_23_2018.csv", header = T, sep = ",",  
#                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## best fit models for each species
#ordsums <- ordsums[ordsums$kappa >0 & ordsums$kappaP <0.05,]

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_5_7_2018.csv", header = T, sep = ",",  
                  quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
#ordsums <- ordsums[ordsums$rarest.bin>=10 & ordsums$kappaP <0.05 & ordsums$kappa >0,]

#ordsums <- ordsums[ordsums$rule5broke == 1,]
#summary(ordsums$kappa)
spp <- ordsums$species.code


formulae <- data.frame(cbind(ordsums$species.code,ordsums$formu), stringsAsFactors = F)
str(formulae)
colnames(formulae) <- c("species","formula")

#spp <- formulae$species[order(formulae$species)]
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
#proj4string(bio) <-  "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
## environmental variables
bioD <- as.data.frame(bio)
str(bioD)
summary(bioD)
bioD[,1:4] <- bioD[,1:4]/10


#setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/") 
## location of maxent ranges
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/")
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN")
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN")
ordsums$class1_5_7 <- -99
ordsums$class2_5_7 <- -99
ordsums$class3_5_7 <- -99
#ordsums$AUC <- -99
ordsums$medianClass_5_7 <- -99
ordsums$mode_5_7 <- -99
ordsums$vr_5_7 <- -99

for (i in 15:length(spp)){
  
 species <- edd[edd$species == spp[i],]
 species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3))
 
 #species[,7:10] <- species[,7:8]/10
 
 M <- clm(formula(formulae$formula[formulae$species==spp[i]]), data = species)
 #AIC(M)
 
 range <- raster(paste0("Binary_asciis/",spp[i],".asc"))
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
 #plot(abun)
 
 ordsums$class1_5_7[i] <- length(p$abundance[p$abundance_presence == 1 & !is.na(p$abundance_presence)])/length(p$abundance[!is.na(p$presence) & p$presence==1])
 ordsums$class2_5_7[i] <- length(p$abundance[p$abundance_presence == 2 & !is.na(p$abundance_presence)])/length(p$abundance[!is.na(p$presence) & p$presence==1])
 ordsums$class3_5_7[i] <- length(p$abundance[p$abundance_presence == 3 & !is.na(p$abundance_presence)])/length(p$abundance[!is.na(p$presence) & p$presence==1])
 
 
 writeRaster(abun, paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/",spp[i],".asc"),overwrite=T)
 writeRaster(hi_abun,paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/",spp[i],".asc"),overwrite=T)
 # #  
 # #  get.auc <- read.table(paste0("FINAL_OUTPUT/",spp[i],"/maxentResults.csv"), header = T, sep = ",", quote= "\"", 
 # #                        comment.char= "", stringsAsFactors = F, strip.white = T)  
 # #  ordsums$AUC[i] <- get.auc$Test.AUC[11]


  
     ordsums$medianClass_5_7[i] <- median(p$abundance_presence[p$presence ==1 & !is.na(p$presence)])

     
     
     if(ordsums$class1_5_7[i] > ordsums$class2_5_7[i] & ordsums$class1_5_7[i] > ordsums$class3_5_7[i]){
       ordsums$mode_5_7[i] <- 1
       ordsums$vr_5_7[i] <- 1-ordsums$class1_5_7[i]
     }
     
     if(ordsums$class2_5_7[i] > ordsums$class1_5_7[i] & ordsums$class2_5_7[i] > ordsums$class3_5_7[i]){
       ordsums$mode_5_7[i] <- 2
       ordsums$vr_5_7[i] <- 1-ordsums$class2_5_7[i]
     }
     
     if(ordsums$class3_5_7[i] > ordsums$class2_5_7[i] & ordsums$class3_5_7[i] > ordsums$class1_5_7[i]){
       ordsums$mode_5_7[i] <- 3
       ordsums$vr_5_7[i] <- 1-ordsums$class3_5_7[i]
     }
  
     
     ## range size
     ######## calculate range size
     #plot(range)
     proj4string(range) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
     sp <- names(range)
     names(range) <- "layer"
     range <- rasterToPolygons(range, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
     range <- spTransform(range, "+init=epsg:5070") ## Albers
     ordsums$areaSQKM_5_7[ordsums$species.code==sp] <- area(range)*1e-6
     
     ######## calculate area in class 3  
     abun.i <- hi_abun
     proj4string(abun.i) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
     #sp <- names(abun.i)
     names(abun.i) <- "layer"
     abun.i <- rasterToPolygons(abun.i, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
     
     
     if(!is.null(abun.i)){
       abun.i <- spTransform(abun.i, "+init=epsg:5070") ## Albers
       ordsums$hiAbunSQKM_5_7[ordsums$species.code==sp] <- area(abun.i)*1e-6  
     } else {
       ordsums$hiAbunSQKM_5_7[ordsums$species.code==sp] <- 0   
     }
     
     
     
     
     ## infilling
     #species<- raster(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN/",spp[i],".asc"))
     # writeRaster(hi_abun,paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/",spp[i],".asc"),overwrite=T)
     #plot(species)
     #names(hi_abun) <- "layer"
     hi_abunD <- as.data.frame(hi_abun)
     # colnames(speciesD) <- "layer"
     #ncells<- 
     
     
     pts <- species
     coordinates(pts) <- c(5,4)
     ext <- extract(hi_abun, pts) ## extract climate values to points
     #colnames(ext)<- "impact"
     pts<-data.frame(cbind(pts$abundance,ext))
     colnames(pts)<- c("observ","impact")
     ordsums$infilling[i]<- length(pts$observ[pts$observ==3 & pts$impact == 1])/length(speciesD$layer[speciesD$layer==1 &!is.na(speciesD$layer==1 )])
     
     #print(i)
     
     
     
     
     
  
  print(i)
  
}

write.csv(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5epv_up_5_7_2018.csv", row.names=F)
#length(ordsums$species.code[ordsums$AUC > 0.70 & ordsums$kappaP < 0.05 & ordsums$MESSgn10 >0.9])
#hist(ordsums$)
#ordsums$range_size[ordsums$species.code==spp[i]] <- length(rangeD$presence[rangeD$presence == 1 & !is.na(rangeD$presence)])
#ordsums$impact_range_size

#ordsums$prop3 <- NULL
#write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_5_2018.csv", row.names=F)

#### calculate size of full study region for proportional estimates
studyregion <- hi_abun
#names(studyregion) <- "layer"
stuDy <- as.data.frame(studyregion)
stuDy$layer[stuDy$layer == 1 & !is.na(stuDy$layer)] <- 0
studyregion <- raster(nrows=nrow(studyregion), ncols=ncol(studyregion), ext=extent(studyregion), vals=stuDy$layer)
#studyregion$layer[studyregion$layer == 1 & !is.na(studyregion$layer)] <- 0

proj4string(studyregion) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
studyregion <- rasterToPolygons(studyregion, na.rm=T, dissolve=T)
studyregion <- spTransform(studyregion, "+init=epsg:5070")
studyregion <- area(studyregion)*1e-6
#7843649


ordsums$estab_US <- ordsums$areaSQKM_5_7/studyregion
ordsums$impac_US <- ordsums$hiAbunSQKM_5_7/studyregion
ordsums$impac_estab <- ordsums$hiAbunSQKM_5_7/ordsums$areaSQKM_5_7

summary(as.factor(ordsums$USE))
write.csv(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5epv_up_5_7_2018.csv", row.names=F)


#### explore variation in abundance within Range ############

#  #spp <- ("FRAL4", "ROMU", "PUMOL", "CYSC4", "LISI")
#  #species <- as.data.frame(raster(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/FRAL4.asc")))
#  #colnames(species) <- "class"
#  #species$class <- ordered(species$class, levels=c(1,2,3))
#  
#  
#  #ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_6_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  
#  #spp<- ordsums$species.code
#  
#  #ordsums$class1 <- -99
#  #ordsums$class2 <- -99
#  #ordsums$class3 <- -99
#  ordsums$medianClass <- -99
#  ordsums$mode <- -99
#  ordsums$vr <- -99
#  
#  for (i in 1:length(spp)){
#  
#    species <- as.data.frame(raster(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/",spp[i],".asc")))  
#    colnames(species)<-"class"
#    species <- as.data.frame(species[!is.na(species$class) & species$class > 0,])
#    colnames(species)<-"class"
#    
#    #ordsums$class1[i] <- length(species$class[species$class==1])/length(species$class)  
#    #ordsums$class2[i] <- length(species$class[species$class==2])/length(species$class) 
#    #ordsums$class3[i] <- length(species$class[species$class==3])/length(species$class) 
#    ordsums$medianClass[i] <- median(species$class)
#    
#    if(ordsums$class1[i] > ordsums$class2[i] & ordsums$class1[i] > ordsums$class3[i]){
#      ordsums$mode[i] <- 1
#      ordsums$vr[i] <- 1-ordsums$class1[i]
#    }
#     
#    if(ordsums$class2[i] > ordsums$class1[i] & ordsums$class2[i] > ordsums$class3[i]){
#      ordsums$mode[i] <- 2
#      ordsums$vr[i] <- 1-ordsums$class2[i]
#    }
#    
#    if(ordsums$class3[i] > ordsums$class2[i] & ordsums$class3[i] > ordsums$class1[i]){
#      ordsums$mode[i] <- 3
#      ordsums$vr[i] <- 1-ordsums$class3[i]
#    }
#  
#   
#  print(i)
#
#
#}


summary(ordsums$class1)
summary(ordsums$class2)
summary(ordsums$class3)
summary(ordsums$medianClass)
summary(ordsums$mode)
summary(ordsums$vr)
  
hist(ordsums$class1)
hist(ordsums$class2)
hist(ordsums$class3)
hist(ordsums$medianClass)
hist(ordsums$mode)
hist(ordsums$vr, breaks=50)
summary(ordsums$vr)

write.table(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_11_2018.csv", sep=",", row.names=F)

plot(ordsums$hiAbun_prop_range~ordsums$class3)
cor.test(ordsums$hiAbun_prop_range,ordsums$class3,family="spearman") # .9996 correlation

rm(edd,abun,hi_abun,aucs,d,dt,dtest,M,p,range,rangeD,samples,samples1,sp_runs,species,spp_val,spp_val1,tab,test,thresh,train,file,i,j,out_file,samp,samp1,spp,test_n,tp,train_n)
gc()

######################## exploration #####################

hist(ordsums$class_3[1:34], breaks=20)
hist(ordsums$class_2[1:34], breaks=20)
hist(ordsums$class_1[1:34], breaks=20)

summary((ordsums$class_3 + ordsums$class_2 + ordsums$class_1)/ordsums$range_size_all )
hist(ordsums$class_3/ordsums$range_size_all)
summary((ordsums$class_3/(ordsums$class_3 + ordsums$class_2 + ordsums$class_1)))
summary(ordsums$class_3/ordsums$range_size_allvars)

hist(ordsums$class_3/(ordsums$class_3 + ordsums$class_2 + ordsums$class_1))
ordsums$prop3 <- ordsums$class_3/ordsums$range_size_all

hist(ordsums$prop3[ordsums$AUCtest_allvars > .7 & ordsums$ALLVARSgn10 > .9 & ordsums$kappa >0 & ordsums$kappaP < 0.05], breaks=50)
hist(ordsums$kappa[ordsums$prop3 == 1])
hist(ordsums$kappa[ordsums$prop3 == 0])

plot(ordsums$bin3o[ordsums$kappa>.2]/ordsums$no.pts.abun[ordsums$kappa>.2], ordsums$prop3[ordsums$kappa>.2])
length(ordsums$prop3[ordsums$AUCtest_allvars > .7 & ordsums$ALLVARSgn10 > .9 & ordsums$kappa >0.2 & ordsums$kappaP < 0.05] )
hist(ordsums$kappa[ordsums$AUCtest_allvars > .7 & ordsums$ALLVARSgn10 > .9 & ordsums$kappaP < 0.05], breaks=100)
abline(v=.05)
  

hist(ordsums$class_3/ordsums$range_size_all)
hist(ordsums$class_2/ordsums$range_size_all)
hist(ordsums$class_1/ordsums$range_size_all)
summary(ordsums$class_2/ordsums$range_size_all)
median(ordsums$class_2/ordsums$range_size_all)

proportions <- c(ordsums$class_1/ordsums$range_size_all, ordsums$class_2/ordsums$range_size_all, ordsums$class_3/ordsums$range_size_all)
hist(proportions)
summary(proportions)


ordsums2 <- ordsums[ordsums$AUCtest_allvars > .7 & ordsums$ALLVARSgn10 > .9 & ordsums$kappaP < 0.05 & ordsums$kappa > 0,]


hist(ordsums2$class_3)
hist(ordsums2$class_2)
hist(ordsums2$class_1)

summary(ordsums2$class_3)
summary(ordsums2$class_2)
summary(ordsums2$class_1)

  plot(abun)
  plot(hi_abun)
  extent(bio)==extent(range)
  ncell(bio)==ncell(range)
  nrow(bio)==nrow(range)
  ncol(bio)==ncol(range)
  
  head(rangeD)

########################################################
## Create Species Richness Map 
######################################################
#ordsumsUP <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_4_2018.csv", header = T, sep = ",",  
#                        quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
#  
#ordsumsDOWN <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_4_29_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  
#
#
#form.comp <- spp70
#for (i in 1:70){
#  form.comp$UPform[i] <- ordsums$formu[ordsums$species.code == form.comp$x[i]]
#  form.comp$DOform[i] <- ordsumsDOWN$formu[ordsumsDOWN$species.code == form.comp$x[i]]
#  print(i)
#}
#
#write.csv(form.comp[form.comp$DOform != form.comp$UPform,], "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/formComp.csv")
#head(form.comp$x[form.comp$DOform == form.comp$UPform])


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_rounded_down/ALL_ABUN/AIAL.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_rounded_down/HI_ABUN/RAFI.asc"))


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/ALL_ABUN/AIAL.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/ALL_ABUN/LELE10.asc"))


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_rounded_down/ALL_ABUN/EUES.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_rounded_down/ALL_ABUN/CEDI3.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_rounded_down/ALL_ABUN/CIAR4.asc"))


xplot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/EUES.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CEDI3.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CIAR4.asc"))


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/FRAL4.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_rounded_down/ALL_ABUN/FRAL4.asc"))



#ordsums2 <- ordsums[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.95 & ordsums$kappa > 0 & ordsums$kappaP < 0.05,]  
spp59 <- read.csv("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/spplist59.csv", stringsAsFactors = F, header=T)
spp70 <- read.csv("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/spplist70.csv", stringsAsFactors = F, header=T)
filenames <- as.data.frame(spp59) 
filenames <- as.data.frame(spp70$x[spp70$x %in% ordsumsDOWN$species[ordsumsDOWN$kappa >0 & ordsumsDOWN$kappaP<0.05]])   
filenames <- as.data.frame(spp70) 



############################################################
ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5epv_up_5_7_2018.csv", header = T, sep = ",",  
                                                quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
                        
head(ordsums)
ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP < 0.05 & ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > 0.9,]

filenames <- data.frame(ordsums2$species.code)
colnames(filenames) <- "x"
    
library(raster)
library(rgdal)
library(ordinal)
library(reshape2)

### Abundant richness
#filenames <- as.data.frame(ordsums2$species.code)
#filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", filenames$x,".asc")
filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", filenames$x,".asc")


ldf=stack(filenames$files)

overlap=calc(ldf, sum)
plot(overlap)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc"
writeRaster(overlap, filename=out_file,overwrite=T)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/summaries/figures/richness_map_HI_ABUN.pdf",width=11)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()

### Species richness
filenames$files2 <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", filenames$x,".asc")

ldf=stack(filenames$files2)

overlap=calc(ldf, sum)
plot(overlap)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc"
writeRaster(overlap, filename=out_file, overwrite=T)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/summaries/figures/richness_map_FULL.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()






#### 2016 Allen and Bradley binary asciis####

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_16_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  

filenames <- paste0("C:/Users/Localadmin/Google Drive/NSF_GSS_shared/Watch_Lists/range_maps/binary_asciis/sppRichCtyBias_current/sppRichCtyBias_current/sppRichCtyBias_current",ordsums$species.code[ordsums$USE==1],"_SppRichBias_FivePercMTP.asc")
ldf <- stack(filenames)
#plot(ldf)
overlap =calc(ldf,sum)
plot(overlap)
quantile(overlap)

overlapD <- as.data.frame(overlap)
overlapD$layer[overlapD$layer <= 29 & !is.na(overlapD$layer)] <-0
overlapD$layer[overlapD$layer > 29 & !is.na(overlapD$layer)] <-1
overlap <-  raster(nrows=nrow(overlap), ncols=ncol(overlap), ext=extent(overlap), vals=overlapD$layer)
plot(overlap)

full_hotspot <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_FULL.asc")
plot(full_hotspot,add=T)
proj4string(full_hotspot)<-proj4string(overlap)


ncell(full_hotspot)
ncell(overlap)

### ABUNDANCE-WEIGHTED RICHNESS #########
filenames <- as.data.frame(ordsums2$species.code)
filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/", filenames$`ordsums2$species.code`,".asc")

ldf=stack(filenames$files)

overlap=calc(ldf, sum)
plot(overlap)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_ABUN_WEIGHT.asc"
writeRaster(overlap, filename=out_file,overwrite=T)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/summaries/figures/richness_map_ABUN_WEIGHT.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off


current_rich=raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_ABUN_WEIGHT.asc")
quantile(current_rich$richness_map_ABUN_WEIGHT)

#length(rich$richness_map_ABUN_WEIGHT[rich$richness_map_ABUN_WEIGHT > quantile(current_rich$richness_map_ABUN_WEIGHT, 0.75) &
#                                   !is.na(rich$richness_map_ABUN_WEIGHT)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 19%
#length(rich$richness_map_ABUN_WEIGHT[rich$richness_map_HI_ABUN >= quantile(current_rich$richness_map_ABUN_WEIGHT, 0.75) &
#                                   !is.na(rich$richness_map_ABUN_WEIGHT)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 32%


rich <- as.data.frame(current_rich)
rich$richness_map_ABUN_WEIGHT[rich$richness_map_ABUN_WEIGHT <= quantile(rich$richness_map_ABUN_WEIGHT, 0.75, na.rm=T)] <- 0
rich$richness_map_ABUN_WEIGHT[rich$richness_map_ABUN_WEIGHT >  quantile(rich$richness_map_ABUN_WEIGHT, 0.75, na.rm=T)] <- 1

abun_wt_hotspot <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_ABUN_WEIGHT)
writeRaster(abun_wt_hotspot, "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_ABUN_WT.asc", overwrite=T)


comp <- as.data.frame(stack(abun_wt_hotspot,full_hotspot))
head(comp)
colnames(comp) <- c("layer.1", "layer.2")
comp$x[comp$layer.1==0 & comp$layer.2==0 & !is.na(comp$layer.2)] <- "non-hotspot" #neither
comp$x[comp$layer.1==1 & comp$layer.2==0 & !is.na(comp$layer.2)] <- "abundance hotspot" #just abun
comp$x[comp$layer.1==1 & comp$layer.2==1 & !is.na(comp$layer.2)] <- "both hotspots" #both
comp$x[comp$layer.1==0 & comp$layer.2==1 & !is.na(comp$layer.2)] <- "richness hotspot" #just full

length(comp$x[comp$x == "both hotspots" & !is.na(comp$x)])/length(comp$layer.1[comp$layer.1 == 1 & !is.na(comp$layer.1)]) ## percent of abundance hotspot also richness
length(comp$x[comp$x == "both hotspots" & !is.na(comp$x)])/length(comp$layer.1[comp$layer.2 == 1 & !is.na(comp$layer.1)]) ## percent of richness hotspot also abundance
.5*(length(comp$x[comp$x == "both hotspots"& !is.na(comp$x)])/length(comp$layer.1[comp$layer.1 == 1& !is.na(comp$x)]) + length(comp$x[comp$x == "both hotspots"& !is.na(comp$x)])/length(comp$layer.1[comp$layer.2 == 1 & !is.na(comp$x)]))
## 69% overlap


comp$x <- as.factor(comp$x) 

comp_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp$x)
plot(comp_hotspot)
spplot(comp_hotspot)

writeRaster(comp_hotspot,"C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc",overwrite=T)







#  ### repeat for abundance points Maxent
#  filenames = list.files("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/ABUN/",pattern="*.asc", full.names=TRUE)
#  ldf=stack(filenames)
#  
#  overlap2=calc(ldf, sum)
#  
#  out_file="C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.asc"
#  writeRaster(overlap2, filename=out_file)
#  
#  pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.pdf",width=11,height=8.5)
#  colors=rev(heat.colors(100))
#  spplot(overlap2, col.regions=colors)
#  dev.off()

########################################
########### no extrap hi_abun_hotspots
#######################################
ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_4_2018.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

head(ordsums)
ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP < 0.05 & ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > 0.9 & !is.na(ordsums$hiAbunSQKM),]

filenames <- data.frame(ordsums2$species.code)
colnames(filenames) <- "x"

library(raster)
library(rgdal)
library(ordinal)
library(reshape2)

### Abundant richness
#filenames <- as.data.frame(ordsums2$species.code)
#filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", filenames$x,".asc")
filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/", filenames$x,".asc")

b <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/AIAL.asc")

filenames$files[1]=="C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/AIAL.asc"

filenames2 <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/")
filenames2<- filenames2[filenames2 %in% paste0(filenames$x,".asc")]


ldf=stack(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN_no_xtrap/",filenames2))

overlap=calc(ldf, sum)
plot(overlap)

current_rich=overlap
names(current_rich) <- "richness_map_HI_ABUN"
quantile(current_rich$richness_map_HI_ABUN)

rich <- as.data.frame(current_rich)
length(rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN > quantile(current_rich$richness_map_HI_ABUN, 0.75) &
                                   !is.na(rich$richness_map_HI_ABUN)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 19%
length(rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN >= quantile(current_rich$richness_map_HI_ABUN, 0.75) &
                                   !is.na(rich$richness_map_HI_ABUN)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 32%



rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN <= quantile(current_rich$richness_map_HI_ABUN, 0.75)] <- 0
rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN > quantile(current_rich$richness_map_HI_ABUN, 0.75)] <- 1

hi_abun_hotspot <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_HI_ABUN)
#writeRaster(hi_abun_hotspot, "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_HI_ABUN.asc",overwrite=T)
plot(hi_abun_hotspot)




###############################################################################################################
### Hotspot maps
################################################################################################################
#current_rich<- overlap
#names(current_rich) <- "richness_map_HI_ABUN"
current_rich=raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
quantile(current_rich$richness_map_HI_ABUN)

rich <- as.data.frame(current_rich)
length(rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN > quantile(current_rich$richness_map_HI_ABUN, 0.75) &
                                   !is.na(rich$richness_map_HI_ABUN)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 19%
length(rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN >= quantile(current_rich$richness_map_HI_ABUN, 0.75) &
                                   !is.na(rich$richness_map_HI_ABUN)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 32%



rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN <= quantile(current_rich$richness_map_HI_ABUN, 0.75)] <- 0
rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN > quantile(current_rich$richness_map_HI_ABUN, 0.75)] <- 1

hi_abun_hotspot <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_HI_ABUN)
writeRaster(hi_abun_hotspot, "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_HI_ABUN.asc",overwrite=T)
plot(hi_abun_hotspot)

#plot(current_rich, col=rev(heat.colors(100)))
#plot(hi_abun_hotspot, col=rev(heat.colors(2)))

####################################################
current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/richness_map_FULL.asc")
quantile(current_rich$richness_map_FULL)

rich <- as.data.frame(current_rich)
length(rich$richness_map_FULL[rich$richness_map_FULL > quantile(current_rich$richness_map_FULL, 0.75) &
                                   !is.na(rich$richness_map_FULL)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 19%
length(rich$richness_map_FULL[rich$richness_map_FULL >= quantile(current_rich$richness_map_FULL, 0.75) &
                                   !is.na(rich$richness_map_FULL)])/length(bioD$bio_2[!is.na(bioD$bio_2)]) ## 32%



rich$richness_map_FULL[rich$richness_map_FULL <= quantile(current_rich$richness_map_FULL, 0.75)] <- 0
rich$richness_map_FULL[rich$richness_map_FULL > quantile(current_rich$richness_map_FULL, 0.75)] <- 1


full_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_FULL)
plot(full_hotspot)
writeRaster(full_hotspot, "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_FULL.asc", overwrite=T)

####################### other maxent runs##############
#  current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/richness_map_ABUN.asc")
#  quantile(current_rich$richness_map_ABUN)
#  
#  rich <- as.data.frame(current_rich)
#  head(rich)
#  rich$richness_map_ABUN[rich$richness_map_ABUN < 43] <- 0
#  rich$richness_map_ABUN[rich$richness_map_ABUN >= 43] <- 1
#  
#  
#  abun_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_ABUN)
#  plot(abun_hotspot)
#  
#  
#  
#  comp <- as.data.frame(stack(abun_hotspot,full_hotspot))
#  head(comp)
#  comp$x[comp$layer.1==0 & comp$layer.2==0 ] <- 0 #neither
#  comp$x[comp$layer.1==1 & comp$layer.2==0 ] <- 1 #just abun
#  comp$x[comp$layer.1==1 & comp$layer.2==1 ] <- 3 #just full
#  comp$x[comp$layer.1==0 & comp$layer.2==1 ] <- 6 #both
#  
#  length(comp$x[comp$x == 3])/length(comp$layer.1[comp$layer.1 == 1])
#  length(comp$x[comp$x == 3])/length(comp$layer.1[comp$layer.2 == 1])
#  
#  comp_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp$x)
#  ## based on richness maps here, adapt to summed abundance maps
#  plot(comp_hotspot)
#  
# 
#  comp3$dev <- (comp3$prop - mean(comp3$prop,na.rm=T))/sd(comp3$prop, na.rm=T)
#  summary(comp3$dev)
#  hist(comp3$dev)
#  quantile(comp3$dev,na.rm=T)
#  
#  comp3$sum[comp3$prop>=1] <- 2
#  comp3$sum[comp3$prop<1 & comp3$prop>-1] <- 1
#  comp3$sum[comp3$prop<=-1] <- 0
#  
#  comp3$sum2[comp3$dev>=2] <- 2
#  comp3$sum2[comp3$dev<2 & comp3$dev>-2] <- 1
#  comp3$sum2[comp3$dev<=-2] <- 0
#  comp3 <- as.data.frame(stack(hi_abun_rich,full_max_rich))
#  head(comp3)
#  comp3$prop<- comp3$richness_map_HI_ABUN/comp3$richness_map_FULL
#  
#  hist(comp3$prop)
#  mean(comp3$prop, na.rm=T)
#  median(comp3$prop, na.rm=T)
#  sd(comp3$prop, na.rm=T)
#  summary(comp3$prop)
################################################

comp <- as.data.frame(stack(hi_abun_hotspot,full_hotspot))
head(comp)
colnames(comp) <- c("layer.1", "layer.2")

comp$x[comp$layer.1==0 & comp$layer.2==0 & !is.na(comp$layer.2)] <- "non-hotspot" #neither
comp$x[comp$layer.1==1 & comp$layer.2==0 & !is.na(comp$layer.2)] <- "abundance hotspot" #just abun
comp$x[comp$layer.1==1 & comp$layer.2==1 & !is.na(comp$layer.2)] <- "both hotspots" #both
comp$x[comp$layer.1==0 & comp$layer.2==1 & !is.na(comp$layer.2)] <- "richness hotspot" #just full

length(comp$x[comp$x == "both hotspots" & !is.na(comp$x)])/length(comp$layer.1[comp$layer.1 == 1 & !is.na(comp$layer.1)]) ## percent of abundance hotspot also richness
length(comp$x[comp$x == "both hotspots" & !is.na(comp$x)])/length(comp$layer.1[comp$layer.2 == 1 & !is.na(comp$layer.1)]) ## percent of richness hotspot also abundance
.5*(length(comp$x[comp$x == "both hotspots"& !is.na(comp$x)])/length(comp$layer.1[comp$layer.1 == 1& !is.na(comp$x)]) + length(comp$x[comp$x == "both hotspots"& !is.na(comp$x)])/length(comp$layer.1[comp$layer.2 == 1 & !is.na(comp$x)]))

## 47% overlap

#Abun <- 89975
#occ <- 111185
comp$x <- as.factor(comp$x) 

comp_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp$x)
plot(comp_hotspot)
spplot(comp_hotspot)

writeRaster(comp_hotspot,"C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc",overwrite=T)


################## proportion of abundant species
abun_prop=stack("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/richness_map_FULL.asc",
                "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
abun_prop <- as.data.frame(abun_prop)
abun_prop$abun_prop <- abun_prop$richness_map_HI_ABUN/abun_prop$richness_map_FULL

summary(abun_prop$abun_prop)
hist(abun_prop$abun_prop)

#abun_prop$prop_hot[abun_prop$abun_prop >= quantile(abun_prop$abun_prop, .75, na.rm=T) & !is.na(abun_prop$abun_prop)] <- 1
#abun_prop$prop_hot[abun_prop$abun_prop < quantile(abun_prop$abun_prop, .75, na.rm=T) & !is.na(abun_prop$abun_prop)] <- 0



#comp$prop <- comp$layer.1/comp$layer.2
abun_prop <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=abun_prop$abun_prop)#,
                   #raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=abun_prop$prop_hot))

plot(abun_prop)

#pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/summaries/figures/hi_abun_proportion.pdf",width=11,height=8.5)
#colors=rev(heat.colors(100))
#spplot(abun_prop, col.regions=colors, main="Proportion of highly abundant species")
#dev.off()

writeRaster(abun_prop,"C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hi_abun_proportion.asc",overwrite=T)
#writeRaster(abun_prop$layer.2,"C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hi_abun_prop_hot.asc",overwrite=T)
 

###############################
#### Calculate Range Sizes ####
###############################
library(raster)
library(rgdal)
library(rgeos)
#filenames <- as.data.frame(ordsums2$species.code)
#filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", filenames$`ordsums2$species.code`,".asc")
#filenames$files2 <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", filenames$`ordsums2$species.code`,".asc")

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_4_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  
ordsums2 <- ordsums[ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > 0.9 & ordsums$kappa > 0 & ordsums$kappaP < 0.05,]  
head(ordsums)

### Do all 155 species
ranges <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", full.names = T)
abuns  <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN/", full.names = T)

#EPSG <- make_EPSG()
#epsg <- EPSG[EPSG$code == 5070,]
#epsg <- epsg[!is.na(epsg$code),]

#library(rgeos)

for (i in 1:length(ranges)){
  
  ######## calculate range size
  range.i <- raster(ranges[i])
  proj4string(range.i) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  sp <- names(range.i)
  names(range.i) <- "layer"
  range.i <- rasterToPolygons(range.i, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
  range.i <- spTransform(range.i, "+init=epsg:5070") ## Albers
  ordsums$areaSQKM[ordsums$species.code==sp] <- area(range.i)*1e-6
  
  ######## calculate area in class 3  
  abun.i <- raster(abuns[i])
  proj4string(abun.i) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  #sp <- names(abun.i)
  names(abun.i) <- "layer"
  abun.i <- rasterToPolygons(abun.i, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
  
  
  if(!is.null(abun.i)){
  abun.i <- spTransform(abun.i, "+init=epsg:5070") ## Albers
  ordsums$hiAbunSQKM[ordsums$species.code==sp] <- area(abun.i)*1e-6  
  } else {
  ordsums$hiAbunSQKM[ordsums$species.code==sp] <- 0   
  }
  print(i)
  #write.csv(i,"C:/Users/Localadmin/Google Drive/AREAprogress.csv", row.names=F)
}

hist(ordsums$hiAbunSQKM/ordsums$areaSQKM)


#### calculate size of full study region for proportional estimates
studyregion <- raster(ranges[i])
names(studyregion) <- "layer"
stuDy <- as.data.frame(studyregion)
stuDy$layer[stuDy$layer == 1 & !is.na(stuDy$layer)] <- 0
studyregion <- raster(nrows=nrow(studyregion), ncols=ncol(studyregion), ext=extent(studyregion), vals=stuDy$layer)
#studyregion$layer[studyregion$layer == 1 & !is.na(studyregion$layer)] <- 0

proj4string(studyregion) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
studyregion <- rasterToPolygons(studyregion, na.rm=T, dissolve=T)
studyregion <- spTransform(studyregion, "+init=epsg:5070")
studyregion <- area(studyregion)*1e-6

 ordsums$range_prop_area <- ordsums$areaSQKM/studyregion
 ordsums$hiAbun_prop_full <- ordsums$hiAbunSQKM/studyregion
 ordsums$hiAbun_prop_range <- ordsums$hiAbunSQKM/ordsums$areaSQKM
 
 hist(ordsums$range_prop_area)
 hist(ordsums$hiAbun_prop_full)
 hist(ordsums$hiAbun_prop_range)

write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_5_6_2018.csv", row.names=F)
#ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_4_6_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  
#ordsums2 <- ordsums[ordsums$AUC > 0.7 & ordsums$MESSgn10 > 0.9 & ordsums$kappa > 0 & ordsums$kappaP < 0.05,]  

summary(ordsums2$hiAbun_prop_range) 
min(ordsums2$hiAbun_prop_range)
#0-83%, median 9%   ## 3%-20%
med.abun <- ordsums2[ordsums2$hiAbun_prop_range >.07 & ordsums2$hiAbun_prop_range < .12,]



summary(ordsums2$areaSQKM)
summary(ordsums2$range_prop_area)



#################################
####### Ecoregion Calculations
#################################


extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
proj4string(extentShape) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

eco <- readOGR(dsn="C:/Users/Localadmin/Documents/Ecoregion_Shapefiles", layer="NA_CEC_Eco_Level1")
head(eco)

ecoD <- as.data.frame(eco)
ecoD$NA_L1NAME <- as.character(ecoD$NA_L1NAME)
ecoD$NA_L1KEY <- as.character(ecoD$NA_L1KEY)
ecoD$NA_L1CODE <- as.character(ecoD$NA_L1CODE)
ecoD$row <- row.names(ecoD)
head(ecoD)
head(eco)



extentShape <- spTransform(extentShape, proj4string(eco))
plot(extentShape, col="grey")
plot(eco, add=T)

extentShape <- gBuffer(extentShape, byid=TRUE, width=0)
eco <- gBuffer(eco, byid=TRUE, width=0)

eco <- gIntersection(extentShape, eco, byid=T)
plot(eco, add=T, lwd=3, border="blue")

eco.rows <- as.data.frame(row.names(eco), stringsAsFactors=F)
colnames(eco.rows) <- "polygon"
eco.rows$row <-  gsub(pattern="0 ", replacement = "", x=eco.rows$polygon)
#colnames(eco.rows) <- "row"
head(eco.rows)
eco.rows <- merge(eco.rows, ecoD, by="row", all=F)
head(eco.rows)
row.names(eco.rows) <- eco.rows$polygon

#ecoD[ecoD$row %in% ]
eco <- SpatialPolygonsDataFrame(Sr=eco, data=eco.rows)
head(eco)

plot(extentShape)
#unique(as.character(eco$NA_L1KEY))
plot(eco[eco$NA_L1NAME=="NORTH AMERICAN DESERTS",], add=T, col="blue")#check

eco <- unionSpatialPolygons(eco, IDs=eco@data$NA_L1NAME)

IDs <- as.data.frame(row.names(eco), stringsAsFactors=F)
colnames(IDs)<-"IDs"
row.names(IDs) <- IDs$IDs

eco <- SpatialPolygonsDataFrame(Sr=eco, data=IDs)
eco$IDs <- as.character(eco$IDs)
head(eco)
str(eco@data)


plot(extentShape)
##CHECK, all good
#spplot(eco)
plot(eco[eco$IDs=="WATER",], col="blue", add=T)
plot(eco[eco$IDs=="NORTH AMERICAN DESERTS",], col="tan", add=T)
plot(eco[eco$IDs=="MEDITERRANEAN CALIFORNIA",], col="yellow", add=T)
plot(eco[eco$IDs=="SOUTHERN SEMIARID HIGHLANDS",], col="orange", add=T)
plot(eco[eco$IDs=="TEMPERATE SIERRAS",], col="red", add=T)
plot(eco[eco$IDs=="TROPICAL WET FORESTS",], col="purple", add=T)
plot(eco[eco$IDs=="NORTHERN FORESTS" ,], col="dark green", add=T)
plot(eco[eco$IDs=="NORTHWESTERN FORESTED MOUNTAINS",], col="light green", add=T)
plot(eco[eco$IDs=="MARINE WEST COAST FOREST",], col="turquoise", add=T)
plot(eco[eco$IDs=="EASTERN TEMPERATE FORESTS" ,], col="green", add=T)
plot(eco[eco$IDs=="GREAT PLAINS" ,], col="brown", add=T)

eco<-spTransform(eco, "+init=epsg:5070")

eco$area <- area(eco)*(1e-6)
eco$prop <- eco$area/sum(eco$area)


writeOGR(eco, dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco", driver="ESRI Shapefile")
##################################


library(rgdal)
library(raster)
library(rgeos)

eco <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco", stringsAsFactors = F)

#hi_abun_richness <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
#full_richness    <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")

hi_abun_hotspots <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_HI_ABUN.asc")
full_hotspot     <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_FULL.asc")

compare <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc" )
props <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hi_abun_proportion.asc")


proj4string(hi_abun_richness)="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(full_richness   )="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(hi_abun_hotspots)="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(full_hotspot    )="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(props)="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(compare)="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"


####### BY AREA #############
#proj4string(hi_abun_richness) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
######################## repeat for abun (MOVE BELOW)
abun_hotspotPOLY <- rasterToPolygons(hi_abun_hotspots, na.rm=T, dissolve=T)
abun_hotspotPOLY <- spTransform(abun_hotspotPOLY, "+init=epsg:5070")
#abun_hotspotPOLY$area <- area(abun_hotspotPOLY)*1e-6



abun_hotspotPOLY <- gIntersection(eco,abun_hotspotPOLY, byid=T)
gIds <- data.frame(names(abun_hotspotPOLY))
colnames(gIds) <- "gIds"
row.names(gIds) <- gIds$gIds

abun_hotspotPOLY <- SpatialPolygonsDataFrame(Sr=abun_hotspotPOLY, data=gIds)
#abun_hotspotPOLY@data
abun_hotspotPOLY$area <- area(full_hotspotPOLY)*1e-6
#abun_hotspotPOLY@data
#abun_hotspotPOLY2$area <- area(abun_hotspotPOLY2)*1e-6

write.table(abun_hotspotPOLY2@data,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/eco_abun_hot.csv",sep=",")
writeOGR(abun_hotspotPOLY2, dsn="C:/Users/Localadmin/Documents/eco_abun_hot", layer="eco", driver="ESRI Shapefile" )

######################## full
#eco<-spTransform(eco, "+init=epsg:5070")

full_hotspotPOLY <- rasterToPolygons(full_hotspot, na.rm=T, dissolve=T)
full_hotspotPOLY <- spTransform(full_hotspotPOLY, "+init=epsg:5070")
full_hotspotPOLY$area <- area(full_hotspotPOLY)*1e-6

full_hotspotPOLY2 <- gIntersection(eco,full_hotspotPOLY, byid=T)
plot(full_hotspotPOLY2)

gIds <- data.frame(names(full_hotspotPOLY2))
colnames(gIds) <- "gIds"
row.names(gIds) <- gIds$gIds

full_hotspotPOLY <- SpatialPolygonsDataFrame(Sr=full_hotspotPOLY2, data=gIds)
full_hotspotPOLY@data
full_hotspotPOLY$area <- area(full_hotspotPOLY)*1e-6
full_hotspotPOLY@data

write.table(full_hotspotPOLY@data,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/eco_full_hot.csv",sep=",")
writeOGR(full_hotspotPOLY, dsn="C:/Users/Localadmin/Documents/eco_full_hot", layer="eco_full_hot", driver="ESRI Shapefile" )


writeOGR(eco, dsn="C:/Users/Localadmin/Desktop/intersection/eco", layer="eco", driver="ESRI Shapefile")
writeOGR(full_hotspotPOLY, dsn="C:/Users/Localadmin/Desktop/intersection/full", layer="full", driver="ESRI Shapefile")
writeOGR(abun_hotspotPOLY, dsn="C:/Users/Localadmin/Desktop/intersection/abun", layer="abun", driver="ESRI Shapefile")



#################################
plot(full_hotspot)
plot(hi_abun_hotspot)

eco <- spTransform(eco, proj4string(hi_abun_richness))

plot(eco)
plot(eco[10,], add=T,col="blue")
### full hotspot X ecoregion approximation
ri_hot <- extract(full_hotspot$hotspot_FULL, eco, fun=mean, na.rm=TRUE, sp = T)
ab_hot <- extract(hi_abun_hotspots$hotspot_HI_ABUN, eco, fun=mean, na.rm=TRUE, sp = T)


ab_hot@data
write.csv(cbind(ri_hot@data,ab_hot@data),"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/hot_eco.csv", sep=",")
          
##############################
#   #eco$IDs
#   #eco1 <- eco[eco$IDs == "TROPICAL WET FORESTS" | eco$IDs == "MARINE WEST COAST FOREST",]
#   
#   mean.abun <- extract(hi_abun_richness$richness_map_HI_ABUN , eco, fun=mean, na.rm=TRUE, sp = T)
#   mean.rich <- extract(full_richness$richness_map_FULL, eco, fun=mean, na.rm=TRUE, sp = T)
#   
#   
#   stdv.rich <- extract(full_richness$richness_map_FULL, eco, fun=sd, na.rm=TRUE, sp = T)

### everything in proj4string towgs (props 2!!!), ##check layer name!


#write.csv(ordsums, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/infilling.csv")


##### switch to hotspot instead of eco ########


hotPOLY <- rasterToPolygons(compare, na.rm=T, dissolve=T)


mean.abun <- extract(hi_abun_richness$richness_map_HI_ABUN, hotPOLY, fun=mean, na.rm=T, sp = T)
stdv.abun <- extract(hi_abun_richness$richness_map_HI_ABUN, hotPOLY, fun=sd, na.rm=T, sp=T)
mean.rich <- extract(full_richness$richness_map_FULL, hotPOLY, fun=mean, na.rm=T, sp=T)
stdv.rich <- extract(full_richness$richness_map_FULL, hotPOLY, fun=sd, na.rm=T, sp=T)
mean.prop <- extract(props$hi_abun_proportion, hotPOLY, fun=mean, na.rm=T, sp=T)
stdv.prop <- extract(props$hi_abun_proportion, hotPOLY, fun=sd, na.rm=T, sp=T)
max.prop <- extract(props$hi_abun_proportion, hotPOLY, fun=max, na.rm=T, sp=T)
min.prop <- extract(props$hi_abun_proportion, hotPOLY, fun=min, na.rm=T, sp=T)


comparSum <- cbind(mean.abun,stdv.abun,mean.rich,stdv.abun)
write.table(comparSum, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ComparSum.csv")


mean.prop@data
stdv.prop@data
min.prop@data
max.prop@data

plot(hotPOLY[3,], col="black")
#####################################
### mean.abun <- extract(hi_abun_richness$richness_map_HI_ABUN, ri, fun=sd, na.rm=TRUE, sp = T)
### stdv.abun <- extract(hi_abun_richness$richness_map_HI_ABUN, eco, fun=sd, na.rm=TRUE, sp = T)
### mean.prop <- extract(props$hi_abun_proportion, eco, fun=sd, na.rm=TRUE, sp = T)
### stdv.prop <- extract(props$hi_abun_proportion, eco, fun=sd, na.rm=TRUE, sp = T)
### 
### 
### # UNPOUND
### #eco$mean_rich <- mean.rich$richness_map_FULL
### #eco$mean_abun <- mean.abun$richness_map_HI_ABUN
### #eco$mean_prop <- mean.prop$hi_abun_proportion
### #eco$stdv_abun <- stdv.abun$richness_map_HI_ABUN
### #eco$stdv_rich <- stdv.rich$richness_map_FULL
### #eco$stdv_prop <- stdv.prop$hi_abun_proportion
### #eco$cova_abun <- eco$stdv_abun/eco$mean_abun
### #eco$cova_rich <- eco$stdv_rich/eco$mean_rich
### #eco$cova_prop <- eco$stdv_prop/eco$mean_prop

write.table(eco@data,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/eco_richness.csv",sep=",")
writeOGR(eco, dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/eco_richness", layer="eco", driver="ESRI Shapefile" )

#################################
#eco@data
#plot(eco[eco$IDs == "WATER",])
#plot(eco[eco$IDs == "NORTHERN FORESTS",], col="lightgreen", add=T, lwd=1)
#plot(eco[eco$IDs == "EASTERN TEMPERATE FORESTS",], col="green", add=T, lwd=1)
#plot(eco[eco$IDs == "WATER",], col="firebrick2", border="transparent", add=T, lwd=1)

eco@data
rm(abunpts,bio,compare,extentShape,maxentRange,mean.abun,mean.rich,messplot,messpts,occ.pts,occ.ran,ordinal,ordsums,ordsums2,species,speciesSum,stdv.rich,i,spp)
#############hotspot ecoregions




########################
###### correlation of species richness and abundance
########################
hi_abun_hotspots <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_HI_ABUN.asc")
full_hotspot     <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_FULL.asc")

hi_abun_richness <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
full_richness  <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")

correlations <- cbind(as.data.frame(hi_abun_richness), as.data.frame(full_richness))
head(correlations)
correlations<-correlations[!is.na(correlations$richness_map_HI_ABUN),]
hist(correlations$richness_map_HI_ABUN)
hist(correlations$richness_map_FULL)

cor.test(correlations$richness_map_HI_ABUN,correlations$richness_map_FULL, pairwise=T, method="spearman")
## .59 rho, p<0.0001
M <- glm(correlations$richness_map_HI_ABUN~correlations$richness_map_FULL)
coefs1 <- as.data.frame(M$coefficients)
M2 <- glm(correlations$richness_map_HI_ABUN~poly(correlations$richness_map_FULL,2))
summary(M)
summary(M2)
M$aic - M2$aic

M3 <- glm(correlations$richness_map_HI_ABUN~correlations$richness_map_FULL+I(correlations$richness_map_FULL^2))
summary(M3)
coefs<- as.data.frame(M3$coefficients)


x <- seq(1:50)
y <- coefs$`M3$coefficients`[1]*x^2 + coefs$`M3$coefficients`[2]*x + coefs$`M3$coefficients`[1]
y1 <- x*coefs1$`M$coefficients`[2] + coefs1$`M$coefficients`[2]

library(scales)
plot(correlations$richness_map_HI_ABUN~correlations$richness_map_FULL, col=alpha("black", .005), pch=16)
lines(y~x, col="red")
lines(y1~x, col="red")

###############################################
########## EXAMPLE SPECIES + FIGURES
###############################################
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
#install.packages("RColorBrewer")
library(RColorBrewer)


ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP<0.05 & ordsums$AUC > 0.7 & ordsums$MESSgn10 > .9,]
ordsums2[ordsums2$no.pts.full < 1000 & ordsums2$range_prop_area > .20,]

display.brewer.all(n=8, type="all", colorblindFriendly=T)

hi_abun_richness <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
full_richness  <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")
props <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hi_abun_proportion.asc")


plot(extentShape, lwd=5)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), add=T)
plot(extentShape, lwd=5)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), add=T)


extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

abunpts <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv",
                      sep=",", header=T, stringsAsFactors = F)
head(abunpts)
abunpts <- abunpts[abunpts$species=="POCU6",]
coordinates(abunpts) <- c(5,4)

occ.pts <- read.table("C:/Users/Localadmin/Documents/MaxEnt_modeling/species/POCU6.csv", sep=",", stringsAsFactors = F, header=T)
str(occ.pts)
coordinates(occ.pts) <- c(2,3)

occ.ran <- occ.pts[sample(1:nrow(occ.pts), 130),]
#abun.ran <- abunpts[sample(1:nrow(abunpts), 50),]
abun.ran <- occ.ran[sample(1:nrow(occ.ran), 40),]

#occ.ran <- occ.pts[sample(1:nrow(occ.pts), 20)]

plot(full_richness)
plot(extentShape, lwd=3, col="grey92", add=T)
plot(occ.ran,lwd=3,col="grey20",add=T)
plot(abun.ran, pch=19, col="firebrick2", add=T)



maxentRange <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/POCU6.asc")
plot(maxentRange)
names(maxentRange) <- "layer"
maxentRange <- rasterToPolygons(maxentRange, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)

plot(extentShape, col="grey92")
plot(maxentRange, add=T, col="mediumblue",border="mediumblue")
plot(maxentRange, add=T, col="firebrick",border="firebrick")


#ordinal <- readAsciiGrid("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc") 
ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$POCU6
ordinal$layer[ordinal$POCU6==0] <- NA

ordinal$layer2 <- 0
ordinal$layer2[ordinal$POCU6==3] <- 1
ordinal$layer2[ordinal$POCU6 != 3 | is.na(ordinal$POCU6)] <- NA


display.brewer.all(n=8, type="all", colorblindFriendly=T)
plot(extentShape, col="darkslategrey", lwd=2)   ### black works too
plot(ordinal$layer, add=T, col= (brewer.pal(3, name="YlOrRd")))

plot(extentShape, col="grey92")
plot(ordinal$layer2, col=c("firebrick"), add=T)






bio <- stack("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc")
             "C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/POCU6.asc")

bio$bio_2[bio$POCU6 == 0] <- NA
bio$bio_5[bio$POCU6 == 0] <- NA
bio$bio_8[bio$POCU6 == 0] <- NA
bio$bio_15[bio$POCU6 == 0] <- NA
#bio$nlcd_5[bio$POCU6 == 0] <- NA
#bio$nlcd_8[bio$POCU6 == 0] <- NA
#bio$nlcd_7[bio$ALPE4_FULL_all_vars == 0] <- NA
#bio$nlcd_8[bio$ALPE4_FULL_all_vars == 0] <- NA
bio$bio_2 <- bio$bio_2/10
bio$bio_5 <- bio$bio_5/10
bio$bio_8 <- bio$bio_8/10
#bio$bio_6 <- bio$bio_6/10

bio <- dropLayer(bio,5)



messpts <- data.frame(cbind(abunpts$bio_2, abunpts$bio_5, abunpts$bio_8, abunpts$bio_15))#), abunpts$nlcd_5, abunpts$nlcd_8))
head(messpts)
colnames(messpts)<-c("bio_2","bio_5", "bio_8", "bio_15" )#, "nlcd_5", "nlcd_8")

library(dismo)
messplot <- mess(bio,messpts,full=F)
#plot(messplot)

messplot$mess2 <- -99
#messplot$mess2[messplot$mess >= -10 & messplot$mess <= 10] <- 2
messplot$mess2[messplot$mess < -10] <- 1
messplot$mess2[messplot$mess >= -10] <- 3
messplot$mess2[is.na(messplot$mess) | messplot$mess == Inf] <- NA

#summary(as.data.frame(messplot))
display.brewer.all(n=8, type="all", colorblindFriendly=T)


bio <- stack("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc")
             
    
plot(bio$bio_2, col=brewer.pal(8, name="Greys") )     
plot(bio$bio_5, col=brewer.pal(8, name="Greys") )  
plot(bio$bio_8, col=brewer.pal(8, name="Greys") ) 
plot(bio$bio_15, col=brewer.pal(8, name="Greys") ) 

#rev((brewer.pal(3, name="RdYlBu"))))

hi_abun_hotspot <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_HI_ABUN.asc")
full_hotspot <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_FULL.asc")
compare <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc")

eco2 <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco")
eco2 <- spTransform(eco2, "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
crs(eco2) <- NA
plot(eco2)


maxentRange2 <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/PUMOL.asc")
names(maxentRange2) <- "layer"
maxentRange2 <- rasterToPolygons(maxentRange2, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)


maxentRange3 <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/CYSC4.asc")
plot(maxentRange3)
names(maxentRange3) <- "layer"
maxentRange3 <- rasterToPolygons(maxentRange3, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)



ordinal2 <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/PUMOL.asc") 
ordinal2$layer <- 0
ordinal2$layer <- ordinal2$PUMOL
ordinal2$layer[ordinal2$PUMOL==0] <- NA

ordinal2$layer2 <- 0
ordinal2$layer2[ordinal2$PUMOL==3] <- 1
ordinal2$layer2[ordinal2$PUMOL != 3 | is.na(ordinal2$PUMOL)] <- NA

ordinal3 <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CYSC4.asc") 
ordinal3$layer <- 0
ordinal3$layer <- ordinal3$CYSC4
ordinal3$layer[ordinal3$CYSC4==0] <- NA
#3
ordinal3$layer2 <- 0
ordinal3$layer2[ordinal3$CYSC4==3] <- 1
ordinal3$layer2[ordinal3$CYSC4 != 3 | is.na(ordinal3$CYSC4)] <- NA

ordinalLOW <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALPE4.asc")
ordinalLOW$layer <- 0
ordinalLOW$layer <- ordinalLOW$ALPE4
ordinalLOW$layer[ordinalLOW$ALPE4==0 | is.na(ordinalLOW$ALPE4) ] <- NA


ordinalAVG <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc" )
ordinalAVG$layer <- 0
ordinalAVG$layer <- ordinalAVG$LOMA6
ordinalAVG$layer[ordinalAVG$LOMA6==0 | is.na(ordinalAVG$LOMA6) ] <- NA


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ROMU.asc"))


#############################
####### plots for GRC #######
#############################


### abun stack with transparent
### plain field with black points/red points

#plot(occ.ran,lwd=5,col="orangered2", cex=2)
#plot(abun.ran,lwd=5,col="firebrick", cex=2, add=T)
#plot(abun.ran,lwd=5,col="deepskyblue1", cex=2, add=T)

png("C:/Users/Localadmin/Google Drive/GRC_figs/points_occ_BIG.png",width=800,height=500)
plot(full_richness, axes=F, box=F,legend=F)
plot(extentShape, lwd=2, col="grey92", add=T)
plot(occ.pts,lwd=1, col="navyblue",add=T, cex=1)
#plot(abun.ran, pch=16, col="firebrick2", add=T, cex=2.5)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/points_abun_BIG.png",width=800,height=500)
plot(full_richness, axes=F, box=F,legend=F)
plot(extentShape, lwd=2, col="grey92", add=T)
plot(occ.pts,lwd=1,col="navyblue",add=T, cex=1)
plot(abunpts, pch=20, col="firebrick2", add=T, cex=1)
#plot(abun.ran, pch=16, col="firebrick2", add=T, cex=2.5)
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/points.png",width=600,height=400)
plot(full_richness, axes=F, box=F,legend=F)
plot(extentShape, lwd=2, col="grey92", add=T)
plot(occ.ran,lwd=5,col="navyblue",add=T, cex=2)
plot(abun.ran, pch=16, col="firebrick2", add=T, cex=2.5)
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/abun.png",width=600,height=400)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F, cex=2)
plot(extentShape, col="transparent", lwd=2, add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/rich.png",width=600,height=400)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/range.png",width=600,height=400, bg='transparent')
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(maxentRange, add=T, col="mediumblue",border="mediumblue", bg='transparent')
plot(extentShape, col="transparent", add=T, lwd=2, bg='transparent')
#plot(maxentRange, add=T, col="firebrick",border="firebrick")
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/abundanceRange.png",width=600,height=400)
plot(full_richness, axes=F, box=F, legend=F)
plot(extentShape, col="black", lwd=2, add=T)   ### black works too
plot(ordinal$layer, add=T, col= (brewer.pal(8, name="YlOrRd")), legend=F)
plot(extentShape, col="transparent", lwd=1.5, add=T) 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/hiAbunRange.png",width=600,height=400, bg='transparent')
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(ordinal$layer2, col=c("firebrick"), axes=F, box=F, legend=F, add=T, bg='transparent')
plot(extentShape, col="transparent", lwd=2, add=T, bg='transparent') 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/mess.png",width=600,height=400)
plot(full_richness, axes=F, box=F, legend=F)
plot(extentShape, col="grey99",lwd=3, add=T)
plot(messplot$mess2, add=T, col=c("blue","slategray2"),box=F,legend=F, axes=F)
plot(extentShape, col="transparent",lwd=3,add=T)
plot(abun.ran, pch=16, cex=3, col="firebrick2", add=T) 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/hi_hot.png",width=600,height=400)
plot(hi_abun_hotspot, col=c("grey90", "orangered2"), box=F, axes=F, legend=F)
plot(extentShape, col="transparent",lwd=3,add=T)
#plot(eco2, col="transparent",add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/ri_hot.png",width=600,height=400)
plot(full_hotspot, col=c("grey95", "deepskyblue3"), box=F, axes=F, legend=F)
plot(extentShape, col="transparent",lwd=3,add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/compare.png",width=600,height=400)
plot(compare, col=c("orangered2","darkslateblue","grey95", "deepskyblue1"), box=F, axes=F, legend=F)
plot(compare, col=c("transparent",adjustcolor(c("black"),.35), "transparent", "transparent"), add=T, box=F, axes=F, legend=F)
plot(extentShape, col="transparent",lwd=2,add=T)
dev.off()                    #"deepskyblue1","black","grey95", "firebrick1"


png("C:/Users/Localadmin/Google Drive/GRC_figs/abunECO.png",width=600,height=400)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F)
plot(extentShape, col="transparent", lwd=2, add=T)
plot(eco2, col="transparent", add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/richECO.png",width=600,height=400)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), axes=F, box=F)
plot(extentShape, lwd=2, add=T)
plot(eco2, col="transparent", add=T)
dev.off()

#png("C:/Users/Localadmin/Google Drive/GRC_figs/redpoints.png",width=600,height=400)
#plot(full_richness, axes=F, box=F,legend=F)
#plot(extentShape, lwd=2, col="grey92", add=T)
##plot(occ.ran,lwd=3,col="navy",add=T)
#plot(occ.ran, pch=19, col="firebrick2", add=T)
#dev.off()

### E variables
png("C:/Users/Localadmin/Google Drive/GRC_figs/diurnal.png",width=600,height=400, bg='transparent')
plot(bio$bio_2, col=brewer.pal(8, name="Greys") , legend=F, box=F, axes=F, bg='transparent')
plot(extentShape, col="transparent",lwd=2,add=T, bg='transparent')
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/maxtemp.png",width=600,height=400, bg='transparent')
plot(bio$bio_5, col=brewer.pal(8, name="Greys") , legend=F, box=F, axes=F, bg='transparent')
plot(extentShape, col="transparent",lwd=2,add=T, bg='transparent')
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/mintempcold.png",width=600,height=400, bg='transparent')
plot(bio$bio_8, col=brewer.pal(8, name="Greys") , legend=F, box=F, axes=F, bg='transparent')
plot(extentShape, col="transparent",lwd=2,add=T, bg='transparent')
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/precipvar.png",width=600,height=400, bg='transparent')
plot(bio$bio_15, col=brewer.pal(8, name="Greys"), legend=F, box=F, axes=F, bg='transparent') 
plot(extentShape, col="transparent",lwd=2,add=T, bg='transparent')
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/range2.png",width=600,height=400, bg='transparent')
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(maxentRange2, add=T, col="mediumblue",border="mediumblue", bg='transparent')
plot(extentShape, col="transparent", add=T, lwd=2, bg='transparent')
#plot(maxentRange, add=T, col="firebrick",border="firebrick")
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/range3.png",width=600,height=400, bg='transparent')
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(maxentRange3, add=T, col="mediumblue",border="mediumblue", bg='transparent')
plot(extentShape, col="transparent", add=T, lwd=2, bg='transparent')
#plot(maxentRange, add=T, col="firebrick",border="firebrick")
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/hiAbunRange2.png",width=600,height=400, bg='transparent')
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(ordinal2$layer2, col=c("firebrick"), axes=F, box=F, legend=F, add=T, bg='transparent')
plot(extentShape, col="transparent", lwd=2, add=T, bg='transparent') 
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/hiAbunRange3.png",width=600,height=400, bg='transparent')
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(ordinal3$layer2, col=c("firebrick"), axes=F, box=F, legend=F, add=T, bg='transparent')
plot(extentShape, col="transparent", lwd=2, add=T, bg='transparent') 
dev.off()


####### big figs

png("C:/Users/Localadmin/Google Drive/GRC_figs/points_occ_BIG.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F,legend=F, bg="transparent")
plot(extentShape, lwd=2, col="grey92", add=T, bg="transparent")
plot(occ.pts,lwd=1, col="navyblue",add=T, cex=1, bg="transparent")
#plot(abun.ran, pch=16, col="firebrick2", add=T, cex=2.5)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/points_abun_BIG.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F,legend=F, bg="transparent")
plot(extentShape, lwd=2, col="grey92", add=T, bg="transparent")
plot(occ.pts,lwd=1,col="navyblue",add=T, cex=1, bg="transparent")
plot(abunpts, pch=20, col="firebrick2", add=T, cex=1, bg="transparent")
#plot(abun.ran, pch=16, col="firebrick2", add=T, cex=2.5)
dev.off()



png("C:/Users/Localadmin/Google Drive/GRC_figs/abunBIG.png",width=1200,height=700, bg="transparent")
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F, cex=5, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent")
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/richBIG.png",width=1200,height=700, bg="transparent")
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), cex=2, axes=F, box=F, bg="transparent")
plot(extentShape, lwd=2, add=T, bg="transparent")
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/rangeBIG.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F, legend=F, bg='transparent')
plot(extentShape, col="grey92", add=T, bg='transparent')
plot(maxentRange, add=T, col="mediumblue",border="mediumblue", bg='transparent')
plot(extentShape, col="transparent", add=T, lwd=2, bg='transparent')
#plot(maxentRange, add=T, col="firebrick",border="firebrick")
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/abundanceRangeBIG.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="black", lwd=2, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= (brewer.pal(8, name="YlOrRd")), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=1.5, add=T, bg="transparent") 
dev.off()



png("C:/Users/Localadmin/Google Drive/GRC_figs/abundanceRangeLOW.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="black", lwd=2, add=T, bg="transparent")   ### black works too
plot(ordinalLOW$layer, add=T, col= (brewer.pal(8, name="YlOrRd")), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=1.5, add=T, bg="transparent") 
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/abundanceRangeAVG.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="black", lwd=2, add=T, bg="transparent")   ### black works too
plot(ordinalAVG$layer, add=T, col= (brewer.pal(8, name="YlOrRd")), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=1.5, add=T, bg="transparent") 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/messBIG.png",width=1200,height=700, bg="transparent")
plot(full_richness, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99",lwd=3, add=T, bg="transparent")
plot(messplot$mess2, add=T, col=c("blue","slategray2"),box=F,legend=F, axes=F, bg="transparent")
plot(extentShape, col="transparent",lwd=3,add=T, bg="transparent")
plot(abun.ran, pch=16, cex=3, col="firebrick2", add=T, bg="transparent") 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/hi_hotBIG.png",width=1200,height=700, bg="transparent")
plot(hi_abun_hotspot, col=c("grey90", "orangered2"), box=F, axes=F, legend=F, bg="transparent")
plot(extentShape, col="transparent",lwd=3,add=T, bg="transparent")
#plot(eco2, col="transparent",add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/ri_hotBIG.png",width=1200,height=700, bg="transparent")
plot(full_hotspot, col=c("grey95", "deepskyblue3"), box=F, axes=F, legend=F, bg="transparent")
plot(extentShape, col="transparent",lwd=3,add=T, bg="transparent")
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/compareBIG.png",width=1200,height=700, bg="transparent")
plot(compare, col=c("orangered2","darkslateblue","grey95", "deepskyblue1"), box=F, axes=F, legend=F, bg="transparent")
plot(compare, col=c("transparent",adjustcolor(c("black"),.35), "transparent", "transparent"), add=T, box=F, axes=F, legend=F, bg="transparent")
plot(extentShape, col="transparent",lwd=2,add=T, bg="transparent")
dev.off()                    #"deepskyblue1","black","grey95", "firebrick1"



########## 3MT
extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

hi_abun_hotspot <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_HI_ABUN.asc")





png("C:/Users/Localadmin/Google Drive/GRC_figs/3MThi_hot.png",width=500,height=300, bg="transparent")
plot(hi_abun_hotspot, col=c("grey99", "firebrick3"), box=F, axes=F, legend=F, bg="transparent")
#plot(eco2, col="transparent",add=T)
dev.off()


#knotweed at high abundance
ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$POCU6
ordinal$layer[ordinal$POCU6==2] <- 1
ordinal$layer[is.na(ordinal$POCU6)] <- NA

#ordinal$layer2 <- 0
#ordinal$layer2[ordinal$BRPA4==3] <- 1
#ordinal$layer2[ordinal$BRPA4 != 3 | is.na(ordinal$BRPA4)] <- NA
png("C:/Users/Localadmin/Google Drive/GRC_figs/3MThi_knotGREEN.png",width=600,height=400, bg="transparent")
plot(ordinal$layer, col= c("white", "white", "darkgreen"), box=F, axes=F, legend=F )
plot(extentShape, col="transparent",lwd=3,add=T, bg="transparent")
dev.off()

#png("C:/Users/Localadmin/Google Drive/GRC_figs/3MThi_knot.png",width=500,height=300, bg="transparent")
#plot(ordinal$layer, col= c("white", "white", "darkgreen"), box=F, axes=F, legend=F )
#plot(extentShape, col="transparent",lwd=3,add=T, bg="transparent")
#dev.off()






library(rgdal)
library(raster)
library(maptools)
library(rgeos)
#install.packages("RColorBrewer")
library(RColorBrewer)


ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP<0.05 & ordsums$AUC > 0.7 & ordsums$MESSgn10 > .9,]
ordsums2[ordsums2$no.pts.full < 1000 & ordsums2$range_prop_area > .20,]

display.brewer.all(n=8, type="all", colorblindFriendly=T)

compare <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc")
hi_abun_richness <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
full_richness  <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")
props <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hi_abun_proportion.asc")


plot(extentShape, lwd=5)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), add=T)
plot(extentShape, lwd=5)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), add=T)


extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

abunpts <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv",
                      sep=",", header=T, stringsAsFactors = F)
head(abunpts)
abunpts <- abunpts[abunpts$species=="BRPA4",]
coordinates(abunpts) <- c(5,4)

occ.pts <- read.table("C:/Users/Localadmin/Documents/MaxEnt_modeling/species/BRPA4.csv", sep=",", stringsAsFactors = F, header=T)
str(occ.pts)
coordinates(occ.pts) <- c(2,3)

#occ.ran <- occ.pts[sample(1:nrow(occ.pts), 130),]
#abun.ran <- abunpts[sample(1:nrow(abunpts), 50),]
#abun.ran <- occ.ran[sample(1:nrow(occ.ran), 40),]

#occ.ran <- occ.pts[sample(1:nrow(occ.pts), 20)]

plot(full_richness)
plot(extentShape, lwd=3, col="grey92", add=T)
plot(occ.ran,lwd=3,col="grey20",add=T)
plot(abun.ran, pch=19, col="firebrick2", add=T)



maxentRange <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")
plot(maxentRange)
names(maxentRange) <- "layer"
maxentRange <- rasterToPolygons(maxentRange, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)

plot(extentShape, col="grey92")
plot(maxentRange, add=T, col="mediumblue",border="mediumblue")
plot(maxentRange, add=T, col="firebrick",border="firebrick")


#ordinal <- readAsciiGrid("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$BRPA4
ordinal$layer[ordinal$BRPA4==0] <- NA

ordinal$layer2 <- 0
ordinal$layer2[ordinal$BRPA4==3] <- 1
ordinal$layer2[ordinal$BRPA4 != 3 | is.na(ordinal$BRPA4)] <- NA


display.brewer.all(n=8, type="all", colorblindFriendly=T)
plot(extentShape, col="darkslategrey", lwd=2)   ### black works too
plot(ordinal$layer, add=T, col= (brewer.pal(3, name="YlOrRd")))

plot(extentShape, col="grey92")
plot(ordinal$layer2, col=c("firebrick"), add=T)






bio <- stack(#"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc")
  "C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")

#bio$bio_2[bio$BRPA4 == 0] <- NA
bio$bio_8[bio$BRPA4 == 0] <- NA
bio$nlcd_6[bio$BRPA4 == 0] <- NA
#bio$bio_15[bio$BRPA4 == 0] <- NA

#bio$bio_2 <- bio$bio_2/10
#bio$bio_5 <- bio$bio_5/10
bio$bio_8 <- bio$bio_8/10
#bio$bio_6 <- bio$bio_6/10

bio <- dropLayer(bio,3)



messpts <- data.frame(cbind(abunpts$bio_8, abunpts$nlcd_6))#), abunpts$nlcd_5, abunpts$nlcd_8))
head(messpts)
colnames(messpts)<-c("bio_8", "nlcd6" )#, "nlcd_5", "nlcd_8")

library(dismo)
messplot <- mess(bio,messpts,full=F)
#plot(messplot)

messplot$mess2 <- -99
#messplot$mess2[messplot$mess >= -10 & messplot$mess <= 10] <- 2
messplot$mess2[messplot$mess < -10] <- 1
messplot$mess2[messplot$mess >= -10] <- 3
messplot$mess2[is.na(messplot$mess) | messplot$mess == Inf] <- NA



display.brewer.all(n=8, type="all", colorblindFriendly=T)

plot(props, col= (brewer.pal("PuRd",n=9)), box=F, axes=F)
plot(extentShape, col="transparent",lwd=2,add=T, bg="transparent")


plot(full_richness)

###############
### THESIS FIGS NEW
###############

library(rgdal)
library(raster)
library(maptools)
library(rgeos)
#install.packages("RColorBrewer")
library(RColorBrewer)


#ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP<0.05 & ordsums$AUC > 0.7 & ordsums$MESSgn10 > .9,]
#ordsums2[ordsums2$no.pts.full < 1000 & ordsums2$range_prop_area > .20,]

display.brewer.all(n=8, type="all", colorblindFriendly=T)

compare <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc")
hi_abun_richness <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
full_richness  <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")
props <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hi_abun_proportion.asc")

compare <- projectRaster(compare, crs="+init=epsg:5070" )
hi_abun_richness <- projectRaster(hi_abun_richness , crs="+init=epsg:5070" )
full_richness <- projectRaster(full_richness, crs="+init=epsg:5070" )
props <- projectRaster(props, crs="+init=epsg:5070" )

#plot(extentShape, lwd=5)
#plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), add=T)
#plot(extentShape, lwd=5)
#plot(full_richness, col=brewer.pal(8, name="YlGnBu"), add=T)


extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

abunpts <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv",
                      sep=",", header=T, stringsAsFactors = F)
head(abunpts)
abunpts <- abunpts[abunpts$species=="BRPA4",]
coordinates(abunpts) <- c(5,4)

occ.pts <- read.table("C:/Users/Localadmin/Documents/MaxEnt_modeling/species/BRPA4.csv", sep=",", stringsAsFactors = F, header=T)
str(occ.pts)
coordinates(occ.pts) <- c(2,3)


proj4string(extentShape) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(abun.pts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(occ.pts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

extentShape <- spTransform(extentShape, "+init=epsg:5070" )
abun.pts <- spTransform(abun.pts, "+init=epsg:5070" )
occ.pts <- spTransform(occ.pts, "+init=epsg:5070" )
#occ.ran <- occ.pts[sample(1:nrow(occ.pts), 130),]
#abun.ran <- abunpts[sample(1:nrow(abunpts), 50),]
#abun.ran <- occ.ran[sample(1:nrow(occ.ran), 40),]

#occ.ran <- occ.pts[sample(1:nrow(occ.pts), 20)]

#plot(full_richness)
#plot(extentShape, lwd=3, col="grey92", add=T)
#plot(occ.ran,lwd=3,col="grey20",add=T)
#plot(abun.ran, pch=19, col="firebrick2", add=T)



maxentRange <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")
#plot(maxentRange)
names(maxentRange) <- "layer"
maxentRange <- rasterToPolygons(maxentRange, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
proj4string(maxentRange) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
maxentRange <- spTransform(maxentRange,"+init=epsg:5070")

#plot(extentShape, col="grey92")
#plot(maxentRange, add=T, col="mediumblue",border="mediumblue")
#plot(maxentRange, add=T, col="firebrick",border="firebrick")


#ordinal <- readAsciiGrid("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$BRPA4
ordinal$layer[ordinal$BRPA4==0] <- NA

ordinal$layer2 <- 0
ordinal$layer2[ordinal$BRPA4==3] <- 1
ordinal$layer2[ordinal$BRPA4 != 3 | is.na(ordinal$BRPA4)] <- NA

proj4string(ordinal) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ordinal <- projectRaster(ordinal, crs="+init=epsg:5070" )

display.brewer.all(n=8, type="all", colorblindFriendly=T)
#plot(extentShape, col="darkslategrey", lwd=2)   ### black works too
#plot(ordinal$layer, add=T, col= (brewer.pal(3, name="YlOrRd")))
#
#plot(extentShape, col="grey92")
#plot(ordinal$layer2, col=c("firebrick"), add=T)






bio <- stack(#"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
  "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
  #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc")
  "C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")

proj4string(bio) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
bio <- projectRaster(bio, crs="+init=epsg:5070" )

#bio$bio_2[bio$BRPA4 == 0] <- NA
bio$bio_8[bio$BRPA4 == 0] <- NA
bio$nlcd_6[bio$BRPA4 == 0] <- NA
#bio$bio_15[bio$BRPA4 == 0] <- NA

#bio$bio_2 <- bio$bio_2/10
#bio$bio_5 <- bio$bio_5/10
bio$bio_8 <- bio$bio_8/10
#bio$bio_6 <- bio$bio_6/10

bio <- dropLayer(bio,3)



messpts <- data.frame(cbind(abunpts$bio_8, abunpts$nlcd_6))#), abunpts$nlcd_5, abunpts$nlcd_8))
head(messpts)
colnames(messpts)<-c("bio_8", "nlcd6" )#, "nlcd_5", "nlcd_8")

library(dismo)
messplot <- mess(bio,messpts,full=F)
#plot(messplot)

messplot$mess2 <- -99
#messplot$mess2[messplot$mess >= -10 & messplot$mess <= 10] <- 2
messplot$mess2[messplot$mess < -10] <- 1
messplot$mess2[messplot$mess >= -10] <- 3
messplot$mess2[is.na(messplot$mess) | messplot$mess == Inf] <- NA



display.brewer.all(n=8, type="all", colorblindFriendly=T)

#plot(props, col= (brewer.pal("PuRd",n=9)), box=F, axes=F)
#plot(extentShape, col="transparent",lwd=2,add=T, bg="transparent")



#######################################################





png("C:/Users/Localadmin/Google Drive/GRC_figs/points_in_range.png",width=800,height=500)
plot(full_richness, axes=F, box=F,legend=F, col="grey92")
plot(maxentRange, add=T, col="deepskyblue",border="transparent", bg='transparent')
plot(extentShape, lwd=2, col="transparent",add=T)
plot(occ.pts,col="black",add=T, cex=1, lwd=1)
plot(abunpts, pch=16, col="firebrick1", add=T, cex=1)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/abundanceRange.png",width=800,height=500)
plot(full_richness, axes=F, box=F, legend=F)
plot(extentShape, col="black", lwd=2, add=T)   ### black works too
plot(ordinal$layer, add=T, col= (brewer.pal(8, name="YlOrRd")), legend=F)
plot(extentShape, col="transparent", lwd=1.5, add=T) 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/mess.png",width=800,height=500)
plot(full_richness, axes=F, box=F, legend=F)
plot(extentShape, col="grey92",lwd=1, add=T)
plot(messplot$mess2, add=T, col=c("blue","slategray2"),box=F,legend=F, axes=F)
plot(extentShape, col="transparent",lwd=3,add=T)
#plot(abun.ran, pch=16, cex=3, col="firebrick2", add=T) 
dev.off()


png("C:/Users/Localadmin/Google Drive/GRC_figs/compare.png",width=800,height=500)
plot(compare, col=c("orangered2","darkslateblue","grey95", "deepskyblue1"), box=F, axes=F, legend=F)
plot(compare, col=c("transparent",adjustcolor(c("black"),.35), "transparent", "transparent"), add=T, box=F, axes=F, legend=F)
plot(extentShape, col="transparent",lwd=2,add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/compareECO.png",width=800,height=500)
plot(compare, col=c("orangered2","darkslateblue","grey95", "deepskyblue1"), box=F, axes=F, legend=F)
plot(compare, col=c("transparent",adjustcolor(c("black"),.35), "transparent", "transparent"), add=T, box=F, axes=F, legend=F)
plot(extentShape, col="transparent",lwd=2,add=T)
plot(eco, lwd=1,col="transparent", add=T)  
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/abun.png",width=800,height=500)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F, cex=2)
plot(extentShape, col="transparent", lwd=2, add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/abunECO.png",width=800,height=500)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F, cex=2)
plot(extentShape, col="transparent", lwd=2, add=T)
plot(eco, col="transparent", add=T) 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/rich.png",width=800,height=500)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/richECO.png",width=800,height=500)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
plot(eco, col="transparent", add=T) 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/props.png",width=800,height=500)
plot(props, col=brewer.pal(9, name="PuRd"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
#plot(eco, col="transparent", add=T) 
dev.off()

png("C:/Users/Localadmin/Google Drive/GRC_figs/propsECO.png",width=800,height=500)
plot(props, col=brewer.pal(8, name="PuRd"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
plot(eco, col="transparent", add=T) 
dev.off()



##############################
































##################
##### THESIS FIGS THIS ONE
##################

library(raster)

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_8_2018.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

png("C:/Users/Localadmin/GRC_figs/infilling.png", width=700, height=300)
hist(ordsums$infilling[ordsums$AUC>.7 &ordsums$kappaP<0.05 & ordsums$kappa >0 & ordsums$MESSgn10 >10],
     main="", ylab="Number of species", xlab="proportion of occupied cells",pch=2,cex.lab=1.5,col="lightblue")
abline(h=0)
dev.off()


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/FRAL4.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc"))


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/NERE.asc")) ###
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CAGL11.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LYMI.asc"))


plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ELAN.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ARMI2.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/TRTE.asc"))  ###
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/TRTE.asc"))



plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BETH.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/PUMOL.asc"))###
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALJU.asc")) ####

plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CYOF.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ONAC.asc")) ###
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CEOR7.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CHMA2.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALPE4.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/MEAZ.asc"))
plot(raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ELUM.asc"))


library(rgdal)
library(raster)
library(maptools)
library(rgeos)
#install.packages("RColorBrewer")
library(RColorBrewer)


ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP<0.05 & ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > .9,]
ordsums2[ordsums2$no.pts.full < 1000 & ordsums2$range_prop_area > .20,]

display.brewer.all(n=8, type="all", colorblindFriendly=T)

compare <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc")
hi_abun_richness <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
full_richness  <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")
props <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hi_abun_proportion.asc")


#plot(extentShape, lwd=5)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), add=T)
plot(extentShape, lwd=5)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), add=T)


extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

abunpts <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv",
                      sep=",", header=T, stringsAsFactors = F)
head(abunpts)
abunpts <- abunpts[abunpts$species=="BRPA4",]
coordinates(abunpts) <- c(5,4)

occ.pts <- read.table("C:/Users/Localadmin/Documents/MaxEnt_modeling/species/BRPA4.csv", sep=",", stringsAsFactors = F, header=T)
str(occ.pts)
coordinates(occ.pts) <- c(2,3)

occ.ran <- occ.pts[sample(1:nrow(occ.pts), 130),]
abun.ran <- abunpts[sample(1:nrow(abunpts), 50),]
#abun.ran <- occ.ran[sample(1:nrow(occ.ran), 40),]

#occ.ran <- occ.pts[sample(1:nrow(occ.pts), 20)]

plot(full_richness)
plot(extentShape, lwd=3, col="grey92", add=T)
plot(occ.ran,lwd=3,col="grey20",add=T)
plot(abun.ran, pch=19, col="firebrick2", add=T)



maxentRange <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")
plot(maxentRange)
names(maxentRange) <- "layer"
maxentRange <- rasterToPolygons(maxentRange, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)

plot(extentShape, col="grey92")
plot(maxentRange, add=T, col="mediumblue",border="mediumblue")
plot(maxentRange, add=T, col="firebrick",border="firebrick")


#ordinal <- readAsciiGrid("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$BRPA4
ordinal$layer[ordinal$BRPA4==0] <- NA

ordinal$layer2 <- 0
ordinal$layer2[ordinal$BRPA4==3] <- 1
ordinal$layer2[ordinal$BRPA4 != 3 | is.na(ordinal$BRPA4)] <- NA


display.brewer.all(n=8, type="all", colorblindFriendly=T)
plot(extentShape, col="darkslategrey", lwd=2)   ### black works too
plot(ordinal$layer, add=T, col= (brewer.pal(3, name="YlOrRd")))

plot(extentShape, col="grey92")
plot(ordinal$layer2, col=c("firebrick"), add=T)






bio <- stack(#"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
            #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc")
             "C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")

#bio$bio_2[bio$BRPA4 == 0] <- NA
bio$bio_8[bio$BRPA4 == 0] <- NA
bio$nlcd_6[bio$BRPA4 == 0] <- NA
#bio$bio_15[bio$BRPA4 == 0] <- NA

#bio$bio_2 <- bio$bio_2/10
#bio$bio_5 <- bio$bio_5/10
bio$bio_8 <- bio$bio_8/10
#bio$bio_6 <- bio$bio_6/10

bio <- dropLayer(bio,3)



messpts <- data.frame(cbind(abunpts$bio_8, abunpts$nlcd_6))#), abunpts$nlcd_5, abunpts$nlcd_8))
head(messpts)
colnames(messpts)<-c("bio_8", "nlcd6" )#, "nlcd_5", "nlcd_8")

library(dismo)
messplot <- mess(bio,messpts,full=F)
#plot(messplot)

messplot$mess2 <- -99
#messplot$mess2[messplot$mess >= -10 & messplot$mess <= 10] <- 2
messplot$mess2[messplot$mess < -10] <- 1
messplot$mess2[messplot$mess >= -10] <- 3
messplot$mess2[is.na(messplot$mess) | messplot$mess == Inf] <- NA



display.brewer.all(n=8, type="all", colorblindFriendly=T)

plot(props, col= (brewer.pal("PuRd",n=9)), box=F, axes=F)
plot(extentShape, col="transparent",lwd=2,add=T, bg="transparent")


plot(full_richness)

proj4string(full_richness) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(ordinal)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(messplot)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(compare)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(hi_abun_richness)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

full_richness <- projectRaster(full_richness, crs="+init=epsg:5070")
ordinal <- projectRaster(ordinal, crs="+init=epsg:5070")
messplot <- projectRaster(messplot, crs="+init=epsg:5070")
compare <- projectRaster(compare, crs="+init=epsg:5070")
hi_abun_richness <- projectRaster(hi_abun_richness, crs="+init=epsg:5070")
proj4string(props) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
props<- projectRaster(props,crs="+init=epsg:5070")

eco <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco", stringsAsFactors = F)

proj4string(occ.pts) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(abunpts) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(extentShape) <-"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(eco) <- "+init=epsg:5070"

occ.pts <- spTransform(occ.pts,"+init=epsg:5070")
abunpts <- spTransform(abunpts,"+init=epsg:5070")
extentShape <- spTransform(extentShape,"+init=epsg:5070")
#eco <- spTransform(eco,"+init=epsg:5070")
proj4string(maxentRange) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
maxentRange <- spTransform(maxentRange,"+init=epsg:5070")

#install.packages('GISTools')
#install.packages('maps')
library(GISTools)
library(maps)


###############
### THESIS FIGS YES THIS ONE
###############

options(scipen = 999)



png("C:/Users/Localadmin/Google Drive/figs/points_in_range.png",width=800,height=500)
plot(full_richness, axes=F, box=F,legend=F, col="grey92")
plot(maxentRange, add=T,col="deepskyblue",border="transparent", bg='transparent')
plot(extentShape, lwd=2, col="transparent",add=T)
plot(occ.pts,col="black",add=T, cex=1, lwd=1)
plot(abunpts, pch=16, col="firebrick1", add=T, cex=1)
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

png("C:/Users/Localadmin/Google Drive/figs/abundanceRange.png",width=800,height=500)
plot(full_richness, axes=F, box=F, legend=F)
plot(extentShape, col="grey99", lwd=2, add=T)   ### black works too
plot(ordinal$layer, add=T, col= c("yellow","orange","firebrick3"), legend=F)
plot(extentShape, col="transparent", lwd=1.5, add=T) 
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

png("C:/Users/Localadmin/Google Drive/figs/mess.png",width=800,height=500)
plot(full_richness, axes=F, box=F, legend=F)
plot(extentShape, col="grey99",lwd=1, add=T)
plot(messplot$mess2, add=T, col=c("red","slategray1"),box=F,legend=F, axes=F)
plot(extentShape, col="transparent",lwd=3,add=T)
#plot(abun.ran, pch=16, cex=3, col="firebrick2", add=T) 
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()


png("C:/Users/Localadmin/Google Drive/figs/compare.png",width=800,height=500)
plot(compare, col=c("orangered2","darkslateblue","grey98", "deepskyblue1"), box=F, axes=F, legend=F)
plot(compare, col=c("transparent",adjustcolor(c("black"),.15), "transparent", "transparent"), add=T, box=F, axes=F, legend=F)
#plot(extentShape, col="transparent",lwd=2,add=T)
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
plot(eco,add=T,col="transparent",lwd=2)
dev.off()

#png("C:/Users/Localadmin/Google Drive/GRC_figs/compareECO.png",width=800,height=500)
#plot(compare, col=c("orangered2","darkslateblue","grey95", "deepskyblue1"), box=F, axes=F, legend=F)
#plot(compare, col=c("transparent",adjustcolor(c("black"),.35), "transparent", "transparent"), add=T, box=F, axes=F, legend=F)
#plot(extentShape, col="transparent",lwd=2,add=T)
#plot(eco, lwd=1,col="transparent", add=T)  
#dev.off()

#png("C:/Users/Localadmin/Google Drive/GRC_figs/abun.png",width=800,height=500)
#plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F, cex=2)
#plot(extentShape, col="transparent", lwd=2, add=T)
#map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
#dev.off()

png("C:/Users/Localadmin/Google Drive/figs/abunECO.png",width=800,height=500)
plot(hi_abun_richness, col=brewer.pal(8, name="YlOrRd"), axes=F, box=F, cex=2)
plot(extentShape, col="transparent", lwd=2, add=T)
plot(eco, col="transparent", add=T) 
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

#png("C:/Users/Localadmin/Google Drive/GRC_figs/rich.png",width=800,height=500)
#plot(full_richness, col=brewer.pal(8, name="YlGnBu"), cex=2, axes=F, box=F)
#plot(extentShape, lwd=2, add=T)
#map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
#dev.off()

png("C:/Users/Localadmin/Google Drive/figs/richECO.png",width=800,height=500)
plot(full_richness, col=brewer.pal(8, name="YlGnBu"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
plot(eco, col="transparent", add=T) 
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

#png("C:/Users/Localadmin/Google Drive/GRC_figs/props.png",width=800,height=500)
#plot(props, col=brewer.pal(9, name="PuRd"), cex=2, axes=F, box=F)
#plot(extentShape, lwd=2, add=T)
##plot(eco, col="transparent", add=T) 
#dev.off()

png("C:/Users/Localadmin/Google Drive/figs/propsECO.png",width=800,height=500)
plot(props, col=brewer.pal(8, name="PuRd"), cex=2, axes=F, box=F)
plot(extentShape, lwd=2, add=T)
plot(eco, col="transparent", add=T) 
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()
#### props


### RANGES: NERE, LOMA6, TRTE
### ABUN: PUMOL,ALJU,LOMA6,OMAC
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(RColorBrewer)
library(GISTools)
library(maps)




extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

#proj4string(extentShape) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
extentShape <- spTransform(extentShape, "+init=epsg:5070")
#plot(extentShape)


ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/PUMOL.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_21_PUMOL.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALJU.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_21_ALJU.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_21_LOMA6.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ONAC.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_21_ONAC.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()


#############################################

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/NERE.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
ordinal$layer[is.na(ordinal$sp==0)] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")
plot(ordinal)

png("C:/Users/Localadmin/Google Drive/figs/Range_5_21_NERE.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("deepskyblue3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/LOMA6.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
ordinal$layer[is.na(ordinal$sp==0)] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")
plot(ordinal)

png("C:/Users/Localadmin/Google Drive/figs/Range_5_21_LOMA6.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("deepskyblue3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/TRTE.asc") 
names(ordinal) <- "sp"
ordinal$layer <- 0
ordinal$layer <- ordinal$sp
ordinal$layer[ordinal$sp==0] <- NA
ordinal$layer[is.na(ordinal$sp==0)] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")
plot(ordinal)

png("C:/Users/Localadmin/Google Drive/figs/Range_5_21_TRTE.png",width=800,height=500, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=0.1, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("deepskyblue3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=2, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

#### redoing some presentation figs for Thesis ###########
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(RColorBrewer)
library(GISTools)
library(maps)

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_8_2018.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

ordsums2 <- ordsums[ordsums$kappaP < 0.05 & ordsums$kappa > 0 & ordsums$MESS_abun2range_SEL > .9 & ordsums$AUC > .7,]
hist(ordsums2$impac_estab)
summary(ordsums2$impac_estab)
ordsums2$species.code[ordsums2$impac_estab < 0.03 & ordsums2$impac_estab > 0  & ordsums2$kappa > 0.1]
ordsums2$species.code[ordsums2$impac_estab > 0.5 & ordsums2$kappa > 0.1]
ordsums2$species.code[ordsums2$impac_estab > 0.09 & ordsums2$impac_estab < 0.161 & ordsums2$kappa > 0.1]



plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/FRAL4.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CHMA2.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/VETH.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEFU.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/EUAL13.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/SOHA.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/TUFA.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/EUAL13.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LYJA.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ROMU.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALJU.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ARMI2.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/SEPU7.asc")) #NO
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ELAN.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEMA3.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMO2.asc")) ####
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc"))

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc"))

ordsums2[ordsums2$species.code=="POCU6",]

##### I NEED: 
## unprojected POCU6
## projected POCU6
## projected low
## projected med


##
ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$POCU6
ordinal$layer[ordinal$POCU6==0] <- NA

extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_14_POCU6_un.png",width=800,height=500)
plot(ordinal$layer, axes=F, box=F, legend=F)
plot(extentShape, col="grey99", lwd=1.5, add=T)   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod1","orange","firebrick3"), legend=F)
plot(extentShape, col="transparent", lwd=2, add=T) 
#map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()


proj4string(extentShape) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
extentShape <- spTransform(extentShape, "+init=epsg:5070")
plot(extentShape)

proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")
options(scipen = 999)


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_14_POCU.png",width=1600,height=1000)
plot(ordinal, axes=F, box=F, legend=F)
plot(extentShape, col="grey99", lwd=1.5, add=T)   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod1","orange","firebrick3"), legend=F)
plot(extentShape, col="transparent", lwd=1.5, add=T) 
map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_18_POCU.png",width=1600,height=1000, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=1.5, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=1.3, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEFU.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$HEFU
ordinal$layer[ordinal$HEFU==0] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")

png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_14_HEFU.png",width=1600,height=1000)
plot(ordinal, axes=F, box=F, legend=F)
plot(extentShape, col="grey99", lwd=1.5, add=T)   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod1","orange","firebrick3"), legend=F)
plot(extentShape, col="transparent", lwd=1.3, add=T) 
map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_18_HEFU.png",width=1600,height=1000, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=1.5, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=1.3, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()




ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$LOMA6
ordinal$layer[ordinal$LOMA6==0] <- NA
proj4string(ordinal) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal$layer, crs="+init=epsg:5070")

png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_14_LOMA6.png",width=1600,height=1000)
plot(ordinal, axes=F, box=F, legend=F)
plot(extentShape, col="grey99", lwd=1.5, add=T)   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod1","orange","firebrick3"), legend=F)
plot(extentShape, col="transparent", lwd=1.3, add=T) 
map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()


png("C:/Users/Localadmin/Google Drive/figs/abundanceRange_5_18_LOMA.png",width=1600,height=1000, bg="transparent")
plot(ordinal, axes=F, box=F, legend=F, bg="transparent")
plot(extentShape, col="grey99", lwd=1.5, add=T, bg="transparent")   ### black works too
plot(ordinal$layer, add=T, col= c("lightgoldenrod2","lightgoldenrod2","firebrick3"), legend=F, bg="transparent")
plot(extentShape, col="transparent", lwd=1.3, add=T, bg="transparent") 
#map.scale(x=-2.7e+06,y=800000,ratio=F,relwidth=0.2,cex=1.5,bg="transparent")
#north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
dev.off()

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
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc")
plot(bio)
plot(bio$nlcd_5)
##############################





########### look at some popular nh invasives ###############
abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)
abun_all <- data.frame(cbind(abun_all$species,abun_all$abundance,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- c("PLANT_CODE", "ABUNDANCE", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")
abun_all <- abun_all[abun_all$PLANT_CODE %in% ordsums2$species.code,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL)
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)

coordinates(abun_all) <- c(3,4)





plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/FRAL4.asc")) 
plot(abun_all[abun_all$PLANT_CODE == "FRAL4" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")
plot(abun_all[abun_all$PLANT_CODE == "FRAL4" & abun_all$ABUNDANCE == 2,], add= T, pch=20, col="pink")
plot(abun_all[abun_all$PLANT_CODE == "FRAL4" & abun_all$ABUNDANCE == 1,], add= T, pch=20, col="white")

read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
           comment.char= "", stringsAsFactors = F, strip.white = T)





plot(readAsciiGrid("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/EUAL13.asc"))
plot(abun_all[abun_all$PLANT_CODE == "EUAL13" & abun_all$ABUNDANCE == 1,], add= T, pch=19, col="red")
eual13 <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)

eual13 <- eual13[eual13$species== "EUAL13",]
eual13$abundance <- ordered(eual13$abundance, levels=c(1,2,3))

ordsums[ordsums$species.code=="EUAL13" | ordsums$species.code == "FRAL4",]
library(ordinal)
M <- clm(abundance~poly(bio_5, 2) + bio_12 + poly(bio_15, 2) + nlcd_3 + nlcd_4 + poly(nlcd_5, 2) + poly(nlcd_6, 2), data=eual13, Hess=T) 
summary(M)

plot(bio$bio_12)
plot(bio$nlcd_4,add=T)
plot(bio)


library(maptools)
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc"))
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOMA6.asc"))

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALPE4.asc"))



plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOJA.asc"))
plot(abun_all[abun_all$PLANT_CODE == "LOJA" & abun_all$ABUNDANCE == 3,], add= T, pch=19, col="red")
plot(abun_all[abun_all$PLANT_CODE == "LOJA" & abun_all$ABUNDANCE == 2,], add= T, pch=19, col="pink")
plot(abun_all[abun_all$PLANT_CODE == "LOJA" & abun_all$ABUNDANCE == 1,], add= T, pch=19, col="white")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/PUMOL.asc"))

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ROMU.asc"), col=rev(heat.colors(100)))  
plot(abun_all[abun_all$PLANT_CODE == "ROMU" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="purple")
plot(abun_all[abun_all$PLANT_CODE == "ROMU" & abun_all$ABUNDANCE == 2,], add= T, pch=20, col="blue")
plot(abun_all[abun_all$PLANT_CODE == "ROMU" & abun_all$ABUNDANCE == 1,], add= T, pch=20, col="green")


plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/RHCA3.asc"))
plot(abun_all[abun_all$PLANT_CODE == "LOJA" & abun_all$ABUNDANCE == 3,], add= T, pch=19, col="red")


plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ACPL.asc"))  
plot(abun_all[abun_all$PLANT_CODE == "ACPL" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")


plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BETH.asc"))
plot(abun_all[abun_all$PLANT_CODE == "BETH" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEHE.asc")) 
plot(abun_all[abun_all$PLANT_CODE == "HEHE" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/MIVI.asc")) 
plot(abun_all[abun_all$PLANT_CODE == "MIVI" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc"))
plot(abun_all[abun_all$PLANT_CODE == "POCU6" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEMA3.asc"))
plot(abun_all[abun_all$PLANT_CODE == "HEMA3" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ELUM.asc"))
plot(abun_all[abun_all$PLANT_CODE == "ELUM" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/CEOR7.asc"))
plot(abun_all[abun_all$PLANT_CODE == "CEOR7" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")

############## old plotting ###############
#library(rgdal)
#library(sp)
#library(maptools)
#library(raster)


species <- abun_all[abun_all$PLANT_CODE=="POCU6",]
plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc"))
plot(abun_all[abun_all$PLANT_CODE == "POCU6" & abun_all$ABUNDANCE == 3,], add= T, pch=20, col="red")
plot(abun_all[abun_all$PLANT_CODE == "POCU6" & abun_all$ABUNDANCE == 2,], add= T, pch=20, col="grey")
plot(abun_all[abun_all$PLANT_CODE == "POCU6" & abun_all$ABUNDANCE == 1,], add= T, pch=20, col="white")


ordsums[ordsums$species.code=="POCU6",]




bio <- stack("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_6.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_12.asc",
             "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_15.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_3.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_4.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_6.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_7.asc",
             #"C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_8.asc")
             "C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/POCU6.asc")

bio$bio_2[bio$POCU6 == 0] <- NA
bio$bio_5[bio$POCU6 == 0] <- NA
bio$bio_8[bio$POCU6 == 0] <- NA
bio$bio_15[bio$POCU6 == 0] <- NA
#bio$nlcd_5[bio$POCU6 == 0] <- NA
#bio$nlcd_8[bio$POCU6 == 0] <- NA
#bio$nlcd_7[bio$ALPE4_FULL_all_vars == 0] <- NA
#bio$nlcd_8[bio$ALPE4_FULL_all_vars == 0] <- NA
bio$bio_2 <- bio$bio_2/10
bio$bio_5 <- bio$bio_5/10
bio$bio_8 <- bio$bio_8/10
#bio$bio_6 <- bio$bio_6/10

bio <- dropLayer(bio,5)




abunpts <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv",
                      sep=",", header=T, stringsAsFactors = F)
head(abunpts)
abunpts <- abunpts[abunpts$species=="POCU6",]
coordinates(abunpts) <- c(5,4)

occ.pts <- read.table("C:/Users/Localadmin/Documents/MaxEnt_modeling/species/POCU6.csv", sep=",", stringsAsFactors = F, header=T)
str(occ.pts)
coordinates(occ.pts) <- c(2,3)




extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")

occ.ran <- occ.pts[sample(1:nrow(occ.pts), 20)]

plot(extentShape, col="grey92")
plot(occ.pts[sample(1:nrow(occ.pts), 150),],lwd=3,col="grey20",add=T)
plot(abunpts[sample(1:nrow(abunpts), 50),], pch=19, col="firebrick2", add=T)


maxentRange <- readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/POCU6.asc")
plot(maxentRange)
maxentRange <- maxentRange[maxentRange$POCU6.asc == 1, ]
plot(maxentRange)



maxentRange <- data.frame(maxentRange)
#maxentRangeD$suitability[maxentRangeD$ALJU_FULL_all_vars.asc == 1] <- 1
#maxentRangeD$suitability[maxentRangeD$ALJU_FULL_all_vars.asc == 0] <- 0
#maxentRangeD$suitability <- as.factor(maxentRangeD$suitability)
## data.frame version of raster
rast <- as.data.frame(raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc"))
#rast$suitable <- -99
rast$suitable[!is.na(rast$bio_2)] <- maxentRange$POCU6_FULL_all_vars.asc
#rast$suitable <- as.factor(rast$suitable)
maxentRange <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=rast$suitable)
maxentRange <- as.data.frame(maxentRange)
maxentRange$suitable[maxentRange$layer == 0] <- "unsuitable"
maxentRange$suitable[maxentRange$layer == 1] <- "suitable"
maxentRange$suitable <- as.factor(maxentRange$suitable)
maxentRange <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=maxentRange$suitable)
#plot(maxentRange)
spplot(maxentRange)


maxentLogit <- readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/FULL_all_vars/FINAL_OUTPUT_FULL_all_vars/POCU6/POCU6_avg.asc") 
ordinal <- readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc") 
#ordinal <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=ordinal$ALJU.asC)
ordinal <- as.data.frame(ordinal)
vals <- ordinal$POCU6.asc
rast <- as.data.frame(raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_2.asc"))
#rast$suitable <- -99
#rast$abundance <- -99
rast$abundance[!is.na(rast$bio_2)] <- vals
#rast$suitable <- as.factor(rast$suitable)
ordinal <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=rast$abundance)
ordinal <- as.data.frame(ordinal)
#ordinal$abundance[ordinal$layer == 3] <- "high"
#ordinal$abundance[ordinal$layer == 2] <- "medium"
#ordinal$abundance[ordinal$layer == 1] <- "low"
#ordinal$abundance[ordinal$layer == 0] <- "unsuitable"
#ordinal$abundance <- as.factor(ordinal$abundance)
ordinal <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=ordinal$layer)
plot(ordinal)
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALPE4.asc"))
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/EUAL13.asc")) #___
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/LOJA.asc"))   #___
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/PUMOL.asc"))
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ROMU.asc"))   #__
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/RHCA3.asc"))
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ALPE4.asc"))  #*****
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/ACPL.asc"))   #__
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BETH.asc"))   #__
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEHE.asc")) 
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/MIVI.asc")) 
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/POCU6.asc")) #****
#plot(readAsciiGrid("file:///C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/HEMA3.asc"))

messpts <- data.frame(cbind(abunpts$bio_2, abunpts$bio_5, abunpts$bio_6, abunpts$bio_8, abunpts$nlcd_5, abunpts$nlcd_8))
head(messpts)
colnames(messpts)<-c("bio_2","bio_5", "bio_6", "bio_8", "nlcd_5", "nlcd_8")

library(dismo)
messplot <- mess(bio,messpts,full=F)
plot(messplot)
#messvals <- as.data.frame(messplot$rmess)
#str(messvals)
#messvals <- messvals[messvals$rmess != Inf,]
#hist(messvals)
messplot <- stack(messplot, bio$bio_2, "C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc")

messD <- as.data.frame(messplot)
head(messD)
#messD$rmess[messD$rmess == Inf] <- NA #rep(NA, length(messD$rmess[messD$rmess == Inf]))
messD$mess2[!is.na(messD$bio_5) & is.na(messD$bio_2)] <- "unsuitable (no prediction)"
messD$mess2[messD$mess > 10 & messD$mess < Inf] <- "interpolation habitat"
messD$mess2[messD$mess <= 10 & messD$mess >= -10] <- "marginal habitat"
messD$mess2[messD$mess < -10] <- "extrapolation habitat"
messD$mess2<- as.factor(messD$mess2)
messplot2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=messD$mess2)
#messD$mess2[messD$rmess == Inf & is.na(messD$layer.1)] <- "absent"
#plot(messplot2)
#plot(ordinal)

################## fix vars
library(ordinal)
library(reshape2)
abunpts$abundance <- ordered(as.factor(abunpts$abundance), levels=c(1,2,3))
M <- clm(abundance~poly(bio_2,2) + bio_5 + poly(bio_8,2) + poly(bio_12,2) + poly(nlcd_5,2)+ 
      poly(nlcd_6,2) + poly(nlcd_7,2) + poly(nlcd_8,2), data=abunpts, Hess=T)
p <- predict(M,as.data.frame(bio),type="prob")
pd <- data.frame(p)
head(pd)
pd


#p <- melt(p)
prob3 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=pd$fit.3)

#extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
#spplot(extentShape, col="grey")


#spplot(prob3, add=T)

#pts=list("sp.points", occ.pts, pch=19, col="dark grey")
#spplot(maxentLogit, col.regions=rev(terrain.colors(100))) +
# spplot(occ.pts, "PLANT_CODE", col.regions=c("dark grey"), colorkey=TRUE, sp.layout=list(pts)) #pts=list("sp.points", occ.pts, pch=19, col="darkorange3")

#colors=c("light green","black")
pts=list("sp.points", occ.pts, pch=19, col="darkorange4")
pts2=list("sp.points", abunpts, pch=19, col="darkorange1")
spplot(maxentRange, col.regions=c("light green","dark grey")) + 
  spplot(occ.pts, "PLANT_CODE", col.regions=c("darkorange4"), colorkey=TRUE, sp.layout=list(pts)) + 
     spplot(abunpts, "species", col.regions=c("darkorange1"), colorkey=TRUE, sp.layout=list(pts2))  #+
#spplot(abunpts, add=T, col="maroon")

#colors <-
spplot(messplot2, col.regions=c( "orange", "light green", "light yellow", "black") ) #+
  #spplot(abunpts, "species", col.regions=c("darkorange3"), colorkey=TRUE, sp.layout=list(pts) ) 


#spplot(prob3, col.regions=rev(terrain.colors(100)))

#plot(extentShape, col= "dark grey")
#prob3[prob3$layer==0] <- 
#plot(prob3, add=T)



#abunpts2 <- abunpr
#pts=list("sp.points", abunpts, pch=19, col=c("light grey", "dark grey", "black"))
spplot(ordinal$layer, col.regions=c( "dark orange", "light green", "gold", "dark grey"), main="Predicted abundance") #+ 
  #spplot(abunpts, "abundance", col.regions=c("light grey", "dark grey", "black"), colorkey=TRUE, sp.layout=list(pts) ) 




pts=list("sp.points", occ.pts, pch=20, col="saddlebrown")
pts2=list("sp.points", abunpts, pch=".", col="orange1")
#pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/rangeALJU.pdf",width=11,height=8.5)
spplot(maxentRange, col.regions=c("forestgreen","dark grey"), main= "Binary Range") + 
  spplot(occ.pts, "PLANT_CODE", pch=20, col.regions=c("saddlebrown"), colorkey=TRUE, sp.layout=list(pts)) + 
  spplot(abunpts, "species", pch=20, col.regions=c("orange1"), colorkey=TRUE, sp.layout=list(pts2)) 
dev.off()

#pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/messALJU.pdf",width=11,height=8.5)
spplot(ordinal$layer, col.regions=c( "dark orange", "light green", "gold", "dark grey"), main="Predicted abundance") #+ 
dev.off()

#pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/abunALJU.pdf",width=11,height=8.5)
spplot(messplot2, col.regions=c("orange","light green","light yellow","dark grey"), main="MESS analysis")#( "dark grey", "orange", "light green", "light yellow") )
dev.off()


plot(hi_abun_rich)

writeRaster(hi_abun_rich,"C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/abun_richness_map_FULL_all_vars.asc")
writeRaster(hi_abun_hotspot,"C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/abun_hotspot_map_FULL_all_vars.asc")
writeRaster(full_hotspot,"C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/rich_hotspot_map_FULL_all_vars.asc")



############### GRC figs


## blank map with points            CHECK
## layers                           CHECK
## rangemap (use pixelated?)        CHECK


## range map                        CHECK
## range map + all pts              CHECK
## range map + abun pts             CHECK

## mess + abun pts                  CHECK
## ordinal                          CHECK

## richness map
## hotspot map (rich)
## abun richness map
## hotspot map (abun)
## hotspot compare



maxentRange$layer2 <- maxentRange$layer
#maxentRange$layer2[maxentRange$layer2==1] <- 1
maxentRange$layer2[maxentRange$layer2==2] <- NA


png("C:/Users/Localadmin/Desktop/GRC_figs/range.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(maxentRange$layer2, col="darkolivegreen2", add=T)
dev.off()
## EXPORT

png("C:/Users/Localadmin/Desktop/GRC_figs/range_occpts.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(maxentRange$layer2, col="darkolivegreen2", add=T)
plot(occ.pts, add=T, pch=20, col="navyblue")
dev.off()
## EXPORT

png("C:/Users/Localadmin/Desktop/GRC_figs/range_abunpts.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(maxentRange$layer2, col="darkolivegreen2", add=T)
plot(abunpts, add=T, pch=20, col="blue")
dev.off()
## EXPORT

png("C:/Users/Localadmin/Desktop/GRC_figs/range_abun_occpts.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(maxentRange$layer2, col="darkolivegreen2", add=T)
plot(occ.pts, add=T, pch=20, col="navyblue")
plot(abunpts, add=T, pch=20, col="blue")
dev.off()
## EXPORT



png("C:/Users/Localadmin/Desktop/GRC_figs/range_abun_occpts2.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(maxentRange$layer2, col="darkolivegreen2", add=T)
plot(occ.pts, add=T, pch=20, col="navyblue")
plot(abunpts, add=T, pch=20, col="firebrick3")
dev.off()
## EXPORT


ordinal$class1 <- 0
ordinal$class2 <- 0
ordinal$class3 <- 0
ordinal$class1[ordinal$layer == 1] <- 1
ordinal$class2[ordinal$layer == 2] <- 1
ordinal$class3[ordinal$layer == 3] <- 1
plot(ordinal$class3)
ordinal$class1[ordinal$class1==0]<- NA
ordinal$class2[ordinal$class2==0]<- NA
ordinal$class3[ordinal$class3==0]<- NA

png("C:/Users/Localadmin/Desktop/GRC_figs/abun.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(ordinal$class1, add=T, col="darkolivegreen3")
plot(ordinal$class2, add=T, col="yellow")
#plot(ordinal$class3, add=T, col="firebrick3")
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/abun3.png", width=1200, height=675, bg="transparent")
plot(ordinal$class3, col="firebrick3")
dev.off()

#rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
#heat.colors(n, alpha = 1)
#terrain.colors(n, alpha = 1)
#topo.colors(n, alpha = 1)
#cm.colors(n, alpha = 1)

messplot2$interp <- -99
messplot2$marg <- -99
messplot2$extrap <- -99

messplot2$interp[messplot2$layer == 2] <- 1
messplot2$marg[messplot2$layer == 3]  <- 1
messplot2$extrap[messplot2$layer == 1]  <- 1

messplot2$interp[messplot2$interp == -99]<- NA
messplot2$marg[messplot2$marg== -99]  <- NA
messplot2$extrap[messplot2$extrap == -99]<- NA

plot(messplot2$extrap)
#messplot2$interp[messplot2$layer == 4 | is.na(messplot2$layer)]<- NA
#messplot2$marg[messplot2$layer == 4| is.na(messplot2$layer)]  <- NA
#messplot2$extrap[messplot2$layer == 4| is.na(messplot2$layer)]<- NA

png("C:/Users/Localadmin/Desktop/GRC_figs/mess.png", width=1200, height=675, bg="transparent")
plot(maxentRange$layer, col="grey60")
plot(messplot2$interp, add=T, col= "steelblue1")
plot(messplot2$marg,   add=T, col= "steelblue1")
#plot(messplot2$extrap, add=T, col= "mediumvioletred")
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/mess2.png", width=1200, height=675, bg="transparent")
plot(messplot2$extrap, col= "violetred1")
dev.off()

#bio5 <- raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc")
#bio8 <- raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc")
#nlcd5 <- raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc")

#plot(bio5)
#plot(bio8)
#plot(nlcd5)

png("C:/Users/Localadmin/Desktop/GRC_figs/env1.png", width=1200, height=675, bg="transparent")
plot(raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_5.asc"))
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/env2.png", width=1200, height=675, bg="transparent")
plot(raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/bio_8.asc"))
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/env3.png", width=1200, height=675, bg="transparent")
plot(raster("C:/Users/Localadmin/Documents/Environmental_Data_2_8_2018/nlcd_5.asc"))
dev.off()

shp <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer="us_shape")

png("C:/Users/Localadmin/Desktop/GRC_figs/pts1.png", width=1200, height=675, bg="transparent")
plot(shp, col="light grey")
plot(occ.pts, col="navyblue", add=T, pch=20)
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/abunrich.png", width=1200, height=675, bg="transparent")
plot(current_rich, col=rev(heat.colors(100)))
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/speciesrich.png", width=1200, height=675, bg="transparent")
plot(current_rich, col=rev(heat.colors(100)))
dev.off()

full_hotspot$hotspot <- -1
full_hotspot$hotspot[full_hotspot$layer == 1] <- 1
full_hotspot$hotspot[full_hotspot$hotspot == -1] <- NA

hi_abun_hotspot$hotspot <- -1
hi_abun_hotspot$hotspot[hi_abun_hotspot$layer == 1] <- 1
hi_abun_hotspot$hotspot[hi_abun_hotspot$hotspot == -1] <- NA


png("C:/Users/Localadmin/Desktop/GRC_figs/abunhot.png", width=1200, height=675, bg="transparent")
plot(shp, col="grey80")
plot(hi_abun_hotspot$hotspot, col="firebrick3", legend=F, add=T)
dev.off()

png("C:/Users/Localadmin/Desktop/GRC_figs/richhot.png", width=1200, height=675, bg="transparent")
plot(shp, col="grey40")
plot(full_hotspot$hotspot, col="yellow", legend=F, add=T)
dev.off()

comp <- as.data.frame(stack(hi_abun_hotspot,full_hotspot))
head(comp)
#comp$x[comp$layer.1==0 & comp$layer.2==0 ] <- 1 #neither
comp$abun[comp$layer.1==1 & comp$layer.2==0 ] <- 2 #just abun
comp$both[comp$layer.1==1 & comp$layer.2==1 ] <- 3 #both
comp$full[comp$layer.1==0 & comp$layer.2==1 ] <- 4 #just full

abun <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=comp$abun)
both <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=comp$both)
full <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=comp$full)

png("C:/Users/Localadmin/Desktop/GRC_figs/compare.png", width=1200, height=675, bg="transparent")
plot(shp, col="grey40")
plot(abun, col="firebrick2", add=T, legend=F)
plot(both, col="orange", add=T, legend=F)
plot(full, col="yellow", add=T, legend=F)
dev.off()

















############################
### leave-one-out analysis

library(maptools)
library(raster)
library(rgdal)
library(colorRamps)
library(ordinal)
library(reshape2)

#ordsums <- ordsums[ordsums$rule5broke == 1,]
#summary(ordsums$kappa)
#ordsums <- ordsums[ordsums$AUC > 0.7 & ordsums$kappa >0 & ordsums$kappaP <0.05 & ordsums$MESSgn10 >.9,]

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_8_2018b.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
ordsums <- ordsums[order(-ordsums$USE, ordsums$species.code),]



spp <- ordsums$species.code[ordsums$AUC > 0.7 & ordsums$kappa >0 & ordsums$kappaP <0.05 & ordsums$MESS_abun2range_SEL >.9]


for ( i in 1:length(spp)) {

  if(i == 1){
    
    all <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", spp[2:length(spp)],".asc")  
    abun <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", spp[2:length(spp)],".asc")  
  }
  
  if(i == length(spp)){
    all <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", spp[1:length(spp)-1],".asc") 
    abun <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", spp[1:length(spp)-1],".asc")  
  }
  
  if(i != 1 & i != length(spp)){
    all <- c(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", spp[1:(i-1)],".asc"),
             paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", spp[(i+1):length(spp)],".asc"))
    
    abun <- c(paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", spp[1:(i-1)],".asc"),
               paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", spp[(i+1):length(spp)],".asc"))
  }
 
  all  = stack(all)
  all  = calc(all, sum)
  
  
  abun = stack(abun)
  abun  = calc(abun, sum)
  
  
  all <- as.data.frame(all)
  all$layer[all$layer <= quantile(all$layer, .75, na.rm=T) & !is.na(all$layer)] <- 0 
  all$layer[all$layer > quantile(all$layer, .75, na.rm=T) & !is.na(all$layer)] <- 1
  
  
  abun <- as.data.frame(abun)
  abun$layer2[abun$layer <= quantile(abun$layer, .75, na.rm=T) & !is.na(abun$layer)] <- 0
  abun$layer2[abun$layer > quantile(abun$layer, .75, na.rm=T) & !is.na(abun$layer)] <- 1

  comp <- cbind(all,abun$layer2)
  ordsums$l1o[i] <- length(comp$layer[comp$'abun$layer2'==1 & comp$layer==1 & !is.na(comp$layer)])/length(comp$layer[!is.na(comp$layer) & comp$layer==1])  
  ordsums$l1o2[i] <- length(comp$layer[comp$'abun$layer2'==1 & comp$layer==1 & !is.na(comp$layer)])/length(comp$layer[!is.na(comp$layer) & comp$'abun$layer2'==1])  
  ordsums$l1o3[i] <- length(comp$layer[comp$'abun$layer2'==1 & comp$layer==1 & !is.na(comp$layer)])
  
  
  
  print(i)
}


write.table(ordsums,"C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_31_2018.csv",row.names=F,sep=",")
hist(ordsums$l1o)
hist(ordsums$l1o2)
hist(ordsums$l1o3)

cor(ordsums$l1o,ordsums$l1o3)
cor(ordsums$l1o2,ordsums$l1o3)
plot(ordsums$l1o2,ordsums$l1o3)

mean(ordsums$l1o2[ordsums$USE==1],na.rm=T)
sd(ordsums$l1o2[ordsums$USE==1],na.rm=T)

hist(ordsums$l1o2*100, col="lightblue", breaks=20, xlim=c(35,60), ylim=c(0,25), ylab="Number of Species", xlab="Percent overlap (%)", main="Leave-one-out hotspot comparisons")
abline(v=0.50, lwd=4)
ordsums$species.code[ordsums$l1o2 < .40]
ordsums$species.code[ordsums$l1o2 > .55]
#ordsums$species.code[order(ordsums$l1o2)][1:5]
#ordsums$species.code[rev(order(ordsums$l1o2))][1:5]
ordsums$l1o2[ordsums$species.code == "COMA2"]
ordsums$l1o2[ordsums$species.code == "HYPE"]

ordsums$kappa[order(ordsums$l1o2)][1:5]
ordsums$kappa[rev(order(ordsums$l1o2))][1:5]
head(ordsums)
ordsums[ordsums$species.code=="FRAL4",]





#################

#### species trait composition

library(raster)
library(RColorBrewer)
display.brewer.all(n=8, type="all", colorblindFriendly=T)



ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_8_2018b.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
ordsums <- ordsums[ordsums$USE==1,]

unique(ordsums$habit)
unique(ordsums$duration_final)

woody <- ordsums$species.code[ordsums$habit == "tree" | ordsums$habit == "shrub_subshrub"]
p.herb <- ordsums$species.code[(ordsums$habit == "forb_herb" | ordsums$habit == "grass") & ordsums$duration_final == "Perennial" ]
np.herb <- ordsums$species.code[(ordsums$habit == "forb_herb" | ordsums$habit == "grass") & ordsums$duration_final != "Perennial" ]

#per  <- ordsums$species.code[ordsums$duration_final == "Perennial"]
#np   <- ordsums$species.code[ordsums$duration_final != "Perennial"]
#summary(as.factor(ordsums$habit[ordsums$duration_final != "Perennial"]))
#& ordsums$duration_final == "Perennial",]
#non.per <- ordsums[ordsums$duration_final != "Perennial",]

rich <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc")
abun <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
riHOT <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_FULL.asc")
abHOT <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_HI_ABUN.asc")
coHOT <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_compare.asc")


richD <- as.data.frame(rich)
abunD <- as.data.frame(abun)


#filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", woody,".asc")
filenames <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", woody,".asc")

ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)
overlapD <- as.data.frame(overlap)

vals <- overlapD$layer/richD$richness_map_FULL

woody_rich_prop <- raster(nrows=nrow(rich),ncols=ncol(rich),ext=extent(rich),vals=vals)
#plot(woody_rich_prop, col=brewer.pal("YlOrRd",n=9))
#plot(riHOT, add=T, col=c("grey","transparent") )

plot(woody_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of establishers that are woody")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )




######################################################
filenames <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", woody,".asc")

ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)
overlapD <- as.data.frame(overlap)

vals <- overlapD$layer/abunD$richness_map_HI_ABUN

woody_abun_prop <- raster(nrows=nrow(rich),ncols=ncol(rich),ext=extent(rich),vals=vals)
#plot(woody_abun_prop, col=brewer.pal("YlOrRd",n=9))
#plot(abHOT, add=T, col=c("grey","transparent") )

plot(woody_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of impacters that are woody")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F)





##################################################################
filenames <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", per.herb,".asc")

ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)
overlapD <- as.data.frame(overlap)

vals <- overlapD$layer/richD$richness_map_FULL

perHerb_rich_prop <- raster(nrows=nrow(rich),ncols=ncol(rich),ext=extent(rich),vals=vals)
#plot(woody_rich_prop, col=brewer.pal("YlOrRd",n=9))
#plot(riHOT, add=T, col=c("grey","transparent") )

plot(perHerb_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of establishers that are perennial herbs")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )




######################################################
filenames <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", per.herb,".asc")

ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)
overlapD <- as.data.frame(overlap)

vals <- overlapD$layer/abunD$richness_map_HI_ABUN

perHerb_abun_prop <- raster(nrows=nrow(rich),ncols=ncol(rich),ext=extent(rich),vals=vals)
#plot(woody_abun_prop, col=brewer.pal("YlOrRd",n=9))
#plot(abHOT, add=T, col=c("grey","transparent") )

plot(perHerb_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of impacters that are perennial herbs")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F)




####################################################################################
filenames <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", np.herb,".asc")

ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)
overlapD <- as.data.frame(overlap)

vals <- overlapD$layer/richD$richness_map_FULL

npHerb_rich_prop <- raster(nrows=nrow(rich),ncols=ncol(rich),ext=extent(rich),vals=vals)
#plot(woody_rich_prop, col=brewer.pal("YlOrRd",n=9))
#plot(riHOT, add=T, col=c("grey","transparent") )

plot(npHerb_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of establishers that are nonperennial herbs")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )




######################################################
filenames <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", np.herb,".asc")

ldf=stack(filenames)

overlap=calc(ldf, sum)
plot(overlap)
overlapD <- as.data.frame(overlap)

vals <- overlapD$layer/abunD$richness_map_HI_ABUN

npHerb_abun_prop <- raster(nrows=nrow(rich),ncols=ncol(rich),ext=extent(rich),vals=vals)
#plot(woody_abun_prop, col=brewer.pal("YlOrRd",n=9))
#plot(abHOT, add=T, col=c("grey","transparent") )

plot(npHerb_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of impacters that are nonperennial herbs")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F)








eco <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco", stringsAsFactors = F)
plot(eco)

#eco <- spTransform(eco, "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
#plot(eco)

proj4string(woody_rich_prop) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(woody_abun_prop) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(perHerb_rich_prop) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(perHerb_abun_prop) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(npHerb_rich_prop) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(npHerb_abun_prop) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

woody_rich_prop <- projectRaster(woody_rich_prop, crs="+init=epsg:5070")
woody_abun_prop <- projectRaster(woody_abun_prop, crs="+init=epsg:5070")
perHerb_rich_prop <- projectRaster(perHerb_rich_prop, crs="+init=epsg:5070")
perHerb_abun_prop <- projectRaster(perHerb_abun_prop, crs="+init=epsg:5070")
npHerb_rich_prop <- projectRaster(npHerb_rich_prop, crs="+init=epsg:5070")
npHerb_abun_prop <- projectRaster(npHerb_abun_prop, crs="+init=epsg:5070")




############### ecoregion plots

pdf(file="C:/Users/Localadmin/Google Drive/figs/traits_eco.pdf",width=15, height=7.5)
par(mfrow=c(2,3))

plot(woody_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of woody plants", ylab="Establishers", mar=c(0,0,0,0) )
plot(eco, add=T, col=c("transparent"), axes=F, lwd=0.5)

plot(perHerb_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of perennial herbs")
plot(eco, add=T, col=c("transparent"), axes=F, lwd=0.5 )

plot(npHerb_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of nonperennial herbs")
plot(eco, add=T, col=c("transparent"), axes=F, lwd=0.5 )

plot(woody_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, ylab= "Impacters")
plot(eco, add=T, col=c("transparent"), axes=F, lwd=0.5 )


plot(perHerb_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F)
plot(eco, add=T, col=c("transparent"), axes=F, lwd=0.5 )


plot(npHerb_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F)
plot(eco, add=T, col=c("transparent"), axes=F, lwd=0.5 )

dev.off()

proj4string(coHOT) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
coHOT <- projectRaster(coHOT, crs="+init=epsg:5070")


############### hotspot plots
pdf(file="C:/Users/Localadmin/Google Drive/figs/traits_hots.pdf",width=15, height=7.5)
par(mfrow=c(2,3))

plot(woody_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of woody plants", ylab="Establishers", mar=c(0,0,0,0) )
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )

plot(perHerb_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of perennial herbs")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )

plot(npHerb_rich_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, main= "Proportion of nonperennial herbs")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )

plot(woody_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F, ylab= "Impacters")
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )

plot(perHerb_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F)
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )

plot(npHerb_abun_prop, col=brewer.pal("YlOrRd",n=9), axes=F, box=F)
plot(coHOT, add=T, col=c("transparent", "transparent", "grey", "transparent"), axes=F, box=F, legend=F )


dev.off()

#BLOOP