#### Modeling abundance for terrestrial invasive plants in CONUS
###  Mitch O'Neill, began work on this code in 2017
##   last modified: 12/12/2020 - added more comments to clarify code


#### files needed
## 1. environmental layers - Environmental_Data_02_08_2018 folder
## 3. bias file 
## 4. Table with each species 
## 5. occurence + abundance dataset
## 6. MaxEnt must be installed on your machine

## define filepath of MaxEnt on computer
maxent.location='C:/Users/Localadmin/Documents/MaxEnt/maxent.jar'


## define directory of environmental data for fitting
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs",sep="")

## species list from "ordsums"
sp.list <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums.csv", header = T, sep = ",", 
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

spp_bias_log = paste0("java -jar ",maxent.location, 
                      " nowarnings noprefixes -E responsecurves jackknife outputformat=logistic removeduplicates noaskoverwrite replicates=10 noproduct nothreshold nohinge writeplotdata noautofeature biastype=3")
## settings for bias surface model: linear + quadratic terms, thinning, logistic, 10 replicates, bias surface

for(i in 1:length(sp.list$code)){

  ## define point location samples (all species in one .csv), must be set up with three columns (species, lon, lat, in that order)
  samples=paste("samplesfile=", paste("C:/Users/Localadmin/Documents/Maxent_modeling/species", 
                                      paste(sp.list$code[i], "csv", sep="."), sep="/"), sep="")
 
  ## create a directory for MaxEnt output (a folder with the name of the species) and set it as the output directory for MaxEnt
  directory<- paste("C:/Users/Localadmin/Documents/MaxEnt_modeling/FINAL_OUTPUT", sp.list$code[i], sep="/")
  dir.create(directory)
  output=paste("outputdirectory=", directory, sep="") 
  
  ## define the bias file
  bias=paste("biasfile=","C:/Users/Localadmin/Documents/MaxEnt_modeling/BIAS_OUTPUT_FULL/BIAS_avgFULL.asc",sep="")
  
  ## call model to run, including the settings, output directory, environmental data, invasive species data, and bias
  system(paste(spp_bias_log,output,environmental,samples,bias, "autorun"))
  
  
  print(i) ## show progress within loop
  
}

## clear up memory
rm(sp.list,bias,directory,environmental,i,maxent.location,output,samples,spp_bias_log)
gc()

####################################################################################
## Find 95% threshold for each species AND get training/test sample sizes
###################################################################################
## maxent output represents relative occurence rate (ROR)
## we want a binary map, so we need to know what ROR value to use to define "suitable" vs "unsuitable"
## all maxent surfaces must sum to 1, so the threshold value can vary from species to species
## therefore we need to calculate the threshold by species (below)

## need a species list again
spp=sp.list ## maybe find and replace and have just one sp.list later


#create dataframe to hold 95% MTP logistic threshold values for each run
spp_val=data.frame("Species",15,0,"CtySppRich_Bias",100000, 100000,stringsAsFactors=F) 
colnames(spp_val)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n","train_n")



## Set directory to maxent output
setwd("C:/Users/Localadmin/Documents/MaxEnt_modeling/FINAL_OUTPUT/")


for (i in spp){        #loop through each species
  for (j in 0:9){      #loop for each of the iterations (ten replicate model runs) for each species
    
    #read file for iteration j point-level model output 
    file=paste(i,"/",i,"_",j,"_samplePredictions.csv",sep="")     
    d=read.table(file, sep=",",header=T)
   
     d$run=j ## what is this?
     
    dt=d[d$Test.or.train == "train",]  #subset iteration data to training data only
    dtest=d[d$Test.or.train == "test",]  #subset iteration data to testing data only
    
   
    
    #get training sample size for this iteration 
    samp=dim(dt)
    train_n=samp[1]
    
    #get test sample size for this iteration
    samp1=dim(dtest)
    test_n=samp1[1]  
    
    tp=round(0.95*samp[1], 0) #number of points for 95% threshold in each iteration
    
    dt=dt[order(dt$Logistic.prediction, decreasing=T),] #order the logistic output largest to smallest
    
    #create dataframe with species, run, logistic value of the point that bounds the 95% percentile, bias model, test and train n
    sp_runs=data.frame(i, j, dt$Logistic.prediction[tp], "CtySppRich_Bias", test_n, train_n) 
    colnames(sp_runs)=c("Species", "Run", "X95_MTP_thresh_logistic","Bias", "test_n", "train_n")
    
    spp_val=rbind(spp_val,sp_runs) 

    
  }
  ### get AUC scores
  aucs=read.table(paste(i,"/maxentResults.csv",sep=""), header = T, sep = ",",  
                  quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
  
  ### record on summaries file
  ordsums$AUC[ordsums$species.code == i] <- aucs$Test.AUC[11]

  
  print(i) #track progress
}


## trim junk row from spp_val
spp_val1=spp_val[spp_val$Species!="Species",]


## calculate mean threshold per species per range from spp_val1
thresh=as.data.frame(tapply(spp_val1$X95_MTP_thresh_logistic, list(spp_val1$Species, spp_val1$Bias), mean))
thresh$Species=rownames(thresh) #retrieve species names from rownames


## make table of species with their number of test and training points ### do I need this???
test=as.data.frame(tapply(spp_val1$test_n, list(spp_val1$Species, spp_val1$Bias), mean))
test$Species=rownames(test)
colnames(test)=c("test_n","Species")

train=as.data.frame(tapply(spp_val1$train_n, list(spp_val1$Species, spp_val1$Bias), mean))
train$Species=rownames(train)
colnames(train)=c("train_n","Species")

samples=cbind(test, train)

samples1=subset(samples, select=c("Species","test_n","train_n"))
samples1$total_n=samples1$test_n + samples1$train_n

## save table of thresholds and table of sample sizes
write.table(thresh, file= "C:/Users/Localadmin/Documents/Maxent_modeling/summaries/thresholds.csv" ,sep=",",row.names=F)
write.table(samples1,file= "C:/Users/Localadmin/Documents/Maxent_modeling/summaries/sample_sizes.csv" ,sep=",",row.names=F) 

#########################################################################
### Create binary maps for each species
#########################################################################

###### create binary rasters based on 95% MTP

for (i in unique(thresh$Species)){  
  
  ## read Ascii file of average logisitic from MaxEnt for species i
  file=paste(i,"/",i,"_avg.asc",sep="")
  d1=readAsciiGrid(file) 
  
  ## where logistic is greater than or equal to 95% threshold, set to 1 (present)  
  d1@data[d1@data >= thresh$CtySppRich_Bias[thresh$Species==i]]=1 
  
  ## where logistic is below the 95% threshold, set to 1 (not present) 
  d1@data[d1@data < thresh$CtySppRich_Bias[thresh$Species==i]]=0
  
  ## write out new binary ascii file, named by species
  out_file=paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",i,".asc",sep="")
  write.asciigrid(d1, out_file)
  
  rm(d1) #garbage cleaning
  gc()
  print(i) #track progress
}

### ^^^^^ based on MaxEnt_and_Hotspots.R lines 1-303, after ordinalSDM, pick up at line 561 of maxent and hotspots


####################################################################################
##### Constructing ordinal models abundance as a response to environmental variables


######################################
###### create a list of formulae with

## list candidate terms (all 12 climate and landcover variables + quadratics)
var.list <- c("bio_2", "poly(bio_2, 2)",
              "bio_5", "poly(bio_5, 2)",
              "bio_6", "poly(bio_6, 2)",
              "bio_8", "poly(bio_8, 2)",
              "bio_12", "poly(bio_12, 2)",
              "bio_15", "poly(bio_15, 2)",
              "nlcd_3", "poly(nlcd_3, 2)",
              "nlcd_4", "poly(nlcd_4, 2)",
              "nlcd_5", "poly(nlcd_5, 2)",
              "nlcd_6", "poly(nlcd_6, 2)",
              "nlcd_7", "poly(nlcd_7, 2)",
              "nlcd_8", "poly(nlcd_8, 2)") 

var.sets <- c("variable sets") ## initialize list of variable sets

## loop through all possible numbers of variables
## 12 instead of 24 because poly(variable, 2) includes the linear and quadratic, so including 
## variable + poly(variable, 2) would duplicate the linear term

for(i in 1:12){
  combos.i <-(combn(var.list, i, simplify=F)) ## all combinations for that # of vars
  
  for (j in 1:length(combos.i)){ ## go through each combination with that # of vars
    combo.j <- as.character(collapse(combos.i[[j]], sep= " + ")) 
    ## collapse the combination of vars into a single char string (easier to deal with)
    
    ## As long as the combination does not duplicate the 1st order term; (i.e. bio_5 + poly(bio_5, 2))
    if(!(grepl("bio_2 ", combo.j) & grepl("poly(bio_2", combo.j, fixed=T)) &
       !(grepl("bio_5 ", combo.j) & grepl("poly(bio_5", combo.j, fixed=T)) &
       !(grepl("bio_6 ", combo.j) & grepl("poly(bio_6", combo.j, fixed=T)) &
       !(grepl("bio_8 ", combo.j) & grepl("poly(bio_8", combo.j, fixed=T)) &
       !(grepl("bio_12 ", combo.j) & grepl("poly(bio_12", combo.j, fixed=T)) &
       !(grepl("bio_15 ", combo.j) & grepl("poly(bio_15", combo.j, fixed=T)) &
       !(grepl("nlcd_3 ", combo.j) & grepl("poly(nlcd_3", combo.j, fixed=T)) &
       !(grepl("nlcd_4 ", combo.j) & grepl("poly(nlcd_4", combo.j, fixed=T)) &
       !(grepl("nlcd_5 ", combo.j) & grepl("poly(nlcd_5", combo.j, fixed=T)) &
       !(grepl("nlcd_6 ", combo.j) & grepl("poly(nlcd_6", combo.j, fixed=T)) &
       !(grepl("nlcd_7 ", combo.j) & grepl("poly(nlcd_7", combo.j, fixed=T)) &
       !(grepl("nlcd_8 ", combo.j) & grepl("poly(nlcd_8", combo.j, fixed=T))){
    var.sets <- c(var.sets,combo.j)} ## add it to the list of variable sets
    #print(j)
  }
  
  rm(combos.i)
  gc() ## garbage cleaning
  print(i)
}

formulae <- paste("abundance", var.sets[2:length(var.sets)], sep = " ~ ")
## add response part of formulaes, excluding the initial row

rm(var.sets, combo.j, i, j, var.list) ## garbage cleaning

colnames(formulae) <- c("formula")

## make a column for each variable to record which variables are present for each formula
## this will make sorting easier later on ***
formulae$bio_2 <- sapply(formulae$formula, grepl, pattern="bio_2", fixed=T)
formulae$bio_5 <- sapply(formulae$formula, grepl, pattern="bio_5", fixed=T)
formulae$bio_6 <- sapply(formulae$formula, grepl, pattern="bio_6", fixed=T)
formulae$bio_8 <- sapply(formulae$formula, grepl, pattern="bio_8", fixed=T)
formulae$bio_12 <- sapply(formulae$formula, grepl, pattern="bio_12", fixed=T)
formulae$bio_15 <- sapply(formulae$formula, grepl, pattern="bio_15", fixed=T)
formulae$nlcd_3 <- sapply(formulae$formula, grepl, pattern="nlcd_3", fixed=T)
formulae$nlcd_4 <- sapply(formulae$formula, grepl, pattern="nlcd_4", fixed=T)
formulae$nlcd_5 <- sapply(formulae$formula, grepl, pattern="nlcd_5", fixed=T)
formulae$nlcd_6 <- sapply(formulae$formula, grepl, pattern="nlcd_6", fixed=T)
formulae$nlcd_7 <- sapply(formulae$formula, grepl, pattern="nlcd_7", fixed=T)
formulae$nlcd_8 <- sapply(formulae$formula, grepl, pattern="nlcd_8", fixed=T)

## calculate number of terms (no.vars) and number of quadratic terms (no.poly) ***
formulae$no.vars <- sapply(formulae$formula, str_count, pattern="\\+")
formulae$no.poly <- sapply(formulae$formula, str_count, pattern="poly")
formulae$no.vars <- formulae$no.vars + formulae$no.poly + 1

######################

## load abundance dataset with extracted environmental values
edd <- read.table("edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

## load species list / model summay file
ordsums <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/species_list_4_27b.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)

## set up new columns to assign values in loop below
ordsums$pS <- -99           ## p-val for slope
ordsums$pS10 <- -99
ordsums$pI <- -99           ## p-val for intercept
ordsums$pI10 <- -99
ordsums$c.up <- -99         ## number of concave-up response curves
ordsums$c.down <- -99       ## number of concave-down response curves?***
ordsums$aic <- -99          ## AIC value of selected model
ordsums$no.terms  <- -99    ## number of terms in the selected model
ordsums$no.vars <- -99      ## number of environmental variables in selected model
ordsums$formu <- "FILLER"   ## formula of the selected model
ordsums$null <- -99         ## null model kappa?***
ordsums$kappa <- -99        ## Cohen's kappa value of selected model
ordsums$kappaP <- -99       ## p-value of the kappa
ordsums$bin1o <- -99        ## number of observations in bin 1 (low abundance)
ordsums$bin2o <- -99        ## number of observations in bin 2 (medium abundance)
ordsums$bin3o <- -99        ## number of observations in bin 3 (high abundance)
ordsums$bin1p <- -99        ## number of observations predicted into bin 1 (to calculate kappa) ***
ordsums$bin2p <- -99        ## number of observations predicted into bin 2 (to calculate kappa) ***
ordsums$bin3p <- -99        ## number of observations predicted into bin1 (to calculate kappa)  ***


head(ordsums)

for(s in 1:155){ ## Loop through species list

  species <- edd[edd$species == ordsums$species.code[s],]## subset to species of the iteration
  species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3))
  
  ########### ASSESS COLINEARITY OF VARIABLES ############
  rhos <- as.data.frame(cor(species[,7:18], use="pairwise", method="spearman")) 
  
  ## construct data.frame of matrix of pairwise spearman's rho values for environmental variables
  var.pairs <- combn(c(1:12), 2, simplify=T) ## All possible pairs selected from 12 variables
  
  ##initialize master data frameof pairwise correlations with template
  var1 <- "TEMPLATE" ## one variable from pairwise correlation
  var2 <- "TEMPLATE" ## other variable from pairwise correlation
  rho <- -2 ## rho value from correlation
  corr.pairs <- data.frame(var1,var2,rho, stringsAsFactors = F) ## create master data frame
  
  for (i in 1:NCOL(var.pairs)) { ## loop through list of pairs of the six variables
    
    rho <- abs(rhos[var.pairs[1,i], var.pairs[2,i]]) 
    ## absolute value of spearman's rho at the row and column corresponding to the pair in the iteration
    var1 <- var.pairs[1,i]#row.names(rhos[var.pairs[1,i],]) ## name of the first variable in the pair
    var2 <- var.pairs[2,i]#(rhos[var.pairs[2,i],]) ## name of the second variable in the pair
    
    corr.pairs.i <- data.frame(var1,var2,rho, stringsAsFactors = F) ## concatenate data from the interation
    corr.pairs <- rbind(corr.pairs, corr.pairs.i) ## append data from iteraction to the master data frame

  }
  
  corr.pairs <- corr.pairs[corr.pairs$var1 != "TEMPLATE",] ## remove template row
  corr.pairs <- corr.pairs[corr.pairs$rho >= .6, ] ## select all pairwise correlations where collinearity is likely 
  

  ######################################
  ##### CLM (cumulative log odds models)
  ######################################
  
  
  formulae.s <- formulae[formulae$no.vars <= ordsums$max.vars.5epv[s],]
  ## remove formulas with large amounts of variables for the number of data points (see paper)
  
  
  ## initialize a dataframe to record key information from each model (we will be looping for all variable combinations)
  i <- -1 ## keep track of place in formula list
  edf <- -1 ## number of variables
  aic <- -1 ## Akaike's info criterion
  
  model.sel <- data.frame(i, edf, aic, stringsAsFactors = F)    
  ## concatenate into master data.frame template
  
  if(length(corr.pairs$var1)>0) { ## if there are pairs of variables with high correlations, these 
                                  ## variable combinations must be filtered out
    
    for (i in 1:length(formulae.s$formula)){ ##loop through all formulae in the formula list
      
      acceptable <- TRUE ## by default, assume formula is acceptable
      
      for (c in 1:length(corr.pairs$var1)){ ## loop through all variable pairs with high correlation
        
        if(formulae.s[i,(as.numeric(corr.pairs$var1[c])+1)] == T & ## if the formula contains both of
           formulae.s[i,(as.numeric(corr.pairs$var2[c])+1)] == T){ ## the variables, then it is this
          acceptable <- FALSE                                      ## combo variables is not acceptable
        }
      }
      
      if(acceptable) ## if the formula is acceptable (no highly correlated variables), then proceed with model construction
      
        M <- clm(formula(formulae.s$formula[i]), data = species) ## construct model using data from the species of this iteration of the loop
      edf <- length(M$coefficients) ## record effective degrees of freedom
      aic <- AIC(M) ## extract aic score
      model.sel.i <- data.frame(i, edf, aic, stringsAsFactors = F)
      model.sel <- rbind(model.sel, model.sel.i) ## concatenate the data from this model to the table with the data for all models for the species of this iteration
    }
    
    #***** there was an extra bracket - verify that high collinearity didnt happen? *****
  
  } else { ## if there are no variables with high correlations, then all models can be constructed
    
    for (i in 1:length(formulae.s$formula)){ ## for each combination of variables:
      M <- clm(formula(formulae.s$formula[i]), data = species) ## construct model using data from the species of this iteration of the loop
      edf <- length(M$coefficients) ## record effective degrees of freedom
      aic <- AIC(M) ## extract aic score
      model.sel.i <- data.frame(i, edf, aic, stringsAsFactors = F)
      model.sel <- rbind(model.sel, model.sel.i) ## concatenate the data from this model to the table with the data for all models for the species of this iteration
    }
  }

  ###################  AIC comparison  #####################
  model.sel <- model.sel[model.sel$i != -1,] ## remove template row
  model.sel$dAIC <- model.sel$aic - min(model.sel$aic) ## calculate delta AIC
  model.sel2 <- model.sel[model.sel$dAIC <= 2,] ## select the model with the best score, 
  ## as well as models with scores that are not significantly worse
  
  rm(model.sel, edf, i) ## remove unneeded objects
  
  model.sel2 <- model.sel2[model.sel2$edf == min(model.sel2$edf),] 
  ## of the models with the best scores (<= min + 2), select the simplest models (lowest edf, lowest number of terms)
  
  model.sel2 <- model.sel2[model.sel2$dAIC == min(model.sel2$dAIC),]
  ## select the model with the best AIC score if there are multiple models
  ## with similar scores and the same complexity
  
  M <- clm(formula(formulae.s$formula[model.sel2$i]), data = species)
  #M <- polr(formula(formulae.s[3]), data = species, Hess = T)
  ## reconstruct the best model, using the formula number (i) to re-access it
  
  ## extract summary of the model to access model-fit values
  Msum <- coef(summary(M))
  
  ####################################
  ##########************
  MsumS <- data.frame(Msum[3:nrow(Msum),]) ## 1 back to 2
  ## Slopes
  ordsums$pS[s] <- length(MsumS$Estimate[MsumS$Pr...z..>0.05])
  ordsums$pS10[s] <- length(MsumS$Estimate[MsumS$Pr...z..>0.1])
  
  ## Intercepts
  MsumI <- data.frame(Msum[1:2,]) ## 2 back to 3
  ordsums$pI[s] <- length(MsumI$Estimate[MsumI$Pr...z..>0.05])
  ordsums$pI10[s] <- length(MsumI$Estimate[MsumI$Pr...z..>0.1])
  
  MsumS$name <- as.character(row.names(MsumS))
  MsumS$poly <- sapply(MsumS$name, grepl, pattern="2)2")
  
  ######################################
  
  ordsums$c.up[s] <- length(MsumS$Estimate[MsumS$Estimate > 0 & MsumS$poly == T])   ## count concave up quadratics
  ordsums$c.down[s] <- length(MsumS$Estimate[MsumS$Estimate < 0 & MsumS$poly == T]) ## count concave down quadratics
  
  ordsums$aic[s] <- AIC(M)          ## extract adjusted R-square value of model #adjDsq <-  Dsquared(M,adjust=T) ## extract deviance explained by model
  ordsums$no.terms[s] <- NROW(Msum) ## extract the number of terms
  ordsums$no.vars[s] <- ordsums$no.terms[s] - 2 - str_count(formulae.s$formula[model.sel2$i], "poly") ## extract number or environmental variables
  ordsums$formu[s] <- formulae.s$formula[model.sel2$i]  ## extract formula
  
  null <- clm(abundance~1,data=species, Hess=T) ## null model
  
  if(AIC(null) > (AIC(M) + 2)){ordsums$null[s] <- 0} else {ordsums$null[s] <- 1}  
  ## if the AIC of the selected model is significantly lower than the null model, then the "null" = 0
  ## if not, then "null"=1, i.e. a flag for species where the selected model was no better than random

  ## Use selected model to predict the most likely abundance level
  ## based on the environmental conditions at each observation
  pm <- predict(M, species, type="class")
  pm <- melt(pm)
  pm <- pm[,1]

  pm <- data.frame(cbind(species$abundance,pm))
  colnames(pm) <- c("observed", "predicted") 
  ## dataframe with observed and predicted value for each observation
  
  kappa <- kappa2(pm) ## compare observed and precicted abundance using Cohens kappa
  ordsums$kappaP[s] <- kappa$p.value ## extract p-value of kappa
  ordsums$kappa[s] <- kappa$value    ## extract kappa value
  
  ## record the number of samples predicted and observed at each abundance class ***
  ordsums$bin1o[s] <- length(pm$observed[pm$observed == 1])
  ordsums$bin2o[s] <- length(pm$observed[pm$observed == 2])
  ordsums$bin3o[s] <- length(pm$observed[pm$observed == 3])
  ordsums$bin1p[s] <- length(pm$observed[pm$predicted == 1])
  ordsums$bin2p[s] <- length(pm$observed[pm$predicted == 2])
  ordsums$bin3p[s] <- length(pm$observed[pm$predicted == 3])
  
  
  print(s) ## keep track of which species we are on
}


######### Now predict abundance within establishment ranges predicted by MaxEnt
library(maptools) #### *** confirm which packages are needed
library(raster)   #### *** did I add packages earlier for stuff too?
library(rgdal)    #### *** make sure I didn't need
library(colorRamps)
library(ordinal)
library(reshape2)
library(rgeos)


edd <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",",  
                  quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)
## eddmaps abundance data (thinned) with environmental data *** INSERT CLEAN FILE HERE

## list of species
spp <- ordsums$species.code ### *** maybe just use spp list in the next object made

## list of selected model formulae for each species
formulae <- data.frame(cbind(ordsums$species.code,ordsums$formu), stringsAsFactors = F)
colnames(formulae) <- c("species","formula")

edd <- edd[edd$species %in% spp,] ### *** maybe just use formulae$spp
## restrict dataset only to species with ordinal models constructed

## load in rasters for environmental variables *** insert here... maybe I can give it a proj4string
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


proj4string(bio) <-  "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

## environmental variables as a dataframe to speed up computing
bioD <- as.data.frame(bio)
str(bioD)
summary(bioD)
bioD[,1:4] <- bioD[,1:4]/10 ## divide bio_2 and bio_8 so that values are interpretable

## create directory for mapped abundance
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/")
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN") ## maps of abundance within range
dir.create("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN")  ## maps of just cells predicted into high abundance (class 3)

## add fields to dataset to record the number of cells in each class as well as sumary stats *** is this needed? remove date?
ordsums$class1_5_7 <- -99
ordsums$class2_5_7 <- -99
ordsums$class3_5_7 <- -99
ordsums$medianClass_5_7 <- -99
ordsums$mode_5_7 <- -99
ordsums$vr_5_7 <- -99


for (i in 1:length(spp)){ ## loop through species *** spp?
  
  species <- edd[edd$species == spp[i],] ## subset dataset to the species of the iteration
  species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3)) ## convert abundance to ordinal factor
  
  M <- clm(formula(formulae$formula[formulae$species==spp[i]]), data = species) 
  ## construct model using the observations and the selected formula
  
  range <- raster(paste0("Binary_asciis/",spp[i],".asc")) ## binary rangemap from MaxEnt for the species
  rangeD <- as.data.frame(range) ## dataframe to cut down on computing time
  colnames(rangeD) <- c("presence")
  
  p <- predict(M, bioD, type="class") ## predict the abundance using the model and environmental conditions of each cell
  p <- melt(p)
  p <- p[,1]
  
  p <- data.frame(cbind(p, rangeD$presence)) ## concatenate binary presence/absence and abundance predictions
  colnames(p) <- c("abundance","presence")
  p$abundance_presence <- p$abundance*p$presence ## multiply abundance times presence;
  ## absent = 0, so abundance will be 0, present=1, so abundance will be unchanged
  
  p$hi_abundance <- p$abundance_presence/3 ## add field for whether the cell is high abundance (3/3= 1)
  p$hi_abundance[p$hi_abundance < 1] <- 0  ## not high abundance = 0
  
  ## convert datafram back into raster
  abun <- raster(nrows=nrow(range),ncols=ncol(range),ext=extent(range),crs=proj4string(bio),vals=p$abundance_presence)
  hi_abun <- raster(nrows=nrow(range),ncols=ncol(range),ext=extent(range),crs=proj4string(bio),vals=p$hi_abundance)
  
  ## record proportion of cells in each abundance class within the establishment range
  ordsums$class1_5_7[i] <- length(p$abundance[p$abundance_presence == 1 & !is.na(p$abundance_presence)])/length(p$abundance[!is.na(p$presence) & p$presence==1])
  ordsums$class2_5_7[i] <- length(p$abundance[p$abundance_presence == 2 & !is.na(p$abundance_presence)])/length(p$abundance[!is.na(p$presence) & p$presence==1])
  ordsums$class3_5_7[i] <- length(p$abundance[p$abundance_presence == 3 & !is.na(p$abundance_presence)])/length(p$abundance[!is.na(p$presence) & p$presence==1])
  
  ## write out rasters
  writeRaster(abun, paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/",spp[i],".asc"),overwrite=T)
  writeRaster(hi_abun,paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/",spp[i],".asc"),overwrite=T)

  ## calculate median class and other summary stats
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
  
  ######## calculate range size (project into Alber's equal area projection)
  proj4string(range) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  sp <- names(range)
  names(range) <- "layer"
  range <- rasterToPolygons(range, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
  range <- spTransform(range, "+init=epsg:5070") ## Albers
  ordsums$areaSQKM_5_7[ordsums$species.code==sp] <- area(range)*1e-6
  
  ######## calculate high-abundance range size 
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

    hi_abunD <- as.data.frame(hi_abun)
  
    
  ## calculate infilling *** do we need this?
  pts <- species
  coordinates(pts) <- c(5,4)
  ext <- extract(hi_abun, pts) ## extract climate values to points
  pts<-data.frame(cbind(pts$abundance,ext))
  colnames(pts)<- c("observ","impact")
  ordsums$infilling[i]<- length(pts$observ[pts$observ==3 & pts$impact == 1])/length(speciesD$layer[speciesD$layer==1 &!is.na(speciesD$layer==1 )])

}

### calculate size of full study region for proportional estimates of establishment ranges
studyregion <- hi_abun ## take most recent high abundance map
stuDy <- as.data.frame(studyregion)
stuDy$layer[stuDy$layer == 1 & !is.na(stuDy$layer)] <- 0 ## assign abundance to 0, so all cells in study region are 0
studyregion <- raster(nrows=nrow(studyregion), ncols=ncol(studyregion), ext=extent(studyregion), vals=stuDy$layer)

proj4string(studyregion) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
studyregion <- rasterToPolygons(studyregion, na.rm=T, dissolve=T)
studyregion <- spTransform(studyregion, "+init=epsg:5070") ## Albers
studyregion <- area(studyregion)*1e-6 ## calculate study region area

#### calculate proportions of the study region made up by the establishment and abundance
ordsums$estab_US <- ordsums$areaSQKM_5_7/studyregion
ordsums$impac_US <- ordsums$hiAbunSQKM_5_7/studyregion
ordsums$impac_estab <- ordsums$hiAbunSQKM_5_7/ordsums$areaSQKM_5_7


################################################
################ MESS ANALYSIS #################
################################################

library(rgdal)
library(raster)
library(dismo) ## for mess()


ordsums <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/old_ordsums/ordsums_4_19_2018_updated.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

## Construct a table of the best fit models for each species
## Create new field for each variable - whether not included in each model
formulae <- data.frame(cbind(ordsums$species.code,ordsums$formu), stringsAsFactors = F)
str(formulae)
colnames(formulae) <- c("species","formula")
formulae$bio_2 <- sapply(formulae$formula, grepl, pattern="bio_2", fixed=T)
formulae$bio_5 <- sapply(formulae$formula, grepl, pattern="bio_5", fixed=T)
formulae$bio_6 <- sapply(formulae$formula, grepl, pattern="bio_6", fixed=T)
formulae$bio_8 <- sapply(formulae$formula, grepl, pattern="bio_8", fixed=T)
formulae$bio_12 <- sapply(formulae$formula, grepl, pattern="bio_12", fixed=T)
formulae$bio_15 <- sapply(formulae$formula, grepl, pattern="bio_15", fixed=T)
formulae$nlcd_3 <- sapply(formulae$formula, grepl, pattern="nlcd_3", fixed=T)
formulae$nlcd_4 <- sapply(formulae$formula, grepl, pattern="nlcd_4", fixed=T)
formulae$nlcd_5 <- sapply(formulae$formula, grepl, pattern="nlcd_5", fixed=T)
formulae$nlcd_6 <- sapply(formulae$formula, grepl, pattern="nlcd_6", fixed=T)
formulae$nlcd_7 <- sapply(formulae$formula, grepl, pattern="nlcd_7", fixed=T)
formulae$nlcd_8 <- sapply(formulae$formula, grepl, pattern="nlcd_8", fixed=T)
head(formulae)

## create list of species
spp <- formulae$species

## load environmental variables
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

## adjust temperature values, which are served multiplied by 10
bio$bio_2 <- bio$bio_2/10
bio$bio_5 <- bio$bio_5/10
bio$bio_6 <- bio$bio_6/10
bio$bio_8 <- bio$bio_8/10

cellIDs <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=seq(1:ncell(bio)))
## give a unique ID ("cell ID") to each raster cell

bio <- stack(bio,cellIDs) # "add the cell ID's to the stack of environmental variables

occ.list <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/species", full=T)
## list of occurence datasets for each species
  
## read in the dataset for each species, and concatenate into a single file of all species data
for (i in 1:length(occ.list)) {
  occ_i <- read.table(occ.list[i], sep=",",header=T,stringsAsFactors = F) 
  occ_full <- rbind(occ_full,occ_i)
  print(i)
}

head(occ_full)
coordinates(occ_full)<-c(2,3)
ext <- extract(bio$layer, occ_full)
occ_full <- cbind(occ_full,ext)
occ_full <- data.frame(occ_full)
occ_full$optional<-NULL
colnames(occ_full) <- c("PLANT_CODE","CELLID", "LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")


## for each species thin occurence data to 1 record per grid cell
for(i in 1:length(unique(occ_full$PLANT_CODE))){ ##loop through all species
  occ_i <- occ_full[occ_full$PLANT_CODE == spp[i],] ## subset occurence dataset to species
  cellids <- unique(occ_i$CELLID) ## collect all of the CELLID's for that species - all cells were species occurs
  for (j in 1:length(cellids)){ ## thin to one point per cell
    occ_i_j <- occ_i[occ_i$CELLID == cellids[j],][1,]
    occ_thinned <- rbind(occ_thinned, occ_i_j)
  }
  
  print(i)  
}


coordinates(occ_thinned)  <- c(3,4) ## designate lat,lon as coordinates
ext <- extract(bio, occ_thinned)  ## extract environmental variables at each point
occ_thinned <-cbind(occ_thinned,ext) ## add extracted values to each record


abun_all <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)
## call in abundance dataset, subset to relevant columns for mess analysis - species and location
abun_all <- data.frame(cbind(abun_all$species,abun_all$longitude,abun_all$latitude), stringsAsFactors = F)
colnames(abun_all) <- c("PLANT_CODE","LONGITUDE_DECIMAL", "LATITUDE_DECIMAL") ## fix column names
abun_all <- abun_all[abun_all$PLANT_CODE %in% spp,]
abun_all$LONGITUDE_DECIMAL <- as.numeric(abun_all$LONGITUDE_DECIMAL) ##make sure lat,lon are numeric
abun_all$LATITUDE_DECIMAL <- as.numeric(abun_all$LATITUDE_DECIMAL)

coordinates(abun_all) <- c(2,3) ## designate lat,lon as coordinates
ext <- extract(bio, abun_all) ## extract environmental values to points
abun_all <- data.frame(cbind(abun_all,ext))
abun_all$optional <- NULL

bioD <- as.data.frame(bio) ## make dataframe verstion of environmental variables
summary(bioD)

## Make a new column to store the results of the MESS analysis
## For the environmental variables included in the best fit model
## Compare the environmental space sampled by the abundance points compared to the full predicted establishment range
ordsums$MESS_abun2range_SEL <- -99


# loop throught all species
for (i in 1:length(spp)){

    ## subset abundance data to the species of the iteration
    abun <- abun_all[abun_all$PLANT_CODE==spp[i],]
    abun_vars <- abun[,2:(NCOL(abun)-3)] ## make table of environmental variable values for all abun point for that spp
    
    ## subset occurence dataset (thinned to one obs per species per cell) to spp of iteration
    occ <- data.frame(occ_thinned[occ_thinned$PLANT_CODE==spp[i],])
    occ_vars <- occ[,3:(NCOL(occ)-4)] ## make table of environ veriables for the occurence dataset
    
    range <- raster(paste("C:/Users/Localadmin/Documents/Maxent_modeling/Binary_asciis/",spp[i],".asc",sep=""))
    names(range) <- "suitability"
    
    #################### get rid of cells outside study area (set to NA) #####
    bio_occ <- bio
    bio_occ$bio_2[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$bio_5[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$bio_6[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$bio_8[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$bio_12[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$bio_15[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$nlcd_3[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$nlcd_4[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$nlcd_5[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$nlcd_6[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$nlcd_7[!(bio_occ$layer %in% occ$layer)] <- NA
    bio_occ$nlcd_8[!(bio_occ$layer %in% occ$layer)] <- NA
    ####################
    
    ## remove the cells predicted unsuitable from analysis (set to NA)
    range <- stack(range,bio)
    range$bio_2[range$suitability==0] <- NA
    range$bio_5[range$suitability==0] <- NA
    range$bio_6[range$suitability==0] <- NA
    range$bio_8[range$suitability==0] <- NA
    range$bio_12[range$suitability==0] <- NA
    range$bio_15[range$suitability==0] <- NA
    range$nlcd_3[range$suitability==0] <- NA
    range$nlcd_4[range$suitability==0] <- NA
    range$nlcd_5[range$suitability==0] <- NA
    range$nlcd_6[range$suitability==0] <- NA
    range$nlcd_7[range$suitability==0] <- NA
    range$nlcd_8[range$suitability==0] <- NA
    range <- dropLayer(range,1)
    range <- dropLayer(range,13)
    
    
    ## prep for loop below
    range <- as.data.frame(range)
    bio_occ <- as.data.frame(bio_occ)
    range2 <- cellIDs
    occ2 <- cellIDs
    
    
    ## for each environmental variable - check whether it is included in the best fit ordinal model for the species.
    ## if it's included in the best fit model, include it in the MESS analysis - otherwise, remove the variable
    if(formulae$bio_2[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_2)
      names(bio_2) <- "bio_2"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_2)
      names(bio_2) <- "bio_2"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$bio_2 <- NULL
      occ_vars$bio_2 <- NULL
    }
    
    if(formulae$bio_5[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_5)
      names(bio_2) <- "bio_5"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_5)
      names(bio_2) <- "bio_5"
      occ2<- stack(occ2, bio_2)
      
    } else {
      abun_vars$bio_5 <- NULL
      occ_vars$bio_5 <- NULL
    }
    
    if(formulae$bio_6[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_6)
      names(bio_2) <- "bio_6"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_6)
      names(bio_2) <- "bio_6"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$bio_6 <- NULL
      occ_vars$bio_6 <- NULL
    }
    
    if(formulae$bio_8[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_8)
      names(bio_2) <- "bio_8"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_8)
      names(bio_2) <- "bio_8"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$bio_8 <- NULL
      occ_vars$bio_8 <- NULL
    }
    
    if(formulae$bio_12[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_12)
      names(bio_2) <- "bio_12"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_12)
      names(bio_2) <- "bio_12"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$bio_12 <- NULL
      occ_vars$bio_12 <- NULL
    }
    
    if(formulae$bio_15[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$bio_15)
      names(bio_2) <- "bio_15"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$bio_15)
      names(bio_2) <- "bio_15"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$bio_15 <- NULL
      occ_vars$bio_15 <- NULL
    }
    
    if(formulae$nlcd_3[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_3)
      names(bio_2) <- "nlcd_3"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_3)
      names(bio_2) <- "nlcd_3"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$nlcd_3 <- NULL
      occ_vars$nlcd_3 <- NULL
    }
    
    if(formulae$nlcd_4[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_4)
      names(bio_2) <- "nlcd_4"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_4)
      names(bio_2) <- "nlcd_4"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$nlcd_4 <- NULL
      occ_vars$nlcd_4 <- NULL
    }
    
    if(formulae$nlcd_5[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_5)
      names(bio_2) <- "nlcd_5"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_5)
      names(bio_2) <- "nlcd_5"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$nlcd_5 <- NULL
      occ_vars$nlcd_5 <- NULL
    }
    
    if(formulae$nlcd_6[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_6)
      names(bio_2) <- "nlcd_6"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_6)
      names(bio_2) <- "nlcd_6"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$nlcd_6 <- NULL
      occ_vars$nlcd_6 <- NULL
    }
    
    if(formulae$nlcd_7[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_7)
      names(bio_2) <- "nlcd_7"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_7)
      names(bio_2) <- "nlcd_7"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$nlcd_7 <- NULL
      occ_vars$nlcd_7 <- NULL
    }
    
    if(formulae$nlcd_8[i]){
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=range$nlcd_8)
      names(bio_2) <- "nlcd_8"
      range2<- stack(range2, bio_2)
      
      bio_2 <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=bio_occ$nlcd_8)
      names(bio_2) <- "nlcd_8"
      occ2<- stack(occ2, bio_2)
    } else {
      abun_vars$nlcd_8 <- NULL
      occ_vars$nlcd_8 <- NULL
    }
    
    range2 <- dropLayer(range2, 1)
    occ2 <- dropLayer(occ2, 1)
    
    abun_vars <- data.frame(abun_vars)
    occ_vars <- data.frame(occ_vars)
    
    if(nlayers(occ2)== 1){
      
      range2b <- raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=as.vector(range2))
      occ2b <-   raster(nrows=nrow(bio), ncols=ncol(bio), ext=extent(bio), vals=as.vector(occ2))

      range2 <- stack(range2,range2b)
      occ2 <- stack(occ2,occ2b)
      abun_vars$layer <- abun_vars[,1]
      occ_vars$layer <- occ_vars[,1]
      
    }
    
    
    ## Calculate the % of establishment range where the mess value is > -10 
    ## (percent where the environmental space is within or not far beyond 
    ## the environmental space represented in the abundance points)
    mess <- mess(range2, abun_vars, full=FALSE) 
    messD <- as.data.frame(mess)
    colnames(messD) <- "mess"
    messD <- messD[!is.na(messD$mess) & messD$mess < Inf,]
    ordsums$MESS_abun2range_SEL[i] <- length(messD[messD > -10])/length(messD)
    
  
  
  print(i)
  
}



ordsums$MESS_abun2range_SEL[ordsums$MESS_abun2range_SEL== -99] <- NA
summary(ordsums)

# Line 1520

################################################
########### HOTSPOT ANALYSIS ###################
################################################

library(raster)
library(rgdal)
library(ordinal)
library(reshape2)


ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5epv_up_5_7_2018.csv", header = T, sep = ",",  
                      quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)

head(ordsums)
ordsums2 <- ordsums[ordsums$kappa > 0 & ordsums$kappaP < 0.05 & ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > 0.9,]
### filter out species for which the best-fit models did not meet th criteria

filenames <- data.frame(ordsums2$species.code) #make list of species code, with intention of converting them to the filenames that can be used to load in the exported rangemaps
colnames(filenames) <- "x"

### Create map of number of overlapping abundance ranges per grid cell
######################################################################
filenames$files <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/HI_ABUN/", filenames$x,".asc")
### Re-generate the filename for each species high abundance range map

ldf=stack(filenames$files) #stack the abundance range maps
overlap=calc(ldf, sum) #calculate the sum of overlapping abundance ranges for each grid cell
plot(overlap)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc"
writeRaster(overlap, filename=out_file,overwrite=T) #Write out map with number of high abundance species in ach grid cell

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/summaries/figures/richness_map_HI_ABUN.pdf",width=11)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()## export pdf of abundance map


### Create map of number of overlapping establishment ranges per grid cell
######################################################################
filenames$files2 <- paste0("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", filenames$x,".asc")
### Re-generate the filename for each species establishment range map
ldf=stack(filenames$files2)
overlap=calc(ldf, sum) #calculate the sum of overlapping establishment ranges for each grid cell
plot(overlap)

out_file="C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_FULL.asc"
writeRaster(overlap, filename=out_file, overwrite=T)

pdf(file="C:/Users/Localadmin/Documents/Maxent_modeling/summaries/figures/richness_map_FULL.pdf",width=11,height=8.5)
colors=rev(heat.colors(100))
spplot(overlap, col.regions=colors)
dev.off()

###############################################################################################################
### Hotspot maps
################################################################################################################

## ABUNDANCE HOTSPOT
## load in map of the total overlapping abundance ranges per grid square
current_rich=raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc") 
quantile(current_rich$richness_map_HI_ABUN)

rich <- as.data.frame(current_rich) # convert to data frame to perform hotspot calculation

## convert to binary hotspot map: 1=abundance range hotspot (greater than 75% of grid cells), 0= non-hotspot
rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN <= quantile(current_rich$richness_map_HI_ABUN, 0.75)] <- 0
rich$richness_map_HI_ABUN[rich$richness_map_HI_ABUN > quantile(current_rich$richness_map_HI_ABUN, 0.75)] <- 1

## convert back into map and export
hi_abun_hotspot <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_HI_ABUN)
writeRaster(hi_abun_hotspot, "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_HI_ABUN.asc",overwrite=T)
plot(hi_abun_hotspot)

########################################################################
## ESTABLISHMENT HOTSPOT
## load in map of the total overlapping abundance ranges per grid square
current_rich=raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/richness_map_FULL.asc")
quantile(current_rich$richness_map_FULL)

rich$richness_map_FULL[rich$richness_map_FULL <= quantile(current_rich$richness_map_FULL, 0.75)] <- 0
rich$richness_map_FULL[rich$richness_map_FULL > quantile(current_rich$richness_map_FULL, 0.75)] <- 1

## convert back into map and export
full_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=rich$richness_map_FULL)
plot(full_hotspot)
writeRaster(full_hotspot, "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/hotspot_FULL.asc", overwrite=T)

#######################################################################
## HOTSPOT COMPARISON
comp <- as.data.frame(stack(hi_abun_hotspot,full_hotspot)) #stack hotspots
head(comp)
colnames(comp) <- c("layer.1", "layer.2")

comp$x[comp$layer.1==0 & comp$layer.2==0 & !is.na(comp$layer.2)] <- "non-hotspot" #neither hotspot
comp$x[comp$layer.1==1 & comp$layer.2==0 & !is.na(comp$layer.2)] <- "abundance hotspot" #just abundance hotspot
comp$x[comp$layer.1==1 & comp$layer.2==1 & !is.na(comp$layer.2)] <- "both hotspots" #both hotspots
comp$x[comp$layer.1==0 & comp$layer.2==1 & !is.na(comp$layer.2)] <- "richness hotspot" #just establishment hotspot

length(comp$x[comp$x == "both hotspots" & !is.na(comp$x)])/length(comp$layer.1[comp$layer.1 == 1 & !is.na(comp$layer.1)]) ## percent of abundance hotspot also richness
length(comp$x[comp$x == "both hotspots" & !is.na(comp$x)])/length(comp$layer.1[comp$layer.2 == 1 & !is.na(comp$layer.1)]) ## percent of richness hotspot also abundance
.5*(length(comp$x[comp$x == "both hotspots"& !is.na(comp$x)])/length(comp$layer.1[comp$layer.1 == 1& !is.na(comp$x)]) + length(comp$x[comp$x == "both hotspots"& !is.na(comp$x)])/length(comp$layer.1[comp$layer.2 == 1 & !is.na(comp$x)]))

## 47% overlap

#Abun <- 89975
#occ <- 11118
comp$x <- as.factor(comp$x) 

comp_hotspot<- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=comp$x)
#turn comparison back into raster
plot(comp_hotspot)
spplot(comp_hotspot)

writeRaster(comp_hotspot,"C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc",overwrite=T)
# write out the hotspot comparison

################## proportion of abundant species
abun_prop=stack("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/richness_map_FULL.asc",
                "C:/Users/Localadmin/Documents/MaxEnt_modeling/summaries/asciis/richness_map_HI_ABUN.asc")
abun_prop <- as.data.frame(abun_prop)
abun_prop$abun_prop <- abun_prop$richness_map_HI_ABUN/abun_prop$richness_map_FULL

summary(abun_prop$abun_prop)
hist(abun_prop$abun_prop)

abun_prop <- raster(nrows=nrow(current_rich), ncols=ncol(current_rich), ext=extent(current_rich), vals=abun_prop$abun_prop)#,


writeRaster(abun_prop,"C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hi_abun_proportion.asc",overwrite=T)

###############################
#### Calculate Range Sizes ####
###############################
library(raster)
library(rgdal)
library(rgeos)

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_4_2018.csv", header = T, sep = ",", quote= "\"", comment.char= "", stringsAsFactors = F, strip.white = T)  
ordsums2 <- ordsums[ordsums$AUC > 0.7 & ordsums$MESS_abun2range_SEL > 0.9 & ordsums$kappa > 0 & ordsums$kappaP < 0.05,]  
head(ordsums)

### Do all 155 species
ranges <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/", full.names = T)
abuns  <- list.files("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL_5epv_round_up/HI_ABUN/", full.names = T)



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

#1382
