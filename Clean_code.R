#### Cleaned code for publication
###  Mitch O'Neill
##   created: 10/20/2019 
#    last modified: 10/20/2019


#### files needed
## 1. environmental layers
## 2. list of species ---- sp.list <- unique(speciesdata$code)
## 3. bias file ---- should I make it in the code
## 4. ordsums ---- should I make it in the code?

## define filepath of MaxEnt on comuputer
maxent.location='C:/Users/Localadmin/Documents/MaxEnt/maxent.jar'


## define directory of environmental data for fitting
environmental=paste("environmentallayers=","C:/Users/Localadmin/Documents/MaxEnt_modeling/envi_ASCIIs",sep="")

## list of species
sp.list <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/spForms3_4_2018.csv", header = T, sep = ",", 
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
#### INSERT FILE HERE - edd with environmental ***
edd <- read.table("edd_w_environmental_NEW.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

#### INSERT FILE HERE - species list; max.vars needs to be available or created
ordsums <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/species_list_4_27b.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)

ordsums$species.code <- ordsums$species
ordsums <- ordsums[order(ordsums$max.vars.5epv),]

ordsums$pS <- -99           ## *****
ordsums$pS10 <- -99
ordsums$pI <- -99
ordsums$pI10 <- -99
ordsums$c.up <- -99         ## number of concave-up response curves?***
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
    } 
  
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
  
  ### ************
  MsumS <- data.frame(Msum[3:nrow(Msum),]) ## 1 back to 2
  ordsums$pS[s] <- length(MsumS$Estimate[MsumS$Pr...z..>0.05])
  ordsums$pS10[s] <- length(MsumS$Estimate[MsumS$Pr...z..>0.1])
  
  MsumI <- data.frame(Msum[1:2,]) ## 2 back to 3
  ordsums$pI[s] <- length(MsumI$Estimate[MsumI$Pr...z..>0.05])
  ordsums$pI10[s] <- length(MsumI$Estimate[MsumI$Pr...z..>0.1])
  
  MsumS$name <- as.character(row.names(MsumS))
  MsumS$poly <- sapply(MsumS$name, grepl, pattern="2)2")
  
  ordsums$c.up[s] <- length(MsumS$Estimate[MsumS$Estimate > 0 & MsumS$poly == T])
  ordsums$c.down[s] <- length(MsumS$Estimate[MsumS$Estimate < 0 & MsumS$poly == T])
  
  ordsums$aic[s] <- AIC(M)          ## extract adjusted R-square value of model #adjDsq <-  Dsquared(M,adjust=T) ## extract deviance explained by model
  #no.pts2 <- nrow(species)
  ordsums$no.terms[s] <- NROW(Msum) ## extract the number of terms
  ordsums$no.vars[s] <- ordsums$no.terms[s] - 2 - str_count(formulae.s$formula[model.sel2$i], "poly")
  ordsums$formu[s] <- formulae.s$formula[model.sel2$i]
  
  null <- clm(abundance~1,data=species, Hess=T)
  
  if(AIC(null) > (AIC(M) + 2)){ordsums$null[s] <- 0} else {ordsums$null[s] <- 1}  
  
  pm <- predict(M, species, type="class")
  pm <- melt(pm)
  pm <- pm[,1]
  #length(pm$)
  
  pm <- data.frame(cbind(species$abundance,pm))
  colnames(pm) <- c("observed", "predicted")
  
  kappa <- kappa2(pm)
  ordsums$kappaP[s] <- kappa$p.value
  ordsums$kappa[s] <- kappa$value
  
  ordsums$bin1o[s] <- length(pm$observed[pm$observed == 1])
  ordsums$bin2o[s] <- length(pm$observed[pm$observed == 2])
  ordsums$bin3o[s] <- length(pm$observed[pm$observed == 3])
  ordsums$bin1p[s] <- length(pm$observed[pm$predicted == 1])
  ordsums$bin2p[s] <- length(pm$observed[pm$predicted == 2])
  ordsums$bin3p[s] <- length(pm$observed[pm$predicted == 3])
  
  write.csv(ordsums, "progressMACHINE1.csv",row.names=F)
  
  #bin1o <- length(pm$observed[pm$observed == 1])
  #bin2o <- length(pm$observed[pm$observed == 2])
  #bin3o <- length(pm$observed[pm$observed == 3])
  #bin1p <- length(pm$predicted[pm$predicted == 1])
  #bin2p <- length(pm$predicted[pm$predicted == 2])
  #bin3p <- length(pm$predicted[pm$predicted == 3])
  
  
  ###############################
  #reg2 <- summary(glm(as.numeric(predicted)~as.numeric(observed), data=pm))$coefficients
  #regI2 <- reg2[1,1]
  #regS2 <- reg2[2,1]
  # one <- c(length(pm$observed[pm$predicted == 1 & pm$observed == 1]),
  #          length(pm$observed[pm$predicted == 1 & pm$observed == 2]),
  #          length(pm$observed[pm$predicted == 1 & pm$observed == 3]),
  #          length(pm$observed[pm$predicted == 1 & pm$observed == 4]))
  # 
  # two <- c(length(pm$observed[pm$predicted == 2 & pm$observed == 1]),
  #          length(pm$observed[pm$predicted == 2 & pm$observed == 2]),
  #          length(pm$observed[pm$predicted == 2 & pm$observed == 3]),
  #          length(pm$observed[pm$predicted == 2 & pm$observed == 4]))
  # 
  # thr <- c(length(pm$observed[pm$predicted == 3 & pm$observed == 1]),
  #          length(pm$observed[pm$predicted == 3 & pm$observed == 2]),
  #          length(pm$observed[pm$predicted == 3 & pm$observed == 3]),
  #          length(pm$observed[pm$predicted == 3 & pm$observed == 4]))  
  # 
  # fou <- c(length(pm$observed[pm$predicted == 4 & pm$observed == 1]),
  #          length(pm$observed[pm$predicted == 4 & pm$observed == 2]),
  #          length(pm$observed[pm$predicted == 4 & pm$observed == 3]),
  #          length(pm$observed[pm$predicted == 4 & pm$observed == 4])) 
  # 
  # confus <- data.frame(one,two,thr,fou)
  ###############################
  #ordsums.i <- data.frame(species.code,aic, no.pts, no.terms, no.vars, formu, c.up, c.down, pS, pS10, pI, pI10,
  #                          null, kappa, kappaP, bin1o, bin2o, bin3o, bin1p, bin2p, bin3p, stringsAsFactors = F)
  
  #ordsums <- rbind(ordsums,ordsums.i) 
  #print(s)
  #print("SPECIES")
  # if statement to remove vmin species}
  print(s)
}
