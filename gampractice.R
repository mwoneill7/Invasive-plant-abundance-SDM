#### Practicing with GAMS
###  Mitch O'Neill
##   created: 9/4/2017
#    last modified: 10/9/2017

#         ## initialize template masterfile to contain summaries of each model for this species
#         i <- -1 ## keep track of place in formula list
#         n <- -1 ## number of variables
#         aic <- -1 ## Akaike's info criterion
#         bios <- "TEMPLATE" ## list of climate covariates
#         nonsig <- -1 ## 1: there is a non-sig variable, 2: all variables are signif
#         collinear <- -1
#         dAIC <- -1
#         species <- "TEMPLATE"
#         
#         model.selX <- data.frame(i, n, aic, bios, nonsig, collinear, dAIC, species, stringsAsFactors = F)
#         ## concatenate into master data.frame
         
## 1. commenting

library(mgcv)     ## for GAMS
library(raster)   ## for Raster climate data
library(rgeos)    ## for extractig climate data to points
library(reshape2) ## for reformatting model predictions
library(combinat) ## for generating pairs of covariates for assessment of collinearity 
#library(rgdal)

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/")

## read in EDDMapS dataset, thinned to 5km climate grid cells (1 point per species per cell)
edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_10acre.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
#head(edd)

## read in Worldclim climate rasters clipped to L48; 6 variables manually selected by assessing correlations
bio <- stack("C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_2.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_5.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_6.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_8.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_12.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_15.asc")
## assign proj4string to that of one of the raw files before clipping
proj4string(bio) <- proj4string(raster("C:/Users/mwone/Documents/geodata/climate_data/current/bio_1"))

## read in list of model formulae with all 63 possible combos of 6 variables (created with help of combinat)
formulae <- read.table("formulae.txt", header = T, sep = "%", quote= "\"", 
                        comment.char= "", stringsAsFactors = F, strip.white = T)
tp = "tp" ## so that the tp referred to in the formulae can be interpreted as thin-plate in mgcv
  
## read in table with species list and useability info
sp.list <- read.table("SpeciesList10_11_2017.csv", header = T, sep = ",", quote= "\"", 
                           comment.char= "", stringsAsFactors = F, strip.white = T)
#head(sp.list)
sp.list <- sp.list$usda.code[sp.list$potential.use_10ac == 1] ## subset to useable species

species.code <- "TEMPLATE" ## initialize master file to contain summaries of GAMs (final model for each species)
Radj <- -99 ## adjusted R-squared of models
devX <- -99 ## deviance explained by models
no.vars <- -99 ## number of variables in models
nonsig <- -99
p2 <- -99 ## significance (p-val) of each variable (2:bio_2, etc.)
p5 <- -99
p6 <- -99
p8 <- -99
p12 <- -99
p15 <- -99
edf2 <- -99 ## effective degrees of freedom for each smooth (edf2:s(bio_2), etc.)
edf5 <- -99
edf6 <- -99
edf8 <- -99
edf12 <- -99
edf15 <- -99

GAMsums <- data.frame(species.code,Radj,devX,no.vars,nonsig,p2,p5,p6,p8,p12,p15,edf2,edf5,edf6,edf8,edf12,edf15, stringsAsFactors=F)
## concatenate into template dataframe

for (s in 1:20){  #1:length(sp.list)){ ## Loop through species list
  
  #species.code <- "VEIN9"
  species.code <- sp.list[s] ## extract the USDA species code for the species of the iteration
  
  species <- edd[edd$usda == species.code,]## subset to species of the iteration
  species <- data.frame(cbind(log(species$abundance), species$latitude, species$longitude)) ## select only variables needed to save computation time
  colnames(species) <- c("abundance", "latitude", "longitude") ## clean up column names

  if (max(species$abundance) > min(species$abundance)){ ## Also, skip over species with no variation in abundance values.
     
    #par(mfrow=c(1,1)) ## clean up plot space from past iterations
    #jpeg(paste("histograms", paste(species.code, "jpg", sep="."), sep="/" )) ## output file name for histogram of abun
    #hist(species$abundance,  breaks=20) # histogram of abundance with 20 breaks
    #dev.off() ## export plot to filepath specified above
    
    coordinates(species) <- c(3,2) ## specify lon/lat
    proj4string(species) <- CRS("+init=epsg:4326") ## assume WGS84/NAD83 for EDDMapS points
    
    species <- spTransform(species, proj4string(bio)) ## project point data to proj4string of climate data
    
    ext <- extract(bio, species) ## extract climate values to points
    species <- cbind(as.data.frame(species), ext) ## append extracted climate data to point data
    species[,4:7] <- species[,4:7]/10 ## convert temp variables back to Celcius (from Celsius*10)
    #head(species) 
  
    #### MODEL SELECTION #####
    ## initialize template masterfile to contain summaries of each model for this species
    i <- -1 ## keep track of place in formula list
    n <- -1 ## number of variables
    aic <- -1 ## Akaike's info criterion
    bios <- "TEMPLATE" ## list of climate covariates
    nonsig <- -1 ## 1: there is a non-sig variable, 2: all variables are signif
    
    model.sel <- data.frame(i, n, aic, bios, nonsig, stringsAsFactors = F)
    ## concatenate into master data.frame
    
    for (i in 1:length(formulae$form)){ ##loop through all formulae in the formula list
     
      if(i != 63 | length(species$abundance) > 23) { 
        ## skip over model w/ 6 params for species w/ <23 points, for coefficients would outumber data.

        M <- gam(formula(formulae$form[i]), method ="GCV.Cp", data = species)
        ## construct model using formula of the iteration
        ## select smoothness using REML
        ## use data for species of this iteration
      
        Msum <- summary(M)
        Msum <- as.data.frame(Msum$s.table) ## extract table with significance of smooth terms
        n <- NROW(Msum) ## number of frows in smooth table = number of smooths = number of variables
        
        nonsig <- 0 ## assume significance by default
        for (j in 1:n){ ## loop through smooths
          if(Msum$'p-value'[j] > 0.01) { ## if any smooth has a p-value greater than .01
            nonsig <- 1 + nonsig ## then note that there is a nonsig variable
          }
        }
      
        aic <- M$aic ## extract aic score
        
        bios <- as.character(M$pred.formula[2]) ## list covariates
       
        model.sel.i <- data.frame(i, n, aic, bios, nonsig, stringsAsFactors = F)#bios, nonsig, stringsAsFactors = F) ## concatenate data from ieration
        model.sel <- rbind(model.sel, model.sel.i) ## append data from iteration to the master data frame
        #print(i) ## keep track of place in loop
      }
    }
    
    model.sel <- model.sel[model.sel$n != -1,] ## remove template row
    #head(model.sel)
    
    ########### ASSESS CORRELATIONS ############
    rhos <- as.data.frame(cor(ext, use="pairwise", method="spearman")) 
    ## construct data.frame of matrix of pairwise spearman's rho values for climatic variables
  
    var.pairs<- combn(c(1:6), 2, simplify=T) ## All possible pairs selected from 6 variables
    
    # initialize master data frameof pairwise correlations with template
    var1 <- "TEMPLATE" ## one variable from pairwise correlation
    var2 <- "TEMPLATE" ## other variable from pairwise correlation
    rho <- -2 ## rho value from correlation
    corr.pairs <- data.frame(var1,var2,rho, stringsAsFactors = F) ## create master data frame
    
    for (i in 1:NCOL(var.pairs)) { ## loop through list of pairs of the six variables
    
    rho <- abs(rhos[var.pairs[1,i], var.pairs[2,i]]) 
    ## absolute value of spearman's rho at the row and column corresponding to the pair in the iteration
    var1 <- row.names(rhos[var.pairs[1,i],]) ## name of the first variable in the pair
    var2 <- row.names(rhos[var.pairs[2,i],]) ## name of the second variable in the pair
    
    corr.pairs.i <- data.frame(var1,var2,rho, stringsAsFactors = F) ## concatenate data from the interation
    corr.pairs <- rbind(corr.pairs, corr.pairs.i) ## append data from iteraction to the master data frame
    
    #print(i) ## keep track of place in loop
  }
  
    corr.pairs <- corr.pairs[corr.pairs$var1 != "TEMPLATE",] ## remove template row
    corr.pairs <- corr.pairs[corr.pairs$rho > .6, ] ## select all pairwise correlations where collinearity is likely 

    model.sel$collinear <- 0 ## assume correlation between 2 vars is below 0.8
    
    if(length(corr.pairs$rho) > 0){ ## if there were any pairs where rho >= 0.8
      for (i in 1:length(corr.pairs$rho)){ ## then loop through those pairs
        for (j in 1:length(model.sel$bios)){ ## for each pair, loop through all of the possible models
          ## and if the two highly correlated variables appear together in a model,
          if(grepl(corr.pairs$var1[i], model.sel$bios[j]) & grepl(corr.pairs$var2[i], model.sel$bios[j])) { 
            model.sel$collinear[j] <- 1} ## then flag the model as a collinearity risk
          #print(j) ## keep track of place in loop w/in loop
        }
        #print(paste("i =", i, sep = " ")) ## keep track of place in loop
      }
    }
    
    rm(rhos, var1, var2, corr.pairs, var.pairs) ## garbage cleaning
  
    model.sel <- model.sel[model.sel$collinear == 0,] ## exclude models with collinearity
    #model.sel <- model.sel[model.sel$nonsig == 0,] ## exclude models with nonsignificant smooths
    
    ##########################################################
    
    #if (length(model.sel$i) > 0) {
      
      ###################  AIC comparison  #####################
      model.sel$dAIC <- model.sel$aic - min(model.sel$aic) ## calculate delta AIC
      
     #    model.selX.i <- model.sel[model.sel$dAIC <= 10,]
      #   model.selX.i$species <- species.code
       #  model.selX <- rbind(model.selX, model.selX.i)
         
      model.sel2 <- model.sel[model.sel$dAIC <= 2,] ## select the model with the best score, 
      ## as well as models with scores that are not significantly worse
      
      model.sel2 <- model.sel2[model.sel2$n == min(model.sel2$n),] 
      ## of the models with the best scores (<= min + 2), select the simplest models
      
      model.sel2 <- model.sel2[model.sel2$dAIC == min(model.sel2$dAIC),]
      ## select the model with the best AIC score if there are multiple models
      ## with similar scores and the same complexity (variable number)
      
      nonsig <- model.sel2$nonsig
      
      #model.selS <- model.sel[model.sel$nonsig == 0,] 
      #model.selS$dAIC <- model.selS$aic - min(model.selS$aic) 
       
      M <- gam(formula(formulae$form[model.sel2$i]), method="GCV.Cp", data = species)
      ## reconstruct the best model, using the formula number (i) to re-create it
     
      Msum <- summary(M) ## extract summary of the model to access model-fit values
      smoothSum <- as.data.frame(Msum$s.table) ## extract table with edf and p-vals for smooths
      
      Radj <- Msum$r.sq ## extract adjusted R-square value of model
      devX <- Msum$dev.expl ## extract deviance explained by model
      no.vars <- NROW(smoothSum) ## extract the number of smooths (number of variables)
      
      ## extract significance (p-vals) and effective degrees of freedom for each variables
      p2 <- -1 ## start by defaulting each to -1, indicating that the variable is not in the model
      p5 <- -1
      p6 <- -1
      p8 <- -1
      p12 <- -1
      p15 <- -1
      edf2 <- -1
      edf5 <- -1
      edf6 <- -1
      edf8 <- -1
      edf12 <- -1
      edf15 <- -1
      
      bios <- as.character(M$pred.formula[2])
      if (grepl("bio_2", bios)) {
        p2 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_2)"]
        edf2 <- smoothSum$edf[row.names(smoothSum) == "s(bio_2)"]
      }
      if (grepl("bio_5", bios)) {
        p5 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_5)"]
        edf5 <- smoothSum$edf[row.names(smoothSum) == "s(bio_5)"]
      }
      if (grepl("bio_6", bios)) {
        p6 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_6)"]
        edf6 <- smoothSum$edf[row.names(smoothSum) == "s(bio_6)"]
      }
      if (grepl("bio_8", bios)) {
        p8 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_8)"]
        edf8 <- smoothSum$edf[row.names(smoothSum) == "s(bio_8)"]
      }
      if (grepl("bio_12", bios)) {
        p12 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_12)"]
        edf12 <- smoothSum$edf[row.names(smoothSum) == "s(bio_12)"]
      }
      if (grepl("bio_15", bios)) {
        p15 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_15)"]
        edf15 <- smoothSum$edf[row.names(smoothSum) == "s(bio_15)"]
      }
      
      
      jpeg(paste("model_pics", paste(species.code, "jpg", sep="."), sep="/"))
      
      if (no.vars > 4){
        par(mfrow=c(2,3)) 
      }
      if (no.vars == 4){
        par(mfrow=c(2,2)) 
      }
      if (no.vars == 3){
        par(mfrow=c(3,1)) 
      }
      if (no.vars == 2){
        par(mfrow=c(2,1)) 
      }
  
      plot(M)
      dev.off()
      
      GAMsums.i <- data.frame(species.code,Radj,devX,no.vars,nonsig,p2,p5,p6,p8,p12,p15,edf2,edf5,edf6,edf8,edf12,edf15, stringsAsFactors=F)
      GAMsums <- rbind(GAMsums,GAMsums.i)

#        ##### SIGNIF VERSION #####     
#      if (model.sel2$nonsig == 1 & length(model.selS$nonsig) > 0) {
#        ## repeat with best model with all signif terms
#        M <- gam(formula(formulae$form[model.selS$i]), method="GCV.Cp", data = species)
#        ## reconstruct the best model, using the formula number (i) to re-create it
#        
#        Msum <- summary(M) ## extract summary of the model to access model-fit values
#        smoothSum <- as.data.frame(Msum$s.table) ## extract table with edf and p-vals for smooths
#        
#        Radj <- Msum$r.sq ## extract adjusted R-square value of model
#        devX <- Msum$dev.expl ## extract deviance explained by model
#        no.vars <- NROW(smoothSum) ## extract the number of smooths (number of variables)
#        
#        ## extract significance (p-vals) and effective degrees of freedom for each variables
#        p2 <- -1 ## start by defaulting each to -1, indicating that the variable is not in the model
#        p5 <- -1
#        p6 <- -1
#        p8 <- -1
#        p12 <- -1
#        p15 <- -1
#        edf2 <- -1
#        edf5 <- -1
#        edf6 <- -1
#        edf8 <- -1
#        edf12 <- -1
#        edf15 <- -1
#        
#        bios <- as.character(M$pred.formula[2])
#        if (grepl("bio_2", bios)) {
#          p2 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_2)"]
#          edf2 <- smoothSum$edf[row.names(smoothSum) == "s(bio_2)"]
#        }
#        if (grepl("bio_5", bios)) {
#          p5 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_5)"]
#          edf5 <- smoothSum$edf[row.names(smoothSum) == "s(bio_5)"]
#        }
#        if (grepl("bio_6", bios)) {
#          p6 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_6)"]
#          edf6 <- smoothSum$edf[row.names(smoothSum) == "s(bio_6)"]
#        }
#        if (grepl("bio_8", bios)) {
#          p8 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_8)"]
#          edf8 <- smoothSum$edf[row.names(smoothSum) == "s(bio_8)"]
#        }
#        if (grepl("bio_12", bios)) {
#          p12 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_12)"]
#          edf12 <- smoothSum$edf[row.names(smoothSum) == "s(bio_12)"]
#        }
#        if (grepl("bio_15", bios)) {
#          p15 <- smoothSum$`p-value`[row.names(smoothSum) == "s(bio_15)"]
#          edf15 <- smoothSum$edf[row.names(smoothSum) == "s(bio_15)"]
#        }
#        
#        
#        jpeg(paste("model_pics", paste(species.code, "jpg", sep="SIG."), sep="/"))
#        
#        if (no.vars > 4){
#          par(mfrow=c(2,3)) 
#        }
#        if (no.vars == 4){
#          par(mfrow=c(2,2)) 
#        }
#        if (no.vars == 3){
#          par(mfrow=c(3,1)) 
#        }
#        if (no.vars == 2){
#          par(mfrow=c(2,1)) 
#        }
#        
#        plot(M)
#        dev.off()
#        
#        GAMsums.i <- data.frame(species.code,Radj,devX,no.vars,nonsig,p2,p5,p6,p8,p12,p15,edf2,edf5,edf6,edf8,edf12,edf15, stringsAsFactors=F)
#        GAMsums <- rbind(GAMsums,GAMsums.i)
#        
#      }
        
#        ############## EXTRAPOLATION ################
#         
#         bio.sp <- as.data.frame(bio)
#         bio.sp[,1:4] <- bio.sp[,1:4]/10
#         
#         ## if bio_2 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
#         if(grepl("bio_2", bios)) {
#     
#     bio_2 <- c(min(species$bio_2), min(species$bio_2 - 10), max(species$bio_2), max(species$bio_2) + 10)
#     ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#     ## to detect trending at min and max
#     
#     input <- data.frame(bio_2) ## change list to data.frame
#     input$bio_5 <- 0 ## set all other variables to constant value
#     input$bio_6 <- 0
#     input$bio_8 <- 0
#     input$bio_12 <- 0
#     input$bio_15 <- 0
#     
#     ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#     p <- predict.gam(M, input)
#     
#     if (p[1]>p[2]) {bio.sp$bio_2[bio.sp$bio_2 < min(species$bio_2)]  <- min(species$bio_2)
#     } ## TRUE: upward trend on left end
#     ## reassign values below min to simulate a flatline
#     
#     if (p[4]>p[3]) {bio.sp$bio_2[bio.sp$bio_2 > max(species$bio_2)]  <- max(species$bio_2) 
#     }## TRUE: upward trend on right end
#     ## reassign values above max to simulate a flatline
#   }
#         
#         ## if bio_5 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
#         if(grepl("bio_5", bios)) {
#     
#     bio_5 <- c(min(species$bio_5), min(species$bio_5 - 10), max(species$bio_5), max(species$bio_5) + 10)
#     ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#     ## to detect trending at min and max
#     
#     input <- data.frame(bio_5) ## change list to data.frame
#     input$bio_2 <- 0 ## set all other variables to constant value
#     input$bio_6 <- 0
#     input$bio_8 <- 0
#     input$bio_12 <- 0
#     input$bio_15 <- 0
#     
#     ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#     p <- predict.gam(M, input)
#     
#     if (p[1]>p[2]) {bio.sp$bio_5[bio.sp$bio_5 < min(species$bio_5)]  <- min(species$bio_5)
#     } ## TRUE: upward trend on left end
#     ## reassign values below min to simulate a flatline
#     
#     if (p[4]>p[3]) {bio.sp$bio_5[bio.sp$bio_5 > max(species$bio_5)]  <- max(species$bio_5) 
#     }## TRUE: upward trend on right end
#     ## reassign values above max to simulate a flatline
#   }
#         
#         ## if bio_6 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
#         if(grepl("bio_6", bios)) {
#     
#     bio_6 <- c(min(species$bio_6), min(species$bio_6 - 10), max(species$bio_6), max(species$bio_6) + 10)
#     ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#     ## to detect trending at min and max
#     
#     input <- data.frame(bio_6) ## change list to data.frame
#     input$bio_2 <- 0 ## set all other variables to constant value
#     input$bio_5 <- 0
#     input$bio_8 <- 0
#     input$bio_12 <- 0
#     input$bio_15 <- 0
#     
#     ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#     p <- predict.gam(M, input)
#     
#     if (p[1]>p[2]) {bio.sp$bio_6[bio.sp$bio_6 < min(species$bio_6)]  <- min(species$bio_6)
#     } ## TRUE: upward trend on left end
#     ## reassign values below min to simulate a flatline
#     
#     if (p[4]>p[3]) {bio.sp$bio_6[bio.sp$bio_6 > max(species$bio_6)]  <- max(species$bio_6) 
#     }## TRUE: upward trend on right end
#     ## reassign values above max to simulate a flatline
#   }
#         
#         ## if bio_8 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
#         if(grepl("bio_8", bios)) {
#     
#     bio_8 <- c(min(species$bio_8), min(species$bio_8 - 10), max(species$bio_8), max(species$bio_8) + 10)
#     ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#     ## to detect trending at min and max
#     
#     input <- data.frame(bio_8) ## change list to data.frame
#     input$bio_2 <- 0 ## set all other variables to constant value
#     input$bio_5 <- 0
#     input$bio_6 <- 0
#     input$bio_12 <- 0
#     input$bio_15 <- 0
#     
#     ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#     p <- predict.gam(M, input)
#     
#     if (p[1]>p[2]) {bio.sp$bio_8[bio.sp$bio_8 < min(species$bio_8)]  <- min(species$bio_8)
#     } ## TRUE: upward trend on left end
#     ## reassign values below min to simulate a flatline
#     
#     if (p[4]>p[3]) {bio.sp$bio_8[bio.sp$bio_8 > max(species$bio_8)]  <- max(species$bio_8) 
#     }## TRUE: upward trend on right end
#     ## reassign values above max to simulate a flatline
#   }
#         
#         ## if bio_12 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
#         if(grepl("bio_12", bios)) {
#     
#     bio_12 <- c(min(species$bio_12), min(species$bio_12 - 10), max(species$bio_12), max(species$bio_12) + 10)
#     ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#     ## to detect trending at min and max
#     
#     input <- data.frame(bio_12) ## change list to data.frame
#     input$bio_2 <- 0 ## set all other variables to constant value
#     input$bio_5 <- 0
#     input$bio_6 <- 0
#     input$bio_8 <- 0
#     input$bio_15 <- 0
#     
#     ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#     p <- predict.gam(M, input)
#     
#     if (p[1]>p[2]) {bio.sp$bio_12[bio.sp$bio_12 < min(species$bio_12)]  <- min(species$bio_12)
#     } ## TRUE: upward trend on left end
#     ## reassign values below min to simulate a flatline
#     
#     if (p[4]>p[3]) {bio.sp$bio_12[bio.sp$bio_12 > max(species$bio_12)]  <- max(species$bio_12) 
#     }## TRUE: upward trend on right end
#     ## reassign values above max to simulate a flatline
#   }
#         
#         ## if bio_15 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
#         if(grepl("bio_15", bios)) {
#     
#     bio_15 <- c(min(species$bio_15), min(species$bio_15 - 10), max(species$bio_15), max(species$bio_15) + 10)
#     ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#     ## to detect trending at min and max
#     
#     input <- data.frame(bio_15) ## change list to data.frame
#     input$bio_2 <- 0 ## set all other variables to constant value
#     input$bio_5 <- 0
#     input$bio_6 <- 0
#     input$bio_8 <- 0
#     input$bio_12 <- 0
#     
#     
#     ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#     p <- predict.gam(M, input)
#     
#     if (p[1]>p[2]) {bio.sp$bio_15[bio.sp$bio_15 < min(species$bio_15)]  <- min(species$bio_15)
#     } ## TRUE: upward trend on left end
#     ## reassign values below min to simulate a flatline
#     
#     if (p[4]>p[3]) {bio.sp$bio_15[bio.sp$bio_15 > max(species$bio_15)]  <- max(species$bio_15) 
#     }## TRUE: upward trend on right end
#     ## reassign values above max to simulate a flatline
#   }
#         
#        ############# PREDICTIONS #############
#         P <- predict.gam(M, bio.sp)
#         P <- melt(P)
#         
#         #P$value[P$value < -5] <- NA
#         predict.raster <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=P$value)
#         writeRaster(predict.raster, filename= paste("currentPreds", paste(species.code, "asc", sep="."), sep="/" ), format="ascii", overwrite=T)
      
    # } else {
    #   ## use -99 to indicate that models were insufficient
    #   Radj <- -99
    #   devX <- -99
    #   no.vars <- -99
    #   p2 <- -99 
    #   p5 <- -99
    #   p6 <- -99
    #   p8 <- -99
    #   p12 <- -99
    #   p15 <- -99
    #   edf2 <- -99
    #   edf5 <- -99
    #   edf6 <- -99
    #   edf8 <- -99
    #   edf12 <- -99
    #   edf15 <- -99
    #   GAMsums.i <- data.frame(species.code,Radj,devX,no.vars,nonsig,p2,p5,p6,p8,p12,p15,edf2,edf5,edf6,edf8,edf12,edf15, stringsAsFactors=F)
    #   GAMsums <- rbind(GAMsums,GAMsums.i)
    # }
    
  } else {
    ## use -99 to indicate that variation in abundance was insufficient
    Radj <- -99
    devX <- -99
    no.vars <- -99
    nonsig <- -99
    p2 <- -99 
    p5 <- -99
    p6 <- -99
    p8 <- -99
    p12 <- -99
    p15 <- -99
    edf2 <- -99
    edf5 <- -99
    edf6 <- -99
    edf8 <- -99
    edf12 <- -99
    edf15 <- -99
    GAMsums.i <- data.frame(species.code,Radj,devX,no.vars,nonsig,p2,p5,p6,p8,p12,p15,edf2,edf5,edf6,edf8,edf12,edf15, stringsAsFactors=F)
    GAMsums <- rbind(GAMsums,GAMsums.i)
  }
   
  #GAMsums.i <- data.frame(species.code,Radj,devX,no.vars,p2,p5,p6,p8,p12,p15,edf2,edf5,edf6,edf8,edf12,edf15, stringsAsFactors=F)
  #GAMsums <- rbind(GAMsums,GAMsums.i)

  print(s)
  
} 

GAMsums <- GAMsums[GAMsums$species.code != "TEMPLATE",]
write.csv(GAMsums, "GamsumsNEW.csv", row.names=F)

#model.selX <- model.selX[model.selX$species != "TEMPLATE",]
#write.csv(model.selX, "10AICB.csv", row.names=F)

#model.sel2 <- model.sel[model.sel$dAIC <=2 & model.sel$insig == 0,]
#model.sel2 <- model.sel2[model.sel2 == min[model.sel2 ]]

#model.sel <- model.sel[model.sel$species != "TEMPLATE",]
#write.csv(model.sel, "vein9models.csv", row.names=F)


#    ####################### VISUALIZATIONS
#    
#    
#    species.code <- "VEIN9" ## extract the USDA species code for the species of the iteration
#    species <- edd[edd$usda == species.code,]## subset to species of the iteration
#    species <- data.frame(cbind(log(species$abundance), species$latitude, species$longitude)) ## select only variables needed to save computation time
#    
#    colnames(species) <- c("abundance", "latitude", "longitude") ## clean up column names
#    coordinates(species) <- c(3,2) ## specify lon/lat
#    proj4string(species) <- CRS("+init=epsg:4326") ## assume WGS84/NAD83 for EDDMapS points
#    
#    species <- spTransform(species, proj4string(bio)) ## project point data to proj4string of climate data
#    
#    ext <- extract(bio, species) ## extract climate values to points
#    species <- cbind(as.data.frame(species), ext) ## append extracted climate data to point data
#    species[,4:7] <- species[,4:7]/10 ## convert temp variables back to Celcius (from Celsius*10)
#    
#    
#    M <- gam(abundance ~ s(bio_2, fx=F, k=5) +
#               #s(bio_6, fx=F, k=5) +
#               s(bio_8, fx=F, k=5) +
#               s(bio_12, fx=F, k=5), method ="REML", data = species)
#    
#    par(mfrow=c(1,3))
#    plot(M)
#    
#    summary(M)
#    anova(M)
#    hist(M$residuals)
#    #bio_2 + bio_5 + bio_8 + bio_12
#    
#    edd2 <- edd[edd$usda %in% sp.list,]
#    hist(log(edd$abundance))
#    
#    
#    P <- predict.gam(M, bio.sp)
#    P <- melt(P)
#    P$value[P$value <= -5] <- NA
#    
#    #P$value[P$value < -5] <- NA
#    predict.raster <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=P$value)
#    exp(-5)
#    library(rgdal)
#    states<- readOGR(dsn="states", layer= "US_states")
#    dev.off()
#    plot(states, col="grey")
#    plot(predict.raster, add=T)
#    plot(species, col="red", pch=19, add=T)
#    
#    
#    par(mfrow=c(2,2))
#    M <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp")  + s(bio_12, fx=F, k=5, bs="tp")+ s(bio_5, fx=F, k=5, bs="tp"), data=species)
#    plot(M)
#    M2 <- gam(abundance ~s(bio_2, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp"), data=species)
#    plot(M2)
#    


quad <- c(1,2,0,1,1,2,0,0,3,3,0,3,1,1,2,0,1,2,2,2)
gam <- c(2,1,1,1,0,2,1,1,1,2,0,2,1,1,0,1,1,1,2,0)

quad <- c(0.166666667,0.4,0,0.25,0.166666667,0.5,0,0,0.6,0.6,0,0.6,0.2,0.25,0.5,0,0.333333333,0.4,0.4,0.4)
gam <- c(1,0.333333333,0.25,0.2,0,0.5,1,1,0.333333333,0.666666667,0,0.666666667,0.333333333,0.25,0,1,0.5,1,0.4,0)
no.vars.gam <- c(2,3,4,5,1,4,1,1,3,3,2,3,3,4,1,1,2,1,5,2)
no.vars.quad <- c(6,5,4,4,6,4,1,1,5,5,3,5,5,4,4,1,3,5,5,5)

cup <- c(0.5,0.4,0,0.25,0,0.25,0,0,0.2,0.4,0.666666667,0.2,0.4,0.25,0.5,2,0,0.2,0.2,0.2)

 wonkiness<- c(3,3,2,3,1,3,1,3,2,3,2,2,3,2,2,3,1,2,3,2)
plot(wonkiness,gam)

cor.test(wonkiness,gam,method="spearman")
plot(wonkiness,cup)
