#### Ordinal SDMs
###  created: 11/10/2017
##   last modified: 12/19/2017
#    Mitch O'Neill
# "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# library(mgcv)     ## for GAMS

library(raster)   ## for Raster climate data
library(rgeos)    ## for extractig climate data to points
library(reshape2) ## for reformatting model predictions
library(combinat) ## for generating pairs of covariates for assessment of collinearity 
library(glue)     ## for collapsing lists of variables into equations
#library(modEvA)
library(stringr)
library(MASS)

## list candidate terms (all six variables + quadratics)
var.list <- c("bio_2", "poly(bio_2, 2)",
              "bio_5", "poly(bio_5, 2)",
              "bio_6", "poly(bio_6, 2)",
              "bio_8", "poly(bio_8, 2)",
              "bio_12", "poly(bio_12, 2)",
              "bio_15", "poly(bio_15, 2)")

var.sets <- c("variable sets") ## initialize list of variable sets

## loop through all possible numbers of variables (half #vars to avoid duplicating 1st order terms)
for(i in 1:(length(var.list)/2)){
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
       !(grepl("bio_15 ", combo.j) & grepl("poly(bio_15", combo.j, fixed=T)) ){
      var.sets <- c(var.sets,combo.j)} ## add it to the list of variable sets
  }
  
  print(i)
}

formulae <- paste("abundance", var.sets[2:length(var.sets)], sep = " ~ ")
## add response part of formulaes, excluding the initial row

rm(var.sets, combo.j, combos.i, i, j, var.list) ## garbage cleaning

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/")
## read in EDDMapS dataset, thinned to 5km climate grid cells (1 point per species per cell)
edd <- read.table("edd_final_12_18_17.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd$abundance <- ceiling(edd$best)
edd$max<-NULL;edd$med<-NULL;edd$best<-NULL;edd$no.pts<-NULL;edd$no.best<-NULL;head(edd)
coordinates(edd) <- c(4,3)
proj4string(edd) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## from transformations during data prep


## read in Worldclim climate rasters clipped to L48; 6 variables manually selected by assessing correlations
bio <- stack("C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_2.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_5.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_6.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_8.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_12.asc",
             "C:/Users/mwone/Documents/geodata/clipped_climate_data/current/bio_15.asc")
## assign proj4string to that of one of the raw files before clipping
proj4string(bio) <- proj4string(raster("C:/Users/mwone/Documents/geodata/climate_data/current/bio_1"))
bioD <- as.data.frame(bio)
bioD[,1:4] <- bioD[,1:4]/10


edd <- spTransform(edd, proj4string(bio)) ## project point data to proj4string of climate data
ext <- extract(bio, edd) ## extract climate values to points
edd <- cbind(as.data.frame(edd), ext) ## append extracted climate data to point data
edd[,(ncol(edd)-5):(ncol(edd)-2)] <- edd[,(ncol(edd)-5):(ncol(edd)-2)]/10
head(edd)

species.code <- "TEMPLATE" ## initialize master file to contain summaries of GAMs (final model for each species)
aic <- -99 ## adjusted R-squared of models
no.pts <- -99
no.terms <- -99
no.vars <- -99
#p05 <- -99
#p10 <- -99
formu <- "TEMPLATE"
c2    <- -999
c2.2  <- -999
c5    <- -999
c5.2  <- -999
c6    <- -999
c6.2  <- -999
c8    <- -999
c8.2  <- -999
c12   <- -999
c12.2 <- -999
c15   <- -999
c15.2 <- -999

ordsums <- data.frame(species.code,aic,no.pts,no.terms,no.vars,formu,c2,c2.2,c5,c5.2,c6,c6.2,c8,c8.2,c12,c12.2,c15,c15.2, stringsAsFactors = F)
#ordsums <- data.frame(species.code,aic,no.pts,no.terms,no.vars,p05s,p10s,p05i,p10i,formu,stringsAsFactors = F)
sp.list <- unique(edd$species)

for(s in 11:length(sp.list)){ ## Loop through species list
  
  species.code <- sp.list[s] ## extract the USDA species code for the species of the iteration
  
  if(species.code != "ACJA" & species.code !="HIAU") {
    species <- edd[edd$species == species.code,]## subset to species of the iteration
    #species <- data.frame(cbind(species$abundance, species$latitude, species$longitude)) ## select only variables needed to save computation time
    species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3,4))
    
    #par(mfrow=c(1,1)) ## clean up plot space from past iterations
    #jpeg(paste("histograms", paste(species.code, "jpg", sep="."), sep="/" )) ## output file name for histogram of abun
    #hist(as.numeric(species$abundance), main=species.code, xlab="abundance rank")# histogram of abundance with 20 breaks
    #dev.off() ## export plot to filepath specified above
    
    ########### ASSESS CORRELATIONS ############
    rhos <- as.data.frame(cor(species[,(ncol(edd)-5):(ncol(edd))], use="pairwise", method="spearman")) 
    ## construct data.frame of matrix of pairwise spearman's rho values for climatic variables
    var.pairs <- combn(c(1:6), 2, simplify=T) ## All possible pairs selected from 6 variables
    
    ##initialize master data frameof pairwise correlations with template
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
    corr.pairs <- corr.pairs[corr.pairs$rho >= .6, ] ## select all pairwise correlations where collinearity is likely 
    
    formulae.s <- data.frame(formulae, stringsAsFactors = F)
    formulae.s$collinear <- 0
    
    if(length(corr.pairs$rho) > 0){ ## if there were any pairs where rho >= 0.8
      for (i in 1:length(corr.pairs$rho)){ ## then loop through those pairs
        #formulae.s <- formulae[1:12] #initialize
        for (j in 2:length(formulae)){ ## for each pair, loop through all of the possible models
          ## and if the two highly correlated variables appear together in a model,
          if(grepl(corr.pairs$var1[i], formulae[j], fixed=T) & grepl(corr.pairs$var2[i], formulae[j], fixed=T)) { 
            formulae.s$collinear[j] <- 1
          }
        } 
      }
    } 
    
    
    formulae.s <- formulae.s$formulae[formulae.s$collinear == 0]
    
    rm(rhos, var1, var2, corr.pairs, var.pairs) ## garbage cleaning
    
    #### MODEL SELECTION #####
    ## initialize template masterfile to contain summaries of each model for this species
    i <- -1 ## keep track of place in formula list
    edf <- -1 ## number of variables
    aic <- -1 ## Akaike's info criterion
    #nonsig <- -1
    #bios <- "TEMPLATE" ## list of climate covariates|#nonsig <- -1 ## 1: there is a non-sig variable, 2: all variables are signif|#conc.up <- -1
    
    model.sel <- data.frame(i, edf, aic, stringsAsFactors = F)    
    #model.sel <- data.frame(i, n, aic, bios, nonsig, conc.up, stringsAsFactors = F)
    ## concatenate into master data.frame
    
    for (i in 1:length(formulae.s)){ ##loop through all formulae in the formula list
      
      M <- polr(formula(formulae.s[i]), data = species, Hess=T)
      # M <- lrm(formula(formulae.s[i]), data = species)|# M <- clm(formula(formulae.s[i]), data=species)|## construct model using formula of the iteration|## use data for species of this iteration
      
      edf <- M$edf ## REPLACE n
      aic <- AIC(M) ## extract aic score
      
      #Msum <- coef(summary(M))
      #p <- pnorm(abs(Msum[, "t value"]), lower.tail = FALSE) * 2
      #Msum <- data.frame(cbind(Msum, "p value" = p))
      #MsumS <- Msum[1:(nrow(Msum)-3),]
      #nonsig <- length(MsumS$p.value[MsumS$p.value > 0.05])
      
      model.sel.i <- data.frame(i, edf, aic, stringsAsFactors = F)
      model.sel <- rbind(model.sel, model.sel.i) ## append data from iteration to the master data frame
      #print(i) ## keep track of place in loop
    }

    #################### PROB ####################
    #   #     model.sel <- model.sel[model.sel$i != -1,]# & model.sel$conc.up == 0,] ## remove template row
    #   #     #head(model.sel,10)
    #   # formulae[c(212,224,228,292,304,372,384,460,468,540,548,620,628,644)]
    #   # probs <- data.frame(formulae.s[c(12,20,22,28,36,44,56,60,72,76,81,82,83,84,88,92,104,108,120,124,137,138,139,140,145,146,147,148)])
    #   probs <- data.frame(formulae.s[c(56,124,128,144,156,160,164,168,196,200,228,232,236,242,244,248,260)])
    #   colnames(probs) <-c("formula")
    #   probs
    #   probs$bio_15 <- 0
    #   for (i in 1:28){
    #     if(grepl("poly(bio_15", probs$formula[i], fixed=T)){
    #       probs$bio_15[i] <- 1
    #     }
    #   }
    #   summary(probs$bio_15)
    #   probs<- probs[probs$bio_15 == 0,]
    #   plot((species$bio_6^2 + species$bio_6)~species$abundance)
    # } ##???
    
    
    ###################  AIC comparison  #####################
    model.sel <- model.sel[model.sel$i != -1,]
    model.sel$dAIC <- model.sel$aic - min(model.sel$aic) ## calculate delta AIC
    #model.sel <- model.sel[order(model.sel$dAIC),]
    model.sel2 <- model.sel[model.sel$dAIC <= 2,] ## select the model with the best score, 
    ## as well as models with scores that are not significantly worse
    
    #if(min(model.sel2$nonsig == 0)){
    #  model.sel2 <- model.sel2[model.sel2$nonsig ==0]
    #}
    rm(model.sel, edf, i)
    
    model.sel2 <- model.sel2[model.sel2$edf == min(model.sel2$edf),] 
    ## of the models with the best scores (<= min + 2), select the simplest models
    
    model.sel2 <- model.sel2[model.sel2$dAIC == min(model.sel2$dAIC),]
    ## select the model with the best AIC score if there are multiple models
    ## with similar scores and the same complexity (variable number)
    
    M <- polr(formula(formulae.s[model.sel2$i]), data = species, Hess = T)
    #M <- polr(formula(formulae.s[3]), data = species, Hess = T)
    ## reconstruct the best model, using the formula number (i) to re-access it|#library(rms)|#M2 <- clm(formula(formulae.s[model.sel2$i]), data=species)|#M2$coefficients|#M2$coefficients|#AIC(M2)|#AIC(M)|#species$abundance <- ordered(as.numeric(as.character(species$abundance)), levels = c("1","2","3"))|#edf <- M$edf ## REPLACE n
    
    ## extract summary of the model to access model-fit values
    Msum <- coef(summary(M))
    p <- pnorm(abs(Msum[, "t value"]), lower.tail = FALSE) * 2
    Msum <- data.frame(cbind(Msum, "p value" = p))
    #plot(Msum)
    #MsumS <- Msum[1:(nrow(Msum)-3),]
    #MsumI <- Msum[(nrow(Msum)-2):nrow(Msum),]
    
    
    aic <- AIC(M)          ## extract adjusted R-square value of model #adjDsq <-  Dsquared(M,adjust=T) ## extract deviance explained by model
    no.pts <- nrow(species)
    no.terms <- NROW(Msum) ## extract the number of terms
    no.vars <- no.terms - 3 - str_count(formulae.s[model.sel2$i], "poly")
    formu <- formulae.s[model.sel2$i]
    #p05s<- length(MsumS$p.value[MsumS$p.value >= 0.05])
    #p10s<- length(MsumS$p.value[MsumS$p.value > 0.1])
    #p05i<- length(MsumI$p.value[MsumI$p.value >= 0.05])
    #p10i<- length(MsumI$p.value[MsumI$p.value > 0.1])
    
    Msum$name <- as.character(row.names(Msum))
    
    c2    <- Msum$Value[Msum$name == "poly(bio_2, 2)1" | Msum$name == "bio_2"]
    c2.2  <- Msum$Value[Msum$name == "poly(bio_2, 2)2"]
    c5    <- Msum$Value[Msum$name == "poly(bio_5, 2)1" | Msum$name == "bio_5"]
    c5.2  <- Msum$Value[Msum$name == "poly(bio_5, 2)2"]
    c6    <- Msum$Value[Msum$name == "poly(bio_6, 2)1" | Msum$name == "bio_6"]
    c6.2  <- Msum$Value[Msum$name == "poly(bio_6, 2)2"]
    c8    <- Msum$Value[Msum$name == "poly(bio_8, 2)1" | Msum$name == "bio_8"]
    c8.2  <- Msum$Value[Msum$name == "poly(bio_8, 2)2"]
    c12   <- Msum$Value[Msum$name == "poly(bio_12, 2)1" | Msum$name == "bio_12"]
    c12.2 <- Msum$Value[Msum$name == "poly(bio_12, 2)2"]
    c15   <- Msum$Value[Msum$name == "poly(bio_15, 2)1" | Msum$name == "bio_15"]
    c15.2 <- Msum$Value[Msum$name == "poly(bio_15, 2)2"]
    if(length(c2) == 0){c2 <- 0}
    if(length(c2.2) == 0){c2.2 <- 0}
    if(length(c5) == 0){c5 <- 0}
    if(length(c5.2) == 0){c5.2 <- 0}
    if(length(c6) == 0){c6 <- 0}
    if(length(c6.2) == 0){c6.2 <- 0}
    if(length(c8) == 0){c8 <- 0}
    if(length(c8.2) == 0){c8.2 <- 0}
    if(length(c12) == 0){c12 <- 0}
    if(length(c12.2) == 0){c12.2 <- 0}
    if(length(c15) == 0){c15 <- 0}
    if(length(c15.2) == 0){c15.2 <- 0}    
    
    #c2,c2.2,c5,c5.2,c6,c6.2,c8,c8.2,c12,c12.2,c15,c15.2
    ordsums.i <- data.frame(species.code,aic,no.pts,no.terms,no.vars,formu,c2,c2.2,c5,c5.2,c6,c6.2,c8,c8.2,c12,c12.2,c15,c15.2, stringsAsFactors = F)
    ordsums <- rbind(ordsums,ordsums.i) 
    print(s)
  }
}

ordsums <- ordsums[ordsums$species.code != "TEMPLATE",]
## 11:55-12:45 for 78
#hist(ordsums$p10s/(ordsums$no.terms-3),breaks=50)
#plot(ordsums$p10s/(ordsums$no.terms-3)~ordsums$no.pts)

write.csv(ordsums,"ordsumsNEW.csv",row.names=F)



for (i in 1:length(ordsums$species.code)){
  linear <- data.frame(c(ordsums$c2[i], ordsums$c5[i], ordsums$c6[i], ordsums$c8[i], ordsums$c12[i], ordsums$c15[i]))
  quad <- data.frame(c(ordsums$c2.2[i], ordsums$c5.2[i], ordsums$c6.2[i], ordsums$c8.2[i], ordsums$c12.2[i], ordsums$c15.2[i]))
  colnames(linear) <- c("coef"); colnames(quad) <- c("coef")
  linear <- data.frame(linear$coef[linear$coef != 0]); quad <- data.frame(quad$coef[quad$coef != 0])
  colnames(linear) <- c("coef"); colnames(quad) <- c("coef")
  
  
  ordsums$no.lin[i]   <- length(linear$coef) - length(quad$coef)
  ordsums$pos.quad[i] <- length(quad$coef[quad$coef > 0])
  ordsums$neg.quad[i] <- length(quad$coef[quad$coef < 0])
    
  print(i)
}

ordsums <- ordsums[2:length(ordsums$species.code),]
sum(ordsums$no.lin)  ## 142
sum(ordsums$pos.quad)## 49
sum(ordsums$neg.quad)## 91



species <- edd[edd$species == "PHRE",]## subset to species of the iteration
#species <- data.frame(cbind(species$abundance, species$latitude, species$longitude)) ## select only variables needed to save computation time
species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3,4))

M <- polr(abundance~bio_6, data=species,Hess=T)
Msum <- data.frame(coef(summary(M)))
bio_6 <- data.frame(seq(from=min(bioD$bio_6, na.rm=T),to=max(bioD$bio_6,na.rm=T),by=0.1))
colnames(bio_6) <- "bio_6"
xbeta <- (Msum$Value[1]*bio_6$bio_6)

logistic_cdf <- function(x) {
  return( 1/(1+exp(-x) ) )
}

p1 <- logistic_cdf( Msum$Value[nrow(Msum)-2] - xbeta )
p2 <- logistic_cdf( Msum$Value[nrow(Msum)-1] - xbeta ) - logistic_cdf( Msum$Value[nrow(Msum)-2] - xbeta )
p3 <- logistic_cdf( Msum$Value[nrow(Msum)] - xbeta ) - logistic_cdf( Msum$Value[nrow(Msum)-1] - xbeta )
p4 <- 1 - logistic_cdf( Msum$Value[nrow(Msum)] - xbeta )

plot(bio_6$bio_6, p1, lwd=4,type='l', ylab='Prob', main="Linear",ylim=c(0,max(p1,p2,p3,p4)), col='blue', xlab='min temp of coldest month')
lines(bio_6$bio_6, p2,lwd=4, col='turquoise')
lines(bio_6$bio_6, p3,lwd=4, col='green')
lines(bio_6$bio_6, p4,lwd=4, col='orange')



############
species <- edd[edd$species == "AEPO",]## subset to species of the iteration
#species <- data.frame(cbind(species$abundance, species$latitude, species$longitude)) ## select only variables needed to save computation time
species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3,4))

M <- polr(abundance ~ bio_12 + I(bio_12^2), data=species,Hess=T)
#Msum <- data.frame(coef(summary(M)))
bio_12 <- data.frame(seq(from=min(bioD$bio_12, na.rm=T),to=max(bioD$bio_12,na.rm=T),by=0.1))
colnames(bio_12) <- "bio_12"
xbeta <- (2.656845e-05*(bio_12$bio_12^2))+ (-5.629731e-02*bio_12$bio_12)

logistic_cdf <- function(x) {
  return( 1/(1+exp(-x) ) )
}

p1 <- logistic_cdf(-31.18570 - xbeta )
p2 <- logistic_cdf(-29.74032 - xbeta ) - logistic_cdf(-31.18570- xbeta )
p3 <- logistic_cdf(-28.27218  - xbeta ) - logistic_cdf(-29.74032 - xbeta )
p4 <- 1 - logistic_cdf(-28.27218  - xbeta )

plot(bio_12$bio_12, p1, lwd=4,type='l', ylab='Prob', main="Quadratic (+)",ylim=c(0,max(p1,p2,p3,p4)), col='blue', xlab='annual precipitation')
lines(bio_12$bio_12, p2,lwd=4, col='turquoise')
lines(bio_12$bio_12, p3,lwd=4, col='green')
lines(bio_12$bio_12, p4,lwd=4, col='orange')


#head(species)
species <- edd[edd$species == "RUPH",]## subset to species of the iteration
#species <- data.frame(cbind(species$abundance, species$latitude, species$longitude)) ## select only variables needed to save computation time
species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3,4))

M <- polr(abundance ~ bio_8 + I(bio_8^2), data=species,Hess=T)
#Msum <- data.frame(coef(summary(M)))
bio_8 <- data.frame(seq(from=min(bioD$bio_8, na.rm=T),to=max(bioD$bio_8,na.rm=T),by=0.1))
colnames(bio_8) <- "bio_8"
xbeta <- (-0.0140012*(bio_8$bio_8^2))+ (.3962326*bio_8$bio_8)

logistic_cdf <- function(x) {
  return( 1/(1+exp(-x) ) )
}

p1 <- logistic_cdf(-.4736391 - xbeta )
p2 <- logistic_cdf(2.1762984 - xbeta ) - logistic_cdf(-.4736391 - xbeta )
p3 <- logistic_cdf(3.3110255   - xbeta ) - logistic_cdf(2.1762984 - xbeta )
p4 <- 1 - logistic_cdf(3.3110255   - xbeta )

plot(bio_8$bio_8, p1, lwd=4,type='l', ylab='Prob', main="Quadratic (-)",ylim=c(0,max(p1,p2,p3,p4)), col='blue', xlab='mean temp of wettest month')
lines(bio_8$bio_8, p2,lwd=4, col='turquoise')
lines(bio_8$bio_8, p3,lwd=4, col='green')
lines(bio_8$bio_8, p4,lwd=4, col='orange')


species <- edd[edd$species == "ECVU",]
coordinates(species) <- c(5,6)
head(species)
proj4string(species) <- proj4string(bio)

M <- polr(abundance ~ poly(bio_2, 2) + bio_5 + poly(bio_8, 2) + bio_12, data=species, Hess= T)

plot(bio$bio_5)
plot(species[species$abundance == 4,], col= "black", pch=20, add=T)
plot(species[species$abundance == 3,], col= "dark gray", pch=20, add=T)
plot(species[species$abundance == 2,], col= "grey", pch=20, add=T)
plot(species[species$abundance == 1,], col= "light grey", pch=20, add=T)

predict

#edd <- read.table("edd_final_12_18_17.csv", header = T, sep = ",", quote= "\"", 
#                  comment.char= "", stringsAsFactors = F, strip.white = T)

traits <- read.table("file:///C:/Users/mwone/Google Drive/NSF_GSS_shared/Hotspots_and_Abundance/Final/master list.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
head(traits)
head(ordsums)

ordsumsM <- merge(ordsums, traits, by.x="species.code", by.y="NewCode", all=F)
head(ordsumsM)
summary(factor(ordsumsM$duration_final))
summary(factor(ordsumsM$habit))


head(species)
p <- predict(M, bioD)
pm <- melt(p)

predict.raster <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=pm$value)
png("C:/Users/mwone/Desktop/myplot.png", width=1000, height=617, bg = "transparent")
plot(predict.raster)
dev.off()

ecvu <- raster("C:/Users/mwone/Google Drive/NSF_GSS_shared/Watch_Lists/range_maps/binary_asciis/sppRichCtyBias_current/sppRichCtyBias_current/sppRichCtyBias_currentECVU_SppRichBias_FivePercMTP.asc")
plot(ecvu)
summary(ecvu$sppRichCtyBias_currentECVU_SppRichBias_FivePercMTP)
proj4string(ecvu) <- proj4string(species)

ecvu2 <- as.data.frame(ecvu)
ecvu2$sppRichCtyBias_currentECVU_SppRichBias_FivePercMTP[ecvu2$sppRichCtyBias_currentECVU_SppRichBias_FivePercMTP==1] <- NA
ecvu2 <- raster(nrows=nrow(ecvu), ncol=ncol(ecvu), ext=extent(ecvu), crs=proj4string(ecvu), vals=ecvu2$sppRichCtyBias_currentECVU_SppRichBias_FivePercMTP)
plot(ecvu2, col="black")

png("C:/Users/mwone/Desktop/myplot.png", width=1000, height=617, bg = "transparent")
plot(ecvu2,bg="transparent",col="black")
dev.off()
#dev.copy (png,'C:/Users/mwone/Desktop/myplot.png', bg = 'transparent')

###########################################################################################


ordsums <- read.table("ordsums.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)
ordsums2 <- read.table("ordsums2.csv", header = T, sep = ",", quote= "\"", 
                       comment.char= "", stringsAsFactors = F, strip.white = T)
summary(ordsums$p05s)
summary(ordsums2$p05s)
hist(ordsums$no.pts[ordsums$p05s/(ordsums$no.terms-3) == 1 & (ordsums$no.terms-3)==1])


par(mfrow=c(2,2))
hist(ordsums$p05s, main="Variables where p>0.05", xlab="#variables/model")
hist(ordsums$p10s, main="Variables where p>0.10", xlab="#variables/model")
hist(ordsums$p05s/(ordsums$no.terms-3),breaks=20, main="Variables where p>0.05", xlab="proportion of variables in model")
hist(ordsums$p10s/(ordsums$no.terms-3),breaks=20, main="Variables where p>0.10", xlab="proportion of variables in model")





species.code <- "LISI" ## extract the USDA species code for the species of the iteration
#if(species.code != "ACJA" & species.code !="HIAU") {
species <- edd[edd$species == species.code,]## subset to species of the iteration
#species <- data.frame(cbind(species$abundance, species$latitude, species$longitude)) ## select only variables needed to save computation time
species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3,4))
M <- polr(abundance ~ poly(bio_8,2), data = species, Hess=T)
M2 <- polr(abundance ~ bio_8 + I(bio_8^2), data = species, Hess=T)
Msum <- summary(M);Msum <- data.frame(Msum$coefficients);Msum
#Msum <- summary(M2);Msum <- data.frame(Msum$coefficients);Msum

coefs <- c(Msum$Value[1],Msum$Value[2])
M$fitted.values
plot(M$deviance)
M$xlevels



bio_8 <- seq(from=min(bioD$bio_8, na.rm=T),to=max(bioD$bio_8,na.rm=T),by=0.1)
bio8 <- data.frame(poly(bio_8,2))
bio8 <- cbind(bio8,bio_8)
colnames(bio8) <- c("one","two","bio_8")
#head(bio8)
#summary(bio8)
#plot(bio8$X1~bio8$bio_8)
#plot(bio8$X2~bio8$bio_8)
#abline(v=-0.0427)
#borders <- c(-2.65,5.95,9.35,23.25,26.65)
#borders <- data.frame(poly(borders,2))
#head(borders)
#bio8

xbeta <- (Msum$Value[2]*(bio_8^2)) + (Msum$Value[1]*bio_8)

logistic_cdf <- function(x) {
  return( 1/(1+exp(-x) ) )
}
p1 <- logistic_cdf( Msum$Value[nrow(Msum)-2] - xbeta )
p2 <- logistic_cdf( Msum$Value[nrow(Msum)-1] - xbeta ) - logistic_cdf( Msum$Value[nrow(Msum)-2] - xbeta )
p3 <- logistic_cdf( Msum$Value[nrow(Msum)] - xbeta ) - logistic_cdf( Msum$Value[nrow(Msum)-1] - xbeta )
p4 <- 1 - logistic_cdf( Msum$Value[nrow(Msum)] - xbeta )

plot(bio8$bio_8, p1, type='l', ylab='Prob', ylim=c(0,max(p1,p2,p3,p4)), col='blue')
lines(bio8$bio_8, p2, col='green')
lines(bio8$bio_8, p3, col='orange')
lines(bio8$bio_8, p4, col='red')

#abline(v=borders$X2[1], lwd=5, col="light grey")
#abline(v=borders$X2[2], lwd=5, col="light grey")
#abline(v=borders$X2[3], lwd=5, col="light grey")
#abline(v=borders$X2[4], lwd=5, col="light grey")
#abline(v=borders$X2[5], lwd=5, col="light grey")

#abline(v=-2.65, lwd=5, col="light grey")
#abline(v=5.95, lwd=5, col="light grey")
#abline(v=9.35, lwd=5, col="light grey")
#abline(v=23.25, lwd=5, col="light grey")
#abline(v=26.65, lwd=5, col="light grey")
abline(v=min(species$bio_8), lwd=1, col="light grey")
abline(v=max(species$bio_8), lwd=1, col="light grey")

bio_8 <- melt(unique(bioD$bio_8))
colnames(bio_8)<-"bio_8"
abun <- cbind(bio_8,melt(predict(M, bio_8)))
abun <- abun[order(abun$bio_8),]

abun$value[93] ## 1
abun$value[94] ## 2
abun$value[179] ## 2
abun$value[180] ## 3
abun$value[213] ## 3
abun$value[214] ## 4
abun$value[352] ## 4
abun$value[353] ## 3
abun$value[386] ## 3
abun$value[387] ## 2

abun$bio_8[93] ## 1
abun$bio_8[94] ## 2
abun$bio_8[179] ## 2
abun$bio_8[180] ## 3
abun$bio_8[213] ## 3
abun$bio_8[214] ## 4
abun$bio_8[352] ## 4
abun$bio_8[353] ## 3
abun$bio_8[386] ## 3
abun$bio_8[387] ## 2



## -2.65; 5.95; 9.35, 23.25,26.65


predict.raster <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=pm$value)
plot(predict.raster)
head(species)
coordinates(species) <- c(5,6)
proj4string(species) <- proj4string(predict.raster)

plot(predict.raster)

#####################################################################################
## 13 mins for 20 spp
#13*(116/20)/60
Msum
bio8 <- seq(from=min(bioD$bio_8, na.rm=T),to=max(bioD$bio_8,na.rm=T),by=1)
xbeta <- bio8*Msum$Value[3]

logistic_cdf <- function(x) {
  return( 1/(1+exp(-x) ) )
}

p1 <- logistic_cdf( Msum$Value[nrow(Msum)-2] - xbeta )
p2 <- logistic_cdf( Msum$Value[nrow(Msum)-1] - xbeta ) - logistic_cdf( Msum$Value[nrow(Msum)-2] - xbeta )
p3 <- logistic_cdf( Msum$Value[nrow(Msum)] - xbeta ) - logistic_cdf( Msum$Value[nrow(Msum)-1] - xbeta )
p4 <- 1 - logistic_cdf( Msum$Value[nrow(Msum)] - xbeta )

plot(bio8, p1, type='l', ylab='Prob', ylim=c(0,max(p1,p2,p3,p4)), col='blue')
lines(bio8, p2, col='green')
lines(bio8, p3, col='orange')
lines(bio8, p4, col='red')
abline(v=min(species$bio_8), lwd=5, col="light grey")
abline(v=max(species$bio_8), lwd=5, col="light grey")
#legend("topleft", lty=1, col=c("blue", "green", "orange", "red" ), 
#       legend=c("rank 1", "rank 2", "rank 3", "rank 4" ))


## Proceed with caution 
#summary(species$bio_5) ###########################################
### extract significance (p-vals) and effective degrees of freedom for each variables
p2 <- -1  ## start by defaulting each to -1, indicating that the variable is not in the model
p2.2 <- -1
p5 <- -1
p5.2 <- -1
p6 <- -1
p6.2 <- -1
p8 <- -1
p8.2 <- -1
p12 <- -1
p12.2 <- -1
p15 <- -1
p15.2 <- -1

#### RECONFIGURE    

#    bios <- as.character(M$formula[3])
#    if (grepl("bio_2", bios) & !grepl("poly(bio_2,", bios, fixed=T)) {
#      p2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_2"]
#    }
#    if (grepl("poly(bio_2,", bios, fixed=T) ) {
#      p2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_2, 2)1"] 
#      p2.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_2, 2)2"]
#    }
#    if (grepl("bio_5", bios) & !grepl("poly(bio_5,", bios, fixed=T)) {
#      p5 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_5"]
#    }
#    if (grepl("poly(bio_5,", bios, fixed=T) ) {
#      p5 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_5, 2)1"] 
#      p5.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_5, 2)2"]
#    }
#    if (grepl("bio_6", bios) & !grepl("poly(bio_6,", bios, fixed=T)) {
#      p6 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_6"]
#    }
#    if (grepl("poly(bio_6,", bios, fixed=T) ) {
#      p6 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_6, 2)1"] 
#      p6.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_6, 2)2"]
#    }
#    if (grepl("bio_8", bios) & !grepl("poly(bio_8,", bios, fixed=T)) {
#      p8 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_8"]
#    }
#    if (grepl("poly(bio_8,", bios, fixed=T) ) {
#      p8 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_8, 2)1"] 
#      p8.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_8, 2)2"]
#    }
#    if (grepl("bio_12", bios) & !grepl("poly(bio_12,", bios, fixed=T)) {
#      p12 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_12"]
#    }
#    if (grepl("poly(bio_12,", bios, fixed=T) ) {
#      p12 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_12, 2)1"] 
#      p12.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_12, 2)2"]
#    }
#    if (grepl("bio_15", bios) & !grepl("poly(bio_15,", bios, fixed=T)) {
#      p15 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_15"]
#    }
#    if (grepl("poly(bio_15,", bios, fixed=T) ) {
#      p15 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_15, 2)1"] 
#      p15.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_15, 2)2"]
#    }

ordSums.i <- data.frame(species.code,aic,adjDsq,no.pts,no.terms,no.vars,lat.ext,lon.ext,habit,duration,bios,p2,p2.2,p5,p5.2,p6,p6.2,p8,p8.2,p12,p12.2,p15,p15.2, stringsAsFactors = F)
ordSums <- rbind(ordSums,ordSums.i)

########### RESPONSE CURVES ##########
#    jpeg(paste("quad_pics2", paste(species.code, "jpg", sep="."), sep="/"))
#    par(pch="." )
#    
#    ylim <- c(quantile(species$abundance,.05)-abs(quantile(species$abundance,.05)), quantile(species$abundance,.95) + abs(quantile(species$abundance,.95)))
#    
#    if (no.vars > 4){
#      par(mfrow=c(2,3)) 
#    }
#    if (no.vars == 4){
#      par(mfrow=c(2,2)) 
#    }
#    if (no.vars == 3){
#      par(mfrow=c(3,1)) 
#    }
#    if (no.vars == 2){
#      par(mfrow=c(2,1)) 
#    }
#    
#    
#    multiplier <- c(0:100)/100
#    
#    if(grepl("bio_2", bios)) {
#      bio_2 <- min(species$bio_2) + multiplier*(max(species$bio_2)- min(species$bio_2))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      
#      input <- data.frame(bio_2) ## change list to data.frame
#      input$bio_5 <- mean(species$bio_5) ## set all other variables to constant value
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      plot(input$value~input$bio_2, main=round(c(p2,p2.2),digits=5), ylim=ylim)
#      lines(input$value~input$bio_2, lwd=3)
#      
#    }
#    
#    if(grepl("bio_5", bios)) {
#      bio_5 <- min(species$bio_5) + multiplier*(max(species$bio_5)- min(species$bio_5))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      
#      input <- data.frame(bio_5) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      plot(input$value~input$bio_5, main=round(c(p5,p5.2),digits=5), ylim=ylim)
#      lines(input$value~input$bio_5, lwd=3)
#    }
#    
#    if(grepl("bio_6", bios)) {
#      bio_6 <- min(species$bio_6) + multiplier*(max(species$bio_6)- min(species$bio_6))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      
#      input <- data.frame(bio_6) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      plot(input$value~input$bio_6, main=round(c(p6,p6.2),digits=5), ylim=ylim)
#      lines(input$value~input$bio_6, lwd=3)
#    }
#    
#    if(grepl("bio_8", bios)) {
#      bio_8 <- min(species$bio_8) + multiplier*(max(species$bio_8) - min(species$bio_8))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      
#      input <- data.frame(bio_8) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      
#      plot(input$value~input$bio_8, main=round(c(p8,p8.2),digits=5), ylim=ylim)
#      lines(input$value~input$bio_8, lwd=3)
#    }
#    
#    if(grepl("bio_12", bios)) {
#      bio_12 <- min(species$bio_12) + multiplier*(max(species$bio_12)- min(species$bio_12))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      
#      input <- data.frame(bio_12) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_15 <- mean(species$bio_15)
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      plot(input$value~input$bio_12, main=round(c(p12,p12.2),digits=5), ylim=ylim)
#      lines(input$value~input$bio_12, lwd=3)
#    }
#    
#    if(grepl("bio_15", bios)) {
#      bio_15 <- min(species$bio_15) + multiplier*(max(species$bio_15)- min(species$bio_15))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      
#      input <- data.frame(bio_15) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_5 <- mean(species$bio_5)
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      plot(input$value~input$bio_15, main=round(c(p15,p15.2),digits=5), ylim=ylim)
#      lines(input$value~input$bio_15, lwd=3)
#    } 
#    dev.off()
#    

#      ###################EXTRAP##########
#    jpeg(paste("quad_pics3", paste(species.code, "jpg", sep="."), sep="/"))
#    par(pch=".")
#    bio.sp <- bioD
#    
#    if (no.vars > 4){
#      par(mfrow=c(2,3)) 
#    }
#    if (no.vars == 4){
#      par(mfrow=c(2,2)) 
#    }
#    if (no.vars == 3){
#      par(mfrow=c(3,1)) 
#    }
#    if (no.vars == 2){
#      par(mfrow=c(2,1)) 
#    }
#    
#    multiplier<- c(0:200)/200
#  
#    if(grepl("bio_2", bios)) {
#      bio_2 <- min(bioD$bio_2,na.rm=T) + multiplier*(max(bioD$bio_2, na.rm=T) - min(bioD$bio_2, na.rm=T))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      bio_2 <- c(min(species$bio_2) - 10, min(species$bio_2), max(species$bio_2), max(species$bio_2) + 10, bio_2)
#      
#      input <- data.frame(bio_2) ## change list to data.frame
#      input$bio_5 <- mean(species$bio_5) ## set all other variables to constant value
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      input$bio_2orig <- input$bio_2
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      
#      if (p[1,]>p[2,]) {
#        bio.sp$bio_2[bio.sp$bio_2 < min(species$bio_2)]  <- min(species$bio_2)
#        input$bio_2[input$bio_2 < min(species$bio_2)]  <- min(species$bio_2)
#      } ## TRUE: upward trend on left end
#      ## reassign values below min to simulate a flatline
#      
#      if (p[4,]>p[3,]) {
#        bio.sp$bio_2[bio.sp$bio_2 > max(species$bio_2)]  <- max(species$bio_2) 
#        input$bio_2[input$bio_2 > max(species$bio_2)]  <- max(species$bio_2)
#      }## TRUE: upward trend on right end
#      ## reassign values above max to simulate a flatline
#      
#      xmin <- min(species$bio_2)
#      xmax <- max(species$bio_2)
#      x5 <- quantile(bioD$bio_2, .05, na.rm=T)
#      x95 <- quantile(bioD$bio_2, .95, na.rm=T) 
#      
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      colnames(p) <- c("abundance")
#      input <- cbind(input, p)
#      
#      input <- input[5:NROW(input),]
#      input <- input[order(input$bio_2),]
#      plot(exp(input$abundance)~input$bio_2orig, xlab="bio_2", ylab="infested area")
#      polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#      polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#      lines(exp(input$abundance)~input$bio_2orig, lwd=3)
#    }
#    
#    if(grepl("bio_5", bios)) {
#      bio_5 <- min(bioD$bio_5,na.rm=T) + multiplier*(max(bioD$bio_5, na.rm=T) - min(bioD$bio_5, na.rm=T))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      bio_5 <- c(min(species$bio_5) - 10, min(species$bio_5), max(species$bio_5), max(species$bio_5) + 10, bio_5)
#      
#      input <- data.frame(bio_5) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      input$bio_5orig <- input$bio_5
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      
#      if (p[1,]>p[2,]) {
#        bio.sp$bio_5[bio.sp$bio_5 < min(species$bio_5)]  <- min(species$bio_5)
#        input$bio_5[input$bio_5 < min(species$bio_5)]  <- min(species$bio_5)
#      } ## TRUE: upward trend on left end
#      ## reassign values below min to simulate a flatline
#      
#      if (p[4,]>p[3,]) {
#        bio.sp$bio_5[bio.sp$bio_5 > max(species$bio_5)]  <- max(species$bio_5) 
#        input$bio_5[input$bio_5 > max(species$bio_5)]  <- max(species$bio_5)
#      }## TRUE: upward trend on right end
#      ## reassign values above max to simulate a flatline
#      
#      xmin <- min(species$bio_5)
#      xmax <- max(species$bio_5)
#      x5 <- quantile(bioD$bio_5, .05, na.rm=T)
#      x95 <- quantile(bioD$bio_5, .95, na.rm=T) 
#      
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      colnames(p) <- c("abundance")
#      input <- cbind(input, p)
#      
#      input <- input[5:NROW(input),]
#      input <- input[order(input$bio_5),]
#      plot(exp(input$abundance)~input$bio_5orig, xlab="bio_5", ylab="infested area")
#      polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#      polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#      lines(exp(input$abundance)~input$bio_5orig, lwd=3)
#    }
#      
#    if(grepl("bio_6", bios)) {
#      bio_6 <- min(bioD$bio_6,na.rm=T) + multiplier*(max(bioD$bio_6, na.rm=T) - min(bioD$bio_6, na.rm=T))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      bio_6 <- c(min(species$bio_6) - 10, min(species$bio_6), max(species$bio_6), max(species$bio_6) + 10, bio_6)
#      
#      input <- data.frame(bio_6) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      input$bio_6orig <- input$bio_6
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      
#      if (p[1,]>p[2,]) {
#        bio.sp$bio_6[bio.sp$bio_6 < min(species$bio_6)]  <- min(species$bio_6)
#        input$bio_6[input$bio_6 < min(species$bio_6)]  <- min(species$bio_6)
#      } ## TRUE: upward trend on left end
#      ## reassign values below min to simulate a flatline
#      
#      if (p[4,]>p[3,]) {
#        bio.sp$bio_6[bio.sp$bio_6 > max(species$bio_6)]  <- max(species$bio_6) 
#        input$bio_6[input$bio_6 > max(species$bio_6)]  <- max(species$bio_6)
#      }## TRUE: upward trend on right end
#      ## reassign values above max to simulate a flatline
#      
#      xmin <- min(species$bio_6)
#      xmax <- max(species$bio_6)
#      x5 <- quantile(bioD$bio_6, .05, na.rm=T)
#      x95 <- quantile(bioD$bio_6, .95, na.rm=T) 
#      
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      colnames(p) <- c("abundance")
#      input <- cbind(input, p)
#      
#      input <- input[5:NROW(input),]
#      input <- input[order(input$bio_6),]
#      plot(exp(input$abundance)~input$bio_6orig, xlab="bio_6", ylab="infested area")
#      polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#      polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#      lines(exp(input$abundance)~input$bio_6orig, lwd=3)
#    }
#    
#    if(grepl("bio_8", bios)) {
#      bio_8 <- min(bioD$bio_8,na.rm=T) + multiplier*(max(bioD$bio_8, na.rm=T) - min(bioD$bio_8, na.rm=T))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      bio_8 <- c(min(species$bio_8) - 10, min(species$bio_8), max(species$bio_8), max(species$bio_8) + 10, bio_8)
#    
#      input <- data.frame(bio_8) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15 <- mean(species$bio_15)
#      input$bio_8orig <- input$bio_8
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#
#      if (p[1,]>p[2,]) {
#        bio.sp$bio_8[bio.sp$bio_8 < min(species$bio_8)]  <- min(species$bio_8)
#        input$bio_8[input$bio_8 < min(species$bio_8)]  <- min(species$bio_8)
#      } ## TRUE: upward trend on left end
#      ## reassign values below min to simulate a flatline
#      
#      if (p[4,]>p[3,]) {
#        bio.sp$bio_8[bio.sp$bio_8 > max(species$bio_8)]  <- max(species$bio_8) 
#        input$bio_8[input$bio_8 > max(species$bio_8)]  <- max(species$bio_8)
#      }## TRUE: upward trend on right end
#      ## reassign values above max to simulate a flatline
#      
#      input <- input[5:NROW(input),]
#      input <- input[order(input$bio_8),]
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      colnames(p) <- c("abundance")
#      input <- cbind(input, p)
#      
#      xmin <- min(species$bio_8)
#      xmax <- max(species$bio_8)
#      x5 <- quantile(bioD$bio_8, .05, na.rm=T)
#      x95 <- quantile(bioD$bio_8, .95, na.rm=T) 
#
#      plot(exp(input$abundance)~input$bio_8orig, xlab="bio_8", ylab="infested area")
#      polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#      polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#      lines(exp(input$abundance)~input$bio_8orig, lwd=3)
#    }
#    
#    if(grepl("bio_12", bios)) {
#      bio_12 <- min(bioD$bio_12,na.rm=T) + multiplier*(max(bioD$bio_12, na.rm=T) - min(bioD$bio_12, na.rm=T))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      bio_12 <- c(min(species$bio_12) - 10, min(species$bio_12), max(species$bio_12), max(species$bio_12) + 10, bio_12)
#      
#      input <- data.frame(bio_12) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_15 <- mean(species$bio_15)
#      input$bio_12orig <- input$bio_12
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      
#      if (p[1,]>p[2,]) {
#        bio.sp$bio_12[bio.sp$bio_12 < min(species$bio_12)]  <- min(species$bio_12)
#        input$bio_12[input$bio_12 < min(species$bio_12)]  <- min(species$bio_12)
#      } ## TRUE: upward trend on left end
#      ## reassign values below min to simulate a flatline
#      
#      if (p[4,]>p[3,]) {
#        bio.sp$bio_12[bio.sp$bio_12 > max(species$bio_12)]  <- max(species$bio_12) 
#        input$bio_12[input$bio_12 > max(species$bio_12)]  <- max(species$bio_12)
#      }## TRUE: upward trend on right end
#      ## reassign values above max to simulate a flatline
#      
#      input <- input[5:NROW(input),]
#      input <- input[order(input$bio_12),]
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      colnames(p) <- c("abundance")
#      input <- cbind(input, p)
#      
#      xmin <- min(species$bio_12)
#      xmax <- max(species$bio_12)
#      x5 <- quantile(bioD$bio_12, .05, na.rm=T)
#      x95 <- quantile(bioD$bio_12, .95, na.rm=T) 
#      
#      plot(exp(input$abundance)~input$bio_12orig, xlab="bio_12", ylab="infested area")
#      polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#      polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#      lines(exp(input$abundance)~input$bio_12orig, lwd=3)
#    }
#    
#    if(grepl("bio_15", bios)) {
#      bio_15 <- min(bioD$bio_15,na.rm=T) + multiplier*(max(bioD$bio_15, na.rm=T) - min(bioD$bio_15, na.rm=T))
#      ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#      ## to detect trending at min and max
#      bio_15 <- c(min(species$bio_15) - 10, min(species$bio_15), max(species$bio_15), max(species$bio_15) + 10, bio_15)
#      
#      input <- data.frame(bio_15) ## change list to data.frame
#      input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#      input$bio_5 <- mean(species$bio_5)
#      input$bio_6 <- mean(species$bio_6)
#      input$bio_8 <- mean(species$bio_8)
#      input$bio_12 <- mean(species$bio_12)
#      input$bio_15orig <- input$bio_15
#      
#      ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      input <- cbind(input, p)
#      
#      if (p[1,]>p[2,]) {
#        bio.sp$bio_15[bio.sp$bio_15 < min(species$bio_15)]  <- min(species$bio_15)
#        input$bio_15[input$bio_15 < min(species$bio_15)]  <- min(species$bio_15)
#      } ## TRUE: upward trend on left end
#      ## reassign values below min to simulate a flatline
#      
#      if (p[4,]>p[3,]) {
#        bio.sp$bio_15[bio.sp$bio_15 > max(species$bio_15)]  <- max(species$bio_15) 
#        input$bio_15[input$bio_15 > max(species$bio_15)]  <- max(species$bio_15)
#      }## TRUE: upward trend on right end
#      ## reassign values above max to simulate a flatline
#      
#      input <- input[5:NROW(input),]
#      input <- input[order(input$bio_15),]
#      p <- predict.glm(M, input)
#      p <- melt(p)
#      colnames(p) <- c("abundance")
#      input <- cbind(input, p)
#      
#      xmin <- min(species$bio_15)
#      xmax <- max(species$bio_15)
#      x5 <- quantile(bioD$bio_15, .05, na.rm=T)
#      x95 <- quantile(bioD$bio_15, .95, na.rm=T) 
#      
#      plot(exp(input$abundance)~input$bio_15orig, xlab="bio_15", ylab="infested area")
#      polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#      polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#      lines(exp(input$abundance)~input$bio_15orig, lwd=3)
#    }
#  
#    dev.off()


#      #### trying again with only signif models ####
#    model.selS <- model.sel[model.sel$nonsig == 0,] 
#    model.selS$dAIC <- model.selS$aic - min(model.selS$aic)
#    if(model.sel2$nonsig == 1 & length(model.selS$nonsig) > 0){
#      
#      model.selS <- model.selS[model.selS$dAIC <= 2,] ## select the model with the best score, 
#      ## as well as models with scores that are not significantly worse
#      
#      model.selS <- model.selS[model.selS$n == min(model.selS$n),] 
#      ## of the models with the best scores (<= min + 2), select the simplest models
#      
#      model.selS <- model.selS[model.selS$dAIC == min(model.selS$dAIC),]
#      ## select the model with the best AIC score if there are multiple models
#      ## with similar scores and the same complexity (variable number)
#      
#      M <- glm(formula(formulae.s[model.selS$i]), data = species)
#      ## reconstruct the best model, using the formula number (i) to re-access it
#      
#      Msum <- summary(M) ## extract summary of the model to access model-fit values
#      Msum2 <- as.data.frame(Msum$coefficients) ## extract table with edf and p-vals for smooths
#      
#      aic <- Msum$aic## extract adjusted R-square value of model
#      adjDsq <-  Dsquared(M,adjust=T) ## extract deviance explained by model
#      no.terms <- NROW(Msum2) ## extract the number of terms
#      no.vars <- no.terms - 1 - str_count(formulae.s[model.selS$i], "poly")
#      ## number of variables (terms - intercept - second orders, because 1st orders accted for)
#      
#      ### extract significance (p-vals) and effective degrees of freedom for each variables
#      p2 <- -1  ## start by defaulting each to -1, indicating that the variable is not in the model
#      p2.2 <- -1
#      p5 <- -1
#      p5.2 <- -1
#      p6 <- -1
#      p6.2 <- -1
#      p8 <- -1
#      p8.2 <- -1
#      p12 <- -1
#      p12.2 <- -1
#      p15 <- -1
#      p15.2 <- -1
#      
#      bios <- as.character(M$formula[3])
#      if (grepl("bio_2 ", bios)) {
#        p2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_2"]
#      }
#      if (grepl("poly(bio_2,", bios, fixed=T) ) {
#        p2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_2, 2)1"] 
#        p2.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_2, 2)2"]
#      }
#      if (grepl("bio_5 ", bios)) {
#        p5 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_5"]
#      }
#      if (grepl("poly(bio_5,", bios, fixed=T) ) {
#        p5 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_5, 2)1"] 
#        p5.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_5, 2)2"]
#      }
#      if (grepl("bio_6 ", bios)) {
#        p6 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_6"]
#      }
#      if (grepl("poly(bio_6,", bios, fixed=T) ) {
#        p6 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_6, 2)1"] 
#        p6.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_6, 2)2"]
#      }
#      if (grepl("bio_8 ", bios)) {
#        p8 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_8"]
#      }
#      if (grepl("poly(bio_8,", bios, fixed=T) ) {
#        p8 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_8, 2)1"] 
#        p8.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_8, 2)2"]
#      }
#      if (grepl("bio_12 ", bios)) {
#        p12 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_12"]
#      }
#      if (grepl("poly(bio_12,", bios, fixed=T) ) {
#        p12 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_12, 2)1"] 
#        p12.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_12, 2)2"]
#      }
#      if (grepl("bio_15 ", bios)) {
#        p15 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "bio_15"]
#      }
#      if (grepl("poly(bio_15,", bios, fixed=T) ) {
#        p15 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_15, 2)1"] 
#        p15.2 <- Msum2$`Pr(>|t|)`[row.names(Msum2) == "poly(bio_15, 2)2"]
#      }
#      
#      ordSums.i <- data.frame(species.code,aic,adjDsq,no.terms,no.vars,bios,p2,p2.2,p5,p5.2,p6,p6.2,p8,p8.2,p12,p12.2,p15,p15.2, stringsAsFactors = F)
#      ordSums <- rbind(ordSums,ordSums.i)
#      
#      ############ RESPONSE CURVES ##########
#      jpeg(paste("quad_pics2", paste(species.code, "jpg", sep="."), sep="/"))
#      par(pch=".")
#      
#      if (no.vars > 4){
#        par(mfrow=c(2,3)) 
#      }
#      if (no.vars == 4){
#        par(mfrow=c(2,2)) 
#      }
#      if (no.vars == 3){
#        par(mfrow=c(3,1)) 
#      }
#      if (no.vars == 2){
#        par(mfrow=c(2,1)) 
#      }
#      
#      multiplier <- c(0:100)/100
#      
#      if(grepl("bio_2", bios)) {
#        bio_2 <- min(species$bio_2) + multiplier*(max(species$bio_2)- min(species$bio_2))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        
#        input <- data.frame(bio_2) ## change list to data.frame
#        input$bio_5 <- mean(species$bio_5) ## set all other variables to constant value
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        plot(input$value~input$bio_2)
#        lines(input$value~input$bio_2, lwd=3)
#        
#      }
#      
#      if(grepl("bio_5", bios)) {
#        bio_5 <- min(species$bio_5) + multiplier*(max(species$bio_5)- min(species$bio_5))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        
#        input <- data.frame(bio_5) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        plot(input$value~input$bio_5)
#        lines(input$value~input$bio_5, lwd=3)
#      }
#      
#      if(grepl("bio_6", bios)) {
#        bio_6 <- min(species$bio_6) + multiplier*(max(species$bio_6)- min(species$bio_6))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        
#        input <- data.frame(bio_6) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        plot(input$value~input$bio_6)
#        lines(input$value~input$bio_6, lwd=3)
#      }
#      
#      if(grepl("bio_8", bios)) {
#        bio_8 <- min(species$bio_8) + multiplier*(max(species$bio_8) - min(species$bio_8))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        
#        input <- data.frame(bio_8) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        plot(input$value~input$bio_8)
#        lines(input$value~input$bio_8, lwd=3)
#      }
#      
#      if(grepl("bio_12", bios)) {
#        bio_12 <- min(species$bio_12) + multiplier*(max(species$bio_12)- min(species$bio_12))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        
#        input <- data.frame(bio_12) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_15 <- mean(species$bio_15)
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        plot(input$value~input$bio_12)
#        lines(input$value~input$bio_12, lwd=3)
#      }
#      
#      if(grepl("bio_15", bios)) {
#        bio_15 <- min(species$bio_15) + multiplier*(max(species$bio_15)- min(species$bio_15))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        
#        input <- data.frame(bio_15) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_5 <- mean(species$bio_5)
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        plot(input$value~input$bio_15)
#        lines(input$value~input$bio_15, lwd=3)
#      } 
#      dev.off()
#      
#      ##EXTRAP
#      jpeg(paste("quad_pics3", paste(species.code, "jpg", sep="SIG."), sep="/"))
#      par(pch=".")
#      bio.sp <- bioD
#      
#      if (no.vars > 4){
#        par(mfrow=c(2,3)) 
#      }
#      if (no.vars == 4){
#        par(mfrow=c(2,2)) 
#      }
#      if (no.vars == 3){
#        par(mfrow=c(3,1)) 
#      }
#      if (no.vars == 2){
#        par(mfrow=c(2,1)) 
#      }
#      
#      multiplier<- c(0:200)/200
#      
#      if(grepl("bio_2", bios)) {
#        bio_2 <- min(bioD$bio_2,na.rm=T) + multiplier*(max(bioD$bio_2, na.rm=T) - min(bioD$bio_2, na.rm=T))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        bio_2 <- c(min(species$bio_2) - 10, min(species$bio_2), max(species$bio_2), max(species$bio_2) + 10, bio_2)
#        
#        input <- data.frame(bio_2) ## change list to data.frame
#        input$bio_5 <- mean(species$bio_5) ## set all other variables to constant value
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        input$bio_2orig <- input$bio_2
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        if (p[1,]>p[2,]) {
#          bio.sp$bio_2[bio.sp$bio_2 < min(species$bio_2)]  <- min(species$bio_2)
#          input$bio_2[input$bio_2 < min(species$bio_2)]  <- min(species$bio_2)
#        } ## TRUE: upward trend on left end
#        ## reassign values below min to simulate a flatline
#        
#        if (p[4,]>p[3,]) {
#          bio.sp$bio_2[bio.sp$bio_2 > max(species$bio_2)]  <- max(species$bio_2) 
#          input$bio_2[input$bio_2 > max(species$bio_2)]  <- max(species$bio_2)
#        }## TRUE: upward trend on right end
#        ## reassign values above max to simulate a flatline
#        
#        xmin <- min(species$bio_2)
#        xmax <- max(species$bio_2)
#        x5 <- quantile(bioD$bio_2, .05, na.rm=T)
#        x95 <- quantile(bioD$bio_2, .95, na.rm=T) 
#        
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        colnames(p) <- c("abundance")
#        input <- cbind(input, p)
#        
#        input <- input[5:NROW(input),]
#        input <- input[order(input$bio_2),]
#        plot(exp(input$abundance)~input$bio_2orig, xlab="bio_2", ylab="infested area")
#        polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#        polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#        lines(exp(input$abundance)~input$bio_2orig, lwd=3)
#      }
#      
#      if(grepl("bio_5", bios)) {
#        bio_5 <- min(bioD$bio_5,na.rm=T) + multiplier*(max(bioD$bio_5, na.rm=T) - min(bioD$bio_5, na.rm=T))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        bio_5 <- c(min(species$bio_5) - 10, min(species$bio_5), max(species$bio_5), max(species$bio_5) + 10, bio_5)
#        
#        input <- data.frame(bio_5) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        input$bio_5orig <- input$bio_5
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        if (p[1,]>p[2,]) {
#          bio.sp$bio_5[bio.sp$bio_5 < min(species$bio_5)]  <- min(species$bio_5)
#          input$bio_5[input$bio_5 < min(species$bio_5)]  <- min(species$bio_5)
#        } ## TRUE: upward trend on left end
#        ## reassign values below min to simulate a flatline
#        
#        if (p[4,]>p[3,]) {
#          bio.sp$bio_5[bio.sp$bio_5 > max(species$bio_5)]  <- max(species$bio_5) 
#          input$bio_5[input$bio_5 > max(species$bio_5)]  <- max(species$bio_5)
#        }## TRUE: upward trend on right end
#        ## reassign values above max to simulate a flatline
#        
#        xmin <- min(species$bio_5)
#        xmax <- max(species$bio_5)
#        x5 <- quantile(bioD$bio_5, .05, na.rm=T)
#        x95 <- quantile(bioD$bio_5, .95, na.rm=T) 
#        
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        colnames(p) <- c("abundance")
#        input <- cbind(input, p)
#        
#        input <- input[5:NROW(input),]
#        input <- input[order(input$bio_5),]
#        plot(exp(input$abundance)~input$bio_5orig, xlab="bio_5", ylab="infested area")
#        polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#        polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#        lines(exp(input$abundance)~input$bio_5orig, lwd=3)
#      }
#      
#      if(grepl("bio_6", bios)) {
#        bio_6 <- min(bioD$bio_6,na.rm=T) + multiplier*(max(bioD$bio_6, na.rm=T) - min(bioD$bio_6, na.rm=T))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        bio_6 <- c(min(species$bio_6) - 10, min(species$bio_6), max(species$bio_6), max(species$bio_6) + 10, bio_6)
#        
#        input <- data.frame(bio_6) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        input$bio_6orig <- input$bio_6
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        if (p[1,]>p[2,]) {
#          bio.sp$bio_6[bio.sp$bio_6 < min(species$bio_6)]  <- min(species$bio_6)
#          input$bio_6[input$bio_6 < min(species$bio_6)]  <- min(species$bio_6)
#        } ## TRUE: upward trend on left end
#        ## reassign values below min to simulate a flatline
#        
#        if (p[4,]>p[3,]) {
#          bio.sp$bio_6[bio.sp$bio_6 > max(species$bio_6)]  <- max(species$bio_6) 
#          input$bio_6[input$bio_6 > max(species$bio_6)]  <- max(species$bio_6)
#        }## TRUE: upward trend on right end
#        ## reassign values above max to simulate a flatline
#        
#        xmin <- min(species$bio_6)
#        xmax <- max(species$bio_6)
#        x5 <- quantile(bioD$bio_6, .05, na.rm=T)
#        x95 <- quantile(bioD$bio_6, .95, na.rm=T) 
#        
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        colnames(p) <- c("abundance")
#        input <- cbind(input, p)
#        
#        input <- input[5:NROW(input),]
#        input <- input[order(input$bio_6),]
#        plot(exp(input$abundance)~input$bio_6orig, xlab="bio_6", ylab="infested area")
#        polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#        polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#        lines(exp(input$abundance)~input$bio_6orig, lwd=3)
#      }
#      
#      if(grepl("bio_8", bios)) {
#        bio_8 <- min(bioD$bio_8,na.rm=T) + multiplier*(max(bioD$bio_8, na.rm=T) - min(bioD$bio_8, na.rm=T))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        bio_8 <- c(min(species$bio_8) - 10, min(species$bio_8), max(species$bio_8), max(species$bio_8) + 10, bio_8)
#        
#        input <- data.frame(bio_8) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15 <- mean(species$bio_15)
#        input$bio_8orig <- input$bio_8
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        if (p[1,]>p[2,]) {
#          bio.sp$bio_8[bio.sp$bio_8 < min(species$bio_8)]  <- min(species$bio_8)
#          input$bio_8[input$bio_8 < min(species$bio_8)]  <- min(species$bio_8)
#        } ## TRUE: upward trend on left end
#        ## reassign values below min to simulate a flatline
#        
#        if (p[4,]>p[3,]) {
#          bio.sp$bio_8[bio.sp$bio_8 > max(species$bio_8)]  <- max(species$bio_8) 
#          input$bio_8[input$bio_8 > max(species$bio_8)]  <- max(species$bio_8)
#        }## TRUE: upward trend on right end
#        ## reassign values above max to simulate a flatline
#        
#        input <- input[5:NROW(input),]
#        input <- input[order(input$bio_8),]
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        colnames(p) <- c("abundance")
#        input <- cbind(input, p)
#        
#        xmin <- min(species$bio_8)
#        xmax <- max(species$bio_8)
#        x5 <- quantile(bioD$bio_8, .05, na.rm=T)
#        x95 <- quantile(bioD$bio_8, .95, na.rm=T) 
#        
#        plot(exp(input$abundance)~input$bio_8orig, xlab="bio_8", ylab="infested area")
#        polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#        polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#        lines(exp(input$abundance)~input$bio_8orig, lwd=3)
#      }
#      
#      if(grepl("bio_12", bios)) {
#        bio_12 <- min(bioD$bio_12,na.rm=T) + multiplier*(max(bioD$bio_12, na.rm=T) - min(bioD$bio_12, na.rm=T))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        bio_12 <- c(min(species$bio_12) - 10, min(species$bio_12), max(species$bio_12), max(species$bio_12) + 10, bio_12)
#        
#        input <- data.frame(bio_12) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_15 <- mean(species$bio_15)
#        input$bio_12orig <- input$bio_12
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        if (p[1,]>p[2,]) {
#          bio.sp$bio_12[bio.sp$bio_12 < min(species$bio_12)]  <- min(species$bio_12)
#          input$bio_12[input$bio_12 < min(species$bio_12)]  <- min(species$bio_12)
#        } ## TRUE: upward trend on left end
#        ## reassign values below min to simulate a flatline
#        
#        if (p[4,]>p[3,]) {
#          bio.sp$bio_12[bio.sp$bio_12 > max(species$bio_12)]  <- max(species$bio_12) 
#          input$bio_12[input$bio_12 > max(species$bio_12)]  <- max(species$bio_12)
#        }## TRUE: upward trend on right end
#        ## reassign values above max to simulate a flatline
#        
#        input <- input[5:NROW(input),]
#        input <- input[order(input$bio_12),]
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        colnames(p) <- c("abundance")
#        input <- cbind(input, p)
#        
#        xmin <- min(species$bio_12)
#        xmax <- max(species$bio_12)
#        x5 <- quantile(bioD$bio_12, .05, na.rm=T)
#        x95 <- quantile(bioD$bio_12, .95, na.rm=T) 
#        
#        plot(exp(input$abundance)~input$bio_12orig, xlab="bio_12", ylab="infested area")
#        polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#        polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#        lines(exp(input$abundance)~input$bio_12orig, lwd=3)
#      }
#      
#      if(grepl("bio_15", bios)) {
#        bio_15 <- min(bioD$bio_15,na.rm=T) + multiplier*(max(bioD$bio_15, na.rm=T) - min(bioD$bio_15, na.rm=T))
#        ## list of minimum value, one below minimum value, maximum value, and one above maximum value
#        ## to detect trending at min and max
#        bio_15 <- c(min(species$bio_15) - 10, min(species$bio_15), max(species$bio_15), max(species$bio_15) + 10, bio_15)
#        
#        input <- data.frame(bio_15) ## change list to data.frame
#        input$bio_2 <- mean(species$bio_2) ## set all other variables to constant value
#        input$bio_5 <- mean(species$bio_5)
#        input$bio_6 <- mean(species$bio_6)
#        input$bio_8 <- mean(species$bio_8)
#        input$bio_12 <- mean(species$bio_12)
#        input$bio_15orig <- input$bio_15
#        
#        ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        input <- cbind(input, p)
#        
#        if (p[1,]>p[2,]) {
#          bio.sp$bio_15[bio.sp$bio_15 < min(species$bio_15)]  <- min(species$bio_15)
#          input$bio_15[input$bio_15 < min(species$bio_15)]  <- min(species$bio_15)
#        } ## TRUE: upward trend on left end
#        ## reassign values below min to simulate a flatline
#        
#        if (p[4,]>p[3,]) {
#          bio.sp$bio_15[bio.sp$bio_15 > max(species$bio_15)]  <- max(species$bio_15) 
#          input$bio_15[input$bio_15 > max(species$bio_15)]  <- max(species$bio_15)
#        }## TRUE: upward trend on right end
#        ## reassign values above max to simulate a flatline
#        
#        input <- input[5:NROW(input),]
#        input <- input[order(input$bio_15),]
#        p <- predict.glm(M, input)
#        p <- melt(p)
#        colnames(p) <- c("abundance")
#        input <- cbind(input, p)
#        
#        xmin <- min(species$bio_15)
#        xmax <- max(species$bio_15)
#        x5 <- quantile(bioD$bio_15, .05, na.rm=T)
#        x95 <- quantile(bioD$bio_15, .95, na.rm=T) 
#        
#        plot(exp(input$abundance)~input$bio_15orig, xlab="bio_15", ylab="infested area")
#        polygon(x=c(x5,x5,x95,x95), y=c(-1,1000,1000,-1), col="light grey")
#        polygon(x=c(xmin,xmin,xmax,xmax), y=c(-1,1000,1000,-1), col="dark grey")
#        lines(exp(input$abundance)~input$bio_15orig, lwd=3)
#      }
#      dev.off()
#    }    

#      ############ PREDICTIONS #############
#    P <- predict.glm(M, bio.sp)
#    P <- melt(P)
#    
#    #P$value[P$value < -5] <- NA
#    predict.raster <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=P$value)
#    writeRaster(predict.raster, filename= paste("quadPreds", paste(species.code, "asc", sep="."), sep="/" ), format="ascii", overwrite=T)


#} else {
## use -99 to indicate that models were insufficient
#  aic <- -99
#  adjDsq <- -99
#  no.pts <- -99
#  no.vars <- -99
#  no.terms - -99
#  p2 <- -99 
#  p5 <- -99
#  p6 <- -99
#  p8 <- -99
#  p12 <- -99
#  p15 <- -99
#  p2.2 <- -99 
#  p5.2 <- -99
#  p6.2 <- -99
#  p8.2 <- -99
#  p12.2 <- -99
#  p15.2 <- -99
#  
#  
#  ordSums.i <- data.frame(species.code,aic,adjDsq,no.pts,no.terms,no.vars,lat.ext,lon.ext,habit,duration,bios,p2,p2.2,p5,p5.2,p6,p6.2,p8,p8.2,p12,p12.2,p15,p15.2, stringsAsFactors = F)
#  ordSums <- rbind(ordSums,ordSums.i)
#}

#print(s)

#} 
#warnings()
#ordSums <- ordSums[ordSums$species.code != "TEMPLATE",]
#write.csv(ordSums, "ordSumsNEW.csv", row.names=F)
## 3 min for 10 spp
## 5 min for 20 spp



