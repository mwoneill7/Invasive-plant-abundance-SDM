#### Ordinal SDMs
###  created: 11/10/2017
##   last modified: 12/19/2017
#    Mitch O'Neill

library(raster)   ## for Raster climate data
library(rgeos)    ## for extractig climate data to points
library(reshape2) ## for reformatting model predictions
library(combinat) ## for generating pairs of covariates for assessment of collinearity 
library(glue)     ## for collapsing lists of variables into equations
#library(modEvA)
library(stringr)
library(MASS)
library(irr)

## list candidate terms (all six variables + quadratics)
var.list <- c("bio_2", "poly(bio_2, 2)",
              "bio_5", "poly(bio_5, 2)",
              "bio_6", "poly(bio_6, 2)",
              "bio_8", "poly(bio_8, 2)",
              "bio_12", "poly(bio_12, 2)",
              "bio_15", "poly(bio_15, 2)",
              "nlcd_2", "poly(nlcd_1 ,2)",
              "nlcd_3", "poly(nlcd_3 ,2)",
              "nlcd_4", "poly(nlcd_4 ,2)",
              "nlcd_5", "poly(nlcd_5 ,2)",
              "nlcd_6", "poly(nlcd_6 ,2)",
              "nlcd_7", "poly(nlcd_7 ,2)",
              "nlcd_8", "poly(nlcd_8 ,2)",
              "nlcd_9", "poly(nlcd_9 ,2)")

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
       !(grepl("bio_15 ", combo.j) & grepl("poly(bio_15", combo.j, fixed=T)) &
       !(grepl("nlcd_2 ", combo.j) & grepl("poly(nlcd_2", combo.j, fixed=T)) &
       !(grepl("nlcd_3 ", combo.j) & grepl("poly(nlcd_3", combo.j, fixed=T)) &
       !(grepl("nlcd_4 ", combo.j) & grepl("poly(nlcd_4", combo.j, fixed=T)) &
       !(grepl("nlcd_5 ", combo.j) & grepl("poly(nlcd_5", combo.j, fixed=T)) &
       !(grepl("nlcd_6 ", combo.j) & grepl("poly(nlcd_6", combo.j, fixed=T)) &
       !(grepl("nlcd_7 ", combo.j) & grepl("poly(nlcd_7", combo.j, fixed=T)) &
       !(grepl("nlcd_8 ", combo.j) & grepl("poly(nlcd_8", combo.j, fixed=T)) &
       !(grepl("nlcd_9 ", combo.j) & grepl("poly(nlcd_9", combo.j, fixed=T))){
      var.sets <- c(var.sets,combo.j)} ## add it to the list of variable sets
  }
  
  print(i)
}

formulae <- paste("abundance", var.sets[2:length(var.sets)], sep = " ~ ")
## add response part of formulaes, excluding the initial row

rm(var.sets, combo.j, combos.i, i, j, var.list) ## garbage cleaning

setwd("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/")
## read in EDDMapS dataset, thinned to 5km climate grid cells (1 point per species per cell)
edd <- read.table("final_dataset_1_25_2018.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd$bio_5<-NULL;edd$bio_6<-NULL

#edd <- read.table("C://Users/mwone/Documents/EDDMapS data/eddmaps_thinned_12_21_2017.csv", header = T, sep = ",", quote= "\"", 
#                  comment.char= "", stringsAsFactors = F, strip.white = T)
#edd$abundance <- ceiling(edd$best)
##### combine 1s and 2s
##edd$abundance[edd$abundance==1]<-2
#hist(edd$abundance)
#edd <- edd[edd$species %in% species.list$usda.code,]
#length(unique(edd$species))
head(edd)
#prob <- edd[2932,]
#proj4string(prob) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


#edd$max<-NULL;edd$med<-NULL;edd$best<-NULL;edd$no.pts<-NULL;edd$no.best<-NULL;head(edd)
coordinates(edd) <- c(5,4)
proj4string(edd) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## from transformations during data prep

#hist(edd$abundance)
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
c.up <- -99
c.down <- -99
pS <- -99
pI <- -99
null <- -99
kappa <- -99
kappaP <- -99
kappaW <- -99
regS <- -99
regI <-99
#c2    <- -999
#c2.2  <- -999
#c5    <- -999
#c5.2  <- -999
#c6    <- -999
#c6.2  <- -999
#c8    <- -999
#c8.2  <- -999
#c12   <- -999
#c12.2 <- -999
#c15   <- -999
#c15.2 <- -999

ordsums <- data.frame(species.code,aic,no.pts,no.terms,no.vars,formu,c.up,c.down,pS,pI,null,kappa,kappaP,kappaW,regS,regI, stringsAsFactors = F)
#ordsums <- data.frame(species.code,aic,no.pts,no.terms,no.vars,formu,c2,c2.2,c5,c5.2,c6,c6.2,c8,c8.2,c12,c12.2,c15,c15.2, stringsAsFactors = F)
#ordsums <- data.frame(species.code,aic,no.pts,no.terms,no.vars,p05s,p10s,p05i,p10i,formu,stringsAsFactors = F)
sp.list <- unique(edd$species)

for(s in 1:length(sp.list)){ ## Loop through species list
  
  species.code <- sp.list[s] ## extract the USDA species code for the species of the iteration

    #if(species.code !="HIAU" &
    #   species.code != "CYSC4" &
    #   species.code != "VIVI" &
    #   species.code != "CICA" & species.code != "SEPU7" & species.code != "MESA" &
    #   species.code != "RUAC3") {
    species <- edd[edd$species == species.code,]## subset to species of the iteration
    #species <- data.frame(cbind(species$abundance, species$latitude, species$longitude)) ## select only variables needed to save computation time
    species$abundance <- ordered(as.factor(species$abundance), levels=c(1,2,3))
    
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
      #tryCatch(c=x, call=x, cond=x, expr=x)
      tryCatch({
        
        M <- polr(formula(formulae.s[i]), data = species, Hess=T)
        edf <- M$edf ## REPLACE n
        aic <- AIC(M) ## extract aic score
        model.sel.i <- data.frame(i, edf, aic, stringsAsFactors = F)
        model.sel <- rbind(model.sel, model.sel.i) ## append data from iteration to the master data frame
        
      },error=function(e){cat(species.code,conditionMessage(e), "\n")})
    
    }

    # M <- lrm(formula(formulae.s[i]), data = species)|# M <- clm(formula(formulae.s[i]), data=species)|## construct model using formula of the iteration|## use data for species of this iteration

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
    Msum <- data.frame(cbind(Msum, "p_value" = p))
    
    MsumI <- Msum[(nrow(Msum)-1):nrow(Msum),] ## 1 back to 2
    pI <- length(MsumI$p_value[MsumI$p_value>0.05])
    
    MsumS <- Msum[1:(nrow(Msum)-2),] ## 2 back to 3
    pS <- length(MsumS$p_value[MsumS$p_value>0.05])
    
    MsumS$name <- as.character(row.names(MsumS))
    MsumS$poly <- sapply(MsumS$name, grepl, pattern="2)2")

    c.up <- length(MsumS$Value[MsumS$Value > 0 & MsumS$poly == T])
    c.down <- length(MsumS$Value[MsumS$Value < 0 & MsumS$poly == T])

    aic <- AIC(M)          ## extract adjusted R-square value of model #adjDsq <-  Dsquared(M,adjust=T) ## extract deviance explained by model
    no.pts <- nrow(species)
    no.terms <- NROW(Msum) ## extract the number of terms
    no.vars <- no.terms - 3 - str_count(formulae.s[model.sel2$i], "poly")
    formu <- formulae.s[model.sel2$i]
    #p05s<- length(MsumS$p.value[MsumS$p.value >= 0.05])
    #p10s<- length(MsumS$p.value[MsumS$p.value > 0.1])
    #p05i<- length(MsumI$p.value[MsumI$p.value >= 0.05])
    #p10i<- length(MsumI$p.value[MsumI$p.value > 0.1])

    #############################################################################
    #c2    <- Msum$Value[Msum$name == "poly(bio_2, 2)1" | Msum$name == "bio_2"]
    #c2.2  <- Msum$Value[Msum$name == "poly(bio_2, 2)2"]
    #c5    <- Msum$Value[Msum$name == "poly(bio_5, 2)1" | Msum$name == "bio_5"]
    #c5.2  <- Msum$Value[Msum$name == "poly(bio_5, 2)2"]
    #c6    <- Msum$Value[Msum$name == "poly(bio_6, 2)1" | Msum$name == "bio_6"]
    #c6.2  <- Msum$Value[Msum$name == "poly(bio_6, 2)2"]
    #c8    <- Msum$Value[Msum$name == "poly(bio_8, 2)1" | Msum$name == "bio_8"]
    #c8.2  <- Msum$Value[Msum$name == "poly(bio_8, 2)2"]
    #c12   <- Msum$Value[Msum$name == "poly(bio_12, 2)1" | Msum$name == "bio_12"]
    #c12.2 <- Msum$Value[Msum$name == "poly(bio_12, 2)2"]
    #c15   <- Msum$Value[Msum$name == "poly(bio_15, 2)1" | Msum$name == "bio_15"]
    #c15.2 <- Msum$Value[Msum$name == "poly(bio_15, 2)2"]
    #if(length(c2) == 0){c2 <- 0}
    #if(length(c2.2) == 0){c2.2 <- 0}
    #if(length(c5) == 0){c5 <- 0}
    #if(length(c5.2) == 0){c5.2 <- 0}
    #if(length(c6) == 0){c6 <- 0}
    #if(length(c6.2) == 0){c6.2 <- 0}
    #if(length(c8) == 0){c8 <- 0}
    #if(length(c8.2) == 0){c8.2 <- 0}
    #if(length(c12) == 0){c12 <- 0}
    #if(length(c12.2) == 0){c12.2 <- 0}
    #if(length(c15) == 0){c15 <- 0}
    #if(length(c15.2) == 0){c15.2 <- 0}    
    #############################################################################
    
    ###check against NULL
    null <- polr(abundance~1,data=species, Hess=T)
    if(AIC(null) > (AIC(M) + 2)){null <- 0} else {null <- 1}  
    
    pm <- predict(M, species)
    pm <- melt(pm)
    pm <- cbind(species$abundance,pm)
    colnames(pm) <- c("observed", "predicted")
    #pm$observed <- as.numeric(pm$observed)
    #pm$predicted <- as.numeric(pm$predicted)
    #cor(as.numeric(as.character(pm$observed)),as.numeric(as.character(pm$predicted)),method="spearman")
    #levels(pm$observed); levels(pm$predicted)
    
    kappa <- kappa2(pm)
    kappaP <- kappa$p.value
    kappa <- kappa$value
    
      
    kappaW <- kappa2(pm, "equal")  
    kappaW <- kappaW$value
    
    reg <- summary(glm(as.numeric(predicted)~as.numeric(observed), data=pm))$coefficients
    regS <- reg[1,1]
    regI <- reg[2,1]
    
    ###############################
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
    #ordsums.i <- data.frame(species.code,aic,no.pts,no.terms,no.vars,formu,c2,c2.2,c5,c5.2,c6,c6.2,c8,c8.2,c12,c12.2,c15,c15.2,kappa,kappaW,regS,regI,p, stringsAsFactors = F)
    ordsums.i <- data.frame(species.code,aic,no.pts,no.terms,no.vars,formu,c.up,c.down,pS,pI,null,kappa,kappaP,kappaW,regS,regI, stringsAsFactors = F)
    
    ordsums <- rbind(ordsums,ordsums.i) 
    print(s)
  }
}

ordsums <- ordsums[ordsums$species.code != "TEMPLATE",]


hist(ordsums$kappa)
length(ordsums$species.code[ordsums$kappaP <=0.05 & !is.na(ordsums$kappaP)]) ##108
length(ordsums$species.code[ordsums$kappaP <=0.01 & !is.na(ordsums$kappaP)]) ##100
par(mfrow=c(2,2))
hist(ordsums$kappa[ordsums$kappaP <=0.05 & !is.na(ordsums$kappaP)], breaks=40)
hist(ordsums$kappa[ordsums$kappaP <=0.01 & !is.na(ordsums$kappaP)], breaks=40)
hist(ordsums$kappaW[ordsums$kappaP <=0.05 & !is.na(ordsums$kappaP)], breaks=40)
hist(ordsums$kappaW[ordsums$kappaP <=0.01 & !is.na(ordsums$kappaP)], breaks=40)

#hist(log(ordsums$kappaW[ordsums$kappaP <=0.05 & !is.na(ordsums$kappaP)]))


## 11:55-12:45 for 78
#hist(ordsums$p10s/(ordsums$no.terms-3),breaks=50)
#plot(ordsums$p10s/(ordsums$no.terms-3)~ordsums$no.pts)
hist(ordsums$kappa, main="three bins", xlab="kappa")
hist(ordsums$kappaW)
nrow(ordsums[ordsums$kappaW>.1,]) ##103
nrow(ordsums[ordsums$kappa>.1,])  ##92
nrow(ordsums[ordsums$kappaW>.2,]) ##63
nrow(ordsums[ordsums$kappa>.2,])  ##54
summary()

write.csv(ordsums,"ordsums1_26THREEBINS.csv",row.names=F)

summary(ordsums$kappa)
summary(ordsums$null)
nulls <- ordsums[ordsums$null == 1,] ## 11 with low kappas
plot(ordsums$kappa~ordsums$kappaW)
abline(a=0,b=1)

summary(ordsums$kappa)
summary(ordsums$kappaW)
summary(ordsums$regS)
summary(ordsums$regI)
summary(ordsums$null[ordsums$kappa <= 0])

cor(ordsums$pS,ordsums$kappa) #-0.2
#summary(ordsums.i$c.up)
#summary(ordsums.i$c.down)

ordsums[ordsums$kappa > 0.5,] ##ruar9
hist(ordsums$kappa, xlab="kappa", main="kappa")
hist(ordsums$kappaW)

cor(ordsums$kappa,ordsums$no.pts)
cor(ed.listing$kappa,ed.listing$rarest.bin)

ordsums[ordsums$null ==1,]
ordsums[ordsums$kappa <0 ,]
nrow(ordsums[ordsums$kappa == 0 ,])
#pm2 <- pm
#pm2$observed <- as.numeric(pm2$observed)
#pm2$predicted <- as.numeric(pm2$predicted)
#pm2$observed[pm2$observed==1] <- 2
#pm2$predicted[pm2$observed==1] <- 2
#pm2$observed <- ordered(as.factor(pm2$observed), levels=c(2,3,4))
#pm2$predicted <- ordered(as.factor(pm2$predicted), levels=c(2,3,4))
#kappa2(pm,"equal")
#kappa2(pm)


###########################################################

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




reg <- "regina phyllange"
#####################################
for (i in 1:10) {
tryCatch({
  print(i)
  p <- i+1
  if (i==7) stop("Urgh, the iphone is in the blender !")
}, error=function(e){cat(reg, conditionMessage(e), "\n")})
}

