#### Practicing with GAMS
###  Mitch O'Neill
##   created: 9/4/2017
#    last modified: 9/24/2017

library(mgcv)
library(raster)
library(rgeos)
library(reshape2)
library(combinat)
library(rgdal)

setwd("C:/Users/mwone/Documents/geodata/")

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_09_12_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
head(edd)

species <- edd[edd$usda == "VEIN9", ]#& edd$abundance < 50,]
species <- data.frame(cbind(log(species$abundance), species$latitude, species$longitude))
colnames(species) <- c("abundance", "latitude", "longitude")
head(species)

hist(species$abundance,  breaks=20)

coordinates(species) <- c(3,2)
proj4string(species) <- CRS("+init=epsg:4326")

bio <- stack("clipped_climate_data/current/bio_2.asc",
             "clipped_climate_data/current/bio_5.asc",
             "clipped_climate_data/current/bio_6.asc",
             "clipped_climate_data/current/bio_8.asc",
             "clipped_climate_data/current/bio_12.asc",
             "clipped_climate_data/current/bio_15.asc")

proj4string(bio) <- proj4string(raster("climate_data/current/bio_1"))
species <- spTransform(species, proj4string(bio))

#plot(bio$bio_5)
#plot(species, add=T, pch=19, cex=0.5)
head(species)
ext <- extract(bio, species)
species <- cbind(species, ext)
species[,2:5] <- species[,2:5]/10

head(species)


n=6 ## using all six variables
M <- gam(species$abundance ~ s(bio_2, fx=F, k=5, bs="tp") +
                             s(bio_5, fx=F, k=5, bs="tp") +
                             s(bio_6, fx=F, k=5, bs="tp") +
                             s(bio_8, fx=F, k=5, bs="tp") +
                             s(bio_12, fx=F, k=5, bs="tp") + 
                             s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data=species)

aic <- M$aic
bio_2 <- 1
bio_5 <- 1
bio_6 <- 1
bio_8 <- 1
bio_12 <- 1
bio_15 <- 1
bios <- "bio_2 bio_5 bio_6 bio_8 bio_12 bio_15"
model.sel <- data.frame(n,aic, bio_2, bio_5, bio_6, bio_8, bio_12, bio_15, bios, stringsAsFactors = F)


#library(combinat)
variables <- c(1,2,3,4,5,6) # c(1:length(cand.vars))

for (n in 1:5){ ## eventually 1:length(variables)

  var.sets <- as.data.frame(combn(variables, n, simplify=F))
  #var.sets ## each column is a list of variables
  
  for (i in 1:NCOL(var.sets)){ 
    var.list.i <- var.sets[,i]
   
    if(n==1){
      M <- gam(species$abundance ~ s(ext[,var.list.i[1]], fx=F, k=5, bs="tp"), method="GCV.Cp")
    }
    if(n==2){ M <- gam(species$abundance ~ s(ext[,var.list.i[1]], fx=F, k=5, bs="tp") +
                                           s(ext[,var.list.i[2]], fx=F, k=5, bs="tp"), method="GCV.Cp")
    }
    if(n==3){
      M <- gam(species$abundance ~ s(ext[,var.list.i[1]], fx=F, k=5, bs="tp") +
                                   s(ext[,var.list.i[2]], fx=F, k=5, bs="tp") +
                                   s(ext[,var.list.i[3]], fx=F, k=5, bs="tp"), method="GCV.Cp")
    }
    if(n==4){
      M <- gam(species$abundance ~ s(ext[,var.list.i[1]], fx=F, k=5, bs="tp") +
                                   s(ext[,var.list.i[2]], fx=F, k=5, bs="tp") +
                                   s(ext[,var.list.i[3]], fx=F, k=5, bs="tp") +
                                   s(ext[,var.list.i[4]], fx=F, k=5, bs="tp"), method="GCV.Cp")
    }
    if(n==5){
       M <- gam(species$abundance ~ s(ext[,var.list.i[1]], fx=F, k=5, bs="tp") +
                                    s(ext[,var.list.i[2]], fx=F, k=5, bs="tp") +
                                    s(ext[,var.list.i[3]], fx=F, k=5, bs="tp") +
                                    s(ext[,var.list.i[4]], fx=F, k=5, bs="tp") +
                                    s(ext[,var.list.i[5]], fx=F, k=5, bs="tp"), method="GCV.Cp")
    }
    
    aic <- M$aic
    bios <- paste(as.character(colnames(ext)[var.list.i]), collapse=" ")
    if (grepl("bio_2", bios)) {bio_2 <- 1} else {bio_2 <- 0}
    if (grepl("bio_5", bios)) {bio_5 <- 1} else {bio_5 <- 0}
    if (grepl("bio_6", bios)) {bio_6 <- 1} else {bio_6 <- 0}
    if (grepl("bio_8", bios)) {bio_8 <- 1} else {bio_8 <- 0}
    if (grepl("bio_12", bios)) {bio_12 <- 1} else {bio_12 <- 0}
    if (grepl("bio_15", bios)) {bio_15 <- 1} else {bio_15 <- 0}
    
    model.sel.i <- data.frame(n,aic, bio_2, bio_5, bio_6, bio_8, bio_12, bio_15, bios, stringsAsFactors = F)
    model.sel <- rbind(model.sel, model.sel.i)
    print(i)
  }
  print(paste("n =", n, sep = " "))
}


cor <- cor(ext, use="pairwise", method="spearman")
rhos <- as.data.frame(corr$r)

variables <-c(1:6)
var.pairs<- combn(variables, 2, simplify=T)

var1 <- "TEMPLATE"
var2 <- "TEMPLATE"
rho <- -2
corr.pairs <- data.frame(var1,var2,rho, stringsAsFactors = F)

for (i in 1:NCOL(var.pairs)) {

  var1 <- row.names(rhos[var.pairs[1,i],])
  var2 <- row.names(rhos[var.pairs[2,i],])
  rho <- abs(rhos[var.pairs[1,i], var.pairs[2,i]])
  corr.pairs.i <- data.frame(var1,var2,rho, stringsAsFactors = F)
  corr.pairs <- rbind(corr.pairs, corr.pairs.i)
  
  print(i)
}

corr.pairs <- corr.pairs[corr.pairs$var1 != "TEMPLATE",]
corr.pairs <- corr.pairs[corr.pairs$rho >= .8, ] ## change to 8, but be careful

model.sel$collinear <- 0

if(length(corr.pairs$rho) > 0){
  for (i in 1:length(corr.pairs$rho)){
    for (j in 1:length(model.sel$bios)){
      
      if(grepl(corr.pairs$var1[i], model.sel$bios[j]) & grepl(corr.pairs$var1[i], model.sel$bios[j])) { 
        model.sel$collinear[j] <- 1}
      print(j)
    }
    print(paste("i =", i, sep = " "))
  }
}

model.sel <- model.sel[model.sel$collinear == 0,]
model.sel$dAIC <- model.sel$aic - min(model.sel$aic)

model.sel <- model.sel[model.sel$dAIC <=2,]
model.sel <- model.sel[model.sel$n == min(model.sel$n),]
model.sel <- model.sel[model.sel$dAIC == min(model.sel$dAIC),]

mod.sel <- as.data.frame(t(model.sel[,3:8]), stringsAsFactors=F)
colnames(mod.sel) <- c("mod")
mod.sel$vars <- row.names(mod.sel)
vars <- mod.sel$vars[mod.sel$mod==1]

M1 <- gam(abundance ~ s(bio_2, fx=F, k=3, bs="tp") + s(bio_5, fx=F, k=3, bs="tp") + 
                      s(bio_8, fx=F, k=3, bs="tp") + s(bio_12, fx=F, k=3, bs="tp") +
                      s(bio_15, fx=F, k=3, bs="tp"), method="GCV.Cp", data =species)
par(mfrow=c(2,3))
plot(M1)
summary(M1)


################ all possible models #############################################################################
M2 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M3 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M4 <- gam(abundance ~ s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M5 <- gam(abundance ~ s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M6 <- gam(abundance ~ s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M7 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M8 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M9 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M10 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M11 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M12 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M13 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M14 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M15 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M16 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M17 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M18 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M19 <- gam(abundance ~ s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M20 <- gam(abundance ~ s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M21 <- gam(abundance ~ s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M22 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M23 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M24 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M25 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M26 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M27 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M28 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M29 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M30 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M31 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M32 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M33 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M34 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M35 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M36 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M37 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M38 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M39 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M40 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M41 <- gam(abundance ~ s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M42 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M43 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M44 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M45 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M46 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M47 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M48 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M49 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M50 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M51 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M52 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M53 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M54 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M55 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M56 <- gam(abundance ~ s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)  
M57 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M58 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M59 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M60 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M61 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M62 <- gam(abundance ~ s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)
M63 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") + s(bio_5, fx=F, k=5, bs="tp") + s(bio_6, fx=F, k=5, bs="tp") + s(bio_8, fx=F, k=5, bs="tp") + s(bio_12, fx=F, k=5, bs="tp") + s(bio_15, fx=F, k=5, bs="tp"), method="GCV.Cp", data =species)



##########################################
M2 <- gam(abundance ~ s(bio_2, fx=F, k=5, bs="tp") +
                      s(bio_6, fx=F, k=5, bs="tp") +
                      s(bio_12, fx=F, k=5, bs="tp"),
                      method="GCV.Cp", data=species) 
par(mfrow=c(1,3)); plot(M2)
summary(M2)

bio.sp <- as.data.frame(bio) ## make copy of bioclimatic variables for manipulation
bio.sp[,1:4] <- bio.sp[,1:4]/10

## if bio_2 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
if(grepl("bio_2", as.character(M2$pred.formula)[2])) {
  
  bio_2 <- c(min(species$bio_2), min(species$bio_2 - 1), max(species$bio_2), max(species$bio_2) + 1)
  ## list of minimum value, one below minimum value, maximum value, and one above maximum value
  ## to detect trending at min and max
  
  input <- data.frame(bio_2) ## change list to data.frame
  input$bio_5 <- 0 ## set all other variables to constant value
  input$bio_6 <- 0
  input$bio_8 <- 0
  input$bio_12 <- 0
  input$bio_15 <- 0
  
  ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
  p <- predict.gam(M2, input)
  
  if (p[2]>p[1]) {bio.sp$bio_2[bio.sp$bio_2 < min(species$bio_2)]  <- min(species$bio_2)
  } ## TRUE: upward trend on left end
  ## reassign values below min to simulate a flatline
  
  if (p[4]>p[3]) {bio.sp$bio_2[bio.sp$bio_2 > max(species$bio_2)]  <- max(species$bio_2) 
  }## TRUE: upward trend on right end
  ## reassign values above max to simulate a flatline
}

## if bio_5 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
if(grepl("bio_5", as.character(M2$pred.formula)[2])) {
  
  bio_5 <- c(min(species$bio_5), min(species$bio_5 - 1), max(species$bio_5), max(species$bio_5) + 1)
  ## list of minimum value, one below minimum value, maximum value, and one above maximum value
  ## to detect trending at min and max
  
  input <- data.frame(bio_5) ## change list to data.frame
  input$bio_2 <- 0 ## set all other variables to constant value
  input$bio_6 <- 0
  input$bio_8 <- 0
  input$bio_12 <- 0
  input$bio_15 <- 0
  
  ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
  p <- predict.gam(M2, input)
  
  if (p[2]>p[1]) {bio.sp$bio_5[bio.sp$bio_5 < min(species$bio_5)]  <- min(species$bio_5)
  } ## TRUE: upward trend on left end
  ## reassign values below min to simulate a flatline
  
  if (p[4]>p[3]) {bio.sp$bio_5[bio.sp$bio_5 > max(species$bio_5)]  <- max(species$bio_5) 
  }## TRUE: upward trend on right end
  ## reassign values above max to simulate a flatline
}

## if bio_6 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
if(grepl("bio_6", as.character(M2$pred.formula)[2])) {
  
  bio_6 <- c(min(species$bio_6), min(species$bio_6 - 1), max(species$bio_6), max(species$bio_6) + 1)
  ## list of minimum value, one below minimum value, maximum value, and one above maximum value
  ## to detect trending at min and max
  
  input <- data.frame(bio_6) ## change list to data.frame
  input$bio_2 <- 0 ## set all other variables to constant value
  input$bio_5 <- 0
  input$bio_8 <- 0
  input$bio_12 <- 0
  input$bio_15 <- 0
  
  ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
  p <- predict.gam(M2, input)
  
  if (p[2]>p[1]) {bio.sp$bio_6[bio.sp$bio_6 < min(species$bio_6)]  <- min(species$bio_6)
  } ## TRUE: upward trend on left end
  ## reassign values below min to simulate a flatline
  
  if (p[4]>p[3]) {bio.sp$bio_6[bio.sp$bio_6 > max(species$bio_6)]  <- max(species$bio_6) 
  }## TRUE: upward trend on right end
  ## reassign values above max to simulate a flatline
}

## if bio_8 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
if(grepl("bio_8", as.character(M2$pred.formula)[2])) {
  
  bio_8 <- c(min(species$bio_8), min(species$bio_8 - 1), max(species$bio_8), max(species$bio_8) + 1)
  ## list of minimum value, one below minimum value, maximum value, and one above maximum value
  ## to detect trending at min and max
  
  input <- data.frame(bio_8) ## change list to data.frame
  input$bio_2 <- 0 ## set all other variables to constant value
  input$bio_5 <- 0
  input$bio_6 <- 0
  input$bio_12 <- 0
  input$bio_15 <- 0
  
  ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
  p <- predict.gam(M2, input)
  
  if (p[2]>p[1]) {bio.sp$bio_8[bio.sp$bio_8 < min(species$bio_8)]  <- min(species$bio_8)
  } ## TRUE: upward trend on left end
  ## reassign values below min to simulate a flatline
  
  if (p[4]>p[3]) {bio.sp$bio_8[bio.sp$bio_8 > max(species$bio_8)]  <- max(species$bio_8) 
  }## TRUE: upward trend on right end
  ## reassign values above max to simulate a flatline
}

## if bio_12 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
if(grepl("bio_12", as.character(M2$pred.formula)[2])) {
  
  bio_12 <- c(min(species$bio_12), min(species$bio_12 - 1), max(species$bio_12), max(species$bio_12) + 1)
  ## list of minimum value, one below minimum value, maximum value, and one above maximum value
  ## to detect trending at min and max
  
  input <- data.frame(bio_12) ## change list to data.frame
  input$bio_2 <- 0 ## set all other variables to constant value
  input$bio_5 <- 0
  input$bio_6 <- 0
  input$bio_8 <- 0
  input$bio_15 <- 0
  
  ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
  p <- predict.gam(M2, input)
  
  if (p[2]>p[1]) {bio.sp$bio_12[bio.sp$bio_12 < min(species$bio_12)]  <- min(species$bio_12)
  } ## TRUE: upward trend on left end
  ## reassign values below min to simulate a flatline
  
  if (p[4]>p[3]) {bio.sp$bio_12[bio.sp$bio_12 > max(species$bio_12)]  <- max(species$bio_12) 
  }## TRUE: upward trend on right end
  ## reassign values above max to simulate a flatline
}

## if bio_15 is in the formula, then assess upward trends in extrapolations & simulate flatlines if necessary
if(grepl("bio_15", as.character(M2$pred.formula)[2])) {
  
  bio_15 <- c(min(species$bio_15), min(species$bio_15 - 1), max(species$bio_15), max(species$bio_15) + 1)
  ## list of minimum value, one below minimum value, maximum value, and one above maximum value
  ## to detect trending at min and max
  
  input <- data.frame(bio_15) ## change list to data.frame
  input$bio_5 <- 0 ## set all other variables to constant value
  input$bio_6 <- 0
  input$bio_8 <- 0
  input$bio_12 <- 0
  input$bio_15 <- 0
  
  ## predict using dataframe of max, max+1, min, min-1, w/ all other vars constant
  p <- predict.gam(M2, input)
  
  if (p[2]>p[1]) {bio.sp$bio_15[bio.sp$bio_15 < min(species$bio_15)]  <- min(species$bio_15)
  } ## TRUE: upward trend on left end
  ## reassign values below min to simulate a flatline
  
  if (p[4]>p[3]) {bio.sp$bio_15[bio.sp$bio_15 > max(species$bio_15)]  <- max(species$bio_15) 
  }## TRUE: upward trend on right end
  ## reassign values above max to simulate a flatline
}


P <- predict.gam(M2, bio.sp)
P1 <- melt(P)
P1$value[P1$value < -12] <- NA
P1r <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=P1$value)
par(mfrow=c(1,1))

states<- readOGR(dsn="states", layer= "US_states")
plot(states, col="grey")
plot(P1r, add=T)
plot(species, col="red", pch=19, add=T)





