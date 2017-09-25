#### Practicing with GAMS
###  Mitch O'Neill
##   created: 9/4/2017
#    last modified: 9/24/2017

library(mgcv)
library(raster)
library(rgeos)

setwd("C:/Users/mwone/Documents/geodata/")

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_09_12_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
head(edd)

wisi <- edd[edd$usda == "WISI" & edd$abundance < 50,]
wisi <- data.frame(cbind(log(wisi$abundance), wisi$latitude, wisi$longitude))
colnames(wisi) <- c("abundance", "latitude", "longitude")
head(wisi)

coordinates(wisi) <- c(3,2)
proj4string(wisi) <- CRS("+init=epsg:4326")
plot(wisi)

bio <- stack("clipped_climate_data/current/bio_2.asc",
             "clipped_climate_data/current/bio_5.asc",
             "clipped_climate_data/current/bio_6.asc",
             "clipped_climate_data/current/bio_8.asc",
             "clipped_climate_data/current/bio_12.asc",
             "clipped_climate_data/current/bio_15.asc")

proj4string(bio) <- proj4string(raster("climate_data/current/bio_1"))
wisi <- spTransform(wisi, proj4string(bio))

plot(bio$bio_2)
plot(wisi, add=T, pch=19)


library(psych)
bioM <- as.matrix(bio)
corr <- corr.test(bioM, y = NULL, use = "pairwise", method="spearman", adjust="bonferroni", alpha=.05, ci=F)
corr$r

ext <- extract(bio, wisi)
wisi <- cbind(wisi, ext)
head(wisi)


## gam = generalized additive model
## response ~ smooth function of predictors;
## fx=F, k=-1: amount of smoothing is not fixed, will be optimized using cross-validation
## bs="cr": cubic regression splines, 
## bs="tp": thin-plate regression bs= "cs": cubic regression shrinkage, bs="ts": thin-plate shrinkage



## THIN PLATE
M1 <- gam(wisi$abundance ~ s(wisi$bio_2, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_5, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_6, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_8, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_12, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_15, fx=F, k=-1, bs="tp"), 
                           method="GCV.Cp") 
plot(M1)
summary(M1)
## signif: none??

## THIN PLATE WITH REML
M2 <- gam(wisi$abundance ~ s(wisi$bio_2, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_5, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_6, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_8, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_12, fx=F, k=-1, bs="tp") +
                           s(wisi$bio_15, fx=F, k=-1, bs="tp"), 
                           method="REML") 
plot(M2)
summary(M2)
## signif: none??


## THIN PLATE WITH REML
M3 <- gam(abundance ~ s(bio_2, fx=F, k=-1, bs="ts") +
                           s(bio_5, fx=F, k=-1, bs="ts") +
                           s(bio_6, fx=F, k=-1, bs="ts") +
                           s(bio_8, fx=F, k=-1, bs="ts") +
                           s(bio_12, fx=F, k=-1, bs="ts") +
                           s(bio_15, fx=F, k=-1, bs="ts"), 
                           method="REML", data=wisi) 
plot(M3)
summary(M3)
## signif: bio 5 and bio 6 (max and min temp), everything else drops out?
## no precips were kept -- maybe test with max and or min?


bioD <- as.data.frame(bio, stringsAsFactors=F)
P3 <- predict.gam(M3, newdata=bioD)
head(P3)

library(reshape2)
P3m <- melt(P3)
head(P3m)
P3m$Var1 <- NULL
colnames(P3m) <- "predicted_abun"
head(P3m)


ras <- raster(ncols=bio@ncols, nrows=bio@nrows, ext=bio@extent, crs=bio@crs, vals=P3m$predicted_abun) ## specify extent and proj4string based on bio1
plot(ras)
## that does not look right


## take out extrapolation zone
P3m$predicted_abun[P3m$predicted_abun < min(wisi$abundance) | P3m$predicted_abun > max(wisi$abundance)] <- NA
ras <- raster(ncols=bio@ncols, nrows=bio@nrows, ext=bio@extent, crs=bio@crs, vals=P3m$predicted_abun) ## specify extent and proj4string based on bio1
plot(ras)



