#### Practicing with GAMS
###  Mitch O'Neill
##   created: 9/4/2017
#    last modified: 9/24/2017

library(mgcv)
library(raster)
library(rgeos)
library(reshape2)

setwd("C:/Users/mwone/Documents/geodata/")

edd <- read.table("C:/Users/mwone/Documents/EDDMapS data/eddmaps_thinned_09_12_2017.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
head(edd)
################################################################################
################################################################################

species <- edd[edd$usda == "VEIN9" & edd$abundance < 50,]
species <- data.frame(cbind(log(species$abundance), species$latitude, species$longitude))
colnames(species) <- c("abundance", "latitude", "longitude")
head(species)

coordinates(species) <- c(3,2)
proj4string(species) <- CRS("+init=epsg:4326")

hist(species$abundance)

bio <- stack("clipped_climate_data/current/bio_2.asc",
             "clipped_climate_data/current/bio_5.asc",
             "clipped_climate_data/current/bio_6.asc",
             "clipped_climate_data/current/bio_8.asc",
             "clipped_climate_data/current/bio_11.asc",
             "clipped_climate_data/current/bio_12.asc",
             "clipped_climate_data/current/bio_15.asc")

proj4string(bio) <- proj4string(raster("climate_data/current/bio_1"))
species <- spTransform(species, proj4string(bio))

plot(bio$bio_5)
plot(species, add=T, pch=19, cex=0.5)


ext <- extract(bio, species)
species <- cbind(species, ext)
head(species)
head(ext)

library(psych)
corr <- corr.test(ext, y = NULL, use = "pairwise", method="spearman", adjust="bonferroni", alpha=.05, ci=F)
corr$r


## THIN PLATE WITH REML
M1 <- gam(abundance ~ s(bio_2, fx=F, k=10, bs="tp") +
            s(bio_6, fx=F, k=10, bs="tp") +
            s(bio_12, fx=F, k=10, bs="tp"),
            method="REML", data=species) 
plot(M1)
summary(M1)


P1 <- predict.gam(M1, as.data.frame(bio))

P1 <- melt(P1)

P1r <- raster(nrows=nrow(bio), ncol=ncol(bio), ext=extent(bio), crs=proj4string(bio), vals=P1$value)

plot(P1r)


#############################################################################




##########################################################################
##################### VISUALIZING ADDITIVENESS ###########################

M2 <- gam(abundance ~ s(bio_6, fx=F, k=10, bs="tp") +
                      s(bio_2, fx=F, k=10, bs="tp") , method="REML", data=species)
plot(M2)


bio_6 <- c(min(species$bio_6):max(species$bio_6))bio_2 <- c(92:141)

predictc <- data.frame(bio_6)
predictc$row <- as.numeric(row.names(predictc))
head(predictc)
predictc$bio_2 <- min(species$bio_2) + ((predictc$row - 1)*((max(species$bio_2)-min(species$bio_2))/(max(species$bio_6)-min(species$bio_6))))
#max(predictc$bio_2)

P2c <- predict.gam(M2, predictc)
P2c <- melt(P2c)
plot(P2c$value)

######################################
predict6 <- data.frame(bio_6)
predict6$bio_2 <- 125

P26 <- predict.gam(M2, predict6)
P26 <- melt(P26)

######################################
predict2 <- data.frame(predictc$bio_2)
predict2$bio_6 <- 29
colnames(predict2) <- c("bio_2", "bio_6")

P22 <- predict.gam(M2, predict2)
P22 <- melt(P22)
colnames(P22) <- c("Var2", "value2")

Pfull <- cbind(P26,P22)
Pfull$val3 <- Pfull$value + Pfull$value2
plot(Pfull$val3)

colnames(P2c) <- c("VarC", "val4")
Pfull<- cbind(P2c, Pfull)

par(mfrow=c(1,1))
plot(Pfull$val3~Pfull$val4)
cor.test(Pfull$val3,Pfull$val4)

######################################
#repeat with diff numbers

bio_6 <- c(min(species$bio_6):max(species$bio_6))bio_2 <- c(92:141)

predictc <- data.frame(bio_6)
predictc$row <- as.numeric(row.names(predictc))
head(predictc)
predictc$bio_2 <- min(species$bio_2) + ((predictc$row - 1)*((max(species$bio_2)-min(species$bio_2))/(max(species$bio_6)-min(species$bio_6))))
#max(predictc$bio_2)

P2c <- predict.gam(M2, predictc)
P2c <- melt(P2c)
plot(P2c$value)

######################################
predict6 <- data.frame(bio_6)
predict6$bio_2 <- 100
#plot(M2)
P26 <- predict.gam(M2, predict6)
P26 <- melt(P26)

par(mfrow=c(2,2))
plot(M2)
plot(P26)
######################################
predict2 <- data.frame(predictc$bio_2)
predict2$bio_6 <- 60
colnames(predict2) <- c("bio_2", "bio_6")

P22 <- predict.gam(M2, predict2)
P22 <- melt(P22)
colnames(P22) <- c("Var2", "value2")
plot(P22$value2)

Pfull <- cbind(P26,P22)
Pfull$val3 <- Pfull$value + Pfull$value2
plot(Pfull$val3)

colnames(P2c) <- c("VarC", "val4")
Pfull<- cbind(P2c, Pfull)

par(mfrow=c(1,1))
plot(Pfull$val3~Pfull$val4)
cor.test(Pfull$val3,Pfull$val4)
##########################################################################
#rm(ext,predict2,predict6)#P1,P2,P22,P26,edd,P2b,P2c,Pfull,predict,predictb,predictc,bio_2,bio_6,corr,M1,P1r)

M2 <- gam(abundance ~ s(bio_6, fx=F, k=10, bs="tp") +
            s(bio_2, fx=F, k=10, bs="tp") , method="REML", data=species)
summary(M2)

par(mfrow=c(2,2))
plot(M2)

plot(species$bio_6, species$abundance)
plot(species$bio_2, species$abundance)

## get smooth from model
smooth <- M2$smooth[1]

## knots... NULL
smooth$knots

## basis class... NULL
attr(smooth, "class")

## coefficients; named vector
beta <- M2$coefficients

## model matrix; intercepts?
mat <- predict.gam(M2, type="lpmatrix")

## visualize a basis function.... needs fine tuning
#plot(species$abundance, mat[, "s(x).20"], type = "l", main = "20th basis")

pred.linear <- mat %*% beta ## requires numeric/complex matrix/vector argument
pred.linear


pred.linear <- M2$linear.predictors ## don't know how to use result
plot(pred.linear)

pred.linear <- predict.gam(M2, type="link") ## same output
plot(pred.linear)

str(pred.linear)
head(pred.linear)
class(pred.linear)

?predict.gam

M2$var.summary
summary(M2$weights)
M2$model
M2$coefficients
plot(M2$residuals)
M2$fitted.values

terms <- predict.gam(M2, type="terms")
head(terms)
link <- predict.gam(M2)

################################################################################
#### DETECTING UPWARD TRENDING EXTRAPOLATION AND SIMULATIING A FLATLINE ########
################################################################################

bio.sp <- as.data.frame(bio) ## make copy of bioclimatic variables for manipulation

## if bio_2 is in the formula, then assess upward trends in extrapolations,
## and simulate flatlines if necessary
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