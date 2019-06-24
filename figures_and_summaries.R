### Data summaries, calculations, some figutes
##  code by Mitch O'Neill
#   written: 4/29/2018
#   updated: 5/25/2018


ordsums <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_8_2018b.csv", sep=",", header=T, stringsAsFactors = F)


mean(ordsums$infilling[ordsums$USE==1 & !is.na(ordsums$infilling)])
median(ordsums$infilling[ordsums$USE==1 & !is.na(ordsums$infilling)])
sd(ordsums$infilling[ordsums$USE==1 & !is.na(ordsums$infilling)])
hist(ordsums$infilling[ordsums$USE==1 & !is.na(ordsums$infilling)])
length(ordsums$species.code[ordsums$AUC>.7 & ordsums$MESS_abun2range_SEL < 0.9])
length(ordsums$species.code[ordsums$AUC>.7 & ordsums$MESS_abun2range_SEL > 0.9 & (ordsums$kappaP >0.05 | ordsums$kappa <=0 )])

options(scipen=999)

## not needed in publication
#pdf("C:/Users/Localadmin/Google Drive/figs/kappas.png",width=700,height=400)
#hist(ordsums$kappa[ordsums$USE==1],breaks=20, xlab="Cohen's kappa", pch=2,cex.lab=1.5,col="lightblue", ylab="Number of species", main="")
#dev.off()

## 25.4 mm per inch

pdf("C:/Users/Localadmin/Google Drive/figs/EstImpInf.pdf",width=79/25.4,height=150/25.4)

par(mfrow=c(3,1))
hist(ordsums$areaSQKM_5_7[ordsums$USE==1],breaks=10, ylab="Number of species", xlab=expression("Area" ~ (km^2)), cex.lab=1,col="lightblue", main="",
     xlim=c(0,6000000))
#title("(a)", adj=0, cex.main=1.5)
mtext("(a)", side=3, line=2, adj=-0.25)


hist(ordsums$hiAbunSQKM_5_7[ordsums$USE==1], main="",ylab="Number of species", xlab=expression("Area" ~ (km^2)), cex.lab=1,col="lightblue",
     xlim=c(0,6000000))
#title("(b)", adj=0, cex.main=1.5)
mtext("(b)", side=3, line=2, adj=-0.25)

hist(ordsums$infilling[ordsums$USE==1]*100, main="",xlab="Percent infilling (%)", ylab="Number of species",cex.lab=1,col="lightblue")
mtext("(c)", side=3, line=2, adj=-0.25)
dev.off()



ordsums$hiAbunSQKM_5_7
summary(ordsums)


M <- lm(ordsums$areaSQKM_5_7[ordsums$USE==1]~ordsums$thinned_occ[ordsums$USE==1])
M2 <- lm(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]~ordsums$no.pts[ordsums$USE==1])
summary(M)
summary(M2)



#par(mfrow=c(2,1))
#png("C:/Users/Localadmin/Google Drive/figs/samplesize.png",width=700,height=500)
pdf("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/figs/samplesize.pdf",width=110/25.4,height=85/25.4, bg="transparent")
par(mai=c(.5,.5,.1,1.3))
plot(ordsums$areaSQKM_5_7[ordsums$USE==1]~ordsums$thinned_occ[ordsums$USE==1], pch=20, col="medium blue", cex.lab=0.8, cex.axis= 0.5,
     xlab="Sample size", ylab= expression("Range size"~(km^2)), cex.main=0.001, mgp=c(1.3,.4,0))
points(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]~ordsums$no.pts[ordsums$USE==1], pch="+", cex=1,col="firebrick2")
abline(M$coefficients[1],M$coefficients[2], col="darkblue", lwd=2)
abline(M2$coefficients[1],M2$coefficients[2], col="firebrick3", lwd=2)

points(8300,3200000,pch=19,col="medium blue",xpd=1)
points(8300,2200000,pch="+",col="firebrick2",xpd=1, cex=1.5)

text(8500,3200000, "Establishment", pos=4, xpd=1, cex=.65)
text(8500,3000000, "range", pos=4, xpd=1, cex=.65)
text(8500,2200000, "Impact range", pos=4, xpd=1, cex=.65)

dev.off()




#Possible figure sizes: single column = 79mm, 2/3rd column = 110mm, double column = 168mm, maximum height of figure = 230mm.
M3 <- lm(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]~ordsums$areaSQKM_5_7[ordsums$USE==1])
#bottom, left, top, right

pdf("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/figs/estab_imp2.pdf",width=110/25.4,height=75/25.4, bg="transparent")
#par(mai=c(0.4,0.4,0.1,0.1))#, mar=c(0,0,0,0))
# bottom, left, top, right
par(mar=c(2.5,2.5,1,1))
options(scipen = 999)
plot(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]~ordsums$areaSQKM_5_7[ordsums$USE==1],pch=20, cex=1, cex.axis= 0.5, cex.lab=0.8,
     xlab= expression("Establishment range size"~(km^2)), ylab= expression("Impact range size"~(km^2)), bg="transparent", mgp=c(1.4,.4,0))
abline(M3$coefficients[1],M3$coefficients[2], col="blue", lwd=1) 

summary(M3)
text(3500000,1000000,"y = 0.14x + 41573.62", col="dark blue",pos=4, xpd=1, cex=.7)
text(3500000,870000,"adjusted R  = 0.22", col="dark blue",pos=4, xpd=1, cex=.7)
text(4350000,900000, "2", col="dark blue", pos=4, xpd=1,cex=0.3)



dev.off()
#summary(M3)

#### CORRELATION between full range and abundance range INSTEAD ######
ordsums <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_31_2018.csv", sep=",", header=T, stringsAsFactors = F)
hist(ordsums$hiAbunSQKM_5_7[ordsums$USE==1])
hist(ordsums$areaSQKM_5_7[ordsums$USE==1])

hist(log(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]))
hist(log(ordsums$areaSQKM_5_7[ordsums$USE==1]))
## non-normal distribution, variances are not the same?

ordsums70 <- ordsums[ordsums$USE==1,]
## there are ties in data (can't calculate p-value for spearman)
cor.test(ordsums70$areaSQKM_5_7,ordsums70$hiAbunSQKM_5_7, method="kendall")
## tau = 0.46
## p = 1.905e * 10^-8 = .00000001905

pdf("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/figs/estab_imp2_correlation.pdf",width=110/25.4,height=75/25.4, bg="transparent")
#par(mai=c(0.4,0.4,0.1,0.1))#, mar=c(0,0,0,0))
# bottom, left, top, right
par(mar=c(2.5,2.5,1,1))
options(scipen = 999)
plot(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]~ordsums$areaSQKM_5_7[ordsums$USE==1],pch=20, cex=1, cex.axis= 0.5, cex.lab=0.8,
     xlab= expression("Establishment range size"~(km^2)), ylab= expression("Impact range size"~(km^2)), bg="transparent", mgp=c(1.4,.4,0))

text(3500000,1000000,"Kendall's tau = 0.46", col="dark blue",pos=4, xpd=1, cex=.7)
text(3500000,870000,"p-value = 1.9 x 10", col="dark blue",pos=4, xpd=1, cex=.7)
text(4930000,900000, "-8", col="dark blue", pos=4, xpd=1,cex=0.4)

dev.off()
#plot(ordsums$hiAbunSQKM_5_7[ordsums$USE==1]~ordsums$areaSQKM_5_7[ordsums$USE==1],pch=19, cex.lab=2,
#     xlab= "Establishment range size (square km)", ylab= "Impact range size (square km)", cex.axis=2, bg="transparent")
#points(ordsums$hiAbunSQKM_5_7[ordsums$species.code=="ALJU"]~ordsums$areaSQKM_5_7[ordsums$species.code=="ALJU"], col="orange", pch=16)
#points(ordsums$hiAbunSQKM_5_7[ordsums$species.code=="PUMOL"]~ordsums$areaSQKM_5_7[ordsums$species.code=="PUMOL"], col="red", pch=16)
#points(ordsums$hiAbunSQKM_5_7[ordsums$species.code=="ONAC"]~ordsums$areaSQKM_5_7[ordsums$species.code=="ONAC"], col="green", pch=16)

########################## TRAIT COMPOSITION ######################
all <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/master list.csv" , header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
summary(as.factor(all$L48_present))
summary(as.factor(all$list))
summary(all$Number.of.States)

all <- all[all$L48_present =="YES" & (
  (all$Federal.Noxious== 1 & !is.na(all$Federal.Noxious)) |
    ((all$list == "IPAUS" | all$list == "IPAUS_USDA") & !is.na(all$list)) |
    all$Number.of.States > 0), ]
## 1089
all <- all[all$habit != "aquatic",] #1088
all <- all[!((!is.na(all$HI) | !is.na(all$AK)) & all$Number.of.States==1 & (all$list == "USDA" | all$list == "") & is.na(all$Federal.Noxious==1)),]
## 1078

all$habit[all$NewCode=="TRSE"] <- "tree" ## from USDA PLANTS database
#summary(as.factor(all$habit))
all <- all[ !is.na(all$habit) & all$habit != "",] ## 1045
#summary(as.factor(all$duration_final))
#hist(as.numeric(as.factor(all$duration_final)))
#1052 with trait data
write.table(all, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/all_contig.csv", sep=",", row.names=F)
traits <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/all_contig.csv" , header = T, sep = ",", quote= "\"", 
                     comment.char= "", stringsAsFactors = F, strip.white = T)


traits$common.name[!(traits$NewCode %in% ordsums$species.code[ordsums$USE==1])]
traits$common.name[(traits$NewCode %in% ordsums$species.code[ordsums$USE==0 & (ordsums$habit=="forb_herb" | ordsums$habit =="grass") & ordsums$duration_final !="Perennial" & (ordsums$kappaP >0.05 | ordsums$kappa < 0)])]
## perennial pepperweed
traits$habit[traits$common.name=="poison hemlock, poison-hemlock" ]
traits$common.name[(traits$NewCode %in% ordsums$species.code[ordsums$USE==0 & (ordsums$habit=="forb_herb" | ordsums$habit =="grass") & ordsums$duration_final =="Perennial" & (ordsums$kappaP >0.05 | ordsums$kappa < 0)])]

ordsums <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_31_2018.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)
head(ordsums)
str(ordsums)
summary(ordsums)
hist(ordsums$AUC)
55/64
summary(as.factor(ordsums$duration[ordsums$USE==1]))
summary(as.factor(ordsums$habit[ordsums$USE==1]))/70

#traits <- read.csv("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/all.csv", header=T,stringsAsFactors = F)
ordsums$habit <- ordered(ordsums$habit, levels=c("forb_herb","grass", "vine","tree", "shrub_subshrub", "MULTIPLE", "palms", "cactus"))
traits$habit <- ordered(traits$habit, levels=c("forb_herb","grass", "vine","tree", "shrub_subshrub", "MULTIPLE", "palms", "cactus"))
#summary(as.factor(ordsums$USE))

habits <- table(traits$habit)
habits <- habits/nrow(traits)


habits2 <- table(ordsums$habit[!is.na(ordsums$aic)])
habits2 <- habits2/nrow(ordsums[!is.na(ordsums$aic),])
habits3 <- table(ordsums$habit[!is.na(ordsums$aic) & ordsums$USE==1])
habits3 <- habits3/nrow(ordsums[!is.na(ordsums$aic) & ordsums$USE==1,])

habits<- (rbind(habits,habits2,habits3))
colnames(habits) <- c("forb", "grass", "vine", "tree", "shrub", "multiple", "palm", "cactus")

ordsums$duration_final <- ordsums$duration
ordsums$duration_final[ordsums$duration_final != "Perennial" & ordsums$duration_final != "Annual" & ordsums$duration_final != "Biennial"] <- "Multiple"
traits$duration_final[traits$duration_final != "Perennial" & traits$duration_final != "Annual" & traits$duration_final != "Biennial"] <- "Multiple"


traits$duration_final <-ordered(traits$duration_final, 
                                levels=c("Perennial", "Biennial", "Annual", "Multiple"))
ordsums$duration_final <-ordered(ordsums$duration_final, 
                                 levels=c("Perennial", "Biennial", "Annual", "Multiple"))
durations <- table(traits$duration_final)
durations <- durations/nrow(traits)

durations2 <- table(ordsums$duration_final[!is.na(ordsums$aic)])
durations2 <- durations2/nrow(ordsums[!is.na(ordsums$aic),])
durations3 <- table(ordsums$duration_final[!is.na(ordsums$aic) & ordsums$USE==1])
durations3 <- durations3/nrow(ordsums[!is.na(ordsums$aic) & ordsums$USE==1,])

durations <- rbind(durations,durations2,durations3)
colnames(durations) <- c( "perennial", "biennial", "annual", "multiple")

colnames(habits) <- c("","","","","","","","")
colnames(durations) <- c("","","","" )



pdf("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/figs/traits.pdf",width=110/25.4,height=110/25.4)
par(mfrow=c(2,1), mai=c(0.3,0.7,0.49,0.2))

barplot(habits*100, main="", ylim=c(0,100), #cex=.75,
        xlab="Habit", ylab= " ", cex.axis=.70,
        col=c("white","lightblue","mediumblue"), beside=TRUE, #bty="L", #mar=c(1,1,1,1), oma=c(3,5,0,0),
        legend = c("all species","data-sufficient species", "species post-screening"), 
        args.legend = list(x="topright", bty="n", cex=.8),
        axes=T)
abline(h=0)
text(-5, 125,"(a)",xpd=1, cex=1.2)
text(2.4,-12,"Forb",xpd=1, cex=0.75)
text(6.5,-12,"Grass",xpd=1, cex=0.75)
text(10.6,-12,"Vine",xpd=1, cex=0.75)
text(14.4,-12,"Tree",xpd=1, cex=0.75)
text(18.4,-12,"Shrub",xpd=1, cex=0.75)
text(22.7,-12.5,"Multiple",xpd=1, cex=0.75) #22.9 is too high
text(26.7,-12,"Palm",xpd=1, cex=0.75)
text(30.5,-12,"Cactus",xpd=1, cex=0.75)
mtext("Percent of species (%)", side=2,line=2, at=c(50))


barplot(durations*100, main="" , ylim=c(0,100),
        xlab="Duration", ylab= " ", cex.axis=.70, cex.main=2,
        col=c("white","lightblue","mediumblue"), beside=TRUE, #bty="L", #mar=c(1,1,1,1), oma=c(3,5,0,0),
        legend = c("all species","data-sufficient species", "species post-screening"), 
        args.legend = list(x="topright", bty="n",cex=.8))
abline(h=0)
text(-2,125,"(b)",xpd=1, cex=1.2)
mtext("Percent of species (%)", side=2,line=2, at=c(50))
text(2.5,-12,"Perennial",xpd=1, cex=0.75)
text(6.5,-12,"Biennial",xpd=1, cex=0.75)
text(10.5,-12,"Annual",xpd=1, cex=0.75)
text(14.5,-12,"Multiple",xpd=1, cex=0.75)
dev.off()
## 1045 out of 1078

#######################################


################### model screening ####
ordsums2 <- ordsums[ordsums$AUC >0.7,]
146/155 # 94%
ordsums2 <- ordsums2[ordsums2$MESSgn10 > .90,]
85/146 # 58%
ordsums2 <- ordsums2[ordsums2$kappaP < 0.05 & ordsums2$kappa > 0,]
64/85 # 75%
(146-85)/146
mean(ordsums2$AUC)
sd(ordsums2$AUC)

summary(ordsums2$hiAbunSQKM)
ordsums2[ordsums2$class3 == 0,]
ordsums2[ordsums2$class1 == 0,]
ordsums2[ordsums2$class2 == 0,]


(155-146)/155
(146-85)/155
(85-64)/155

length(ordsums$species.code[ordsums$USE ==1 & ordsums$kappa > 0.2])

26/64


###############infilling series###########


inf<-read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/infilling.csv",  header = T, sep = ",", quote= "\"", 
                comment.char= "", stringsAsFactors = F, strip.white = T)
inf <- inf[inf$kappaP<0.05 & inf$AUC>0.7 & inf$kappa>0 & inf$MESSgn10>.9,]
cbind(inf$species.code,ordsums2$species.code)
ordsums2 <- cbind(ordsums2,inf$infilling)
hist(ordsums2$`inf$infilling`)


png("C:/Users/Localadmin/Google Drive/GRC_figs/infil.png",width=700,height=300)
hist(ordsums2$`inf$infilling`, main="",xlab="Porportion of cells occupied", ylab="Number of species",pch=2,cex.lab=1.5,col="lightblue")
dev.off()


M<-lm(ordsums2$hiAbun_prop_range~ordsums2$areaSQKM)
summary(M)

M<-lm(ordsums2$hiAbunSQKM~ordsums2$areaSQKM)
summary(M)
plot(ordsums2$hiAbunSQKM~ordsums2$areaSQKM)
abline(b=.1421,a= 41150)

ordsums$duration[ordsums$duration=="filler"] <- "Perennial"

M <- lm(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration) + as.factor(ordsums2$habit) + as.factor(ordsums2$duration)*as.factor(ordsums2$habit))
summary(M)

anova <- aov(ordsums2$hiAbunSQKM~as.factor(ordsums2$habit))
anova

anova2 <- aov(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration))
summary(anova2)

plot(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration))
plot(ordsums2$hiAbunSQKM~as.factor(ordsums2$habit))

ordsums2$perenn[ordsums2$duration == "Perennial"] <- "YES"
ordsums2$perenn[ordsums2$duration != "Perennial"] <- "NO"

summary(as.factor(ordsums2$perenn))

anova <- aov(ordsums2$hiAbunSQKM~as.factor(ordsums2$perenn))
summary(anova)




summary(ordsums2$habit)
ordsums2$treegrass[as.character(ordsums2$habit) == "tree"|
                     as.character(ordsums2$habit) == "shrub_subshrub"] <- "YES"
ordsums2$treegrass[as.character(ordsums2$habit) == "grass"|
                     as.character(ordsums2$habit) == "forb_herb"] <- "NO"
#anova <- aov(ordsums2$areaSQKM[!is.na(ordsums2$treegrass)]~as.factor(ordsums2$treegrass[!is.na(ordsums2$treegrass)]))
anova <- aov(ordsums2$areaSQKM[!is.na(ordsums2$treegrass)]~as.factor(ordsums2$treegrass[!is.na(ordsums2$treegrass)]))
summary(anova)
plot(ordsums2$areaSQKM~as.factor(ordsums2$treegrass))


summary(aov(ordsums2$areaSQKM~ordsums2$habit))
summary(aov(ordsums2$areaSQKM~ordsums2$perenn))
plot(ordsums2$areaSQKM~as.factor(ordsums2$perenn))


length(ordsums2$kappa[ordsums2$kappa <0.2])
length(ordsums2$kappa[ordsums2$kappa <0.4 & ordsums2$kappa >0.2])
length(ordsums2$kappa[ordsums2$kappa <0.6 & ordsums2$kappa >0.4])      
length(ordsums2$kappa[ordsums2$kappa <0.8 & ordsums2$kappa >0.6])      
length(ordsums2$kappa[ordsums2$kappa >0.8 & ordsums2$kappa >0.6])
38/64
22/64
2/64  
1/64
38+22+2+1



###################### ECOREGION COMPOSITION #####
eco.hot <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ecoXhot_5_8.csv", sep=",", header=T)



eco.hot <- as.matrix(eco.hot[,2:NCOL(eco.hot)])
head(eco.hot)
colnames(eco.hot) <- c("","","","","","","","","","")
## 900 500


pdf("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/figs/EcoregionComp2.pdf",width=168/25.4,height=200/25.4)
par(mar = c(14, 5, 2, 2))
barplot(eco.hot*100, main=" ", ylim=c(0,100), cex.lab=1.25, #cex=2,
        ylab= "Percent of area in each ecoregion (%)", xlab= "", cex.axis=1, 
        col=c("white","lightblue","firebrick3"), beside=TRUE,#, #bty="L", #mar=c(1,1,1,1), oma=c(3,5,0,0),
        legend = c("Coterminuous U.S.","Establishment hotspots","Abundance hotspots"), 
        args.legend = list(x="topright", bty="n", cex=1.4)#,
        #args.axis = list(side=1, las=2)
)


axis(side=1, tck=0, at=c(-1,40), label=c("", ""))

abline(h=0, lwd=2)
#abline(v=0, lwd=2)

mtext("Eastern Temperate Forests   ", side=1, las=2, at=c(2.35,-40), cex=1.2)
mtext("Great Plains   ", side=1, las=2, at=c(6.35,-40), cex=1.2)
mtext("North American Deserts   ", side=1, las=2, at=c(10.35,-40), cex=1.2)
mtext("NW Forested Mountains   ", side=1, las=2, at=c(14.35,-40), cex=1.2)
mtext("Northern Forests   ", side=1, las=2, at=c(18.35,-40), cex=1.2)
mtext("Mediterranean California   ", side=1, las=2, at=c(22.35,-40), cex=1.2)
mtext("Temperate Sierras   ", side=1, las=2, at=c(26.35,-40), cex=1.2)
mtext("Marine West Coast Forests   ", side=1, las=2, at=c(30.35,-40), cex=1.2)
mtext("Southern Semiarid Highlands   ", side=1, las=2, at=c(34.35,-40), cex=1.2)
mtext("Tropical Wet Forests   ", side=1, las=2, at=c(38.35,-40), cex=1.2)

dev.off()






eco.hot <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/hot_eco.csv", sep=",", header=T)
eco.hot <- rbind(eco.hot$IDs[1:10],eco.hot$hotspot_FULL[1:10],eco.hot$hotspot_HI_ABUN[1:10])
eco.hot <- cbind(eco.hot[,1:2], eco.hot[,5], eco.hot[,7:6], eco.hot[,4], eco.hot[,9], eco.hot[,3], eco.hot[,8], eco.hot[,10])
eco.hot <- eco.hot[2:3,]

eco.hot1 <- c(.32,.291,.184,.11,.041,.021,.013,.010,.008,.002)
eco.hot1b <- c(.8,.045,.035,.025,.12,.001,.0008,.008,.0001,0)
sum(eco.hot1b)
eco.hot1c <- c(.48, .18,.07,.17,.02,.03,.01,.04,.0001,0)
sum(eco.hot1c)

eco.hot2 <- rbind(eco.hot1,eco.hot1b,eco.hot1c)
eco.hot2







####################################
###### making a pretty table
####################################

ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_31_2018.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)

#infilling <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/infilling.csv", header = T, sep = ",", quote= "\"", 
#                         comment.char= "", stringsAsFactors = F, strip.white = T)
#leaveoneout <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_1_2018.csv", header = T, sep = ",", quote= "\"", 
#                          comment.char= "", stringsAsFactors = F, strip.white = T)
#samplesizes<- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/summaries_copy_4_07/sample_sizes.csv", header = T, sep = ",", quote= "\"", 
#                         comment.char= "", stringsAsFactors = F, strip.white = T)


edd <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd$abundance <- ordered(as.factor(edd$abundance), levels=c(1,2,3))



all <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/master list.csv" , header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)

#for (i in 1:155){
#  if(ordsums$AUC[i] > .7 & ordsums$kappa[i] > 0 & ordsums$kappaP[i] <0.05 & ordsums$MESS_abun2range_SEL[i] > .9){
#    ordsums$USE[i] <- 1
#  } else {
#    ordsums$USE[i] <- 0 
#  }
#}

summary(as.factor(ordsums$USE))
#l1o2 <- ordsums$l1o2[!is.na(ordsums$l1o2)]

#ordsums <- rbind(ordsums[ordsums$USE==1,],ordsums[ordsums$USE==0,])
#ordsums <- ordsums[order(ordsums$species.code),]
#ordsums$l1o2[ordsums$USE==1]


tab <- data.frame(ordsums$species.code, stringsAsFactors = F)
head(tab)
colnames(tab)<-"code"

#tab$habit <- "BLOOP"
#tab$duration <- "BLOOP"
tab$species <- "BLOOP"
tab$occ_cells <- -999
tab$abun_cells <- -999
tab$AUC <- " "
tab$mess <- " "
tab$kappa <- " "

tab$b2 <- " "
tab$b5 <- " "
tab$b6 <- " "
tab$b8 <- " "
tab$b12 <-" "
tab$b15 <-" "
tab$n3 <- " "
tab$n4 <- " "
tab$n5 <- " "
tab$n6 <- " "
tab$n7 <- " "
tab$n8 <- " "


tab$rangekm <- " "
tab$rangeP <- " "
tab$impackm <- " "
tab$impacP <- " "
tab$infilling <- " "
tab$l1o <- " "

library(ordinal)
library(glue)

for (i in 1:length(tab$code)){
  tab$species[i] <- all$Accepted.Latin.Name[all$NewCode == tab$code[i]]
  tab$occ_cells[i] <- ordsums$thinned_occ[ordsums$species.code==tab$code[i]]
  tab$abun_cells[i] <- length(edd$species[edd$species== tab$code[i]])
  tab$AUC[i] <-   as.character(round(ordsums$AUC[ordsums$species.code == tab$code[i]],3))
  if(as.numeric(tab$AUC[i])<0.7){
    tab$AUC[i] <- paste0(tab$AUC[i],"^")
  }
  
  tab$kappa[i] <- as.character(round(ordsums$kappa[ordsums$species.code == tab$code[i]],2))
  if(as.numeric(tab$kappa[i]) <0.01){
    tab$kappa[i] <- "<0.01"
  }
  
  
  if(ordsums$kappaP[ordsums$species.code == tab$code[i]] < .1 & !is.na(ordsums$kappaP[ordsums$species.code == tab$code[i]])){
    tab$kappa[i] <- paste0(tab$kappa[i], "*")
    if(ordsums$kappaP[ordsums$species.code==tab$code[i]] < 0.05 & !is.na(ordsums$kappaP[ordsums$species.code == tab$code[i]])){
      tab$kappa[i] <- paste0(tab$kappa[i], "*")
    } 
  } 
  
  if(ordsums$kappaP[ordsums$species.code == tab$code[i]] > .05 | is.na(ordsums$kappaP[ordsums$species.code == tab$code[i]])){
    tab$kappa[i] <- paste0(tab$kappa[i], "^")
  }
  
  
  #  tab$habit[i] <- all$habit[all$NewCode==tab$code[i]]
  #  tab$duration[i] <- all$duration_final[all$NewCode==tab$code[i]]
  
  data <- edd[edd$species==tab$code[i],]
  
  M <- clm(formula(ordsums$formu[ordsums$species.code==tab$code[i]]), data=data, Hess=T)
  msum <- as.data.frame(coef(summary(M)))
  msum<- msum[3:NROW(msum),]
  #vars <- collapse(row.names(msum)[3:NROW(msum)],sep=" ")
  #sapply(formulae$formula, grepl, pattern="bio_2", fixed=T)
  msum$var <- row.names(msum)
  
  msum$b2 <- sapply(msum$var, grepl, pattern="bio_2", fixed=T)
  msum$b5 <- sapply(msum$var, grepl, pattern="bio_5", fixed=T) 
  msum$b6 <- sapply(msum$var, grepl, pattern="bio_6", fixed=T) 
  msum$b8 <- sapply(msum$var, grepl, pattern="bio_8", fixed=T) 
  msum$b12 <- sapply(msum$var, grepl, pattern="bio_12", fixed=T) 
  msum$b15 <- sapply(msum$var, grepl, pattern="bio_15", fixed=T) 
  msum$n3 <- sapply(msum$var, grepl, pattern="nlcd_3", fixed=T) 
  msum$n4 <- sapply(msum$var, grepl, pattern="nlcd_4", fixed=T) 
  msum$n5 <- sapply(msum$var, grepl, pattern="nlcd_5", fixed=T) 
  msum$n6 <- sapply(msum$var, grepl, pattern="nlcd_6", fixed=T) 
  msum$n7 <- sapply(msum$var, grepl, pattern="nlcd_7", fixed=T) 
  msum$n8 <- sapply(msum$var, grepl, pattern="nlcd_8", fixed=T) 
  
  if(NROW(msum[msum$b2==T,]) == 1){
    if(msum$Estimate[msum$b2==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$b2==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$b2==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$b2[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$b2==T,]) == 2){
    if(msum$Estimate[msum$b2==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$b2==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$b2==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$b2[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$b5==T,]) == 1){
    if(msum$Estimate[msum$b5==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$b5==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$b5==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$b5[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$b5==T,]) == 2){
    if(msum$Estimate[msum$b5==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$b5==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$b5==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$b5[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$b6==T,]) == 1){
    if(msum$Estimate[msum$b6==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$b6==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$b6==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$b6[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$b6==T,]) == 2){
    if(msum$Estimate[msum$b6==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$b6==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$b6==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$b6[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$b8==T,]) == 1){
    if(msum$Estimate[msum$b8==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$b8==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$b8==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$b8[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$b8==T,]) == 2){
    if(msum$Estimate[msum$b8==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$b8==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$b8==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$b8[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$b12==T,]) == 1){
    if(msum$Estimate[msum$b12==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$b12==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$b12==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$b12[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$b12==T,]) == 2){
    if(msum$Estimate[msum$b12==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$b12==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$b12==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$b12[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$b15==T,]) == 1){
    if(msum$Estimate[msum$b15==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$b15==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$b15==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$b15[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$b15==T,]) == 2){
    if(msum$Estimate[msum$b15==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$b15==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$b15==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$b15[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$n3==T,]) == 1){
    if(msum$Estimate[msum$n3==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$n3==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$n3==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$n3[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$n3==T,]) == 2){
    if(msum$Estimate[msum$n3==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$n3==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$n3==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$n3[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$n4==T,]) == 1){
    if(msum$Estimate[msum$n4==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$n4==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$n4==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$n4[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$n4==T,]) == 2){
    if(msum$Estimate[msum$n4==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$n4==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$n4==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$n4[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$n5==T,]) == 1){
    if(msum$Estimate[msum$n5==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$n5==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$n5==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$n5[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$n5==T,]) == 2){
    if(msum$Estimate[msum$n5==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$n5==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$n5==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$n5[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$n6==T,]) == 1){
    if(msum$Estimate[msum$n6==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$n6==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$n6==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$n6[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$n6==T,]) == 2){
    if(msum$Estimate[msum$n6==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$n6==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$n6==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$n6[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$n7==T,]) == 1){
    if(msum$Estimate[msum$n7==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$n7==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$n7==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$n7[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$n7==T,]) == 2){
    if(msum$Estimate[msum$n7==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$n7==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$n7==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$n7[i] <- paste0(sign,sig)
  }
  
  if(NROW(msum[msum$n8==T,]) == 1){
    if(msum$Estimate[msum$n8==T]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "L"
    if(msum$`Pr(>|z|)`[msum$n8==T]< 0.10){
      sig <- "L*"
      if(msum$`Pr(>|z|)`[msum$n8==T] <0.05){
        sig <- "L**"
      }
    }
    
    tab$n8[i] <- paste0(sign,sig)
  }
  if(NROW(msum[msum$n8==T,]) == 2){
    if(msum$Estimate[msum$n8==T][2]>0){
      sign <- "+"
    } else{
      sign <- "-"
    }
    
    sig <- "Q"
    if(msum$`Pr(>|z|)`[msum$n8==T][2]< 0.10){
      sig <- "Q*"
      if(msum$`Pr(>|z|)`[msum$n8==T][2] <0.05){
        sig <- "Q**"
      }
    }
    
    tab$n8[i] <- paste0(sign,sig)
  }
  
  
  tab$mess[i]      <- as.character(round((100-ordsums$MESSgn10[ordsums$species.code==tab$code[i]]*100), 1))
  
  if(as.numeric(tab$mess[i]) < 0.1 & as.numeric(tab$mess[i]) > 0) {
    tab$mess[i] <- "<0.1"
  }
  
  
  if( as.numeric(tab$mess[i]) > 10){
    tab$mess[i] <- paste0(tab$mess[i],"^")
  }
  
  
  
  tab$rangekm[i]   <- as.character(round(ordsums$areaSQKM_5_7[ordsums$species.code==tab$code[i]]/1000000, 1))
  if(as.numeric(tab$rangekm[i]) < 0.1 & as.numeric(tab$rangekm[i]) > 0){
    tab$rangekm[i] <- "<0.1"
  }
  
  tab$rangeP[i]    <- as.character(round(ordsums$estab_US[ordsums$species.code==tab$code[i]]*100, 1))
  if(as.numeric(tab$rangeP[i]) < 0.1 & as.numeric(tab$rangeP[i]) > 0){
    tab$rangeP[i] <- "<0.1"
  }
  
  tab$impackm[i]   <- as.character(round(ordsums$hiAbunSQKM_5_7[ordsums$species.code==tab$code[i]]/1000000, 1))
  if(as.numeric(tab$impackm[i]) < 0.1 & as.numeric(tab$impackm[i]) > 0){
    tab$impackm[i] <- "<0.1"
  }
  
  tab$impacP[i]    <- as.character(round(ordsums$impac_estab[ordsums$species.code==tab$code[i]]*100, 1))
  if(as.numeric(tab$impacP[i]) < 0.1 & as.numeric(tab$impacP[i]) > 0){
    tab$impacP[i] <- "<0.1"
  }
  
  tab$infilling[i] <- as.character(round(ordsums$infilling[ordsums$species.code==tab$code[i]]*100, 1))
  if(ordsums$infilling[ordsums$species.code==tab$code[i]]*100 < 0.1 & ordsums$infilling[ordsums$species.code==tab$code[i]]*100 > 0 & !is.na(ordsums$infilling[ordsums$species.code==tab$code[i]])){
    tab$infilling[i] <- "<0.1"
  }
  if(i<71){
    tab$l1o[i]       <- as.character(round(ordsums$l1o2[ordsums$species.code==tab$code[i]]*100, 1))
    
  }
  print(i)
}


str(tab)

tab$species
summary(as.factor(tab$b2[1:70]))
summary(as.factor(tab$b5[1:70]))
summary(as.factor(tab$b6[1:70]))
summary(as.factor(tab$b8[1:70]))
summary(as.factor(tab$b12[1:70]))
summary(as.factor(tab$b15[1:70]))
summary(as.factor(tab$n3[1:70]))
summary(as.factor(tab$n4[1:70]))
summary(as.factor(tab$n5[1:70]))
summary(as.factor(tab$n6[1:70]))
summary(as.factor(tab$n7[1:70]))
summary(as.factor(tab$n8[1:70]))


#tab$infilling[tab$infilling==0]

write.csv(tab, "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/prettyTable_5_31.csv", row.names=F)



summary(ordsums$l1o2)

getwd()
summary(tab[1:64,])
summary(tab$abun_cells/tab$occ_cells)

tab[is.na(tab$infilling),]

tab$habit[tab$code=="TRSE6"] <- "tree"
summary(as.factor(tab$duration))
summary(as.factor(tab$habit))
mean(tab$rangekm*(1/tab$rangeP))
tab$trait_combo <- paste(tab$duration, tab$habit, sep=" ")
summary(as.factor(tab$trait_combo))
summary(as.factor(tab$trait_combo[1:64]))
summary(as.factor(tab$duration[1:64]))
summary(ordsums$USE[1:64])


#######################################################
#### check residuals and correlation b/t prob and class.. anova?
#######################################################
library(ordinal)
library(irr)
library(reshape2)
ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/old_ordsums/model_summaries_4_11_2018.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)



for (i in 1:155){
  if(ordsums$AUC[i] > .7 & ordsums$kappa[i] > 0 & ordsums$kappaP[i] <0.05 & ordsums$MESSgn10[i] > .9){
    ordsums$USE[i] <- 1
  } else {
    ordsums$USE[i] <- 0 
  }
}


edd <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv", header = T, sep = ",", quote= "\"", 
                  comment.char= "", stringsAsFactors = F, strip.white = T)
edd$abundance <- ordered(as.factor(edd$abundance), levels=c(1,2,3))



for(i in 1:length(ordsums$species.code))
  
  
  data <- edd[edd$species==ordsums$species.code[i],]

M <- clm(formula(ordsums$formu[i]), data=data, Hess=T)

pm <- predict(M, data, type="cum.prob")
pm <- melt(pm[1])
pm <- pm[,1]

summary(pm)
#length(pm$)

pm2 <- predict(M, data, type="class")
pm2 <- melt(pm2)
pm2 <- pm2[,1]
cbind(pm,pm2)
#plot(pm$predicted[1:138]~pm2)
plot(pm[1:138]~pm2)
plot(pm~data$abundance, xlab="observed class", ylab="probability of class 3")  

pm <- data.frame(cbind(data$abundance,pm))
colnames(pm) <- c("observed", "predicted")

plot(pm$predicted~as.factor(pm$observed))


kappa <- kappa2(pm)
ordsums$kappaP[s] <- kappa$p.value
ordsums$kappa[s] <- kappa$value



###################################
###################################
varsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/variable_sums.csv", sep=",", header=T)
#eco.hot <- as.data.frame(eco.hot)
#eco.hot <- as.matrix(eco.hot)
varsums <- as.matrix(varsums)
head(varsums)
colnames(varsums)
row.names(varsums) <- varsums[,1]
varsums<- varsums[,2:13]
varsums[,1:12] <- as.numeric(varsums[,1:12])
varsums

#png("C:/Users/Localadmin/Google Drive/GRC_figs/ecoXhot.png",width=900,height=500)
barplot(varsums, main=" ", ylim=c(0,50),#cex=2,
        xlab="Predictor variable", ylab= "Frequency in best-fit models", cex.lab=1.5, cex.axis=1.5, cex.main=2,
        col=c("lightblue","blue","white","firebrick1"), #bty="L", #mar=c(1,1,1,1), oma=c(3,5,0,0),
        legend = c("negative linear","positive linear", "concave-down quadratic", "concave-up quadratic"), 
        args.legend = list(x="topright", bty="n"))
dev.off()

#####################
#### traits analysis
#####################


ordsums <- read.table("C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/model_summaries_5_6_2018.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)

for (i in 1:155){
  if(ordsums$AUC[i] > .7 & ordsums$kappa[i] > 0 & ordsums$kappaP[i] <0.05 & ordsums$MESS_abun2range_SEL[i] > .9){
    ordsums$USE[i] <- 1
  } else {
    ordsums$USE[i] <- 0 
  }
}

summary(as.factor(ordsums$USE))


ordsums2 <- ordsums[ordsums$USE==1,]
summary(as.factor(ordsums2$habit))
summary(as.factor(ordsums2$duration_final))
plot(ordsums2$hiAbunSQKM~as.factor(ordsums2$habit))
plot(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration))
ordsums2$duration2[ordsums2$duration_final=="Perennial"] <- "Perennial"
ordsums2$duration2[ordsums2$duration_final!="Perennial"] <- "Non-Perennial"
plot(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration2))
plot(ordsums2$areaSQKM~as.factor(ordsums2$duration2))


aov1 <- aov(ordsums2$areaSQKM~as.factor(ordsums2$duration2))
aov2 <- aov(ordsums2$areaSQKM~as.factor(ordsums2$habit))

aov3 <- aov(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration2))
aov4 <- aov(ordsums2$hiAbunSQKM~as.factor(ordsums2$habit))

aov5 <- aov(ordsums2$hiAbun_prop_range~as.factor(ordsums2$duration2))
aov6 <- aov(ordsums2$hiAbun_prop_range~as.factor(ordsums2$habit))



#head(ordsums2$duration2)
summary(aov1) ## signif
summary(aov2) ## signif
summary(aov3)
summary(aov4)
summary(aov5)
summary(aov6)

aov2$coefficients


Tuk2 <- TukeyHSD(aov2)
Tuk2

mean(ordsums2$areaSQKM[ordsums2$habit == "MULTIPLE"])
mean(ordsums2$areaSQKM[ordsums2$habit == "forb_herb"])
mean(ordsums2$areaSQKM[ordsums2$habit == "shrub_subshrub"])
mean(ordsums2$areaSQKM[ordsums2$habit == "grass"])
mean(ordsums2$areaSQKM[ordsums2$habit == "tree"])
mean(ordsums2$areaSQKM[ordsums2$habit == "vine"])

ordsums2$habit2 <- ordered(ordsums2$habit, levels=c( "forb_herb", "MULTIPLE", "shrub_subshrub", "grass", "tree", "vine"))
Tuk2
plot(Tuk2)

## signif diffs: tree + forb
##               vine + forb
##               vine + multiple

plot(ordsums2$areaSQKM~as.factor(ordsums2$duration2), ylim=c(0,6e+06),
     xlab="range size (km)", ylab="duration")     

plot(ordsums2$areaSQKM~as.factor(ordsums2$habit2), ylim=c(0,6e+06),
     xlab="range size (km)", ylab="habit")
points(.95,5.5e+06, pch="a", col="black")
points(1.05,5.55e+06, pch="b", col="black")

points(1.88,5.5e+06, pch="a", col="black")
points(1.98,5.55e+06, pch="b", col="black")
points(2.06,5.5e+06, pch="c", col="black")
points(2.16,5.55e+06, pch="d", col="black")

points(2.80,5.5e+06, pch="a", col="black")
points(2.90,5.55e+06, pch="b", col="black")
points(3.0,5.5e+06, pch="c", col="black")
points(3.10,5.55e+06, pch="d", col="black")
points(3.20,5.5e+06, pch="e", col="black")

points(3.80,5.5e+06, pch="a", col="black")
points(3.90,5.55e+06, pch="b", col="black")
points(4.0,5.5e+06, pch="c", col="black")
points(4.10,5.54e+06, pch="d", col="black")
points(4.20,5.5e+06, pch="e", col="black")

points(4.88,5.5e+06, pch="b", col="black")
points(4.98,5.475e+06, pch="c", col="black")
points(5.08,5.5e+06, pch="d", col="black")
points(5.18,5.475e+06, pch="e", col="black")

points(5.98,5.475e+06, pch="c", col="black")
points(6.08,5.5e+06, pch="d", col="black")
points(6.18,5.475e+06, pch="e", col="black")
#points(5.05,3.25e+06, pch="b", col="black")



#Tuk4 <- TukeyHSD(aov4)
#Tuk4 
#summary(aov4)


plot(ordsums2$areaSQKM~as.factor(ordsums2$duration2))
plot(ordsums2$hiAbunSQKM~as.factor(ordsums2$duration2))
plot(ordsums2$hiAbun_prop_range~as.factor(ordsums2$duration2))


citation("dismo")




##############################

library(rgdal)
states <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/states", layer="US_states")
proj4string(states)

states <- spTransform(states,"+init=epsg:5070")

png("C:/Users/Localadmin/Google Drive/figs/states.png", height=1000, width=1600)
plot(states, col="grey95",lwd=2)
dev.off()


png("C:/Users/Localadmin/Google Drive/figs/floridaborder.png", height=1000, width=1600, bg="transparent")
plot(states[as.character(states$NAME)=="Florida",], col="transparent",lwd=10, bg="transparent") 
dev.off()




hist(ordsums$l1o2*100, col="light blue")



##################
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(RColorBrewer)
library(laticeExtra)
library(sp)

compare <- raster("C:/Users/Localadmin/Documents/Maxent_modeling/summaries/asciis/hotspot_compare.asc")
proj4string(compare)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
compare <- projectRaster(compare, crs="+init=epsg:5070")

eco <- readOGR(dsn="C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco", stringsAsFactors = F)
proj4string(eco) <- "+init=epsg:5070"

extentShape = readOGR(dsn = "C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
extentShape <- spTransform(extentShape,"+init=epsg:5070")


#####################################

#D <- (as.data.frame(compare))
#head(D)
#
#D$hotspot_compare2[D$hotspot_compare>.7 & D$hotspot_compare <1.3 & !is.na(D$hotspot_compare) ] <- "Abundance hotspot"
#D$hotspot_compare2[D$hotspot_compare>1.7 & D$hotspot_compare <2.3 & !is.na(D$hotspot_compare) ] <- "Both hotspots"
#D$hotspot_compare2[D$hotspot_compare>3.7 & D$hotspot_compare <4.3 & !is.na(D$hotspot_compare) ] <- "Richness hotspot"
#D$hotspot_compare2[D$hotspot_compare>2.7 & D$hotspot_compare <3.3 & !is.na(D$hotspot_compare) ] <- "Not hotspot"
#
#summary(as.factor(D$hotspot_compare2))
#
#D$hotspot_compare2 <- ordered(D$hotspot_compare2, levels = c("Abundance hotspot", "Both hotspots", "Richness hotspot", "Not hotspot"))
#str(D$hotspot_compare2)
#
#compare <- raster(nrows=nrow(compare),ncols=ncol(compare),ext=extent(compare),vals=D$hotspot_compare2, crs=crs(compare))
#compare$hotspot_compare2<- as.factor(compare$hotspot_compare)
#####################################
#install.packages('plot3D')
library(plot3D)
library(GISTools)
library(maps)


plot(compare, col= c("orangered2","darkslateblue", "grey98","deepskyblue1"), axes=F, box=F, legend=F, #mar=c(5, 0, 4, 10))
     xlim=c(-2.3e+06, 5e+06))

legend(x=2e+06, y=2000000, legend=c("Richness hotspot", "Abundance hotspot", "Both hotspots", "Not hotspot"),
       col=c("deepskyblue","darkslateblue","orangered2", "grey98"),bty="n", xpd=T,
       pch=c(15,15, 15, 15), pt.cex=3)
legend(x=2e+06, y=2000000, legend=c("Richness hotspot", "Abundance hotspot", "Both hotspots", "Not hotspot"),
       col=c("black","black","black", "black"),bty="n", xpd=T,
       pch=c(0,0,0,0), pt.cex=3)


pdf("C:/Users/Localadmin/Google Drive/figs/compare.png",width=168/25.4,height=100/25.4)
spplot(compare$hotspot_compare, col=c("orangered2","darkslateblue","grey98", "deepskyblue1") )

map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.7,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")
plot(eco,add=T,col="transparent",lwd=2)
dev.off()


plot(compare)
bbox(compare)
abline(h=185177)
abline(h=3281388)
abline(v=-2921970, col="blue")
abline(v=-2946360)
plot(extentShape, mar=c(-20,-20,-20,-20))

extent(extentShape)
bbox(extentShape)
plot(compare)
abline(v=-2358781)
abline(v=2258381.8)
abline(h=310447)
abline(h=3166167)


pdf("C:/Users/Localadmin/Google Drive/figs/test2.pdf", width=4.5, height=3)
plot(comp2, mar=c(0,0,0,0), cex=0.01)
#spplot(maxentRange, asp=1, mar=c(0,0,0,0))
dev.off()

comp2 <- crop(compare, extent(extentShape))
plot(comp2)


library(sp)
e <- as(raster::extent(-2200000, 2100000, 250000, 2900000), "SpatialPolygons")
proj4string(e) <- "+init=epsg:5070"
plot(e)
plot(comp2, add=T, xpd=T)
plot()

plot(comp2)


pdf("C:/Users/Localadmin/Google Drive/figs/test3.pdf", width=6, height=9)
par(mfrow=c(3,1))
plot(extentShape)
plot(extentShape)
plot(extentShape)
dev.off()


pdf("C:/Users/Localadmin/Google Drive/figs/test4.pdf", width=3, height=4.5)
par(mfrow=c(3,1), mai=c(0,0,0,0))
plot(extentShape)
plot(extentShape)
plot(extentShape)
dev.off()

pdf("C:/Users/Localadmin/Google Drive/figs/test4.pdf", width=3, height=4.5)
par(mfrow=c(3,1), mai=c(0,0,0,.5))
plot(extentShape)
plot(extentShape)
plot(extentShape)
dev.off()



#plot(bbox(extentShape))
############################
##### Three panel monster
############################
abunpts <- read.table("file:///C:/Users/Localadmin/Google Drive/Invasive-plant-abundance-SDM-files/edd_w_environmental.csv",
                      sep=",", header=T, stringsAsFactors = F)
head(abunpts)
abunpts <- abunpts[abunpts$species=="BRPA4",]
coordinates(abunpts) <- c(5,4)

occ.pts <- read.table("C:/Users/Localadmin/Documents/MaxEnt_modeling/species/BRPA4.csv", sep=",", stringsAsFactors = F, header=T)
str(occ.pts)
coordinates(occ.pts) <- c(2,3)


proj4string(occ.pts) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
proj4string(abunpts) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
occ.pts <- spTransform(occ.pts, "+init=epsg:5070")
abunpts <- spTransform(abunpts, "+init=epsg:5070")

occ.ran <- occ.pts[sample(1:nrow(occ.pts), 130),]
abun.ran <- abunpts[sample(1:nrow(abunpts), 50),]

maxentRange <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/Binary_asciis/BRPA4.asc")
names(maxentRange) <- "layer"
maxentRange <- rasterToPolygons(maxentRange, fun=function(layer){layer > 0}, na.rm=T, dissolve=T)
proj4string(maxentRange) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
maxentRange <- spTransform(maxentRange, "+init=epsg:5070")


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

library(dismo)
messpts <- data.frame(cbind(abunpts$bio_8, abunpts$nlcd_6))#), abunpts$nlcd_5, abunpts$nlcd_8))
head(messpts)
colnames(messpts)<-c("bio_8", "nlcd6" )#, "nlcd_5", "nlcd_8")

messplot <- mess(bio,messpts,full=F)

messplot$mess2 <- -99
messplot$mess2[messplot$mess < -10] <- 1
messplot$mess2[messplot$mess >= -10] <- 3
messplot$mess2[is.na(messplot$mess) | messplot$mess == Inf] <- NA
proj4string(messplot)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
messplot <- projectRaster(messplot, crs="+init=epsg:5070")

ordinal <- raster("C:/Users/Localadmin/Documents/MaxEnt_modeling/ORDINAL/ALL_ABUN/BRPA4.asc") 
ordinal$layer <- 0
ordinal$layer <- ordinal$BRPA4
ordinal$layer[ordinal$BRPA4==0] <- NA

ordinal$layer2 <- 0
ordinal$layer2[ordinal$BRPA4==3] <- 1
ordinal$layer2[ordinal$BRPA4 != 3 | is.na(ordinal$BRPA4)] <- NA
proj4string(ordinal)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ordinal <- projectRaster(ordinal, crs="+init=epsg:5070")

options(scipen=999)


pdf("C:/Users/Localadmin/Google Drive/figs/3panel.pdf",width=168/25.4,height=230/25.4)
par(mfrow=c(3,1), mai=c(0,0,0,1.5))

#par(mfrow=c(3,1))
plot(extentShape, lwd=1.5, col="grey94", asp=1)
plot(maxentRange, add=T,col="deepskyblue",border="transparent", bg='transparent')
plot(extentShape, lwd=1.5, col="transparent",add=T)
plot(occ.pts,col="black",add=T, cex=1, lwd=1)
plot(abunpts, pch=16, col="firebrick1", add=T, cex=1)
legend(x=2e+06, y=2000000, legend=c("Establishment range", "Unsuitable", "Occurrence record", "Abundance record"),
       col=c("deepskyblue","grey94","black", "firebrick1"),bty="o", xpd=T,
       pch=c(15,15, 3, 16), pt.cex=c(3.5,3.5,1.5,1.5), cex=1.5, title="MaxEnt predictions")
legend(x=2e+06, y=2000000, legend=c("Establishment range", "Unsuitable", "Occurrence record", "Abundance record"),
       col=c("black","black","black", "firebrick1"),bty="n", xpd=T,
       pch=c(0,0, 3, 16), pt.cex=c(3.5,3.5,1.5,1.5), cex=1.5, title="MaxEnt predictions")
text(x=-2350000, y=3120000,"(a)", cex=2)


#plot(extentShape, lwd=1.5, col="transparent",asp=1)
plot(extentShape, col="grey96",lwd=1.5)
plot(messplot$mess2, add=T, col=c("firebrick2","slategray2"),box=F,legend=F, axes=F)
plot(extentShape, col="transparent",lwd=1.5,add=T)
text(x=-2350000, y=3120000,"(b)", cex=2)
legend(x=2e+06, y=2000000, legend=c("Analog", "Non-analog", "Unsuitable               "),
       col=c("slategray2","firebrick3","grey96"),bty="o", xpd=T,
       pch=c(15,15, 15), pt.cex=c(3.5,3.5,3.5), cex=1.5, title="MESS results")
legend(x=2e+06, y=2000000, legend=c("Analog", "Non-analog", "Unsuitable               "),
       col=c("black","black","black" ),bty="n", xpd=T,
       pch=c(0,0,0), pt.cex=c(3.5,3.5,3.5), cex=1.5, title="MESS results")


plot(extentShape, lwd=1.5, col="grey99")
plot(ordinal$layer, add=T, col= c("lightgoldenrod","orange","firebrick3"), legend=F)
plot(extentShape, col="transparent", lwd=1.5, add=T) 
legend(x=2e+06, y=2000000, legend=c("High (impact range) ", "Medium", "Low", "Unsuitable"),
       col=c("firebrick2","orange","lightgoldenrod","grey99"),bty="o", xpd=T,
       pch=c(15,15, 15), pt.cex=c(3.5,3.5,3.5), cex=1.5, title="Abundance predictions")
legend(x=2e+06, y=2000000, legend=c("High (impact range) ", "Medium", "Low", "Unsuitable"),
       col=c("black","black","black", "black"),bty="o", xpd=T,
       pch=c(0,0,0,0), pt.cex=c(3.5,3.5,3.5), cex=1.5, title="Abundance predictions")
text(x=-2350000, y=3120000,"(c)", cex=2)
map.scale(x=-2.5e+06,y=800000,ratio=F,relwidth=0.15,cex=0.75,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")

dev.off()


########
ordsums <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_8_2018b.csv", sep=",", header=T, stringsAsFactors = F)
head(ordsums)
splist <- data.frame(ordsums$species.code[ordsums$USE==1])
colnames(splist) <- c("name")

library(rgdal)
library(raster)

eco <- readOGR(dsn="C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco")
#plot(eco)
#eco@data
eco$ID2 <- as.character(eco$IDs)


eco <- eco[eco$ID2=="MARINE WEST COAST FOREST",]
eco <- spTransform(eco, "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
plot(eco)
proj4string(eco)

splist$mwcf <- "filler"

for (i in 1:length(splist$name)){
  sp <- raster(paste0("C:/Users/mwone/Desktop/HI_ABUN/",splist$name[i],".asc"))
  proj4string(sp) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  
  sp <- mask(crop(sp,eco),eco)
  
  sp <- as.data.frame(sp)
  colnames(sp) <- c("sp")
  splist$mwcf[i] <- NROW(sp[sp$sp==1 & !is.na(sp$sp),])/NROW(sp[!is.na(sp$sp),])
  print(i)
}

#sp <- raster(paste0("C:/Users/mwone/Desktop/HI_ABUN/","PUMOL",".asc"))
hist(as.numeric(splist$mwcf))

splist$name[as.numeric(splist$mwcf)>0.4]


plot(sp)

splist

plot(raster("C:/Users/mwone/Desktop/HI_ABUN/AIAL.asc"))
plot(raster("C:/Users/mwone/Desktop/HI_ABUN/ARMI2.asc"))
plot(raster("C:/Users/mwone/Desktop/HI_ABUN/IRPS.asc"))
plot(raster("C:/Users/mwone/Desktop/HI_ABUN/LIVU2.asc")) #*****
plot(raster("C:/Users/mwone/Desktop/HI_ABUN/TAVU.asc")) ## Both NE and NW
plot(raster("C:/Users/mwone/Desktop/HI_ABUN/CEDI3.asc")) 


plot(raster("C:/Users/mwone/Desktop/HI_ABUN/EUES.asc")) 


plot(raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/nlcd/nlcd3.asc"))
plot(raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/nlcd/nlcd4.asc")) # L+ ****
plot(raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/nlcd/nlcd5.asc")) # Q-
plot(raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/nlcd/nlcd6.asc")) # L+ great plains
plot(raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/nlcd/nlcd7.asc")) #
plot(raster("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/nlcd/nlcd8.asc"))




################### HOTSPOR COMPARISON FIGURE
install.packages("GISTools")
install.packages("maps")


library(GISTools)
library(maps)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)

#single column = 79mm, 2/3rd column = 110mm, double column = 168mm, maximum height of figure = 230mm.

extentShape = readOGR(dsn = "C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ArcFiles_2_2_2018/us_shape_2_9_2018", layer = "us_shape")
compare <- raster("C:/Users/mwone/Desktop/asciis/hotspot_compare.asc")
eco <- readOGR(dsn="C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ecoregions", layer="eco", stringsAsFactors = F)
#proj4string(eco) <- "+init=epsg:5070"

proj4string(compare)<- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
compare <- projectRaster(compare, crs="+init=epsg:5070")

options(scipen = 999)

pdf("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/figs/compare.pdf",width=168/25.4,height=3)
par(mai=c(0,0,0,.2)) # c(bottom, left, top, right) 
plot(compare, col=c("orangered2","darkslateblue","grey98", "deepskyblue1"), box=F, axes=F, legend=F)
#plot(compare, col=c("transparent",adjustcolor(c("black"),.15), "transparent", "transparent"), add=T, box=F, axes=F, legend=F)
#plot(extentShape, col="transparent",lwd=2,add=T)
map.scale(x=-2.6e+06,y=800000,ratio=F,relwidth=0.2,cex=0.4,bg="transparent")
north.arrow(xb=-2.4e+06,yb=1200000,len=80000,lab="N",cex=1,bg="transparent")

legend(x=2.0e+06, y=1950000, legend=c("Abundance hotspot", "Establishment hotspot", "Both hotspots", "Non-hotspot"),
       col=c("orangered2","deepskyblue1","darkslateblue", "grey98"),bty="n", xpd=T,
       pch=c(15,15,15,15), pt.cex=c(2.3,2.3,2.3,2.3), cex=.85)
legend(x=2.0e+06, y=1950000, legend=c("Abundance hotspot", "Establishment hotspot", "Both hotspots", "Non-hotspot"),
       col=c("black","black","black", "black"),bty="o", xpd=T,
       pch=c(0,0,0,0), pt.cex=c(2.3,2.3,2.3,2.3), cex=.85)

plot(eco,add=T,col="transparent",lwd=1)
dev.off()


#######################################################################
##### double checking details of this project for clarification for
##### JMA and BB

ordsums <- read.table("C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/ordsums_5_31_2018.csv", header = T, sep = ",", quote= "\"", 
                      comment.char= "", stringsAsFactors = F, strip.white = T)
head(ordsums)

nrow(ordsums[ordsums$rule5broke==1 & ordsums$USE==1,])
nrow(ordsums[ordsums$rule5broke==0 & ordsums$USE==1,])

# Perennial tree

ordsums[ordsums$species.code=="BRPA4",]

median(ordsums$kappa[ordsums$USE==1])
ordsums$kappa[ordsums$species.code=="BRPA4"]

hist(ordsums$impac_estab[ordsums$USE==1])
ordsums$impac_estab[ordsums$species.code=="BRPA4"] #large impact range

median(ordsums$no.pts.abun[ordsums$USE==1])
ordsums$no.pts.abun[ordsums$species.code=="BRPA4"] #few points

median(ordsums$thinned_occ[ordsums$USE==1])
ordsums$thinned_occ[ordsums$species.code=="BRPA4"]

median(ordsums$no.terms[ordsums$USE==1])
ordsums$no.terms[ordsums$species.code=="BRPA4"]

median(ordsums$areaSQKM_5_7[ordsums$USE==1])
ordsums$areaSQKM_5_7[ordsums$species.code=="BRPA4"]

mean(ordsums$range_size[ordsums$USE==1])
ordsums$range_size[ordsums$species.code=="BRPA4"]

median(ordsums$MESSgn10[ordsums$USE==1])
ordsums$MESSgn10[ordsums$species.code=="BRPA4"]
hist(ordsums$MESSgn10[ordsums$USE==1])

ordsums$no.terms[5]
ordsums$no.vars[5]
ordsums$formu[5]

length(ordsums$species.code[ordsums$USE==1 & (ordsums$no.vars > ordsums$no.pts.abun)])
length(ordsums$species.code[ordsums$USE==1 & ((ordsums$no.terms-2) > ordsums$no.vars)])
58/70



#### do all sources of duration/form need to be cited?
traits <- read.table("file:///C:/Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/all_contig.csv" , header = T, sep = ",", quote= "\"", 
                     comment.char= "", stringsAsFactors = F, strip.white = T)

summary(traits$source)

ordsums$source <- "BLANK"

for (i in 65:155){
ordsums$source[i] <- traits$source[traits$NewCode==ordsums$species.coSde[i]]
print(i)
}

ordsums$species.code[64]
head(traits)
unique(traits$source)
unique(ordsums$source)

## maximum abun richness and full richness
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(RColorBrewer)

hi_abun_richness <- raster("C:/Users/mwone/Downloads/richness_map_HI_ABUN.asc")
full_richness  <- raster("C:/Users/mwone/Downloads/richness_map_FULL.asc")
head(hi_abun_richness)
hi_abun <- as.data.frame(hi_abun_richness$richness_map_HI_ABUN)
summary(hi_abun) #max 16
full <- as.data.frame(full_richness$richness_map_FULL)
summary(full) #max 52

  mean(ordsums70$thinned_occ)
  mean(ordsums70$no.pts)
  mean(ordsums70$AUC)
1-mean(ordsums70$MESSgn10)
  mean(ordsums70$kappa)
  mean(ordsums70$kappaP)
  mean(ordsums70$areaSQKM_5_7)/(10^6)
  mean(ordsums70$range_prop_area)
  mean(ordsums70$hiAbunSQKM_5_7)/(10^6)
  mean(ordsums70$hiAbun_prop_range)
  mean(ordsums70$infilling)
  mean(ordsums70$l1o2)

  median(ordsums70$thinned_occ)
  median(ordsums70$no.pts)
  median(ordsums70$AUC)
1-median(ordsums70$MESSgn10)
  median(ordsums70$kappa)
  median(ordsums70$kappaP)
  median(ordsums70$areaSQKM_5_7)/(10^6)
  median(ordsums70$range_prop_area)
  median(ordsums70$hiAbunSQKM_5_7)/(10^6)
  median(ordsums70$hiAbun_prop_range)
  median(ordsums70$infilling)
  median(ordsums70$l1o2)


1-ordsums$MESSgn10[ordsums$species.code=="PARE3"]
1-ordsums70$MESSgn10[ordsums70$species.code=="CHMA2"]

ordsums <- ordsums[order(-ordsums$USE,ordsums$species.code),]
write.table(data.frame(ordsums$species.code,ordsums$MESS_abun2range_SEL),"C://Users/mwone/Google Drive/Invasive-plant-abundance-SDM-files/MESS_for_pretty_table.csv", sep=",", row.names=F)
