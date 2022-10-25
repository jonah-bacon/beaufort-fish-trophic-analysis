#CCSIA Analyses - Jessica R. Glass
#February 23 2019
#Last updated 6/12/2019

#Clear R session
rm(list = ls())

#### INDEX ####
#51 Read in all bulk/CCSIA data
#119 Look at ranges of bulk N and C by site
#187 GT analyses using CCSIA TP (Bradley 2015)
#269 Plot Fork Length versus bulk C and N
#286 SIBER analyses using Bulk Data
#378 SIBER Analysis but with CCSIA_TP (Bradley) instead of bulk N
#467 SIBER with CCSIA_TP (Weighted Means) instead of N
#563 Read in AA-CCSIA full data 
#602 Plotting Source and Trophic AA by locality
#721 AA boxplots by locality
#730 Threonine
#769 Modeling source and trophic AA by locality
#877 Bulk N versus N(phe) Plot and Model
#938 Compare Bluefins and GTs from Providence (Bulk)
#992 SIBER analyses with GTs and BTs from Providence
#1098 SIBER Analyses with GTS and BTs using Weighted Means TP
#1202 Calculate weighted means for trophic and source AAs
#1249 Weighted Means Trophic Position Analyses
#1266 Comparison of all TP estimates, FL and C
#1326 Test differences in TP values using ANOVA
#1387 Test weighted means and C for each locality
#1437 Use TP weighted means for site-comparisons (mean, SD, range)
#1509 GT-Bluefin comparisons using TP weighted means
#1548 Baseline analyses
#1563 basic PCA (separate source and trophic)
#1675 PCA using normalized weighted means
#1764 Plot Sampling Locality Map
#1797 Propogation of Errors

library(devtools)
library(ggplot2)
#devtools::install_github("andrewljackson/SIBER@v2.1.4", build_vignettes = TRUE)
library(SIBER)
#install.packages("rjags")
#install.packages("coda")
library(coda)
library(rjags)


#Set working directory
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA")

#Read in .csv file with all bulk/CCSIA data combined
FullDat <- read.csv("SIA_Samples_Comb_Analysis.csv", head = T)

head(FullDat)
str(FullDat)
colnames(FullDat)

FullDat

#Plot SIA samples by length and locality
plot(FullDat$FL~ FullDat$Locality)

#For both GTs and Bluefins
#plot CCSIA Trophic position (using metrics from Bradley 2015) versus FL
#First remove sample JRG15_259A because of no FL value
RedDat <- FullDat[-c(16),]

plot(RedDat$FL, RedDat$CCSIA_TP,col=RedDat$Group,type="p",pch=20, cex = 1.5, ylab="CCSIATrophic Position")
legend("topleft",legend=as.character(paste(unique(RedDat$Locality))),
       pch=19,col=1:length(unique(RedDat$Group)))

#Test for correlation between FL and TP
FLmod.all <- lm(RedDat$CCSIA_TP~RedDat$FL)
summary(FLmod.all)
#There is a slightly significant correlation between FL and TP p=0.0247

#Add the best fit line to the plot
abline(FLmod.all, lwd = 2, col = "gray")

#Look at correlation between FL and bulk N for Gts
plot(GTDat$FL, GTDat$N,col=GTDat$Group,type="p",pch=20, cex = 1.5, ylab="N")
FLmod.N <- lm(GTDat$N~GTDat$FL)
summary(FLmod.N)
#significantly negative correlation(?)
abline(FLmod.N, lwd = 2, col = "gray")


### For ALL GTs and BTs, separate by site
Moz.Dat <- RedDat[RedDat$Locality=="MOZ",]
FLmod.Moz <- lm(Moz.Dat$CCSIA_TP~Moz.Dat$FL)
summary(FLmod.Moz)
abline(FLmod.Moz, lwd = 2, col = "black", lty = 2)
#p = 0.626

Mahe.Dat <- RedDat[RedDat$Locality=="MAHE",]
FLmod.Mahe <- lm(Mahe.Dat$CCSIA_TP~Mahe.Dat$FL)
summary(FLmod.Mahe)
abline(FLmod.Mahe, lwd = 2, col = "red", lty = 2)
#p = 0.0709

StJo.Dat <- RedDat[RedDat$Locality=="STJO",]
FLmod.StJo <- lm(StJo.Dat$CCSIA_TP~StJo.Dat$FL)
summary(FLmod.StJo)
abline(FLmod.StJo, lwd = 2, col = "green", lty = 2)
#p = 0.287

Prov.Dat <- RedDat[RedDat$Locality=="PROV",]
FLmod.Prov <- lm(Prov.Dat$CCSIA_TP~Prov.Dat$FL)
summary(FLmod.Prov)
abline(FLmod.Prov, lwd = 2, col = "blue", lty = 2)
#p = 0.155

##Subsample GTs
GTDat <- FullDat[FullDat$Species == "GT",]
head(GTDat)
##Subsample Bluefins
BTDat <- FullDat[FullDat$Species == "BT",]

#For JUST GTs, separate by site
GT.Moz.Dat <- GT.Red[GT.Red$Locality=="MOZ",]
GT.Mahe.Dat <- GT.Red[GT.Red$Locality=="MAHE",]
GT.StJo.Dat <- GT.Red[GT.Red$Locality=="STJO",]
GT.Prov.Dat <- GT.Red[GT.Red$Locality=="PROV",]

##look at bulk Nitrogen ranges for GTs
range(GT.Moz.Dat$N) #12.13379 15.35844
sd(GT.Moz.Dat$N) # SD: 0.881
#Range: 3.225
range(GT.Mahe.Dat$N) #Does not work bc of NAs
sd(GT.Mahe.Dat$N, na.rm = T) # sd: 0.637
#13.897 15.859
#Range: 1.962
range(GT.StJo.Dat$N)
sd(GT.StJo.Dat$N) #SD: 0.404
#12.74802 14.32228
#Range: 1.57426
range(GT.Prov.Dat$N)
sd(GT.Prov.Dat$N) #SD: 0.347
#12.33036 13.38122
#Range: 1.05086

#look at TP ranges (CCSIA_SingleAA method- NOT weighted means)
range(GT.Moz.Dat$CCSIA_TP)#does not work bc of NAs
sd(GT.Moz.Dat$CCSIA_TP, na.rm =T)

range(GT.Mahe.Dat$CCSIA_TP)
sd(GT.Mahe.Dat$CCSIA_TP)

range(GT.StJo.Dat$CCSIA_TP)
sd(GT.StJo.Dat$CCSIA_TP)

range(GT.Prov.Dat$CCSIA_TP)
sd(GT.Prov.Dat$CCSIA_TP)


#look at C ranges
range(GT.Moz.Dat$AdjC) 
sd(GT.Moz.Dat$AdjC) #sd: 0.6827247
#-17.98585 -15.50188
#Range: 2.48397
range(GT.Mahe.Dat$AdjC) #does not work bc of NAs
sd(GT.Mahe.Dat$AdjC, na.rm =T)
#-17.49471 -12.32075
#Range: 5.17395
range(GT.StJo.Dat$AdjC)
sd(GT.StJo.Dat$AdjC) #sd: 1.444876
#-15.31208 -10.55126
#Range: 4.76082
range(GT.Prov.Dat$AdjC)
sd(GT.Prov.Dat$AdjC) #sd: 1.398049
#-13.84661 -10.10950
#Range: 3.73711

#look at CCSIA (Bradley TP means
mean(Moz.Dat$CCSIA_TP, na.rm = T) 
#4.246254
mean(Mahe.Dat$CCSIA_TP, na.rm = T) 
#3.830443
mean(StJo.Dat$CCSIA_TP, na.rm = T)
#4.016871
mean(Prov.Dat$CCSIA_TP, na.rm = T)
#3.724368

#Boxplots of TP
boxplot(Mahe.Dat$CCSIA_TP, Moz.Dat$CCSIA_TP, StJo.Dat$CCSIA_TP, Prov.Dat$CCSIA_TP, main = "TP Ranges: 1.Mahe, 2.Moz, 3.StJo, 4.Prov")


####GT Analyses####
#Plot Adjusted C vs Nitrogen
plot(GTDat$AdjC, GTDat$N, col=GTDat$Group,type="p",pch=20, cex = 1.5)
legend("topright",legend=as.character(paste(unique(GTDat$Locality))),
       pch=19,col=1:length(unique(GTDat$Group)))
#add labels to indicate YEAR
with(GTDat, text(GTDat$N~GTDat$AdjC, labels = GTDat$Year, cex = 0.5, pos = 4))

#Remove sample JRG15_259A because of no FL value
GT.Red <- GTDat[-c(16),]
#And remove the Acoustic column because of many NAs
GT.Red <- GTDat[,-c(5)]

#plot CCSIA Trophic position (using metrics from Bradley 2015) versus FL
plot(GT.Red$FL, GT.Red$CCSIA_TP,col=GT.Red$Group,type="p",pch=20, cex = 1.5, ylab="CCSIATrophic Position", main = "GT FL vs CCSIA TP")
legend("topleft",legend=as.character(paste(unique(GT.Red$Locality))),
       pch=19,col=1:length(unique(GT.Red$Group)))

#Test for correlation between FL and TP
FLmod <- lm(GT.Red$CCSIA_TP~GT.Red$FL)
summary(FLmod)
#There is a slightly significant correlation between FL and TP p=0.0707

#Add the best fit line to the plot
abline(FLmod, lwd = 2, col = "gray")

#add labels for dispersal distance
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$FL, labels = GT.Red$DispDist_km, cex = 0.5, pos = 2))
#add labels for acoustic tags
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$FL, labels = GT.Red$Acoustic, cex = 0.5, pos = 4))

#plot CCSIA Trophic position (using metrics from Bradley 2015) versus C
plot(GT.Red$AdjC, GT.Red$CCSIA_TP,col=GT.Red$Group,type="p",pch=20, cex = 1.5, ylab="CCSIATrophic Position", main = "GT CCSIA_TP vs C")
legend("topright",legend=as.character(paste(unique(GT.Red$Locality))),
       pch=19,col=1:length(unique(GT.Red$Group)))
#add labels for dispersal distance
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$AdjC, labels = GT.Red$DispDist_km, cex = 0.7, pos = 3))
#add labels for dispersal direction
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$AdjC, labels = GT.Red$DispDirect, cex = 0.7, pos = 4))

#add labels for sample name
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$AdjC, labels = GT.Red$Sample_ID, cex = 0.4, pos = 3))

#add labels for year
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$AdjC, labels = GT.Red$Year, cex = 0.5, pos = 2))

TP_C.mod <- lm(GT.Red$CCSIA_TP ~ GT.Red$AdjC)
summary(TP_C.mod)
#significant effect of Carbon on TP
abline(TP_C.mod, lty = 2, col = "gray")

#Only plot Mozambique samples (multi-year)
plot(Moz.Dat$AdjC, Moz.Dat$CCSIA_TP, col = as.factor(Moz.Dat$Year), type="p",pch=20, cex = 1.5, ylab="CCSIATrophic Position", main = "Moz GT CCSIA_TP vs C")
with(Moz.Dat, text(Moz.Dat$CCSIA_TP~Moz.Dat$AdjC, labels = Moz.Dat$Year, cex = 0.5, pos = 3))

#Only plot Mahe samples (multi-year)
plot(Mahe.Dat$AdjC, Mahe.Dat$CCSIA_TP, col = as.factor(Mahe.Dat$Year), type="p",pch=20, cex = 1.5, ylab="CCSIATrophic Position", main = "Mahe GT CCSIA_TP vs C")
with(Mahe.Dat, text(Mahe.Dat$CCSIA_TP~Mahe.Dat$AdjC, labels = Mahe.Dat$Year, cex = 0.5, pos = 3))

plot(Mahe.Dat$AdjC, Mahe.Dat$FL, type="p",pch=20, cex = 1.5, ylab="FL", main = "Mahe GT FL vs C")
Mahe.C.mod <- lm(Mahe.Dat$FL ~ Mahe.Dat$AdjC)
summary(Mahe.C.mod)
abline(Mahe.C.mod, lty = 2, col = "gray")
#there is a significant negative correlation between FL and Adj C
#p = 0.00588

#Plot Bulk TP versus FL
plot(GT.Red$FL,GT.Red$BulkTP,col=GT.Red$Group,type="p",pch=20, cex=1.5, ylab="Bulk Trophic Position", main = "GT FL vs Bulk TP")
legend("bottomright",legend=as.character(paste(unique(GT.Red$Locality))),
       pch=19,col=1:length(unique(GT.Red$Group)))


#Plot Bulk TP versus CCSIA TP
plot(GT.Red$BulkTP,GT.Red$CCSIA_TP,col=GT.Red$Group,type="p",pch=20, ylab="CCSIA Trophic Position")
legend("topleft",legend=as.character(paste(unique(GT.Red$Locality))),
       pch=19,col=1:length(unique(GT.Red$Group)))
#add line of y-intercept=0 and slope=1
abline(0,1)
#label by year
with(GT.Red, text(GT.Red$CCSIA_TP~GT.Red$BulkTP, labels = GT.Red$Year, cex = 0.5, pos = 4))
#it looks like TP values using bulk data were underestimated for MAU and PROV, and overestimated for STJo and most of Mahe

#Plot FL vs Adj C
plot(GT.Red$FL~GT.Red$AdjC,col=GT.Red$Group,type="p", cex = 1.5,pch=20, xlab = "Carbon", ylab = "Fork Length")
legend("bottomleft",legend=as.character(paste(unique(GT.Red$Locality))),
       pch=19,col=1:length(unique(GT.Red$Group)))
with(GT.Red, text(GT.Red$FL~GT.Red$AdjC, labels = GT.Red$Sample_ID, cex = 0.5, pos = 4))
Carb.mod <- lm(GT.Red$FL ~ GT.Red$AdjC)
summary(Carb.mod)
#No effect of FL on Carbon
abline(Carb.mod, lty = 2, col = "Gray")

#Plot FL vs Bulk N
plot(GT.Red$N~GT.Red$FL,col=GT.Red$Group,type="p", cex = 1.5,pch=20, ylab = "Bulk Nitrogen", xlab = "Fork Length")
legend("topright",legend=as.character(paste(unique(GT.Red$Locality))),
       pch=19,col=1:length(unique(GT.Red$Group)))
with(GT.Red, text(GT.Red$FL~GT.Red$AdjC, labels = GT.Red$Sample_ID, cex = 0.5, pos = 4))


### Calculate Bayesian Ellipses in Package 'SIBER'
#Link to Tutorial: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html

## For SIBER analyses, remove the single Mtentu sample and the Mauritius samples (last 3 lines)
#GT.Red <- GT.Red[GT.Red$Locality=="MTENTU",]
GT.Red <- GT.Red[-c(57),]

#GT.Red <- GT.Red[GT.Red$Locality=="MAU",]
GT.Red <- GT.Red[-c(55,56),]

SIBER.dat <- read.csv("SIBER_GT_Data.csv")

#Set the seed so results will be the same for subsequent runs
set.seed(1)

GT.SIBER <- createSiberObject(SIBER.dat)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions. - Can change p.interval to represent X % of the data (e.g. p = 0.95 represents 95% of the data) or setting P to null produces  ML standard ellipses (40%)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


par(mfrow=c(1,1))
plotSiberObject(GT.SIBER,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                y.limits = as.vector(c(11, 16)),
                x.limits = as.vector(c(-19,-10)),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

legend("bottomleft",legend=as.character(paste(GT.SIBER$all.groups)),
       pch=19,col=1:length(unique(GT.SIBER$group.names)))
#Levels: 1 (red) = Moz, 2 (black) = Mahe, 3 (blue) = St. Jo, 4 (green) = Prov

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(GT.SIBER, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(GT.SIBER, n = 100, p.interval = 0.80, ci.mean = T,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(GT.SIBER)
print(group.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(GT.SIBER, parms, priors)

#Comparing among groups using the Standard Ellipse Area

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


#### SIBER Analysis but with CCSIA_TP (Bradley) instead of bulk N ####
SIBER.GT.TP <- read.csv("SIBER_GT_TP.csv") 

#Set the seed so results will be the same for subsequent runs
set.seed(1)

GT.SIBER.TP <- createSiberObject(SIBER.GT.TP)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions. - Can change p.interval to represent X % of the data (e.g. p = 0.95 represents 95% of the data) or setting P to null produces  ML standard ellipses (40%)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

dev.off()

par(mfrow=c(1,1))
plotSiberObject(GT.SIBER.TP,
                ax.pad = NULL, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(2,1),
                y.limits = as.vector(c(3, 5)),
                x.limits = as.vector(c(-19,-10)),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = ("CCSIA TP")
)

legend("bottomleft",legend=as.character(paste(unique((GT.SIBER.TP$all.groups)))),
       pch=19,col=1:length(unique(GT.SIBER$group.names)))
#Levels: 1 (red) = Moz, 2 (black) = Mahe, 3 (blue) = St. Jo, 4 (green) = Prov


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(GT.SIBER.TP, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(GT.SIBER.TP, n = 100, p.interval = 0.80, ci.mean = T,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(GT.SIBER.TP)
print(group.ML)
#       1.MOZ     2.MAHE   3.STJO    4.PROV
#TA   1.1539382 2.203265 2.382412 1.3298684
#SEA  0.3757328 1.010640 1.286533 0.7158378
#SEAc 0.3966069 1.094860 1.447349 0.8181003


# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(GT.SIBER.TP, parms, priors)

#Comparing among groups using the Standard Ellipse Area

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


#### SIBER with CCSIA_TP (Weighted Means) instead of N ####
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA")
SIBER.GT.TPWT <- read.csv("SIBER_GT_TPWT.csv") 

#Set the seed so results will be the same for subsequent runs
set.seed(1)

GT.SIBER.TPWT <- createSiberObject(SIBER.GT.TPWT)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions. - Can change p.interval to represent X % of the data (e.g. p = 0.95 represents 95% of the data) or setting P to null produces  ML standard ellipses (40%)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

dev.off()

par(mfrow=c(1,1))
plotSiberObject(GT.SIBER.TPWT,
                ax.pad = NULL, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(2,1),
                y.limits = as.vector(c(3, 5.5)),
                x.limits = as.vector(c(-19,-10)),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = ("CCSIA TP")
                
)

legend("bottomleft",legend=as.character(paste(unique((GT.SIBER.TPWT$all.groups)))),
       pch=19,col=1:length(unique(GT.SIBER.TPWT$group.names)))
#Levels: 1 (red) = Moz, 2 (black) = Mahe, 3 (blue) = St. Jo, 4 (green) = Prov


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(GT.SIBER.TP, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(GT.SIBER.TP, n = 100, p.interval = 0.80, ci.mean = T,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(GT.SIBER.TPWT)
print(group.ML)
        #1.MOZ   2.MAHE   3.STJO   4.PROV
#TA   1.7010546 3.029373 2.324151 1.817983
#SEA  0.5899620 1.391486 1.018150 1.146651
#SEAc 0.6227377 1.507443 1.145419 1.310458


# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(GT.SIBER.TPWT, parms, priors)

#Comparing among groups using the Standard Ellipse Area

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


#Try an ANOVA <- **Need to look into this: how to test statistical differences between niche width
#SEA.B.df.full <- as.data.frame(SEA.B)
#anova(lm(SEA.B.df.full))
#summary(lm(SEA.B.df.full))

#### Read in CCSIA full Data ####
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA")
CCAA.dat <- read.csv("CCSIA_Full_workingdata.csv")

head(CCAA.dat)
str(CCAA.dat)

#Get min and max and summary stats for dataset
summary(CCAA.dat)

#now have 67 samples (BT and GTs)
#Pull out BTs
CCAA.BT <- CCAA.dat[c(14,15,17, 39:44),]
CCAA.BT

#Pull out GTs
CCAA.GT <- CCAA.dat[-c(14,15,17, 39:44),]
colnames(CCAA.GT)

#Remove Mtentu and Mauritius samples and sample from Mahe with NA for bulk N
CCAA.GT.Red <- CCAA.GT[-c(16,17,42,45),]
CCAA.GT.Red$Locality <- droplevels(CCAA.GT.Red$Locality, exclude = c("MAU", "MTENTU"))
CCAA.GT.Red$Locality

dev.off()
summary(CCAA.GT)

#Plot Chikarashi TP with old and new TEF and ß values
plot(CCAA.GT$TrophicPosition ~ CCAA.GT$NewTP_Bradley2015, col = CCAA.GT$Locality, type = "p", pch = 20)
#Look at differences between these values
CCAA.GT$NewTP_Bradley2015 - CCAA.GT$TrophicPosition

#look at boxplots of AAs
boxplot(CCAA.GT[,c(4:17)], main = "GT Amino Acids", names = c("Ala#", "Asx#", "Glx#", "Gly+", "Ile#", "Leu#", "Lys+", "Met+", "Phe+", "Pro#", "Ser+", "Thr*", "Tyr+", "Val#"), xlab = "# trophic  + source *metabolic")
mycol <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue50")
boxplot(CCAA.BT[,c(4:17)], add = T, col = mycol)


#Plot TP versus the difference of Thr and Phe
plot((CCAA.GT$Thr-CCAA.GT$Phe)[1:45], CCAA.GT$NewTP_Bradley2015[1:45], col = CCAA.GT$Locality, type = "p", pch = 20)
legend("topright",legend=unique(CCAA.GT$Locality),
       pch=19,col=unique(CCAA.GT$Locality))


#Source AAs: Gly, Phe, Lys
#Trophic AAs: Ala, Leu, Glu

#Source: Phe
dev.off()
par(mfrow = c(1,1))
plot(CCAA.GT$FL, CCAA.GT$Phe, col = CCAA.GT$Locality, type = "p", pch = 20)
Phe.mod <- lm(CCAA.GT$Phe~CCAA.GT$FL)
summary(Phe.mod)
abline(Phe.mod, lty = 2, col = "gray")
#non-significant negative relationship (p = 0.143)

#Source: Gly
plot(CCAA.GT$FL, CCAA.GT$Gly, col = CCAA.GT$Locality, type = "p", pch = 20)
Gly.mod <- lm(CCAA.GT$Gly~CCAA.GT$FL)
summary(Gly.mod)
abline(Gly.mod, lty = 2, col = "gray")
#significantly negative (p = 0.032)

#Source: Lys
plot(CCAA.GT$FL, CCAA.GT$Lys, col = CCAA.GT$Locality, type = "p", pch = 20)
Lys.mod <- lm(CCAA.GT$Lys~CCAA.GT$FL)
summary(Lys.mod)
abline(Lys.mod, lty = 2, col = "gray")
#negative but non-significant (p = 0.261)


#Trophic: Glu
plot(CCAA.GT$FL, CCAA.GT$Glx, col = CCAA.GT$Locality, type = "p", pch = 20)
Glx.mod <- lm(CCAA.GT$Glx~CCAA.GT$FL)
summary(Glx.mod)
abline(Glx.mod, lty = 2, col = "gray")
#non-significant (p = 0.554)

#Trophic: Ala
plot(CCAA.GT$FL, CCAA.GT$Ala, col = CCAA.GT$Locality, type = "p", pch = 20)
Ala.mod <- lm(CCAA.GT$Ala~CCAA.GT$FL)
summary(Ala.mod)
abline(Ala.mod, lty = 2, col = "gray")
#non-significant (p = 0.377)

#Trophic: Leu
plot(CCAA.GT$FL, CCAA.GT$Leu, col = CCAA.GT$Locality, type = "p", pch = 20)
Leu.mod <- lm(CCAA.GT$Leu~CCAA.GT$FL)
summary(Leu.mod)
abline(Leu.mod, lty = 2, col = "gray")
#non-significant (p = 0.61)

#plot all together
dev.off()
par(mfrow = c(2,3))

plot(CCAA.GT$FL, CCAA.GT$Phe, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Source: Phe", xlab = "FL")
abline(Phe.mod, lty = 2, col = "gray")
with(CCAA.GT, text(CCAA.GT$Phe~CCAA.GT$FL, labels = CCAA.GT$Locality, cex = 0.5, pos = 4))

plot(CCAA.GT$FL, CCAA.GT$Gly, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Source: Gly", xlab = "FL")
abline(Gly.mod, lty = 2, col = "gray")

plot(CCAA.GT$FL, CCAA.GT$Lys, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Source: Lys", xlab = "FL")
abline(Lys.mod, lty = 2, col = "gray")

plot(CCAA.GT$FL, CCAA.GT$Glx, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Trophic: Glx", xlab = "FL")
abline(Glx.mod, lty = 2, col = "gray")

plot(CCAA.GT$FL, CCAA.GT$Ala, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Trophic: Ala", xlab = "FL")
abline(Ala.mod, lty = 2, col = "gray")

plot(CCAA.GT$FL, CCAA.GT$Leu, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Trophic: Leu", xlab = "FL")
abline(Leu.mod, lty = 2, col = "gray")


#Plot Fl vs Carbon - no relationship
par(mfrow = c(1,1))
plot(GT.Red$FL, GT.Red$AdjC,col=GT.Red$Group,type="p", cex = 1.5,pch=20, xlab = "Fork Length", ylab = "Carbon")
C.model.GT <- lm(GT.Red$AdjC ~ GT.Red$FL)
summary(C.model.GT)
abline(C.model.GT, lty = 2, col = "gray")
#no relationship between Carbon and FL (p = 0.9714)

#Plot different N sources vs C
plot(CCAA.GT$AdjC, CCAA.GT$Phe, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Source: Phe", xlab = "C")
Phe.mod.C <- lm(CCAA.GT$Phe ~ CCAA.GT$AdjC)
summary(Phe.mod.C) #significant (p = 0.0407)
abline(Phe.mod.C, lty = 2, col = "gray")
with(CCAA.GT, text(CCAA.GT$Phe~CCAA.GT$AdjC, labels = CCAA.GT$Locality, cex = 0.5, pos = 4))

plot(CCAA.GT$AdjC, CCAA.GT$Gly, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Source: Gly", xlab = "C")
Gly.mod.C <- lm(CCAA.GT$Gly ~ CCAA.GT$AdjC)
summary(Gly.mod.C) #significant (p = 0.005957)
abline(Gly.mod.C, lty = 2, col = "gray")

plot(CCAA.GT$AdjC, CCAA.GT$Lys, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Source: Lys", xlab = "C")
Lys.mod.C <- lm(CCAA.GT$Lys ~ CCAA.GT$AdjC)
summary(Lys.mod.C) #non-significant (p = 0.92232)
abline(Lys.mod.C, lty = 2, col = "gray")

#trophic AAs
plot(CCAA.GT$AdjC, CCAA.GT$Glx, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Trophic: Glx", xlab = "C")
Glx.mod.C <- lm(CCAA.GT$Glx ~ CCAA.GT$AdjC)
summary(Glx.mod.C) #significant (p = 0.00161)
abline(Glx.mod.C, lty = 2, col = "gray")

plot(CCAA.GT$AdjC, CCAA.GT$Ala, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Trophic: Ala", xlab = "C")
Ala.mod.C <- lm(CCAA.GT$Ala ~ CCAA.GT$AdjC)
summary(Ala.mod.C) #significant (p = 0.000138)
abline(Ala.mod.C, lty = 2, col = "gray")

plot(CCAA.GT$AdjC, CCAA.GT$Leu, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Trophic: Leu", xlab = "C")
Leu.mod.C <- lm(CCAA.GT$Leu ~ CCAA.GT$AdjC)
summary(Leu.mod.C) #significant (p < 0.0001)
abline(Leu.mod.C, lty = 2, col = "gray")

## Model Carbon and Source and Trophic AA's by site
CCAA.Moz <- CCAA.dat[CCAA.dat$Locality=="PDO",]
CCAA.Mahe <- CCAA.dat[CCAA.dat$Locality=="MAHE",]
CCAA.StJo <- CCAA.dat[CCAA.dat$Locality=="STJO",]
CCAA.Prov <- CCAA.dat[CCAA.dat$Locality=="PROV",]

#boxplots by locality
boxplot(CCAA.Moz[,c(4:17)], main = "GT Amino Acids")
mycol <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue50")
boxplot(CCAA.Mahe[,c(4:17)], add = T, col = mycol)
mycol2 <- rgb(0, 0, 105, max = 105, alpha = 100, names = "brightblue")
boxplot(CCAA.StJo[,c(4:17)], add = T, col = mycol2)
mycol3 <- rgb(0, 0, 55, max = 155, alpha = 100, names = "gray")
boxplot(CCAA.Prov[,c(4:17)], add = T, col = mycol3)

#### Threonine ####
#plot GT Threonine by fish size
plot(CCAA.GT$Thr, CCAA.GT$FL,col = CCAA.GT$Locality, type = "p", pch = 20, main = "Threonine vs Fork Length: GTs")
legend("topright",legend=unique(CCAA.GT$Locality),
       pch=19,col=unique(CCAA.GT$Locality))
#plot BT threonine by fish size
plot(CCAA.BT$Thr~ CCAA.BT$FL,col = CCAA.BT$Locality, type = "p", pch = 20, main = "Threonine vs Fork Length: BTs")
#Plot of GTs AND BTs
plot(CCAA.dat$Thr, CCAA.dat$FL,col = CCAA.dat$Locality, type = "p", pch = 20, main = "Threonine vs Fork Length: GTs and BTs")
legend("topright",legend=unique(CCAA.dat$Locality),
       pch=19,col=unique(CCAA.dat$Locality))
#with(CCAA.GT, text(CCAA.dat$Thr~ CCAA.dat$FL, labels = CCAA.dat$Row.Labels, cex = 0.5, pos = 1))
FL.Thr.Mod <- lm(CCAA.dat$Thr~ CCAA.dat$FL)
summary(FL.Thr.Mod)
abline(FL.Thr.Mod, col = "gray", lty = 2) 


#Plot Threonine by trophic position
#GTs
plot(CCAA.GT$Thr, CCAA.GT$NewTP_Bradley2015,col = CCAA.GT$Locality, type = "p", pch = 20, main = "Threonine vs TP: GTs")
legend("topright",legend=unique(CCAA.GT$Locality),
       pch=19,col=unique(CCAA.GT$Locality))
#BTs
plot(CCAA.BT$Thr, CCAA.BT$NewTP_Bradley2015,col = CCAA.BT$Locality, type = "p", pch = 20, main = "Threonine vs TP: BTs")
legend("topright",legend=unique(CCAA.BT$Locality),
       pch=19,col=unique(CCAA.BT$Locality))
with(CCAA.BT, text(CCAA.BT$Thr, CCAA.BT$NewTP_Bradley2015, labels = CCAA.BT$FL, cex = 0.5, pos = 1))
#GTs and BTs
plot(CCAA.dat$Thr, CCAA.dat$NewTP_Bradley2015,col = CCAA.dat$Locality, type = "p", pch = 20, main = "Threonine vs TP: GTs and BTs", ylab = "TP (Bradley)")
legend("topright",legend=unique(CCAA.dat$Locality),
       pch=19,col=unique(CCAA.dat$Locality))
TP.Thr.Mod <- lm(CCAA.dat$Thr~ CCAA.dat$NewTP_Bradley2015)
summary(TP.Thr.Mod)
abline(TP.Thr.Mod, col = "gray", lty = 2) 

#Threonine boxplots by locality
boxplot(CCAA.Moz[c(1:16),c(15)], CCAA.Mahe[c(1:14),c(15)], CCAA.StJo[c(1:6),c(15)], CCAA.Prov[,c(15)], main = "Threonine: GTs and BTs", vertical = T, names = c("Moz", "Mahe", "StJo", "Prov"))


#Source AAs: Phe, Gly, Lys
#Trophic AAs: Glx, Ala, Leu

#Phenylalanine by locality
Phe.mod.C.Moz <- lm(CCAA.Moz$Phe ~ CCAA.Moz$AdjC)
summary(Phe.mod.C.Moz) #significant (p = 0.01871)
abline(Phe.mod.C.Moz, lty = 2, col = "cyan")

Phe.mod.C.Mahe <- lm(CCAA.Mahe$Phe ~ CCAA.Mahe$AdjC)
summary(Phe.mod.C.Mahe) #slightly significant (p = 0.096709)
abline(Phe.mod.C.Mahe, lty = 2, col = "red")

Phe.mod.C.StJo <- lm(CCAA.StJo$Phe ~ CCAA.StJo$AdjC)
summary(Phe.mod.C.StJo) #non-significant (p = 0.254)
abline(Phe.mod.C.StJo, lty = 2, col = "yellow")

Phe.mod.C.Prov <- lm(CCAA.Prov$Phe ~ CCAA.Prov$AdjC)
summary(Phe.mod.C.Prov) #non-significant (p = 0.32263)
abline(Phe.mod.C.Prov, lty = 2, col = "pink")

#Gly by locality
Gly.mod.C.Moz <- lm(CCAA.Moz$Gly ~ CCAA.Moz$AdjC)
summary(Gly.mod.C.Moz) #significant (p = 0.0670)
abline(Gly.mod.C.Moz, lty = 2, col = "cyan")

Gly.mod.C.Mahe <- lm(CCAA.Mahe$Gly ~ CCAA.Mahe$AdjC)
summary(Gly.mod.C.Mahe) #significant (p = < 0.001)
abline(Gly.mod.C.Mahe, lty = 2, col = "red")

Gly.mod.C.StJo <- lm(CCAA.StJo$Gly ~ CCAA.StJo$AdjC)
summary(Gly.mod.C.StJo) #non-significant (p = 0.811)
abline(Gly.mod.C.StJo, lty = 2, col = "yellow")

Gly.mod.C.Prov <- lm(CCAA.Prov$Gly ~ CCAA.Prov$AdjC)
summary(Gly.mod.C.Prov) #non-significant (p = 0.726)
abline(Gly.mod.C.Prov, lty = 2, col = "pink")

#Lysine by locality
Lys.mod.C.Moz <- lm(CCAA.Moz$Lys ~ CCAA.Moz$AdjC)
summary(Lys.mod.C.Moz) #non-significant (p = 0.838)
abline(Lys.mod.C.Moz, lty = 2, col = "cyan")

Lys.mod.C.Mahe <- lm(CCAA.Mahe$Lys ~ CCAA.Mahe$AdjC)
summary(Lys.mod.C.Mahe) #significant (p = < 0.0433)
abline(Lys.mod.C.Mahe, lty = 2, col = "red")

Lys.mod.C.StJo <- lm(CCAA.StJo$Lys ~ CCAA.StJo$AdjC)
summary(Lys.mod.C.StJo) #non-significant (p = 0.7861)
abline(Lys.mod.C.StJo, lty = 2, col = "yellow")

Lys.mod.C.Prov <- lm(CCAA.Prov$Lys ~ CCAA.Prov$AdjC)
summary(Lys.mod.C.Prov) #non-significant (p = 0.9854)
abline(Lys.mod.C.Prov, lty = 2, col = "pink")

#Glutamine by locality
Glx.mod.C.Moz <- lm(CCAA.Moz$Glx ~ CCAA.Moz$AdjC)
summary(Glx.mod.C.Moz) #non significant (p = 0.58986)
abline(Glx.mod.C.Moz, lty = 2, col = "cyan")

Glx.mod.C.Mahe <- lm(CCAA.Mahe$Glx ~ CCAA.Mahe$AdjC)
summary(Glx.mod.C.Mahe) #non-significant (p = 0.179)
abline(Glx.mod.C.Mahe, lty = 2, col = "red")

Glx.mod.C.StJo <- lm(CCAA.StJo$Glx ~ CCAA.StJo$AdjC)
summary(Glx.mod.C.StJo) #slightly significant (p = 0.126679)
abline(Glx.mod.C.StJo, lty = 2, col = "yellow")

Glx.mod.C.Prov <- lm(CCAA.Prov$Glx ~ CCAA.Prov$AdjC)
summary(Glx.mod.C.Prov) #slightly significant (p = 0.259)
abline(Glx.mod.C.Prov, lty = 2, col = "pink")

#Alanine by locality
Ala.mod.C.Moz <- lm(CCAA.Moz$Ala ~ CCAA.Moz$AdjC)
summary(Ala.mod.C.Moz) #non significant (p = 0.1859)
abline(Ala.mod.C.Moz, lty = 2, col = "cyan")

Ala.mod.C.Mahe <- lm(CCAA.Mahe$Ala ~ CCAA.Mahe$AdjC)
summary(Ala.mod.C.Mahe) #non significant (p = 0.269)
abline(Ala.mod.C.Mahe, lty = 2, col = "red")

Ala.mod.C.StJo <- lm(CCAA.StJo$Ala ~ CCAA.StJo$AdjC)
summary(Ala.mod.C.StJo) #slightly significant (p = 0.1895)
abline(Ala.mod.C.StJo, lty = 2, col = "yellow")

Ala.mod.C.Prov <- lm(CCAA.Prov$Ala ~ CCAA.Prov$AdjC)
summary(Ala.mod.C.Prov) #significant (p = 0.0511)
abline(Ala.mod.C.Prov, lty = 2, col = "pink")

#Leusine by locality
dev.off()
par(mfrow = c(1,1))
Leu.mod.C.Moz <- lm(CCAA.Moz$Leu ~ CCAA.Moz$AdjC)
summary(Leu.mod.C.Moz) #slightly significant (p = 0.0964)
abline(Leu.mod.C.Moz, lty = 2, col = "cyan")

Leu.mod.C.Mahe <- lm(CCAA.Mahe$Leu ~ CCAA.Mahe$AdjC)
summary(Leu.mod.C.Mahe) #slightly significant (p = 0.0238)
abline(Leu.mod.C.Mahe, lty = 2, col = "red")

Leu.mod.C.StJo <- lm(CCAA.StJo$Leu ~ CCAA.StJo$AdjC)
summary(Leu.mod.C.StJo) #slightly significant (p = 0.067082)
abline(Leu.mod.C.StJo, lty = 2, col = "yellow")

Leu.mod.C.Prov <- lm(CCAA.Prov$Leu ~ CCAA.Prov$AdjC)
summary(Leu.mod.C.Prov) #slightly significant (p = 0.06372)
abline(Leu.mod.C.Prov, lty = 2, col = "pink")


#Fit a linear model of Bulk N vs Nphe (predictor) - see Hetherington et al. 2017
dev.off()
par(mfrow = c(1,1))
plot(CCAA.GT$BulkN ~CCAA.GT$Phe, col = CCAA.GT$Locality, type = "p", pch = 20, main = "Bulk N vs N(Phe)")
bulkN.mod <- lm(CCAA.GT$BulkN ~ CCAA.GT$Phe)
summary(bulkN.mod)
with(CCAA.GT, text(CCAA.GT$Phe, CCAA.GT$BulkN, labels = CCAA.GT$Locality, cex = 0.5, pos = 4))
abline(bulkN.mod, lty = 2, col = "gray")
legend("topleft",legend=unique(CCAA.GT$Locality),
       pch=19,col=unique(CCAA.GT$Locality))
#significant effect of source (Phe) on Bulk N (p = 0.0001018)

plot(CCAA.GT.Red$BulkN ~CCAA.GT.Red$Phe, col = CCAA.GT.Red$Locality, type = "p", pch = 20, main = "Bulk N vs N(Phe)")
bulkN.mod <- lm(CCAA.GT.Red$BulkN ~ CCAA.GT.Red$Phe)
summary(bulkN.mod)
abline(bulkN.mod, lty = 2, col = "gray")
legend("topleft",legend=unique(CCAA.GT.Red$Locality),
       pch=19,col=unique(CCAA.GT.Red$Locality))

CCAA.GT.Moz <- CCAA.GT.Red[CCAA.GT.Red$Locality=="PDO",]
CCAA.GT.Mahe <- CCAA.GT.Red[CCAA.GT.Red$Locality=="MAHE",]
CCAA.GT.StJo <- CCAA.GT.Red[CCAA.GT.Red$Locality=="STJO",]
CCAA.GT.Prov <- CCAA.GT.Red[CCAA.GT.Red$Locality=="PROV",]

bulkN.mod.Moz <- lm(CCAA.GT.Moz$BulkN ~ CCAA.GT.Moz$Phe)
summary(bulkN.mod.Moz)
abline(bulkN.mod.Moz, lty = 2, col = "red")

bulkN.mod.Mahe <- lm(CCAA.GT.Mahe$BulkN ~ CCAA.GT.Mahe$Phe)
summary(bulkN.mod.Mahe)
abline(bulkN.mod.Mahe, lty = 2, col = "black")

bulkN.mod.StJo <- lm(CCAA.GT.StJo$BulkN ~ CCAA.GT.StJo$Phe)
summary(bulkN.mod.StJo)
abline(bulkN.mod.StJo, lty = 2, col = "blue")

bulkN.mod.Prov <- lm(CCAA.GT.Prov$BulkN ~ CCAA.GT.Prov$Phe)
summary(bulkN.mod.Prov)
abline(bulkN.mod.Prov, lty = 2, col = "green")

bulkN.Ancova.mod <- aov(CCAA.GT.Red$BulkN ~ CCAA.GT.Red$Phe * CCAA.GT.Red$Locality)
summary.aov(bulkN.Ancova.mod)


#plot the other source AAs (appear to also have a significant positive relationship)
plot(CCAA.dat$BulkN ~CCAA.dat$Gly, col = CCAA.dat$Locality, type = "p", pch = 20)
plot(CCAA.dat$BulkN ~CCAA.dat$Lys, col = CCAA.dat$Locality, type = "p", pch = 20)
with(CCAA.dat, text(CCAA.dat$Lys, CCAA.dat$BulkN, labels = CCAA.dat$Row.Labels, cex = 0.5, pos = 4))

#Examine BulkN and FL
plot(CCAA.GT$BulkN ~ CCAA.GT$FL, col = CCAA.GT$Locality, type = "p", pch = 20)
bulkN_FL.mod <- lm(CCAA.GT$BulkN ~ CCAA.GT$FL)
summary(bulkN_FL.mod)
#negative correlation (p = 0.0604)
abline(bulkN_FL.mod, lty = 2, col = "gray")

#Compare Bluefin and GTs - Mahe and Providence
dev.off()

CCAA.GT.Moz <- CCAA.GT[CCAA.GT$Locality=="PDO",]
CCAA.GT.Mahe <- CCAA.GT[CCAA.GT$Locality=="MAHE",]
CCAA.GT.StJo <- CCAA.GT[CCAA.GT$Locality=="STJO",]
CCAA.GT.Prov <- CCAA.GT[CCAA.GT$Locality=="PROV",]

CCAA.BT.Prov <- CCAA.BT[CCAA.BT$Locality =="PROV",]

range(CCAA.GT.Moz$NewTP_Bradley2015)
range(CCAA.GT.Mahe$NewTP_Bradley2015)
range(CCAA.GT.StJo$NewTP_Bradley2015)
range(CCAA.GT.Prov$NewTP_Bradley2015)

range(CCAA.GT.Prov$BulkN)
range(CCAA.GT.Prov$AdjC)
range(CCAA.GT.Mahe$BulkN)
range(CCAA.GT.Mahe$AdjC)

dev.off()
#Plot GT and BT from Providence C and N values
plot(CCAA.GT.Prov$AdjC, CCAA.GT.Prov$BulkN, col = 1, type = "p", pch = 22, xlim = c(-15,-10), ylim = c(12,14))
points(CCAA.BT.Prov$AdjC, CCAA.BT.Prov$BulkN, col = 4, type = "p", pch = 21)

#Plot GT and BT from Providence TP and FL values
plot(CCAA.GT.Prov$NewTP_Bradley2015 ~ CCAA.GT.Prov$FL, col = 1, type = "p", pch = 16, xlim = c(58,95), ylim = c(3,4.5), xlab = "FL", ylab = "TP_Bradley", cex = 1.2)
points(CCAA.BT.Prov$NewTP_Bradley2015 ~ CCAA.BT.Prov$FL, col = 4, type = "p", pch = 16, cex = 1.2)
legend("bottomright",legend=unique(SIBER.Red$group),
       pch=19,col=c(1,4))

#Show a boxplot of TP between GTs and BTs from Providence
boxplot(CCAA.GT.Prov$NewTP_Bradley2015, CCAA.BT.Prov$NewTP_Bradley2015)

#Plot GT and BT from Providence TP and C values
plot(CCAA.GT.Prov$NewTP_Bradley2015 ~ CCAA.GT.Prov$AdjC, col = 1, type = "p", pch = 16, xlim = c(-15,-10), ylim = c(3.2,4.4), xlab = "C", ylab = "TP_Bradley", cex = 1.2, main = "Providence BTs and GTs")
CCAA.GT.Prov.C.mod <- lm(CCAA.GT.Prov$NewTP_Bradley2015 ~ CCAA.GT.Prov$AdjC)
summary(CCAA.GT.Prov.C.mod)
#non signficant
abline(CCAA.GT.Prov.C.mod, lty = 2, col = "black")
points(CCAA.BT.Prov$NewTP_Bradley2015 ~ CCAA.BT.Prov$AdjC, col = 4, type = "p", pch = 16, cex = 1.2)
legend("bottomleft",legend=unique(SIBER.Red$group),
       pch=19,col=c(1,4))
CCAA.BT.Prov.C.mod <- lm(CCAA.BT.Prov$NewTP_Bradley2015 ~ CCAA.BT.Prov$AdjC)
summary(CCAA.BT.Prov.C.mod)
abline(CCAA.BT.Prov.C.mod, lty = 2, col = "blue")

#Boxplot of Source and Trophic AAs: GTs and BTs
boxplot(CCAA.GT.Prov[,c(4:17)], main = "Amino Acids")
#Create a transparent color and add BTs
mycol <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue50")
boxplot(CCAA.BT.Prov[,c(4:17)], add = T, col = mycol)


#### Conduct SIBER analyses with BT and GT data ####
#Made separate dataset in Excel
#use bulk N and C
SIBER_BT.dat <- read.csv("SIBER_BTGT_Data.csv")
#Set the seed so results will be the same for subsequent runs
set.seed(1)

BT.SIBER <- createSiberObject(SIBER_BT.dat)
#Due to a low sample size N < 5, drop the Mahe samples (for now) - see warning that pops up

SIBER.Red <- SIBER_BT.dat[SIBER_BT.dat$community == "PROV",]
#Need to change the community names to be different too
comm <- as.vector(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2))
SIBER.Red$community <- comm
Prov.SIBER <- createSiberObject(SIBER.Red)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions. - Can change p.interval to represent X % of the data (e.g. p = 0.95 represents 95% of the data) or setting P to null produces  ML standard ellipses (40%)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


par(mfrow=c(1,1))
plotSiberObject(Prov.SIBER,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

legend("bottomleft",legend=unique(SIBER.Red$group),
       pch=19,col=1:length(unique(SIBER.Red$group)))
#Levels: 1 (red) = Moz, 2 (black) = Mahe, 3 (blue) = St. Jo, 4 (green) = Prov


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(BT.SIBER, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(Prov.SIBER, n = 100, p.interval = 0.80, ci.mean = T,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(Prov.SIBER)
print(group.ML)

#         1.GT     2.BT
#TA   2.845000 2.415000
#SEA  1.499084 1.760055
#SEAc 1.713239 2.200068


# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(Prov.SIBER, parms, priors)

#Comparing among groups using the Standard Ellipse Area

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

anova(group.ML)
class(group.ML)

#Perform a t-test of the siberEllipses posteriors
class(SEA.B)
SEA.B.df <- as.data.frame(SEA.B)
t.test(SEA.B.df$V1, SEA.B.df$V2)
#Supposedley these differences are significnt
#mean of x mean of y 
#1.689600  2.385336 

#### BT GT SIBER comparisons using TPwtmns ####
SIBER_BT.dat <- read.csv("SIBER_BTGT_Data_TPwtmns.csv")
#Set the seed so results will be the same for subsequent runs
set.seed(1)

BT.SIBER <- createSiberObject(SIBER_BT.dat)
#Due to a low sample size N < 5, drop the Mahe samples (for now) - see warning that pops up

SIBER.Red <- SIBER_BT.dat[SIBER_BT.dat$community == "PROV",]
#Need to change the community names to be different too
comm <- as.vector(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2))
SIBER.Red$community <- comm
Prov.SIBER <- createSiberObject(SIBER.Red)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions. - Can change p.interval to represent X % of the data (e.g. p = 0.95 represents 95% of the data) or setting P to null produces  ML standard ellipses (40%)
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


par(mfrow=c(1,1))
plotSiberObject(Prov.SIBER,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                y.limits = as.vector(c(3.5, 4.4)),
                x.limits = as.vector(c(-15,-10)),
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = "Trophic Position (weighted means)"
)

#red = GT, black = BT
legend("bottomleft",legend=unique(SIBER.Red$group),
       pch=19,col=c("Red", "black"))

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(BT.SIBER, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(Prov.SIBER, n = 100, p.interval = 0.80, ci.mean = T,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(Prov.SIBER)
print(group.ML)

#         1.BT     2.GT
#TA   1.817983 1.1004178
#SEA  1.146651 0.9451659
#SEAc 1.310458 1.1814574

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(Prov.SIBER, parms, priors)

#Comparing among groups using the Standard Ellipse Area

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

anova(group.ML)
class(group.ML)

#Perform a t-test of the siberEllipses posteriors
class(SEA.B)
SEA.B.df <- as.data.frame(SEA.B)
t.test(SEA.B.df$V1, SEA.B.df$V2)
#These differences are not significnt
#mean of x   mean of y 
#1.298424    1.286939  

#### Weighted Means ####
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA/CCSIA_Analyses")
WtMns <- read.csv("WeightedMeans.csv")

#Drop Mtentu and Mauritius because of low sample size
WtMns.Red <- WtMns[-c(16:17,45),]
WtMns.Red$Locality <- droplevels(WtMns.Red$Locality, exclude = c("MAU", "MTENTU"))
WtMns.Red$Locality

head(WtMns.Red)

par(mfrow = c(1,1))
boxplot(WtMns.Red$WtMnSrc ~ WtMns.Red$Locality, main = "Weighted Means- Source AA")
boxplot(WtMns.Red$WtMnTrp ~ WtMns.Red$Locality, main = "Weighted Means- Trophic AA")
boxplot(WtMns.Red$DeltaTrpSrc ~ WtMns.Red$Locality, main = "Delta Trophic - Source")

plot(WtMns.Red$WtMnSrc ~ WtMns.Red$Locality)

#Plot weighted means trophic and source vs C
plot(WtMns.Red$WtMnSrc~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), main = "Weighted Means (Source) vs C", xlab="Carbon", ylab = "Weighted Means (source)")
WtMnsSrc.C.mod <- lm(WtMns.Red$WtMnSrc~ WtMns.Red$AdjC)
summary(WtMnsSrc.C.mod) #p = 0.0804
abline(WtMnsSrc.C.mod, lty = 2, col = "gray")

plot(WtMns.Red$WtMnTrp~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), main = "Weighted Means (Trophic) vs C", xlab="Carbon", ylab = "Weighted Means (trophic)")
WtMnsTrp.C.mod <- lm(WtMns.Red$WtMnTrp~ WtMns.Red$AdjC)
summary(WtMnsTrp.C.mod) #p < 0.001
abline(WtMnsTrp.C.mod, lty = 2, col = "gray")

#calculate the Weighted Mean Source value by locality
WtMn.Loc <- by(WtMns.Red$WtMnSrc, WtMns.Red$Locality, FUN = mean, simplify = F)
str(WtMn.Loc)
class(WtMn.Loc)

WtMn.Loc[c(1, 4:6)]
boxplot(WtMn.Loc[c(1, 4:6)])

#calculate the Weighted Mean Delta Trophic - Source by locality
WtMn.delta.Loc <- by(WtMns.Red$DeltaTrpSrc, WtMns.Red$Locality, FUN = mean, simplify = F)
WtMn.delta.Loc
class(WtMn.delta.Loc)

WtMn.delta.Loc[c(1, 4:6)]
boxplot(WtMn.delta.Loc[c(1, 4:6)])

head(WtMns.Red)

#### Weighted Means TP ####
#plot weighted means TP versus TP using Bradley 2015
plot(WtMns.Red$TP_WtMn, WtMns.Red$NewTP_Bradley2015, col = as.factor(WtMns.Red$Locality), main = "Weighted Means TP vs Bradley TP", xlab="Weighted Means TP", ylab = "Bradley TP")
abline(0,1)
#The Weighted Means TP values are almost always larger than  Bradley 2015 values 

#Plot weighted means TP versus TP using AA averages and Bradley 2015 values
plot(WtMns.Red$TP_WtMn, WtMns.Red$TP_Avgs, col = as.factor(WtMns.Red$Locality), main = "Weighted Means TP vs AA_Averages TP", xlab="Weighted Means TP", ylab = "AA_Averages TP")
abline(0,1)
#This is more evenly split

#Plot AA Averages TP versus Single AA TP
plot(WtMns.Red$NewTP_Bradley2015, WtMns.Red$TP_Avgs, col = as.factor(WtMns.Red$Locality), main = "Single AA TP vs AA Averages TP", xlab="Single AA TP", ylab = "AA Averages TP")
abline(0,1)
#AA Averages TP values are generally larger than Single AA TP values


### Plot all TPs vs FL ###
plot(WtMns.Red$TP_WtMn ~ WtMns.Red$FL, col = as.factor(WtMns.Red$Locality), main = "Trophic Positions vs FL", xlab="Fork Length", ylab = "Trophic Positions", pch = 15, cex = 1.3, xlim = c(20,120), ylim = c(2.7, 5.2))
points(WtMns.Red$NewTP_Bradley2015 ~ WtMns.Red$FL, col = as.factor(WtMns.Red$Locality), pch = 12, cex = 1.3)
points(WtMns.Red$TP_Avgs ~ WtMns.Red$FL, col = as.factor(WtMns.Red$Locality), pch = 4, cex = 1.3)
points(WtMns.Red$TrophicPosition ~ WtMns.Red$FL, col = as.factor(WtMns.Red$Locality), pch = 2, cex = 1.3)
legend("topright",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality))
legend("topleft",legend=as.factor(paste(c("TP_WtMn", "TP_SingleAA", "TP_Avgs", "TP_Chikarashi"))),pch=c(15,12,4,2), bty = "n")

### Plot all TPs vs C ###
plot(WtMns.Red$TP_WtMn ~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), main = "Trophic Positions vs C", xlab="AdjC", ylab = "Trophic Positions", pch = 15, cex = 1.3, xlim = c(-18,-10), ylim = c(2.7, 5.2))
points(WtMns.Red$NewTP_Bradley2015 ~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), pch = 12, cex = 1.3)
points(WtMns.Red$TP_Avgs ~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), pch = 4, cex = 1.3)
points(WtMns.Red$TrophicPosition ~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), pch = 2, cex = 1.3)
legend("topright",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality), bty = "n", inset = (c(0.3,0)))
legend("topright",legend=as.factor(paste(c("TP_WtMn", "TP_SingleAA", "TP_Avgs", "TP_Chikarashi"))),pch=c(15,12,4,2), bty = "n")

head(WtMns.Red)
dim(WtMns.Red)
## Calculate TP ranges ##
TP.Ranges <- sapply(WtMns.Red[,c(34,35,40,41),], FUN = mean)
range(WtMns.Red[,c(34,35,40,41),])

#Calculate TP summary stats 

#Weighted means
WtMns.Red$TP_WtMn
mean(WtMns.Red$TP_WtMn)
#4.253129
sd(WtMns.Red$TP_WtMn)
#0.3680935
range(WtMns.Red$TP_WtMn)
#3.530494 5.084529

#Bradley 2015 - Single AA
mean(WtMns.Red$NewTP_Bradley2015)
#4.012501
sd(WtMns.Red$NewTP_Bradley2015)
#0.2998395
range(WtMns.Red$NewTP_Bradley2015)
#3.410614 4.723509

#TP Averages
mean(WtMns.Red$TP_Avgs)
#4.197721
sd(WtMns.Red$TP_Avgs)
#0.2764886
range(WtMns.Red$TP_Avgs)
#3.511111 4.841520

#TP Chikaraishi
mean(WtMns.Red$TrophicPosition)
#3.285691
sd(WtMns.Red$TrophicPosition)
#0.2248796
range(WtMns.Red$TrophicPosition)
#2.834276 3.818947


##Use an ANOVA to compare TP methods
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA/CCSIA_Analyses")
TP.Anova <- read.csv("TP_ANOVA.csv")
head(TP.Anova)

#install.packages("ggpubr")
library("ggpubr")
ggline(TP.Anova, x = "Group", y = "TP", merge = F, group = 1,
       add = c("mean_se", "jitter"), 
       ylab = "Method", xlab = "TP")

# Compute the analysis of variance
TP.aov <- aov(TP ~ Group, data = TP.Anova)
# Summary of the analysis
summary(TP.aov) #There are significant differences between groups

#Pairwise comparison between groups
TukeyHSD(TP.aov)
#significant differences between all methods except WtMns and Avgs


#Look at effect of fish size (FL)
plot(WtMns.Red$WtMnSrc~ WtMns.Red$FL, col = as.factor(WtMns.Red$Locality), main = "WeightedMeans(SourceAAs) vs FL")
WtMnSrc.Mod <- lm(WtMns.Red$WtMnSrc~ WtMns.Red$FL)
summary(WtMnSrc.Mod)
abline(WtMnSrc.Mod, lty = 2, col = "gray")
legend("bottomleft",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality))
#Slightly significant (p = 0.0635) effect of FL on WtMnSource

#add year
with(WtMns.Red, text(WtMns.Red$FL, WtMns.Red$WtMnSrc, labels = WtMns.Red$Year, cex = 0.5, pos = 4))

plot(WtMns.Red$WtMnTrp~ WtMns.Red$FL, col = as.factor(WtMns.Red$Locality), main = "WeightedMeans(TrophicAAs) vs FL")
WtMnTrp.Mod <- lm(WtMns.Red$WtMnTrp~ WtMns.Red$FL)
summary(WtMnTrp.Mod)
abline(WtMnTrp.Mod, lty = 2, col = "gray")
legend("bottomleft",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality))
#no effect of FL on WtMnTrophic

## This effect could potentially be driving the differences seen in Mahe compared to the other sites, since Mahe contains mostly juveniles < 60 cm (all but 3 samples)

## Look at effects of Adjusted C
par(mfrow = c(1,1))
plot(WtMns.Red$WtMnSrc~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), main = "WeightedMeans(SourceAAs) vs C")
WtMnSrc.Mod.C <- lm(WtMns.Red$WtMnSrc~ WtMns.Red$AdjC)
summary(WtMnSrc.Mod.C)
abline(WtMnSrc.Mod.C, lty = 2, col = "gray")
legend("topright",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality))
#slightly positive effect of C on source AA weighted means (p = 0.08036)

plot(WtMns.Red$WtMnTrp~ WtMns.Red$AdjC, col = as.factor(WtMns.Red$Locality), main = "WeightedMeans(TrophicAAs) vs C")
WtMnTrp.Mod.C <- lm(WtMns.Red$WtMnTrp~ WtMns.Red$AdjC)
summary(WtMnTrp.Mod.C)
abline(WtMnTrp.Mod.C, lty = 2, col = "gray")
legend("topright",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality))
#significant p < 0.001

### Test weighted means and C for each locality ###
WtMns.Moz <- WtMns.Red[c(WtMns.Red$Locality == "PDO"),]
WtMns.Mahe <- WtMns.Red[c(WtMns.Red$Locality == "MAHE"),]
WtMns.StJo <- WtMns.Red[c(WtMns.Red$Locality == "STJO"),]
WtMns.Prov <- WtMns.Red[c(WtMns.Red$Locality == "PROV"),]

#test source
#mozambique
WtMnSrc.Mod.C.Moz <- lm(WtMns.Moz$WtMnSrc~ WtMns.Moz$AdjC)
summary(WtMnSrc.Mod.C.Moz) #non significant (p = 0.499)
abline(WtMnSrc.Mod.C.Moz, lty = 2, col = "red")

plot(WtMns.Moz$WtMnSrc)

#mahe
WtMnSrc.Mod.C.Mahe <- lm(WtMns.Mahe$WtMnSrc~ WtMns.Mahe$AdjC)
summary(WtMnSrc.Mod.C.Mahe) #non significant (p = 0.1388)
abline(WtMnSrc.Mod.C.Mahe, lty = 2, col = "black")

#StJo
WtMnSrc.Mod.C.StJo <- lm(WtMns.StJo$WtMnSrc~ WtMns.StJo$AdjC)
summary(WtMnSrc.Mod.C.StJo) #non significant (p = 0.671)
abline(WtMnSrc.Mod.C.StJo, lty = 2, col = "blue")

#Providence
WtMnSrc.Mod.C.Prov <- lm(WtMns.Prov$WtMnSrc~ WtMns.Prov$AdjC)
summary(WtMnSrc.Mod.C.Prov) #non significant (p = 0.4403)
abline(WtMnSrc.Mod.C.Prov, lty = 2, col = "green")

#Test Trophic
#mozambique
WtMnTrp.Mod.C.Moz <- lm(WtMns.Moz$WtMnTrp~ WtMns.Moz$AdjC)
summary(WtMnTrp.Mod.C.Moz) #non-significant (p = 0.195)
abline(WtMnTrp.Mod.C.Moz, lty = 2, col = "red")

#mahe
WtMnTrp.Mod.C.Mahe <- lm(WtMns.Mahe$WtMnTrp~ WtMns.Mahe$AdjC)
summary(WtMnTrp.Mod.C.Mahe) #slightly significant (p = 0.0852)
abline(WtMnTrp.Mod.C.Mahe, lty = 2, col = "black")

#stjo
WtMnTrp.Mod.C.StJo <- lm(WtMns.StJo$WtMnTrp~ WtMns.StJo$AdjC)
summary(WtMnTrp.Mod.C.StJo) #significant (p = 0.0315)
abline(WtMnTrp.Mod.C.StJo, lty = 2, col = "blue")

#providence
WtMnTrp.Mod.C.Prov <- lm(WtMns.Prov$WtMnTrp~ WtMns.Prov$AdjC)
summary(WtMnTrp.Mod.C.Prov) #non-significant (p = 0.692)
abline(WtMnTrp.Mod.C.Prov, lty = 2, col = "green")

#### Use Weighted Means TP for analyses ####
#Moz
WtMns.Moz$TP_WtMn
mean(WtMns.Moz$TP_WtMn)
#mean: 4.531495
range(WtMns.Moz$TP_WtMn)
#4.036340 5.084529
#1.048189
sd(WtMns.Moz$TP_WtMn)
#SD: 0.2900845

#Mahe
WtMns.Mahe$TP_WtMn
mean(WtMns.Mahe$TP_WtMn)
#mean: 4.079328
range(WtMns.Mahe$TP_WtMn)
#range:3.530494 4.827766
#1.297272
sd(WtMns.Mahe$TP_WtMn)
#SD: 0.3761083

#St Jo
WtMns.StJo$TP_WtMn
mean(WtMns.StJo$TP_WtMn)
#mean: 4.14974
range(WtMns.StJo$TP_WtMn)
#range: 3.792249 4.643574
# 0.851325
sd(WtMns.StJo$TP_WtMn)
#SD:0.231134

#Providence
WtMns.Prov$TP_WtMn
mean(WtMns.Prov$TP_WtMn)
#mean: 4.039085
range(WtMns.Prov$TP_WtMn)
#range: 3.571320 4.343097
#0.771777
sd(WtMns.Prov$TP_WtMn)
#SD: 0.273632

#Plot TP vs FL
plot(WtMns.Red$FL, WtMns.Red$TP_WtMn,col=WtMns.Red$Locality,type="p",pch=20, cex = 1.5, ylab="TP WtMeans", main = "GT FL vs CCSIA Wt Mns TP")
legend("topleft",legend=as.factor(paste(unique(WtMns.Red$Locality))),
       pch=19,col=unique(WtMns.Red$Locality))

TPwm.FL.mod <- lm(WtMns.Red$TP_WtMn ~ WtMns.Red$FL)
summary(TPwm.FL.mod)

abline(TPwm.FL.mod, col = "gray", pch = 2, lty = 2)

#Measure TP weighted means vs FL by site
TPwm.FL.mod.Moz <- lm(WtMns.Moz$TP_WtMn ~ WtMns.Moz$FL)
summary(TPwm.FL.mod.Moz)
#p= 0.768, Radj = -0.05324
abline(TPwm.FL.mod.Moz, col = "red", pch = 2, lty = 2)

TPwm.FL.mod.Mahe <- lm(WtMns.Mahe$TP_WtMn ~ WtMns.Mahe$FL)
summary(TPwm.FL.mod.Mahe)
#p= 0.0543, Radj = 0.1988
abline(TPwm.FL.mod.Mahe, col = "black", pch = 2, lty = 2)

TPwm.FL.mod.StJo <- lm(WtMns.StJo$TP_WtMn ~ WtMns.StJo$FL)
summary(TPwm.FL.mod.StJo)
#p= 0.2118, Radj = 0.08555
abline(TPwm.FL.mod.StJo, col = "blue", pch = 2, lty = 2)

TPwm.FL.mod.Prov <- lm(WtMns.Prov$TP_WtMn ~ WtMns.Prov$FL)
summary(TPwm.FL.mod.Prov)
#p= 0.7935, Radj = -0.1309 
abline(TPwm.FL.mod.Prov, col = "green", pch = 2, lty = 2)

#### Compare BTs and GTs - weighted means ####

#GT TP from Prov
WtMns.Prov$TP_WtMn

mean(WtMns.Prov$TP_WtMn)
#4.039085
sd(WtMns.Prov$TP_WtMn)
#0.273632
range(WtMns.Prov$TP_WtMn)
#3.571320 4.343097

#BT TP from Prov
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA/CCSIA_Analyses")
WtMns.BT <- read.csv("WeightedMeans_BT.csv")

WtMns.BT.Prov <- WtMns.BT[WtMns.BT$Locality == "PROV",]

mean(WtMns.BT.Prov$TP_WtMn)
#3.957378
sd(WtMns.BT.Prov$TP_WtMn)
#0.2812004
range(WtMns.BT.Prov$TP_WtMn)
#3.519625 4.252643

#Plot BT and GT TP vs Fork Length
plot(WtMns.Prov$TP_WtMn~ WtMns.Prov$FL, xlim = c(55,100), ylim = c(3.5, 4.4), pch = 19, main = "GT and BT TPwtmns vs FL")
points(WtMns.BT.Prov$TP_WtMns ~ WtMns.BT.Prov$FL, col = "blue", pch = 19)
legend("bottomright",legend= c("GT","BT"), pch=19,col= c("black", "blue"))

#Plot BT and GT TP vs C
plot(WtMns.Prov$TP_WtMn~ WtMns.Prov$AdjC, ylim = c(3.5, 4.4), xlim = c(-15, -10), pch = 19, main = "GT and BT TPwtmns vs C")
points(WtMns.BT.Prov$TP_WtMns ~ WtMns.BT.Prov$AdjC, col = "blue", pch = 19)
legend("bottomleft",legend= c("GT","BT"), pch=19,col= c("black", "blue"))
#Add FL labels
with(WtMns.BT.Prov, text(WtMns.BT.Prov$TP_WtMns ~ WtMns.BT.Prov$AdjC, labels = WtMns.BT.Prov$FL, cex = 0.5, pos = 3))
with(WtMns.Prov, text(WtMns.Prov$TP_WtMn ~ WtMns.Prov$AdjC, labels = WtMns.Prov$FL, cex = 0.5, pos = 3))


#### Baselines ####
baselines <- read.csv("baselines.csv")
baselines

#Look at boxplot of AAs
boxplot(baselines[,c(4:17)], main = "Baseline Amino Acids", names = c("Ala#", "Asx#", "Glx#", "Gly+", "Ile#", "Leu#", "Lys+", "Met+", "Phe+", "Pro#", "Ser+", "Thr*", "Tyr+", "Val#"), xlab = "# trophic  + source *metabolic")
#really high values coming from Mauritius B1 (barnacle)

#plot TP versus C (NOT using Bradley 2015 TEF values)
plot(baselines$Trophic.Position~ baselines$X13C, pch = 15, cex = 1.3, col = as.factor(baselines$Locality), main = "TP versus C for Baselines")
legend("topleft",legend=as.factor(paste(unique(baselines$Locality))),
       pch=19,col=unique(baselines$Locality))
with(baselines, text(baselines$Trophic.Position~baselines$X13C, labels = baselines$Row.Labels, cex = 0.5, pos = 1))


### Clustering and PCA ####
#separate source and trophic (and get rid of last line of NAs)
CCAA.GT.Src <- CCAA.GT[c(1:57),c(1,2,3,7,10,12,18)]
CCAA.GT.Tro <- CCAA.GT[c(1:57),c(1,2,3,4,6,9,18)]

## Apply PCA to source AAs ##
GT.pca.src = stats :: prcomp(CCAA.GT.Src[,c(4:6)])
GT.pca.src

#loadings
GT.pca.src$rotation
#scores (PCs)
GT.pca.src$x

#basic plot
plot(GT.pca.src$x[,1], GT.pca.src$x[,2])

# load ggplot2
library(ggplot2)
#install.packages("ggfortify")
library(ggfortify)
library(grid)
library(gridExtra)
library(cluster)

#Create a new dataframe with sample information
df_out <- as.data.frame(GT.pca.src$x)
df_out$samples <- CCAA.GT.Src$Row.Labels
df_out$locality <- CCAA.GT.Src$Locality
df_out$carbon <- CCAA.GT.Src$AdjC
df_out$FL <- CCAA.GT.Src$FL
head(df_out)

#plot with colors by locality
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=locality ))
p<-p+geom_point()
p

#add the percent variance explained by PC1 and PC2
percentage <- round(GT.pca.src$sdev / sum(GT.pca.src$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=locality ))
p<-p+geom_point()+ xlab(percentage[1]) + ylab(percentage[2]) + ggtitle("Source AAs")
p

#plot the features that contribute to the classification
df_out_r <- as.data.frame(GT.pca.src$rotation)
df_out_r$feature <- row.names(df_out_r)
df_out_r

p<-ggplot(df_out_r,aes(x=PC1,y=PC2,label=feature,color=feature ))
p<-p+geom_point() + geom_text(size=3)
p

## Now Apply PCA to trophic AAs 
GT.pca.tro = stats :: prcomp(CCAA.GT.Tro[,c(4:6)])
GT.pca.tro

#loadings
GT.pca.tro$rotation

#Create a new dataframe with sample information
df_out.tro <- as.data.frame(GT.pca.tro$x)
df_out.tro$samples <- CCAA.GT.Tro$Row.Labels
df_out.tro$locality <- CCAA.GT.Tro$Locality
df_out.tro$carbon <- CCAA.GT.Tro$AdjC
df_out.tro$FL <- CCAA.GT.Tro$FL
head(df_out.tro)

p.tro<-ggplot(df_out.tro,aes(x=PC1,y=PC2,color=locality ))
p.tro<-p.tro+geom_point()
p.tro

#add the percent variance explained by PC1 and PC2
percentage.tro <- round(GT.pca.tro$sdev / sum(GT.pca.tro$sdev) * 100, 2)
percentage.tro <- paste( colnames(df_out.tro), "(", paste( as.character(percentage.tro), "%", ")", sep="") )

p.tro<-ggplot(df_out.tro,aes(x=PC1,y=PC2,color=locality ))
p.tro<-p.tro+geom_point() + xlab(percentage.tro[1]) + ylab(percentage.tro[2]) + ggtitle("Trophic AAs")
p.tro


#plot the features that contribute to the classification
df_out_r.tro <- as.data.frame(GT.pca.tro$rotation)
df_out_r.tro$feature <- row.names(df_out_r.tro)

df_out_r.tro

p.tro<-ggplot(df_out_r.tro,aes(x=PC1,y=PC2,label=feature,color=feature ))
p.tro<-p.tro+geom_point() + geom_text(size=3)
p.tro


### PCA for all data (source and trophic) ###
head(CCAA.GT)
CCAA.GT.red <- CCAA.GT[-c(16:17,45),]
head(CCAA.GT.red)

WtMns.Red <- WtMns.Red[,-c(6)]
head(WtMns.Red)
#normalize: scale: divide by SD
#install.packages("BBmisc")
library("BBmisc")
dim(WtMns.Red)
#normalize without FL (b/c of NA value)
WtMns.Norm <- normalize(WtMns.Red[-c(40),-c(3)], method = "scale", range = c(0, 1), margin = 1L)
head(WtMns.Norm)

dim(WtMns.Norm)


####  PCA:  normalized weighted means ####
#apply PCA to: AdjC, BulkN, WtMnTrp, WtMnSrc, and TP_WtMn #
GT.pca = stats :: prcomp(WtMns.Norm[,c(31:32,35:36,39)])
GT.pca


#loadings
GT.pca$rotation
#scores (PCs)
GT.pca$x

#basic plot
plot(GT.pca$x[,1], GT.pca$x[,2])

#Create a new dataframe with sample information
df.GT <- as.data.frame(GT.pca$x)
df.GT$samples <- WtMns.Norm$Row.Labels
df.GT$locality <- WtMns.Norm$Locality
head(df.GT)

#plot with colors by locality
p.GT<-ggplot(df.GT,aes(x=PC1,y=PC2,color=locality ))
p.GT<-p.GT+geom_point()
p.GT

#add the percent variance explained by PC1 and PC2
percentage.GT <- round(GT.pca$sdev / sum(GT.pca$sdev) * 100, 2)
percentage.GT <- paste( colnames(df.GT), "(", paste( as.character(percentage.GT), "%", ")", sep="") )

p.GT<-ggplot(df.GT,aes(x=PC1,y=PC2,color=locality ))
p.GT<-p.GT+geom_point()+ xlab(percentage.GT[1]) + ylab(percentage.GT[2]) + ggtitle("GT Normalized Values")
p.GT

#plot the features that contribute to the classification
df.GT.r <- as.data.frame(GT.pca$rotation)
df.GT.r$feature <- row.names(df.GT.r)
df.GT.r

p.GT.f<-ggplot(df.GT.r,aes(x=PC1,y=PC2,label=feature,color=feature )) 
p.GT.f<-p.GT.f+geom_point() + geom_text(size=3) 
p.GT.f


##Plot all in one using autoplot
##But PC1 and PC2 percentages are different - check this!!
row.names(WtMns.Norm) <- WtMns.Norm$Row.Labels
autoplot(prcomp(WtMns.Norm[,c(31:32,35:36,39)]), data = WtMns.Norm, colour = 'Locality', label = TRUE,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = 'black')

prcomp(WtMns.Norm[,c(31:32,35:36,39)])
GT.pca$rotation
GT.pca$x
summary(GT.pca)
#proportion of variance: PC1: 45.58%; PC2: 42.32%; PC3: 6.394%, PC4: 5.715%
#87.89% of variance explained by the first two PCs


#Plot without point labels
autoplot(prcomp(WtMns.Norm[,c(31:32,35:36,39)]), data = WtMns.Norm, colour = 'Locality', label = F,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = 'black')


#add a 40% confidence interval ellipse
p.GT.e <- p.GT+geom_point()+ xlab(percentage.GT[1]) + ylab(percentage.GT[2]) + theme_classic() + ggtitle("GT Normalized Values") + stat_ellipse(geom = "path", position = "identity", type = "t", level = 0.40, show.legend = NA, inherit.aes = TRUE) 

+ stat_ellipse(type = "t", level = 0.95, linetype = 2)


# Determine number of clusters
GT.Red.CCSIA.TP <- GT.Red$CCSIA_TP
GT.Red.CCSIA.TP <- na.exclude(GT.Red.CCSIA.TP)
class(GT.Red.CCSIA.TP)
GT.Red.CCSIA.TP <- as.matrix(GT.Red.CCSIA.TP)
wss <- (nrow(GT.Red.CCSIA.TP)-1)*sum(apply(GT.Red.CCSIA.TP,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(GT.Red$CCSIA_TP, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)


#### Plot Sampling Locality Map ####
#install.packages("dismo")
library(dismo)
library(rgdal)
library(maps)
#install.packages("mapdata")
library(mapdata)

#Map Figure 1
dev.off()
map("world2", xlim = c(0,110), ylim = c(-50,40), fill = T, col = "light gray", lty = 0)
map.scale(0, -45, relwidth = 0.20, metric = TRUE, ratio = F)
map.axes(cex.axis = 0.8)

#x = longitude, y = latitude
dev.off()
map("world2", lwd = 1.5, xlim = c(25,60), ylim = c(-40,0), fill = T, col = "light gray")
SIA.sites <- read.csv("SIA_SampleSites.csv", header = T)
points(SIA.sites$Long, SIA.sites$Lat, pch = 19, cex = 1.5, col = unique(as.factor(SIA.sites$Site)))
map.axes(cex.axis = 0.8)
map.scale(26, -38, relwidth = 0.20, metric = TRUE, ratio = F)


#Map Figure 2
dev.off()
plot.new()
par(mfrow = c(1,2))
map("world2", lwd = 1.5, xlim = c(0,90), ylim = c(-35,25), fill = T, col = "light gray")
points(SIA.sites$Long, SIA.sites$Lat, pch = 19, cex = 1.5, col = unique(as.factor(SIA.sites$Site)))
map.scale(60, -30, relwidth = 0.20, metric = TRUE, ratio = F)
map.axes(cex.axis = 0.8)


#### Propogation of Errors - using PropogationErrors CSV file ####
setwd("~/Documents/Yale_Grad/Stable_Isotope_Analysis/CCSIA/CCSIA_Analyses")
PropErr<- read.csv("PropogationErrors.csv")

#Drop Mtentu and Mauritius because of low sample size
PropErr.Red <- PropErr[-c(16:17,45),]
PropErr.Red$Locality <- droplevels(PropErr.Red$Locality, exclude = c("MAU", "MTENTU"))
PropErr.Red$Locality

head(PropErr.Red)

#Beta = 3.6 sd = 0.5
#TEF = 5.7 sd = 0.3 

#Variance using simplified equation with partial derivatives
((PropErr.Red$WtMnTrp.VAR + PropErr.Red$WtMnSrc.VAR + (0.5^2)) / 5.7^2 ) + ((PropErr.Red$WtMnTrp - PropErr.Red$WtMnSrc - 3.6 )/5.7^2)^2 * 0.3^2

#SD using simplified equation with partial derivatives
sqrt((PropErr.Red$WtMnTrp.VAR + (PropErr.Red$WtMnSrc.VAR) + (0.5^2)) / 5.7^2 + ((PropErr.Red$WtMnTrp - PropErr.Red$WtMnSrc - 3.6 )/5.7^2)^2 * 0.3^2)

