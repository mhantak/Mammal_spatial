#Mammalian body size is determined by interactions between climate, urbanization, and ecological history traits
#Body mass and head-body length datasets are post initial data filtering - see 'data filtering' within the methods

setwd("/Users/Maggie/Dropbox/Mammal_data/Mammal_spatial_Comm_Bio/GitHub/")

library(plyr); library(dplyr)
library(lme4); library(lmerTest); library(MuMIn)
library(effects); library(ggeffects)
library(ggplot2)
library(tidyr); library(broom); library(tibble); library(purrr)
library(car)
library(MASS)
library(optimx)
library(performance)
library(gridExtra)
library(ggpubr)
library(rr2)

##TRAIT DATA
Mammal.traits <- read.csv("Mammal_traits_extended.csv", header = TRUE, stringsAsFactors = FALSE)
str(Mammal.traits)

Mammal.traits$binomial2 <- as.factor(Mammal.traits$binomial2)
Mammal.traits$hibernation_binary <- as.factor(Mammal.traits$hibernation_binary)
Mammal.traits$buffered_binary <- as.factor(Mammal.traits$buffered_binary)
Mammal.traits$habitat_buffer_three <- as.factor(Mammal.traits$habitat_buffer_three)
Mammal.traits$diurnal_nocturnal <- as.factor(Mammal.traits$diurnal_nocturnal)
Mammal.traits$hibernation <- as.factor(Mammal.traits$hibernation)

#bin sizes -- avg body mass < 500 = small; >500 = large
Large.BM <- Mammal.traits[ which(Mammal.traits$average_body_mass > 500), ] #500
Large.BM["mean_BM_binary"] <- "Large.BM"
Small.BM <- Mammal.traits[ which(Mammal.traits$average_body_mass < 500), ] #500
Small.BM["mean_BM_binary"] <- "Small.BM"
Mammal.traits <- rbind(Large.BM, Small.BM)

#bin sizes -- avg HB length < 200 = small; >200 = large
Large.HBL <- Mammal.traits[ which(Mammal.traits$average_HB_length > 200), ] #200
Large.HBL["mean_HBL_binary"] <- "Large.HBL"
Small.HBL <- Mammal.traits[ which(Mammal.traits$average_HB_length < 200), ] #200
Small.HBL["mean_HBL_binary"] <- "Small.HBL"
Mammal.traits <- rbind(Large.HBL, Small.HBL)

Mammal.traits$mean_BM_binary <- as.factor(Mammal.traits$mean_BM_binary)
Mammal.traits$mean_HBL_binary <- as.factor(Mammal.traits$mean_HBL_binary)
str(Mammal.traits)

##########################################################################################

##Spatial data
Mammal.BM <- read.csv("Mammal_spatial_BM_extended.csv", header = TRUE, stringsAsFactors = FALSE)
str(Mammal.BM)
plyr::count(Mammal.BM$binomial2)

## decade cut off 
rm.dec <-(Mammal.BM$decade > 1939)
Mammal.BM <- Mammal.BM[rm.dec,]

#Remove records that are missing data
Mammal.BM <- Mammal.BM %>% drop_na(NA_L1NAME) #ecoregion
Mammal.BM <- Mammal.BM %>% drop_na(MAT_1yr) #Mean Annual Temp.
Mammal.BM <- Mammal.BM %>% drop_na(pop_density_10km2) #human pop. density
Mammal.BM <- Mammal.BM %>% drop_na(MAP_1yr) #Mean annual precip. 
Mammal.BM <- Mammal.BM %>% drop_na(ave_2_yr_ppt) 
Mammal.BM <- Mammal.BM %>% drop_na(ave_2_yr_tmean) 
Mammal.BM <- Mammal.BM %>% drop_na(ave_5_yr_ppt) 
Mammal.BM <- Mammal.BM %>% drop_na(ave_5_yr_tmean) 

#Remove TROPICAL WET FORESTS & WATER - not enough data
Mam.BM.spatial <- Mammal.BM %>% filter(NA_L1NAME != 'TROPICAL WET FORESTS') %>% droplevels()
Mam.BM.spatial <- Mam.BM.spatial %>% filter(NA_L1NAME != 'WATER') %>% droplevels()

table(Mam.BM.spatial$NA_L1NAME)
plyr::count(Mam.BM.spatial$binomial2)

hist(Mam.BM.spatial$MAP_1yr)
hist(Mam.BM.spatial$MAT_1yr)
hist(Mam.BM.spatial$pop_density_10km2)
hist(Mam.BM.spatial$ave_2_yr_ppt)
hist(Mam.BM.spatial$ave_2_yr_tmean)
hist(Mam.BM.spatial$ave_5_yr_ppt)
hist(Mam.BM.spatial$ave_5_yr_tmean)

#log MAT and log10 pop density 
MAP.log <- log(Mam.BM.spatial$MAP_1yr) 
Mam.BM.spatial <- cbind(Mam.BM.spatial, MAP.log)
MAP.log_2yr <- log(Mam.BM.spatial$ave_2_yr_ppt) 
Mam.BM.spatial <- cbind(Mam.BM.spatial, MAP.log_2yr)
MAP.log_5yr <- log(Mam.BM.spatial$ave_5_yr_ppt) 
Mam.BM.spatial <- cbind(Mam.BM.spatial, MAP.log_5yr)

######log 10 transform human population density then scale
pop.den.log <- log10(Mam.BM.spatial$pop_density_10km2 + 1) # +1 because there are zero values 
Mam.BM.spatial <- cbind(Mam.BM.spatial, pop.den.log)

hist(Mam.BM.spatial$MAP.log)
hist(Mam.BM.spatial$MAP.log_2yr)
hist(Mam.BM.spatial$MAP.log_5yr)
hist(Mam.BM.spatial$pop.den.log)

###log10 trasform HBL length and body mass
BM.log10 <- log10(Mam.BM.spatial$X1st_body_mass) 
Mam.BM.spatial <- cbind(Mam.BM.spatial, BM.log10)
##
HBL.log10 <- log10(Mam.BM.spatial$HB.length) 
Mam.BM.spatial <- cbind(Mam.BM.spatial, HBL.log10)

####
Mam.BM.spatial$binomial2 <- as.factor(Mam.BM.spatial$binomial2)
Mam.BM.spatial$season <- as.factor(Mam.BM.spatial$season)
Mam.BM.spatial$source <- as.factor(Mam.BM.spatial$source)
Mam.BM.spatial$NA_L1NAME <- as.factor(Mam.BM.spatial$NA_L1NAME)
Mam.BM.spatial$sex <- as.factor(Mam.BM.spatial$sex)

#scaling and centering
Mam.BM.spatial <- transform(Mam.BM.spatial, pop.den.log=scale(pop.den.log), MAP.log=scale(MAP.log), MAT_1yr=scale(MAT_1yr), MAP_1yr=scale(MAP_1yr), pop_density_10km2=scale(pop_density_10km2), 
                            ave_2_yr_tmean=scale(ave_2_yr_tmean), MAP.log_2yr=scale(MAP.log_2yr), ave_5_yr_tmean=scale(ave_5_yr_tmean), MAP.log_5yr=scale(MAP.log_5yr))
str(Mam.BM.spatial)

Mam.BM.spatial$decade2 <- Mam.BM.spatial$decade
str(Mam.BM.spatial)
table(Mam.BM.spatial$decade2)

#change date to numbers for model interpretation
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "1940"] <- "0"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "1950"] <- "1"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "1960"] <- "2"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "1970"] <- "3"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "1980"] <- "4"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "1990"] <- "5"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "2000"] <- "6"
Mam.BM.spatial$decade2[Mam.BM.spatial$decade2 == "2010"] <- "7"

str(Mam.BM.spatial)
Mam.BM.spatial$decade2 <- as.numeric(Mam.BM.spatial$decade2)
Mam.BM.spatial$decade <- as.numeric(Mam.BM.spatial$decade)
plyr::count(Mam.BM.spatial$binomial2)

############################################################################################################################################
#spatial HB Length 

Mammal.HBL <- read.csv("Mammal_spatial_HBL_extended.csv", header = TRUE, stringsAsFactors = FALSE)
str(Mammal.HBL)

# decade cut off 
rm.dec2 <-(Mammal.HBL$decade > 1939)
Mammal.HBL <- Mammal.HBL[rm.dec2,]

#Remove records that are missing data
Mammal.HBL <- Mammal.HBL %>% drop_na(NA_L1NAME)
Mammal.HBL <- Mammal.HBL %>% drop_na(MAT_1yr)
Mammal.HBL <- Mammal.HBL %>% drop_na(pop_density_10km2)
Mammal.HBL <- Mammal.HBL %>% drop_na(MAP_1yr)
Mammal.HBL <- Mammal.HBL %>% drop_na(ave_2_yr_ppt) 
Mammal.HBL <- Mammal.HBL %>% drop_na(ave_2_yr_tmean) 
Mammal.HBL <- Mammal.HBL %>% drop_na(ave_5_yr_ppt) 
Mammal.HBL <- Mammal.HBL %>% drop_na(ave_5_yr_tmean) 

table(Mammal.HBL$NA_L1NAME)

#Remove TROPICAL WET FORESTS & WATER - not enough data 
Mam.HBL.sp.spatial <- Mammal.HBL %>% filter(NA_L1NAME != 'TROPICAL WET FORESTS') %>% droplevels()
Mam.HBL.sp.spatial <- Mam.HBL.sp.spatial %>% filter(NA_L1NAME != 'WATER') %>% droplevels()

hist(Mam.HBL.sp.spatial$MAP_1yr)
hist(Mam.HBL.sp.spatial$MAT_1yr)
hist(Mam.HBL.sp.spatial$ave_2_yr_ppt)
hist(Mam.HBL.sp.spatial$ave_2_yr_tmean)
hist(Mam.HBL.sp.spatial$pop_density_10km2)
hist(Mam.HBL.sp.spatial$ave_5_yr_ppt)
hist(Mam.HBL.sp.spatial$ave_5_yr_tmean)

#log MAT and log10 pop density 
MAP.log <- log(Mam.HBL.sp.spatial$MAP_1yr) 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, MAP.log)
MAP.log_2yr <- log(Mam.HBL.sp.spatial$ave_2_yr_ppt) 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, MAP.log_2yr)
MAP.log_5yr <- log(Mam.HBL.sp.spatial$ave_5_yr_ppt) 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, MAP.log_5yr)
######log 10 transform human population density then scale
pop.den.log <- log10(Mam.HBL.sp.spatial$pop_density_10km2 + 1) # +1 because there are zero values 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, pop.den.log)

hist(Mam.HBL.sp.spatial$MAP.log)
hist(Mam.HBL.sp.spatial$MAP.log_2yr)
hist(Mam.HBL.sp.spatial$MAP.log_5yr)
hist(Mam.HBL.sp.spatial$pop.den.log)

###log10 trasform HBL length and body mass
BM.log10 <- log10(Mam.HBL.sp.spatial$X1st_body_mass) 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, BM.log10)
##
HBL.log10 <- log10(Mam.HBL.sp.spatial$HB.length) 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, HBL.log10)

####
Mam.HBL.sp.spatial$binomial2 <- as.factor(Mam.HBL.sp.spatial$binomial2)
Mam.HBL.sp.spatial$season <- as.factor(Mam.HBL.sp.spatial$season)
Mam.HBL.sp.spatial$source <- as.factor(Mam.HBL.sp.spatial$source)
Mam.HBL.sp.spatial$NA_L1NAME <- as.factor(Mam.HBL.sp.spatial$NA_L1NAME)
Mam.HBL.sp.spatial$sex <- as.factor(Mam.HBL.sp.spatial$sex)

#scaling and centering
Mam.HBL.sp.spatial <- transform(Mam.HBL.sp.spatial, pop.den.log=scale(pop.den.log), MAP.log=scale(MAP.log), MAT_1yr=scale(MAT_1yr), MAP_1yr=scale(MAP_1yr), pop_density_10km2=scale(pop_density_10km2), 
                                ave_2_yr_tmean=scale(ave_2_yr_tmean), MAP.log_2yr=scale(MAP.log_2yr), ave_5_yr_tmean=scale(ave_5_yr_tmean), MAP.log_5yr=scale(MAP.log_5yr))
str(Mam.HBL.sp.spatial)

Mam.HBL.sp.spatial$decade2 <- Mam.HBL.sp.spatial$decade
str(Mam.HBL.sp.spatial)
table(Mam.HBL.sp.spatial$decade2)

#change date to numbers for model interpretation
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "1940"] <- "0"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "1950"] <- "1"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "1960"] <- "2"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "1970"] <- "3"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "1980"] <- "4"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "1990"] <- "5"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "2000"] <- "6"
Mam.HBL.sp.spatial$decade2[Mam.HBL.sp.spatial$decade2 == "2010"] <- "7"

str(Mam.HBL.sp.spatial)
Mam.HBL.sp.spatial$decade2 <- as.numeric(Mam.HBL.sp.spatial$decade2)
Mam.HBL.sp.spatial$decade <- as.numeric(Mam.HBL.sp.spatial$decade)
plyr::count(Mam.HBL.sp.spatial$binomial2)

################################################################################################
################################################################################################
##TRAITS WITH SPATIAL BM AND HB LENGTH DATASETS 

#Combine traits with spatial body mass data 
BM.traits.climate.sp <- merge(Mammal.traits, Mam.BM.spatial)
str(BM.traits.climate.sp)
plyr::count(BM.traits.climate.sp$binomial2)

#Combine traits with spatial HB length data 
HBL.traits.climate.sp <- merge(Mammal.traits, Mam.HBL.sp.spatial)
str(HBL.traits.climate.sp)
plyr::count(HBL.traits.climate.sp$binomial2)

################################################################################################
################################################################################################

##phylogenetic generalized linear mixed models 

library(ape)
library(geiger)
library(phytools)

####Mammal tree
mam.tree <- read.tree("mammal2.tre")
str(mam.tree)

##body mass
BM.spatial.species <- read.csv("Mam.BM.spatial.species.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

#checking match
matching.files <- name.check(mam.tree, BM.spatial.species) 
matching.files$data_not_tree

#re-upload
BM.spatial.species <- read.csv("Mam.BM.spatial.species.csv")
plyr::count(BM.spatial.species$binomial2)
plyr::count(BM.traits.climate.sp$binomial2)

pruned.tree.BM.spatial <-drop.tip(mam.tree, mam.tree$tip.label[-na.omit(match(BM.spatial.species[,1],mam.tree$tip.label))])
plyr::count(pruned.tree.BM.spatial$tip.label)
str(pruned.tree.BM.spatial)
plot(pruned.tree.BM.spatial)
plotTree(pruned.tree.BM.spatial,ftype="i",lwd=1,fsize=0.3,type="fan",part=1)

##PGLMM
library(phyr)

######top spatial BM model 
bm.sp.pglmm <- pglmm(BM.log10 ~ MAT_1yr + MAP.log + season + sex + pop.den.log +
                       hibernation + habitat_buffer_three + diurnal_nocturnal + mean_BM_binary +
                       MAT_1yr:mean_BM_binary + 
                       MAT_1yr:pop.den.log + MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal + 
                       (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                     data=BM.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.BM.spatial), bayes = TRUE)

library(rr2)
summary(bm.sp.pglmm)
R2(bm.sp.pglmm)
R2_pred(bm.sp.pglmm)

####################
##plotting interactions 

inla_bm <-  bm.sp.pglmm$inla.model
bm_fe <- inla_bm$summary.fixed %>% 
  dplyr::add_rownames(var = "vars")

bm_int <- filter(bm_fe, vars == "(Intercept)")$mean
bm_MAT <- filter(bm_fe, vars == "MAT_1yr")$mean
bm_popden <- filter(bm_fe, vars == "pop.den.log")$mean
bm_hib.hib <- filter(bm_fe, vars == "hibernationHIB")$mean
bm_hib.none <- filter(bm_fe, vars == "hibernationNONE")$mean
bm_buffN <- filter(bm_fe, vars == "habitat_buffer_threeN")$mean
bm_buffO <- filter(bm_fe, vars == "habitat_buffer_threeO")$mean
bm_actD <- filter(bm_fe, vars == "diurnal_nocturnaldiurnal")$mean
bm_actN <- filter(bm_fe, vars == "diurnal_nocturnalnocturnal")$mean
bm_meanBM <- filter(bm_fe, vars == "mean_BM_binarySmall.BM")$mean
bm_MAT_meanBM <- filter(bm_fe, vars == "MAT_1yr:mean_BM_binarySmall.BM")$mean
bm_popden_meanBM <- filter(bm_fe, vars == "pop.den.log:mean_BM_binarySmall.BM")$mean
bm_MAT_popden <- filter(bm_fe, vars == "MAT_1yr:pop.den.log")$mean
bm_MAT_hib <- filter(bm_fe, vars == "MAT_1yr:hibernationHIB")$mean
bm_MAT_none <- filter(bm_fe, vars == "MAT_1yr:hibernationNONE")$mean
bm_MAT_buffN <- filter(bm_fe, vars == "MAT_1yr:habitat_buffer_threeN")$mean
bm_MAT_buffO <- filter(bm_fe, vars == "MAT_1yr:habitat_buffer_threeO")$mean
bm_popden_actD <- filter(bm_fe, vars == "pop.den.log:diurnal_nocturnaldiurnal")$mean
bm_popden_actN <- filter(bm_fe, vars == "pop.den.log:diurnal_nocturnalnocturnal")$mean

##################
## MAT x Mean BM
MAT.x <- BM.traits.climate.sp$MAT_1yr

popden_high <- 1
popden_mid <- 0
popden_low <- -1

bm.low.popden = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_popden * popden_low) + 
  (bm_MAT_popden * popden_low * MAT.x)

bm.mid.popden = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_popden * popden_mid) + 
  (bm_MAT_popden * popden_mid* MAT.x)

bm.high.popden = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_popden * popden_high) + 
  (bm_MAT_popden * popden_high * MAT.x)

mat.bm1 <- c(MAT.x, MAT.x, MAT.x)
bm_pop <- c(bm.high.popden, bm.mid.popden, bm.low.popden)
popz <- c(rep('high', 140499), rep('mid', 140499), rep('low', 140499))

pop_mat <- data.frame(MAT.p = mat.bm1, BM.p = bm_pop, Human.pop.den = as.factor(popz)) 

tpop <- ggplot() + 
  geom_line(pop_mat, mapping = aes(x = mat.bm1, y = BM.p, color = Human.pop.den), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Body mass (log10)") + 
  scale_color_manual(labels = c("High", "Low", "Mid"), values = c("turquoise4", "grey5", "violetred2")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Human Population Density (log10)") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

tpop

########################################

## MAT x Hibernation 
MAT.x <- BM.traits.climate.sp$MAT_1yr

bm.mat.hibN = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_hib.none) + 
  (bm_MAT_none * MAT.x)

bm.mat.hibHIB = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_hib.hib) + 
  (bm_MAT_hib * MAT.x)

bm.mat.hibtor = bm_int + 
  (bm_MAT * MAT.x)

MATz1 <- c(MAT.x, MAT.x, MAT.x)
BM_MAT_hibernation <- c(bm.mat.hibN, bm.mat.hibHIB, bm.mat.hibtor)
hib1.1 <- c(rep("None", 140499), rep("Hibernator", 140499), rep("Torpor", 140499))

hib_MAT1 <- data.frame(MAT.h = MATz1, BM.h = BM_MAT_hibernation, 
                       Hibernation = as.factor(hib1.1)) 

h.m.bm <- ggplot() + 
  geom_line(hib_MAT1, mapping = aes(x = MAT.h, y = BM.h, color = Hibernation), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Body Mass (log10)") + 
  scale_color_manual(values = c("purple", "black", "orange")) +
  labs(color = "Hibernation") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Hibernation") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

h.m.bm
##############################

#MAT x buffer

bm.mat.buffN = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_buffN) + 
  (bm_MAT_buffN * MAT.x)

bm.mat.buffO = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_buffO) + 
  (bm_MAT_buffO * MAT.x)

bm.mat.buffF = bm_int + 
  (bm_MAT * MAT.x)

mat1 <- c(MAT.x, MAT.x, MAT.x)
bm_buffered <- c(bm.mat.buffN, bm.mat.buffO, bm.mat.buffF)
buffered_categories1 <- c(rep("None", 140499), rep("Obligate", 140499), 
                          rep("Facultative", 140499))

bm_mat_buffered <- data.frame(MAT.b = mat1, BM.b = bm_buffered, Buffered = as.factor(buffered_categories1)) 

bm.buff.mat <- ggplot() + 
  geom_line(bm_mat_buffered, mapping = aes(x = MAT.b , y = BM.b, color = Buffered), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Body Mass (log10)") + 
  scale_colour_viridis_d() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Habitat Buffering") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

bm.buff.mat

#######################
## MAT x mean body size  
MAT.x <- BM.traits.climate.sp$MAT_1yr

bm.MAT.bm.s = bm_int + 
  (bm_MAT * MAT.x) + 
  (bm_meanBM) + 
  (bm_MAT_meanBM * MAT.x)

bm.MAT.bm.l = bm_int + 
  (bm_MAT * MAT.x) 

MATz1 <- c(MAT.x, MAT.x)
BM_MAT_mean.bm <- c(bm.MAT.bm.s, bm.MAT.bm.l)
mean.bm1.1 <- c(rep("Small", 140499), rep("Large", 140499))

mean_bm_MAT1 <- data.frame(MAT.m = MATz1, BM.mean = BM_MAT_mean.bm, 
                           Mean.BM = as.factor(mean.bm1.1)) 

bm.mean.bm <- ggplot() + 
  geom_line(mean_bm_MAT1, mapping = aes(x = MAT.m, y = BM.mean, color = Mean.BM), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Body Mass (log10)") + 
  scale_color_manual(values = c("turquoise3", "lightsalmon")) +
  labs(color = "Mean binned BM") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Mean Body Mass (binned)") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

bm.mean.bm

##############################
#pop den x diurnal/nocturnal

pop.x <- BM.traits.climate.sp$pop.den.log

bm.popden.D = bm_int + 
  (bm_popden * pop.x) + 
  (bm_actD) + 
  (bm_popden_actD * pop.x)

bm.popden.N = bm_int + 
  (bm_popden * pop.x) + 
  (bm_actN) + 
  (bm_popden_actN * pop.x)

bm.popden.B = bm_int + 
  (bm_popden * pop.x)

popden1 <- c(pop.x, pop.x, pop.x)
bm_activity <- c(bm.popden.D, bm.popden.N, bm.popden.B)
activity_categories1 <- c(rep("Diurnal", 140499), rep("Nocturnal", 140499), 
                          rep("Both", 140499))

bm_opden_activity <- data.frame(Popden = popden1, BMAct = bm_activity, Activity = as.factor(activity_categories1)) 

bm.act.popden <- ggplot() + 
  geom_line(bm_opden_activity, mapping = aes(x = Popden , y = BMAct, color = Activity), size = 1.05) +
  labs(x = "Human Population Density (log10)", y = "Body Mass (log10)") + 
  scale_color_manual(values = c("gray62", "steelblue2", "gold1")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Activity Time") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

bm.act.popden

#########################################################################################
#body mass plots
#grid.arrange(tpop, h.m.bm, bm.buff.mat, bm.mean.bm, bm.act.popden, ncol = 2)
ggarrange(tpop, h.m.bm, bm.buff.mat, bm.act.popden, bm.mean.bm, labels = c("A", "B", "C", 'D', "E"), ncol = 2, nrow = 3)

############################################################################################################################
############################################################################################################################

##Head-body length
HBL.spatial.species <- read.csv("Mam.HBL.spatial.species.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#checking match
matching.files2 <- name.check(mam.tree, HBL.spatial.species) 
matching.files2$data_not_tree

#re-upload
HBL.spatial.species <- read.csv("Mam.HBL.spatial.species.csv")
plyr::count(HBL.spatial.species$binomial2)
plyr::count(HBL.traits.climate.sp$binomial2)

pruned.tree.HBL.spatial <-drop.tip(mam.tree, mam.tree$tip.label[-na.omit(match(HBL.spatial.species[,1],mam.tree$tip.label))])
str(pruned.tree.HBL.spatial)
plyr::count(pruned.tree.HBL.spatial$tip.label)
plot(pruned.tree.HBL.spatial)
plotTree(pruned.tree.HBL.spatial,ftype="i",lwd=1,fsize=0.3,type="fan",part=1)

str(HBL.traits.climate.sp)

#####################

######top spatial HBL model 
hbl.sp.pglmm2 <- pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                         hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                         MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                         MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                         (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                       data=HBL.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), bayes = TRUE)

R2(hbl.sp.pglmm2)

#plot interactions 
inla_hbl <-  hbl.sp.pglmm2$inla.model
hbl_fe <- inla_hbl$summary.fixed %>% 
  dplyr::add_rownames(var = "vars")

hbl_int <- filter(hbl_fe, vars == "(Intercept)")$mean
hbl_MAT <- filter(hbl_fe, vars == "MAT_1yr")$mean
hbl_popden <- filter(hbl_fe, vars == "pop.den.log")$mean
hbl_hib.hib <- filter(hbl_fe, vars == "hibernationHIB")$mean
hbl_hib.none <- filter(hbl_fe, vars == "hibernationNONE")$mean
hbl_buffN <- filter(hbl_fe, vars == "habitat_buffer_threeN")$mean
hbl_buffO <- filter(hbl_fe, vars == "habitat_buffer_threeO")$mean
hbl_actD <- filter(hbl_fe, vars == "diurnal_nocturnaldiurnal")$mean
hbl_actN <- filter(hbl_fe, vars == "diurnal_nocturnalnocturnal")$mean
hbl_meanHBL <- filter(hbl_fe, vars == "mean_HBL_binarySmall.HBL")$mean
hbl_MAT_meanHBL <- filter(hbl_fe, vars == "MAT_1yr:mean_HBL_binarySmall.HBL")$mean
hbl_popden_meanHBL <- filter(hbl_fe, vars == "pop.den.log:mean_HBL_binarySmall.HBL")$mean
hbl_MAT_hib <- filter(hbl_fe, vars == "MAT_1yr:hibernationHIB")$mean
hbl_MAT_none <- filter(hbl_fe, vars == "MAT_1yr:hibernationNONE")$mean
hbl_MAT_buffN <- filter(hbl_fe, vars == "MAT_1yr:habitat_buffer_threeN")$mean
hbl_MAT_buffO <- filter(hbl_fe, vars == "MAT_1yr:habitat_buffer_threeO")$mean
hbl_popden_actD <- filter(hbl_fe, vars == "pop.den.log:diurnal_nocturnaldiurnal")$mean
hbl_popden_actN <- filter(hbl_fe, vars == "pop.den.log:diurnal_nocturnalnocturnal")$mean

#################
## MAT x Hibernation 
MAT.x2 <- HBL.traits.climate.sp$MAT_1yr

hbl.mat.hibN = hbl_int + 
  (hbl_MAT * MAT.x2) + 
  (hbl_hib.none) + 
  (hbl_MAT_none * MAT.x2)

hbl.mat.hibHIB = hbl_int + 
  (hbl_MAT * MAT.x2) + 
  (hbl_hib.hib) + 
  (hbl_MAT_hib * MAT.x2)

hbl.mat.hibtor = hbl_int + 
  (hbl_MAT * MAT.x2)

MATz <- c(MAT.x2, MAT.x2, MAT.x2)
HBL_MAT_hibernation <- c(hbl.mat.hibN, hbl.mat.hibHIB, hbl.mat.hibtor)
hib1 <- c(rep("None", 92724), rep("Hibernator", 92724), rep("Torpor", 92724))

hib_MAT <- data.frame(MAT.h2 = MATz, HBL.h2 = HBL_MAT_hibernation, 
                      Hibernation2 = as.factor(hib1)) 

h.m.hbl <- ggplot() + 
  geom_line(hib_MAT, mapping = aes(x = MAT.h2, y = HBL.h2, color = Hibernation2), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Head-body Length (log10)") + 
  scale_color_manual(values = c("purple", "black", "orange")) +
  labs(color = "Hibernation") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Hibernation") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

h.m.hbl

#########################
#MAT x buffer

hbl.mat.buffN = hbl_int + 
  (hbl_MAT * MAT.x2) + 
  (hbl_buffN) + 
  (hbl_MAT_buffN * MAT.x2)

hbl.mat.buffO = hbl_int + 
  (hbl_MAT * MAT.x2) + 
  (hbl_buffO) + 
  (hbl_MAT_buffO * MAT.x2)

hbl.mat.buffF = hbl_int + 
  (hbl_MAT * MAT.x2)

mat2 <- c(MAT.x2, MAT.x2, MAT.x2)
hbl_buffered <- c(hbl.mat.buffN, hbl.mat.buffO, hbl.mat.buffF)
buffered_categories <- c(rep("None", 92724), rep("Obligate", 92724), 
                         rep("Facultative", 92724))

hbl_mat_buffered <- data.frame(MAT.b2 = mat2, HBL.b2 = hbl_buffered, Buffered2 = as.factor(buffered_categories)) 

hbl.buff.mat <- ggplot() + 
  geom_line(hbl_mat_buffered, mapping = aes(x = MAT.b2 , y = HBL.b2, color = Buffered2), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Head-body Length (log10)") + 
  scale_color_viridis_d() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Habitat Buffering") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

hbl.buff.mat

###################
## MAT x mean HBL  
MAT.x2 <- HBL.traits.climate.sp$MAT_1yr

hbl.MAT.bm.s = hbl_int + 
  (hbl_MAT * MAT.x2) + 
  (hbl_meanHBL) + 
  (hbl_MAT_meanHBL * MAT.x2)

hbl.MAT.bm.l = hbl_int + 
  (hbl_MAT * MAT.x2) 

MAT.hbl <- c(MAT.x2, MAT.x2)
HBL_MAT_mean.hbl <- c(hbl.MAT.bm.s, hbl.MAT.bm.l)
mean.hbl1.1 <- c(rep("Small", 92724), rep("Large", 92724))

mean_hbl_MAT1 <- data.frame(MAT.m2 = MAT.hbl, HBL.mean = HBL_MAT_mean.hbl, 
                            Mean.HBL = as.factor(mean.hbl1.1)) 

hbl.mean.hbl <- ggplot() + 
  geom_line(mean_hbl_MAT1, mapping = aes(x = MAT.m2, y = HBL.mean, color = Mean.HBL), size = 1.05) +
  labs(x = "Mean Annual Temperature", y = "Head-body Length (log10)") + 
  scale_color_manual(values = c("turquoise3", "lightsalmon")) +
  labs(color = "Mean binned HBL") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Mean Head-body Length (binned)") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

hbl.mean.hbl

#########################################
#pop den x diurnal/nocturnal

pop.x2 <- HBL.traits.climate.sp$pop.den.log

hbl.popden.D = hbl_int + 
  (hbl_popden * pop.x2) + 
  (hbl_actD) + 
  (hbl_popden_actD * pop.x2)

hbl.popden.N = hbl_int + 
  (hbl_popden * pop.x2) + 
  (hbl_actN) + 
  (hbl_popden_actN * pop.x2)

hbl.popden.B = hbl_int + 
  (hbl_popden * pop.x2)

popden2 <- c(pop.x2, pop.x2, pop.x2)
hbl_activity <- c(hbl.popden.D, hbl.popden.N, hbl.popden.B)
activity_categories2 <- c(rep("Diurnal", 92724), rep("Nocturnal", 92724), 
                          rep("Both", 92724))

hbl_opden_activity <- data.frame(Popden2 = popden2, HBLAct = hbl_activity, Activity2 = as.factor(activity_categories2)) 

hbl.act.popden <- ggplot() + 
  geom_line(hbl_opden_activity, mapping = aes(x = Popden2 , y = HBLAct, color = Activity2), size = 1.05) +
  labs(x = "Human Population Density (log10)", y = "Head-body Length (log10)") + 
  scale_color_manual(values = c("gray62", "steelblue2", "gold1")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Activity Time") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

hbl.act.popden

################################
## pop den x mean HBL  
pop.x2 <- HBL.traits.climate.sp$pop.den.log

hbl.pop.hbl.s = hbl_int + 
  (hbl_popden * pop.x2) + 
  (hbl_meanHBL) + 
  (hbl_popden_meanHBL * pop.x2)

hbl.pop.hbl.l = hbl_int + 
  (hbl_popden * pop.x2) 

pop.hbl <- c(pop.x2, pop.x2)
HBL_pop_mean.hbl <- c(hbl.pop.hbl.s, hbl.pop.hbl.l)
mean.hbl1.2 <- c(rep("Small", 92724), rep("Large", 92724))

mean_hbl_pop1 <- data.frame(POP.m2 = pop.hbl, HBL.mean2 = HBL_pop_mean.hbl, 
                            Mean.HBL2 = as.factor(mean.hbl1.2)) 

hbl.mean.pop <- ggplot() + 
  geom_line(mean_hbl_pop1, mapping = aes(x = POP.m2, y = HBL.mean2, color = Mean.HBL2), size = 1.05) +
  labs(x = "Human Population Density (log10)", y = "Head-body Length (log10)") + 
  scale_color_manual(values = c("red", "gray")) +
  labs(color = "Mean binned HBL") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5)) +
  theme(legend.position="top", legend.justification='left') + 
  theme(legend.title = element_text(size = 12), legend.text=element_text(size=12)) +
  labs(col="Mean Head-body Length (binned)") +
  theme(axis.title.y = element_text(size = rel(1.2))) + theme(axis.title.x = element_text(size = rel(1.2)))

hbl.mean.pop

##################################################
#grid.arrange(h.m.hbl, hbl.buff.mat, hbl.mean.hbl, hbl.act.popden, hbl.mean.pop, ncol = 2)
ggarrange(h.m.hbl, hbl.buff.mat, hbl.act.popden, hbl.mean.hbl, hbl.mean.pop, labels = c("A", "B", "C", 'D', "E"), ncol = 2, nrow = 3)
########################################################################################################################################
########################################################################################################################################

#2 year 

#Body mass

bm.sp.pglmm.2yr <- pglmm(BM.log10 ~ ave_2_yr_tmean + MAP.log_2yr + season + sex + pop.den.log +
                           hibernation + habitat_buffer_three + diurnal_nocturnal + mean_BM_binary +
                           ave_2_yr_tmean:mean_BM_binary + 
                           ave_2_yr_tmean:pop.den.log + ave_2_yr_tmean:hibernation + ave_2_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal + 
                           (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                         data=BM.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.BM.spatial), bayes = TRUE)

#HB length 

hbl.sp.pglmm2.2yr <- pglmm(HBL.log10 ~ ave_2_yr_tmean + MAP.log_2yr + season + pop.den.log + 
                             hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                             ave_2_yr_tmean:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                             ave_2_yr_tmean:hibernation + ave_2_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                             (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                           data=HBL.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), bayes = TRUE)


####################################################################################

#5 year 

#Body mass

bm.sp.pglmm.5yr <- pglmm(BM.log10 ~ ave_5_yr_tmean + MAP.log_5yr + season + sex + pop.den.log +
                           hibernation + habitat_buffer_three + diurnal_nocturnal + mean_BM_binary +
                           ave_5_yr_tmean:mean_BM_binary + 
                           ave_5_yr_tmean:pop.den.log + ave_5_yr_tmean:hibernation + ave_5_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal + 
                           (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                         data=BM.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.BM.spatial), bayes = TRUE)


#HB length 

hbl.sp.pglmm2.5yr <- pglmm(HBL.log10 ~ ave_5_yr_tmean + MAP.log_5yr + season + pop.den.log + 
                             hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                             ave_5_yr_tmean:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                             ave_5_yr_tmean:hibernation + ave_5_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                             (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                           data=HBL.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), bayes = TRUE)

