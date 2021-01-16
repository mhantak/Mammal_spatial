#Mammalian body size is determined by interactions between climate, urbanization, and life history traits
#Body mass and head-body length datasets are post inital data filtering - see 'data filtering' within the methods

setwd("/Users/Maggie/Dropbox/Mammal_data/Mammal_temporal/Github/")

library(plyr); library(dplyr)
library(lme4); library(lmerTest); library(MuMIn)
library(effects); library(ggeffects)
library(ggplot2)
library(tidyr); library(broom); library(tibble); library(purrr)
library(car)
library(MASS)
library(optimx)

##TRAIT DATA
Mammal.traits <- read.csv("Traits_all.csv", header = TRUE, stringsAsFactors = FALSE)
#str(Mammal.traits)

Mammal.traits$binomial2 <- as.factor(Mammal.traits$binomial2)
Mammal.traits$hibernation_binary <- as.factor(Mammal.traits$hibernation_binary)
Mammal.traits$buffered_binary <- as.factor(Mammal.traits$buffered_binary)
Mammal.traits$habitat_buffer <- as.factor(Mammal.traits$habitat_buffer)
Mammal.traits$diurnal_nocturnal <- as.factor(Mammal.traits$diurnal_nocturnal)

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

##########################################################################################

##Spatial data

#BODY MASS
Mammal.BM <- read.csv("Mammal_temporal_BM.csv", header = TRUE, stringsAsFactors = FALSE)
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

#Remove TROPICAL WET FORESTS & WATER - not enough data
Mam.BM.spatial <- Mammal.BM %>% filter(NA_L1NAME != 'TROPICAL WET FORESTS') %>% droplevels()
Mam.BM.spatial <- Mam.BM.spatial %>% filter(NA_L1NAME != 'WATER') %>% droplevels()

table(Mam.BM.spatial$NA_L1NAME)
plyr::count(Mam.BM.spatial$binomial2)

#log MAT and log10 pop density 
MAP.log <- log(Mam.BM.spatial$MAP_1yr) 
Mam.BM.spatial <- cbind(Mam.BM.spatial, MAP.log)
######log 10 transform human population density then scale
pop.den.log <- log10(Mam.BM.spatial$pop_density_10km2 + 1) # +1 because there are zero values 
Mam.BM.spatial <- cbind(Mam.BM.spatial, pop.den.log)

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
Mam.BM.spatial <- transform(Mam.BM.spatial, pop.den.log=scale(pop.den.log), MAP.log=scale(MAP.log), MAT_1yr=scale(MAT_1yr), MAP_1yr=scale(MAP_1yr), pop_density_10km2=scale(pop_density_10km2))
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

############################################################################################################################################
#spatial HB Length 

Mammal.HBL <- read.csv("Mammal_temporal_HBL.csv", header = TRUE, stringsAsFactors = FALSE)
str(Mammal.HBL)

# decade cut off 
rm.dec2 <-(Mammal.HBL$decade > 1939)
Mammal.HBL <- Mammal.HBL[rm.dec2,]

#Remove records that are missing data
Mammal.HBL <- Mammal.HBL %>% drop_na(NA_L1NAME)
Mammal.HBL <- Mammal.HBL %>% drop_na(MAT_1yr)
Mammal.HBL <- Mammal.HBL %>% drop_na(pop_density_10km2)
Mammal.HBL <- Mammal.HBL %>% drop_na(MAP_1yr)

table(Mammal.HBL$NA_L1NAME)

#Remove TROPICAL WET FORESTS & WATER - not enough data 
Mam.HBL.sp.spatial <- Mammal.HBL %>% filter(NA_L1NAME != 'TROPICAL WET FORESTS') %>% droplevels()
Mam.HBL.sp.spatial <- Mam.HBL.sp.spatial %>% filter(NA_L1NAME != 'WATER') %>% droplevels()

#log MAT and log10 pop density 
MAP.log <- log(Mam.HBL.sp.spatial$MAP_1yr) 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, MAP.log)
######log 10 transform human population density then scale
pop.den.log <- log10(Mam.HBL.sp.spatial$pop_density_10km2 + 1) # +1 because there are zero values 
Mam.HBL.sp.spatial <- cbind(Mam.HBL.sp.spatial, pop.den.log)

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
Mam.HBL.sp.spatial <- transform(Mam.HBL.sp.spatial, pop.den.log=scale(pop.den.log), MAP.log=scale(MAP.log), MAT_1yr=scale(MAT_1yr), MAP_1yr=scale(MAP_1yr), pop_density_10km2=scale(pop_density_10km2))
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

################################################################################################
##TRAITS WITH SPATIAL BM AND HB LENGTH DATASETS 

#Combine traits with spatial body mass data 
BM.traits.climate.sp <- merge(Mammal.traits, Mam.BM.spatial)
str(BM.traits.climate.sp)
plyr::count(BM.traits.climate.sp$binomial2)
BM.traits.climate.sp <- BM.traits.climate.sp %>% filter(binomial2 != 'Gulo_gulo') %>% droplevels() #Needed to make sure species are dropped from trait dataset

#Combine traits with spatial HB length data 
HBL.traits.climate.sp <- merge(Mammal.traits, Mam.HBL.sp.spatial)
str(HBL.traits.climate.sp)
plyr::count(HBL.traits.climate.sp$binomial2)
HBL.traits.climate.sp <- HBL.traits.climate.sp %>% filter(binomial2 != 'Gulo_gulo') %>% droplevels() #Needed to make sure species are dropped from trait dataset


################################################################################################
##Modeling 

#Full model 
bm.sp.mod.no.s <- lmer(BM.log10 ~ MAT_1yr + MAP.log + season + sex + pop.den.log + 
                         hibernation_binary + buffered_binary + diurnal_nocturnal + mean_BM_binary +
                         MAT_1yr:mean_BM_binary + pop.den.log:mean_BM_binary +
                         MAT_1yr:pop.den.log + MAT_1yr:hibernation_binary + MAT_1yr:buffered_binary + pop.den.log:diurnal_nocturnal +
                         (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=BM.traits.climate.sp,
                       control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(bm.sp.mod.no.s)

#step-wise selection 
step.mod.bm.sp01.no.s <- step(bm.sp.mod.no.s)
top.bm.sp.mod.no.s <- get_model(step.mod.bm.sp01.no.s)
summary(top.bm.sp.mod.no.s)
#plot(allEffects(top.bm.sp.mod.no.s))
MuMIn::r.squaredGLMM(top.bm.sp.mod.no.s)


######################
#HB LENGTH 

###Full model
hbl.sp.mod.no.s <- lmer(HBL.log10 ~ MAT_1yr + MAP.log + season + sex + pop.den.log + 
                          hibernation_binary + buffered_binary + diurnal_nocturnal + mean_HBL_binary +
                          MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                          MAT_1yr:pop.den.log + MAT_1yr:hibernation_binary + MAT_1yr:buffered_binary + pop.den.log:diurnal_nocturnal +
                          (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=HBL.traits.climate.sp,
                        control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(hbl.sp.mod.no.s)

#step-wise selection 
step.mod.hbl.sp01.no.s <- step(hbl.sp.mod.no.s)
top.hbl.sp.mod1.no.s <- get_model(step.mod.hbl.sp01.no.s)
summary(top.hbl.sp.mod1.no.s)
MuMIn::r.squaredGLMM(top.hbl.sp.mod1.no.s)

################################################################################################################
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
matching.files <- name.check(mam.tree, BM.spatial.species) # all are fine, so re-upload dataset, without header or row.names
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
bm.sp.pglmm.no.s <- pglmm(BM.log ~ MAT_1yr + MAP.log + season + sex + pop.den.log +
                            hibernation_binary + buffered_binary + diurnal_nocturnal + mean_BM_binary +
                            MAT_1yr:mean_BM_binary + 
                            MAT_1yr:pop.den.log + MAT_1yr:hibernation_binary + MAT_1yr:buffered_binary + pop.den.log:diurnal_nocturnal + 
                            (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                          data=BM.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.BM.spatial), bayes = TRUE)



##Head-body length
##spatial - HBL species list 
HBL.spatial.species <- read.csv("Mam.HBL.spatial.species.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#checking match
matching.files2 <- name.check(mam.tree, HBL.spatial.species) # all are fine, so re-upload dataset, without header or row.names
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
hbl.sp.pglmm2.no.s <- pglmm(HBL.log ~ MAT_1yr + MAP.log + season + pop.den.log + 
                              hibernation_binary +  buffered_binary + diurnal_nocturnal + mean_HBL_binary +
                              MAT_1yr:mean_BM_binary + pop.den.log:mean_BM_binary +
                              MAT_1yr:hibernation_binary + MAT_1yr:buffered_binary + pop.den.log:diurnal_nocturnal +
                              (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                            data=HBL.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), bayes = TRUE)







