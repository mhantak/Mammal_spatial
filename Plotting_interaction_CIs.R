setwd("/Users/Maggie/Dropbox/Mammal_data/Mammal_spatial_Comm_Bio/GitHub/")

library(tidyverse)
library(gridExtra)
library(phyr)
library(INLA)

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

#HBL.traits.climate.sp$id_cells <- seq.int(nrow(HBL.traits.climate.sp))
str(HBL.traits.climate.sp)

################################################################################################
################################################################################################
# read in phylogeny
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

##############
#HB length 
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

################################################

# set up spatial correlation covariance matrix
## grab coordinates for each unique cell
#cell_id <- dplyr::distinct(HBL.traits.climate.sp, id_cells, .keep_all = T) %>% 
#  dplyr::select(id_cells, decimallongitude, decimallatitude)
#x.coord <- cell_id$decimallongitude
#y.coord <- cell_id$decimallatitude

#nsite <- length(x.coord)

# generate matrix
#Dist <- matrix(0, nrow = nsite, ncol = nsite)
#for (i in 1:nsite)
#  for (j in 1:nsite)
#    Dist[i, j] <-
#  ((x.coord[i] - x.coord[j]) ^ 2 + (y.coord[i] - y.coord[j]) ^ 2) ^ .5

#range <- 14446371
#sd.space <- 30000 # sd of b0_site
#V.space <- sd.space ^ 2 * exp(-Dist / range)
#rownames(V.space) <- 1:nsite
#colnames(V.space) <- 1:nsite
#V.space <- V.space / max(V.space)
#V.space <- V.space / exp(determinant(V.space)$modulus[1] / nsite)
#row.names(V.space) <- cell_id$id_number
#colnames(V.space) <- cell_id$id_number


#Model to use
#hbl.sp.pglmm2 <- pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
#                         hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
#                         MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
#                         MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
#                         (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
#                       data=HBL.traits.climate.sp, cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), bayes = TRUE)


# interaction between temperature (x-axis) and categorical data (hibernation).
plot_pglmm_inter_3a <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_ag <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep("DT"), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ag <- rbind(mdf, new_inter_term_ag)
  
  new_inter_term.fw <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep("HIB"), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.fw <- rbind(mdf, new_inter_term.fw)
  
  new_inter_term.ug <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep("NONE"), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ug <- rbind(mdf, new_inter_term.ug)
  
  # fit models
  m.ag <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.ag, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  m.fw <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.fw, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  m.ug <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.ug, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag <- m.ag$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Torpor",
           inter_eff2_vals = 0)
  a.fw <- m.fw$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Hibernator",
           inter_eff2_vals = -1)
  a.ug <- m.ug$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "None",
           inter_eff2_vals = 1)
  
  a.total <- rbind(a.ag, a.fw, a.ug)
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    #geom_point(HBL.traits.climate.sp, mapping = aes(x = MAT_1yr, y = HBL.log10),
    #           alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Mean Annual Temperature", y = "Head-body Length (log10)", fill = "Hibernation", color = "Hibernation") + 
    scale_color_manual(values = c("purple", "black", "orange")) +
    scale_fill_manual(values = c("purple", "black", "orange")) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="top", legend.justification='left') 
  
  return(p)
  
}

## run function 
fig3a <- plot_pglmm_inter_3a(model_df = HBL.traits.climate.sp, pred_vals = -3:3)
fig3a

#devtools::install_github("tidyverse/tidyverse")
#install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
###########################################################################

#Habitat buffer
str(HBL.traits.climate.sp)
table(HBL.traits.climate.sp$habitat_buffer_three)

plot_pglmm_inter_3b <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_ag2 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep("F"), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ag2 <- rbind(mdf, new_inter_term_ag2)
  
  new_inter_term.fw2 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep("N"), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.fw2 <- rbind(mdf, new_inter_term.fw2)
  
  new_inter_term.ug2 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep("O"), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ug2 <- rbind(mdf, new_inter_term.ug2)
  
  # fit models
  m.ag2 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.ag2, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  m.fw2 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.fw2, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  m.ug2 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.ug2, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag2 <- m.ag2$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Facultative",
           inter_eff2_vals = 0)
  a.fw2 <- m.fw2$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "None",
           inter_eff2_vals = -1)
  a.ug2 <- m.ug2$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Obligate",
           inter_eff2_vals = 1)
  
  a.total2 <- rbind(a.ag2, a.fw2, a.ug2)
  
  p2 <-  ggplot(a.total2, mapping = aes(x = inter_eff1, y = mean)) +
    #geom_point(HBL.traits.climate.sp, mapping = aes(x = MAT_1yr, y = HBL.log10),
    #           alpha = 0.08) + 
    #geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
    #                          fill = inter_eff2), 
    #            alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Mean Annual Temperature", y = "Head-body Length (log10)", fill = "Habitat Buffering", color = "Habitat Buffering") + 
    scale_color_manual(values = c("yellow1", "cyan4", "purple4")) +
    scale_fill_manual(values = c("yellow1", "cyan4", "purple4")) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="top", legend.justification='left') 
  
  return(p2)
  
}

## run function 
fig3b <- plot_pglmm_inter_3b(model_df = HBL.traits.climate.sp, pred_vals = -3:3)
fig3b


###########################################################################

#mean head-body length - temp 
str(HBL.traits.climate.sp)
table(HBL.traits.climate.sp$mean_HBL_binary)

plot_pglmm_inter_3d <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_ag3 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                   habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                   reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                   mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep("Large.HBL"), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                   X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                   decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                   sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                   HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                   month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                   MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                   pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                   current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                   pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ag3 <- rbind(mdf, new_inter_term_ag3)
  
  new_inter_term.fw3 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                   habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                   reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                   mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep("Small.HBL"), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                   X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                   decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                   sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                   HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                   month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                   MAT_1yr= pred_vals, MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                   pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                   current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                   pop.den.log= rep(0, length(pred_vals)), BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.fw3 <- rbind(mdf, new_inter_term.fw3)
  
  # fit models
  m.ag3 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                    hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                    MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                    MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                    (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                  data = pred.df.ag3, 
                  cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                  bayes = TRUE)
  
  m.fw3 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                    hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                    MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                    MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                    (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                  data = pred.df.fw3, 
                  cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                  bayes = TRUE)

  # make dataframe of predicted results 
  a.ag3 <- m.ag3$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Large",
           inter_eff2_vals = 0)
  a.fw3 <- m.fw3$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Small",
           inter_eff2_vals = -1)
  
  a.total3 <- rbind(a.ag3, a.fw3)
  
  p3 <-  ggplot(a.total3, mapping = aes(x = inter_eff1, y = mean)) +
    #geom_point(HBL.traits.climate.sp, mapping = aes(x = MAT_1yr, y = HBL.log10),
    #           alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Mean Annual Temperature", y = "Head-body Length (log10)", fill = "Mean Head-body Length (binned)", color = "Mean Head-body Length (binned)") + 
    scale_color_manual(values = c("yellow1", "cyan4", "purple4")) +
    scale_fill_manual(values = c("yellow1", "cyan4", "purple4")) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="top", legend.justification='left') 
  
  return(p3)
  
}

## run function 
fig3d <- plot_pglmm_inter_3d(model_df = HBL.traits.climate.sp, pred_vals = -3:3)
fig3d

##########################################################################################

#mean head-body length - pop den 
str(HBL.traits.climate.sp)
table(HBL.traits.climate.sp$mean_HBL_binary)

plot_pglmm_inter_3e <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_ag4 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                   habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                   reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                   mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep("Large.HBL"), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                   X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                   decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                   sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                   HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                   month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                   MAT_1yr= rep(0, length(pred_vals)), MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                   pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                   current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                   pop.den.log= pred_vals, BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ag4 <- rbind(mdf, new_inter_term_ag4)
  
  new_inter_term.fw4 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA, length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                   habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep(NA, length(pred_vals)),
                                   reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                   mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep("Small.HBL"), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                   X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                   decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                   sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                   HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                   month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                   MAT_1yr= rep(0, length(pred_vals)), MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                   pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                   ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                   current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                   pop.den.log= pred_vals, BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.fw4 <- rbind(mdf, new_inter_term.fw4)
  
  # fit models
  m.ag4 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                    hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                    MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                    MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                    (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                  data = pred.df.ag4, 
                  cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                  bayes = TRUE)
  
  m.fw4 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                    hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                    MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                    MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                    (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                  data = pred.df.fw4, 
                  cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                  bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag4 <- m.ag4$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Large",
           inter_eff2_vals = 0)
  a.fw4 <- m.fw4$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Small",
           inter_eff2_vals = -1)
  
  a.total4 <- rbind(a.ag4, a.fw4)
  
  p4 <-  ggplot(a.total4, mapping = aes(x = inter_eff1, y = mean)) +
    #geom_point(HBL.traits.climate.sp, mapping = aes(x = MAT_1yr, y = HBL.log10),
    #           alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Human Population Density (log10)", y = "Head-body Length (log10)", fill = "Mean Head-body Length (binned)", color = "Mean Head-body Length (binned)") + 
    scale_color_manual(values = c("yellow1", "cyan4", "purple4")) +
    scale_fill_manual(values = c("yellow1", "cyan4", "purple4")) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="top", legend.justification='left') 
  
  return(p4)
  
}

## run function 
fig3e <- plot_pglmm_inter_3e(model_df = HBL.traits.climate.sp, pred_vals = -3:3)
fig3e

##############################################################################
#activity time - pop den 
str(HBL.traits.climate.sp)
table(HBL.traits.climate.sp$diurnal_nocturnal)

plot_pglmm_inter_3c <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_ag5 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA,length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep("both"),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= rep(0, length(pred_vals)), MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= pred_vals, BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ag5 <- rbind(mdf, new_inter_term_ag5)
  
  new_inter_term.fw5 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA,length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep("diurnal"),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= rep(0, length(pred_vals)), MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= pred_vals, BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.fw5 <- rbind(mdf, new_inter_term.fw5)
  
  new_inter_term.ug5 <- data.frame(binomial2 =  rep(NA, length(pred_vals)), average_body_mass = rep(NA, length(pred_vals)), hibernation = rep(NA,length(pred_vals)), hibernation_binary = rep(NA,length(pred_vals)),
                                  habitat_buffer = rep(NA, length(pred_vals)), buffered_binary = rep(NA, length(pred_vals)), habitat_buffer_three = rep(NA, length(pred_vals)), diurnal_nocturnal = rep("nocturnal"),
                                  reproductive_rate =  rep(NA, length(pred_vals)), reprdocutive_coding =  rep(NA, length(pred_vals)), mating_system =  rep(NA, length(pred_vals)), average_HB_length =  rep(NA, length(pred_vals)),
                                  mean_BM_binary =  rep(NA, length(pred_vals)), mean_HBL_binary = rep(NA, length(pred_vals)), id_number= rep(NA, length(pred_vals)), X = rep(NA,length(pred_vals)),
                                  X1st_body_mass = rep(NA, length(pred_vals)), catalognumber = rep(NA, length(pred_vals)), collectioncode = rep(NA, length(pred_vals)), decimallatitude = rep(NA, length(pred_vals)), 
                                  decimallongitude = rep(NA, length(pred_vals)), dynamicproperties = rep(NA, length(pred_vals)), institutioncode = rep(NA, length(pred_vals)), lifestage= rep(NA, length(pred_vals)), 
                                  sex = rep(NA, length(pred_vals)), X1st_tail_length= rep(NA, length(pred_vals)), X1st_total_length =  rep(NA, length(pred_vals)), HB.length= rep(NA, length(pred_vals)), 
                                  HBL.log= rep(NA, length(pred_vals)), TotL.log= rep(NA, length(pred_vals)), BM.log = rep(NA, length(pred_vals)), year= rep(NA, length(pred_vals)), 
                                  month= rep(NA, length(pred_vals)), day= rep(NA, length(pred_vals)), season= rep(NA, length(pred_vals)), source= rep(NA, length(pred_vals)), 
                                  MAT_1yr= rep(0, length(pred_vals)), MAP_1yr= rep(NA, length(pred_vals)), decade= rep(NA, length(pred_vals)), pop_density_1km2= rep(NA, length(pred_vals)), 
                                  pop_density_4km2= rep(NA, length(pred_vals)), pop_density_10km2= rep(NA, length(pred_vals)), NA_L1NAME= rep(NA, length(pred_vals)), ave_5_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_5_yr_tmean= rep(NA, length(pred_vals)), ave_4_yr_ppt= rep(NA, length(pred_vals)), ave_4_yr_tmean= rep(NA, length(pred_vals)), ave_3_yr_ppt= rep(NA, length(pred_vals)), 
                                  ave_3_yr_tmean= rep(NA, length(pred_vals)), ave_2_yr_ppt= rep(NA, length(pred_vals)), ave_2_yr_tmean= rep(NA, length(pred_vals)), current_yr_ppt= rep(NA, length(pred_vals)), 
                                  current_yr_tmean= rep(NA, length(pred_vals)), MAP.log= rep(0, length(pred_vals)), MAP.log_2yr= rep(NA, length(pred_vals)), MAP.log_5yr= rep(NA, length(pred_vals)), 
                                  pop.den.log= pred_vals, BM.log10= rep(NA, length(pred_vals)), HBL.log10= rep(NA, length(pred_vals)), decade2= rep(NA, length(pred_vals)))
  
  pred.df.ug5 <- rbind(mdf, new_inter_term.ug5)
  
  # fit models
  m.ag5 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.ag5, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  m.fw5 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.fw5, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  m.ug5 <-  pglmm(HBL.log10 ~ MAT_1yr + MAP.log + season + pop.den.log + 
                   hibernation +  habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                   MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                   MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                   (1 | NA_L1NAME) + (1 | decade2) + (1 | binomial2__), 
                 data = pred.df.ug5, 
                 cov_ranef = list(binomial2 = pruned.tree.HBL.spatial), 
                 bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag5 <- m.ag5$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Both",
           inter_eff2_vals = 0)
  a.fw5 <- m.fw5$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Diurnal",
           inter_eff2_vals = -1)
  a.ug5 <- m.ug5$inla.model$summary.linear.predictor[92725:92731,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Nocturnal",
           inter_eff2_vals = 1)
  
  a.total5 <- rbind(a.ag5, a.fw5, a.ug5)
  
  p5 <-  ggplot(a.total5, mapping = aes(x = inter_eff1, y = mean)) +
    #geom_point(HBL.traits.climate.sp, mapping = aes(x = MAT_1yr, y = HBL.log10),
    #           alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Human Population Density (log10)", y = "Head-body Length (log10)", fill = "Activity Time", color = "Activity Time") + 
    scale_color_manual(values = c("purple", "black", "orange")) +
    scale_fill_manual(values = c("purple", "black", "orange")) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="top", legend.justification='left') 
  
  return(p5)
  
}

## run function 
fig3c <- plot_pglmm_inter_3c(model_df = HBL.traits.climate.sp, pred_vals = -3:3)
fig3c


