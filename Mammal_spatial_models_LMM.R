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

#library(PerformanceAnalytics)
#test_vars1 <- c("MAT_1yr", "MAP.log", "pop.den.log")
#corr.bm1 <- Mam.BM.spatial[test_vars1]
#chart.Correlation(corr.bm1, histogram=TRUE, pch=19) 

#library(corrplot)
#library(RColorBrewer)
#corrplot(cor(corr.bm1),
#         method = "number",
#         type = "upper", tl.col="black")

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

#test_vars1 <- c("MAT_1yr", "MAP.log", "pop.den.log")
#corr.hbl1 <- Mam.HBL.sp.spatial[test_vars1]
#chart.Correlation(corr.hbl1, histogram=TRUE, pch=19) 

#corrplot(cor(corr.hbl1),
#         method = "number",
#         type = "upper", tl.col="black")


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
##Modeling 

#library(PerformanceAnalytics)
#test_vars <- c("MAT_1yr", "MAP.log", "pop.den.log")
#corr.bm <- BM.traits.climate.sp[test_vars]
#chart.Correlation(corr.bm, histogram=TRUE, pch=19) 

##############
#1 year 

#Body mass
str(BM.traits.climate.sp)
#body mass climate 1 year
bm.sp.mod.no.s2 <- lmer(BM.log10 ~ MAT_1yr + MAP.log + season + sex + pop.den.log + 
                          hibernation + habitat_buffer_three + diurnal_nocturnal + mean_BM_binary +
                          MAT_1yr:mean_BM_binary + pop.den.log:mean_BM_binary +
                          MAT_1yr:pop.den.log + MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                          (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=BM.traits.climate.sp,
                        control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(bm.sp.mod.no.s2)

#step-wise selection 
step.mod.bm.sp01.no.s2 <- step(bm.sp.mod.no.s2)
top.bm.sp.mod.no.s2 <- get_model(step.mod.bm.sp01.no.s2)
summary(top.bm.sp.mod.no.s2)

MuMIn::r.squaredGLMM(top.bm.sp.mod.no.s2)
library(car)
vif(top.bm.sp.mod.no.s2)
library(performance)
check_collinearity(top.bm.sp.mod.no.s2)

###################################################
#head body length 
#test_vars <- c("MAT_1yr", "MAP.log", "pop.den.log")
#corr.hbl <- HBL.traits.climate.sp[test_vars]
#chart.Correlation(corr.hbl, histogram=TRUE, pch=19) 

#1 year
hbl.sp.mod.no.s2 <- lmer(HBL.log10 ~ MAT_1yr + MAP.log + season + sex + pop.den.log + 
                           hibernation + habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                           MAT_1yr:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                           MAT_1yr:pop.den.log + MAT_1yr:hibernation + MAT_1yr:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                           (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=HBL.traits.climate.sp,
                         control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(hbl.sp.mod.no.s2)

#step-wise selection 
step.mod.hbl.sp01.no.s2 <- step(hbl.sp.mod.no.s2)
top.hbl.sp.mod1.no.s2 <- get_model(step.mod.hbl.sp01.no.s2)
summary(top.hbl.sp.mod1.no.s2)
MuMIn::r.squaredGLMM(top.hbl.sp.mod1.no.s2)

################################################################################################

#2 year
#Modeling collection year + 1 year before (ave_2_year_tmean or ppt)

#body mass
bm.sp.mod.no.s2.2yr <- lmer(BM.log10 ~ ave_2_yr_tmean + MAP.log_2yr + season + sex + pop.den.log + 
                              hibernation + habitat_buffer_three + diurnal_nocturnal + mean_BM_binary +
                              ave_2_yr_tmean:mean_BM_binary + pop.den.log:mean_BM_binary +
                              ave_2_yr_tmean:pop.den.log + ave_2_yr_tmean:hibernation + ave_2_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                              (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=BM.traits.climate.sp,
                            control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(bm.sp.mod.no.s2.2yr)

#step-wise selection 
step.mod.bm.sp01.no.s2.2yr <- step(bm.sp.mod.no.s2.2yr)
top.bm.sp.mod.no.s2.2yr <- get_model(step.mod.bm.sp01.no.s2.2yr)
summary(top.bm.sp.mod.no.s2.2yr)

MuMIn::r.squaredGLMM(top.bm.sp.mod.no.s2.2yr)

############################################################
#2 year HBL

hbl.sp.mod.no.s2.2yr <- lmer(HBL.log10 ~ ave_2_yr_tmean + MAP.log_2yr + season + sex + pop.den.log + 
                               hibernation + habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                               ave_2_yr_tmean:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                               ave_2_yr_tmean:pop.den.log + ave_2_yr_tmean:hibernation + ave_2_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                               (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=HBL.traits.climate.sp,
                             control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(hbl.sp.mod.no.s2.2yr)

#step-wise selection 
step.mod.hbl.sp01.no.s2.2yr <- step(hbl.sp.mod.no.s2.2yr)
top.hbl.sp.mod1.no.s2.2yr <- get_model(step.mod.hbl.sp01.no.s2.2yr)
summary(top.hbl.sp.mod1.no.s2.2yr)
MuMIn::r.squaredGLMM(top.hbl.sp.mod1.no.s2.2yr)

################################################################################################

#5 year
#Modeling collection year + 5 year before (ave_5_year_tmean or ppt)

#body mass
bm.sp.mod.no.s2.5yr <- lmer(BM.log10 ~ ave_5_yr_tmean + MAP.log_5yr + season + sex + pop.den.log + 
                              hibernation + habitat_buffer_three + diurnal_nocturnal + mean_BM_binary +
                              ave_5_yr_tmean:mean_BM_binary + pop.den.log:mean_BM_binary +
                              ave_5_yr_tmean:pop.den.log + ave_5_yr_tmean:hibernation + ave_5_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                              (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=BM.traits.climate.sp,
                            control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(bm.sp.mod.no.s2.5yr)

#step-wise selection 
step.mod.bm.sp01.no.s2.5yr <- step(bm.sp.mod.no.s2.5yr)
top.bm.sp.mod.no.s2.5yr <- get_model(step.mod.bm.sp01.no.s2.5yr)
summary(top.bm.sp.mod.no.s2.5yr)

MuMIn::r.squaredGLMM(top.bm.sp.mod.no.s2.5yr)

###################################

#5 year HBL

hbl.sp.mod.no.s2.5yr <- lmer(HBL.log10 ~ ave_5_yr_tmean + MAP.log_5yr + season + sex + pop.den.log + 
                               hibernation + habitat_buffer_three + diurnal_nocturnal + mean_HBL_binary +
                               ave_5_yr_tmean:mean_HBL_binary + pop.den.log:mean_HBL_binary +
                               ave_5_yr_tmean:pop.den.log + ave_5_yr_tmean:hibernation + ave_5_yr_tmean:habitat_buffer_three + pop.den.log:diurnal_nocturnal +
                               (1 | NA_L1NAME) + (1 | binomial2) + (1 | decade2), data=HBL.traits.climate.sp,
                             control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(hbl.sp.mod.no.s2.5yr)

#step-wise selection 
step.mod.hbl.sp01.no.s2.5yr <- step(hbl.sp.mod.no.s2.5yr)
top.hbl.sp.mod1.no.s2.5yr <- get_model(step.mod.hbl.sp01.no.s2.5yr)
summary(top.hbl.sp.mod1.no.s2.5yr)
MuMIn::r.squaredGLMM(top.hbl.sp.mod1.no.s2.5yr)

