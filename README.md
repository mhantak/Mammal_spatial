# Mammal_spatial
Data and code for "Mammalian body size is determined by interactions between climate, urbanization, and life history traits"

#Mammal_spatial_models.R
Annotated R code.

#Mammal_temporal_BM.csv - Body mass data set 
X1st_body_mass: mammal body mass (g)
catalognumber: museum catalog number
collectioncode: museum collection code
decimallatitude: latitude for specimen record 
decimallongitude: longtitude for specimen record 
dynamicproperties: extra record information
institutioncode: museum 
lifestage: record life stage
sex: sex of record 
X1st_tail_length: tail length of record
X1st_total_length: total length of record
HB.length: total length minus tail length of each record
HBL.log: log head-body length
TotL.log: log total length 
BM.log: log body mass
binomial2: scientific name of record
year: year of record collection 
month: month of record collection 
day: day of record collection 
season: season of record collection 
source: the source the record was obtained from
MAT_1yr: Mean Annual Temperature
MAP_1yr: Mean Annual Precipitation 
decade: decade of record collection
pop_density_1km2: human population density 1km2
pop_density_4km2: human population density 4km2
pop_density_10km2: human population density 10km2
NA_L1NAME: Ecoregion

#Mammal_temporal_HBL.csv - Head-body length set
X1st_body_mass: mammal body mass (g)
catalognumber: museum catalog number
collectioncode: museum collection code
decimallatitude: latitude for specimen record 
decimallongitude: longtitude for specimen record 
dynamicproperties: extra record information
institutioncode: museum 
lifestage: record life stage
sex: sex of record 
X1st_tail_length: tail length of record
X1st_total_length: total length of record
HB.length: total length minus tail length of each record
HBL.log: log head-body length
TotL.log: log total length 
BM.log: log body mass
binomial2: scientific name of record
year: year of record collection 
month: month of record collection 
day: day of record collection 
season: season of record collection 
source: the source the record was obtained from
MAT_1yr: Mean Annual Temperature
MAP_1yr: Mean Annual Precipitation 
decade: decade of record collection
pop_density_1km2: human population density 1km2
pop_density_4km2: human population density 4km2
pop_density_10km2: human population density 10km2
NA_L1NAME: Ecoregion

#Traits_all.csv - Life history data for each species
binomial2: scientific name of record
average_body_mass: species avergae body mass
hibernation: hibernation information for each species 
hibernation_binary: does the species hibernate (yes/no)
habitat_buffer: habitat buffering information for each species 
buffered_binary: does the species use habitat buffereing (yes/no)
diurnal_nocturnal: is the species diurnal, nocturnal, or both diurnal/nocturnal 
reproductive_rate: species reproductive rate 
reprdocutive_coding: simpilified reproductive rate data 
mating_system: mating type of species 
average_HB_length: species average head-body length 

Mammal2.tre 
Global mammal consensus phylogeny from Upham et al. (2019)

