packs <- c("sp","raster","stringr","rgeos","rgdal","grid","gdalUtils","dplyr","utils","car","reshape2","ggplot2","ggrepel","data.table","stats","readr","sf","ncdf4","readxl","gridExtra","akima","rasterVis")
lapply(packs, require, character.only = TRUE)

setwd("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/SNAPS_NFC_data")

# coordinate reference systems required #
#EMEP50 <- CRS("+init=epsg:9829")
#EMEP50 <- CRS('+ellps=sphere +a=127.4 +e=0 +proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110') # EMEP european grid in FRAME
BNG <- CRS("+init=epsg:27700")  # NAEI data in British National Grid
LL <- CRS("+init=epsg:4326") # For EMEP European data

##########################################################################
#### The following functions create (as per README.md)                ####
####                                                                  ####
####   1. FRAME-EUROPE input file from EMEP data + QAQC               ####
####   2. Re-gridded european concentrations to FRAME-UK boundary     ####
####   3a. FRAME-UK input file for point sources (UK & Eire) + QAQC   ####
####   3b. FRAME-UK input file for diffuse sources (UK & Eire) + QAQC ####
####   4. QAQC + Evaluation of FRAME-UK outputs (with plots)          ####
####   5. A 3 year rolling mean for NH3 concentration                 ####
####                                                                  ####
##########################################################################

##########################################################################################################
##~~ DATA must have been processed via the UK_emissions_model to flat csv/raster files, from NAEI etc ~~##
##########################################################################################################


##### input variables: #######
pollutants <- c("nh3", "sox", "nox")  # "nh3", "nox", "sox"
year <- 2018 # emissions year to run
naei.map.year <- 2018
date <- Sys.Date()
author <- "CEH"
run.name <- "NFC18_R1-4"


##########################################################################
#### 1. EMISSIONS INPUT FILES FOR FRAME-EUROPE
source("./NFC_FRAME_processing/src/FRAME_EUROPE_inputs.R")

# runs by species
frame.input.europe()

#################################
#### RUN FRAME-EU on NEMESIS ####
#################################

##########################################################################
#### 2. RE-GRID FRAME-EUROPE OUTPUT TO 1km FRAME-UK BOUNDARY CONDITIONS
source("./NFC_FRAME_processing/src/regrid_Europe.R")

regrid.frame.europe()

##########################################################################
#### 3a. CREATE POINT EMISSIONS INPUT FILE FOR FRAME-UK
## 18/11/20 : point sources in ireland removed, all on gridded surfaces
source("./NFC_FRAME_processing/src/NFC_point_emissions.R")

NFC.point.source.inputs(uk.write = T, eire.write = F, uk.eire.write = T, run.name = run.name)

##########################################################################
#### 3a. CREATE DIFFUSE EMISSIONS INPUT FILE FOR FRAME-UK
## 18/11/20 : Eire gridded surfaces are now all Eire data (no points)
source("./NFC_FRAME_processing/src/NFC_area_emissions.R")

NFC.area.source.inputs() # UK, all (UK & Eire), Eire + # 1, 5 or 'both'

#################################
#### RUN FRAME-UK on NEMESIS #### ~! NOW RUN ON ROGUE (14/10/20) !~
#################################

##########################################################################
#### 4. PRODUCE N-dep OUPUT SURFACES
source("./NFC_FRAME_processing/src/Ndep_surfaces.R")

produce.ndep.surfaces()

##########################################################################
#### 5. EVALUATION OF NH3 CONCENTRATION OUTPUT + QAQC (AFTER running FRAME)
source("./NFC_FRAME_processing/src/FRAME_evaluation.R")

evaluate.frame.output()

##########################################################################
#### 6. NH3 CONC TABLE, 3-YEAR NH3 MEAN + CRITICAL LEVELS INPUT FILE
source("./NFC_FRAME_processing/src/Critical_levels.R")

single.year.nh3.concs()
produce.3.year.nh3.mean()
process.NH3.for.CLe()




