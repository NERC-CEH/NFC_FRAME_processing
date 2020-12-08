single.year.nh3.concs <- function(){
  
  ##############################################################################
  ####  Function to create a single year NH3 concentration file for NFC     ####
  ####           from FRAME-UK 1km model, ug m3                             ####
  ##############################################################################
  
  print(paste0(Sys.time(),": Creating single year NH3 concentration file from FRAME for NFC ",year,"..."))
  print(paste0(Sys.time(),":      Reading ",year," data and writing..."))
  
  # latest NH3 concentration files
  latest.nh3.uncal <- raster(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output/NH3_gas_conc_non.calib_NFC_",year,"_ugm3_1km_",run.name,".tif"))
  latest.nh3.cal <- raster(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output/NH3_gas_conc_calib_NFC_",year,"_ugm3_1km_",run.name,".tif"))
  
  # read in cells file
  cells <- fread("./NFC_conc_cells.csv")
  
  ## turn into a data table, match onto the the cells file, write out to single year conc file
  dtr.uncal <- as.data.table(as.data.frame(latest.nh3.uncal, centroids=TRUE,xy=TRUE,na.rm=FALSE))
  dtr.uncal.cells <- dtr.uncal[cells, on = c("x","y")]
  
  dtr.cal <- as.data.table(as.data.frame(latest.nh3.cal, centroids=TRUE,xy=TRUE,na.rm=FALSE))
  dtr.cal.cells <- dtr.cal[cells, on = c("x","y")]
  
  # join the unclaibrated and calibrated data together
  dtr.cells <- dtr.uncal.cells[dtr.cal.cells, on = c("x","y")]
  names(dtr.cells) <- c("ie_OS","in_OS","air_NH3.uncalib","air_NH3.calib")
  
  # checks QA
  
  if(identical(nrow(dtr.uncal.cells),nrow(dtr.cal.cells),nrow(dtr.cells))==T){
    NULL
  }else{
    stop(print("ERROR: NUMBER OF ROWS FROM CALIBRATED AND UNCALIBRATED DATA JOIN HAS ALTERED - CHECK"))
  }
  
  if(nrow(dtr.cells[air_NH3.calib > air_NH3.uncalib]) > 0){
    stop(print("ERROR: CALIBRATED CELL VALUES BIGGER THAN UNCALIBRATED - CHECK"))
  }else{NULL}
  
  # write
  write.csv(dtr.cells,paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/FRAME_NFC_",year,"_NH3_concs_ugm3.csv"), row.names = F)
  
  
  print(paste0(Sys.time(),":      COMPLETE."))
  
}
# end of function



produce.3.year.nh3.mean <- function(){
  
  ##############################################################################
  ####   Function to create a 3-year rolling average, with the latest NFC   ####
  ####   year as the 3rd time step (i.e. mean is for year -1)               ####
  ##############################################################################
  
  print(paste0(Sys.time(),": Creating 3-year rolling mean for NH3 from FRAME for NFC: ",year-2," to ",year,"..."))
  
  print(paste0(Sys.time(),":      Reading ",year-2," to ",year," data..."))
  
  # reading data for nominated year plus previous 2
  nh3.yr <- fread(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/FRAME_NFC_",year,"_NH3_concs_ugm3.csv"))
  nh3.yr.min1 <- fread(paste0("./NAEI_NFC_",year-1,"_forFRAME/FRAME_UK/FRAME_NFC_",year-1,"_NH3_concs_ugm3.csv"))
  nh3.yr.min2 <- fread(paste0("./NAEI_NFC_",year-2,"_forFRAME/FRAME_UK/FRAME_NFC_",year-2,"_NH3_concs_ugm3.csv"))
  
  names(nh3.yr)[3:4] <- paste0(names(nh3.yr)[3:4],".",year)
  names(nh3.yr.min1)[3:4] <- paste0(names(nh3.yr.min1)[3:4],".",year-1)
  names(nh3.yr.min2)[3:4] <- paste0(names(nh3.yr.min2)[3:4],".",year-2)
  
  setkey(nh3.yr, ie_OS, in_OS)
  setkey(nh3.yr.min1, ie_OS, in_OS)
  setkey(nh3.yr.min2, ie_OS, in_OS)
  
  # join together the calibrated data & remove uncalibrated
  nh3.3yr <- nh3.yr.min2[nh3.yr.min1,][nh3.yr,]
  
  if(nrow(nh3.3yr) != nrow(nh3.yr)){
    stop(print("ERROR: JOINING OF SINGLE YEARS HAS LOST ROWS OF DATA - CHECK"))
  }else{}
  
  nh3.3yr[,grep("uncalib",names(nh3.3yr)) := NULL]
  
  # make a mean of the calibrated data
  mean.cols <- names(nh3.3yr)[grep("calib",names(nh3.3yr))]
  
  nh3.3yr[, avg.nh3 := rowMeans(.SD), .SDcols = mean.cols ]
  
  # drop and rename columns
  cols.keep <- c("ie_OS","in_OS","avg.nh3")
  nh3.3yr <- nh3.3yr[,..cols.keep]
  names(nh3.3yr) <- c("easting","northing","conc_NH3")
  
  print(paste0(Sys.time(),":      Writing ",year-2," to ",year," data..."))
  
  # write out 3-year mean
  write.csv(nh3.3yr,paste0("./NAEI_NFC_",year,"_forFRAME/CLe/FRAME_NFC_",year-2,"_to_",year,"_mean_calib_NH3_concs_ugm3.csv"), row.names = F)
  
  # rasterize and write also
  r <- rasterFromXYZ(nh3.3yr)
  writeRaster(r, paste0("./NAEI_NFC_",year,"_forFRAME/CLe/FRAME_NFC_",year-2,"_to_",year,"_mean_calib_NH3_concs_ugm3.tif"), overwrite=T)
  
  print(paste0(Sys.time(),":      COMPLETE."))
  
}
# end of function



process.NH3.for.CLe <- function(){
  
  ###############################################################################
  ####    Create a Critical Levels input file from the 3-year mean data      ####
  ###############################################################################
  
  print(paste0(Sys.time(),": Creating Critical Levels NH3 input file for ",year-2," to ",year,"..."))
  
  # extent for Critical Loads modelling
  crit.ext <- extent(0,700000,0,1200000)
  
  # read in the data
  conc.dat <- fread(paste0("./NAEI_NFC_",year,"_forFRAME/CLe/FRAME_NFC_",year-2,"_to_",year,"_mean_calib_NH3_concs_ugm3.csv")) 
  
  print(paste0(Sys.time(),":      Writing CLe input file for NH3 concs ",year-2," to ",year," data..."))
  
  # create header info
  line1 <- paste0("NFC NH3 3-year concentration data from the FRAME-UK 1km model")
  line2 <- paste0("Calibrated by median bias per year")
  line3 <- paste0("Simulation Date/Time: ",Sys.Date())
  line4 <- paste0("Units: Air concentration ug/m3")
  line5 <- paste0("Model Version Number: 9-15-0")
  line6 <- paste0("Authors: NFC - CEH")
  line7 <- paste0("Start year: ",year-2)
  line8 <- paste0("End year: ", year)
  line9 <- "Coordinate_grid_reference: BNG"
  line10 <- ""
  header <- c(line1,line2,line3,line4,line5,line6,line7,line8,line9,line10)
  
  # write out blank file and populate
  file.create(paste0("./NAEI_NFC_",year,"_forFRAME/CLe/NFC_",year-2,"_to_",year,"_mean_NH3_conc_CLe_1km.csv"),overwrite=T)
  newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/CLe/NFC_",year-2,"_to_",year,"_mean_NH3_conc_CLe_1km.csv")
  zz <- file(newfile,"w+")
  writeLines(header, newfile)
  suppressWarnings(write.table(conc.dat,newfile,append=TRUE,row.names=FALSE,sep=","))
  close(zz)

  print(paste0(Sys.time(),":      COMPLETE."))
  
  
}
# end of function
