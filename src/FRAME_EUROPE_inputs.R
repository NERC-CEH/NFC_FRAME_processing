
frame.input.europe <- function(){
  # go through each nominated pollutant to process
  for(species in pollutants){
  
  print(paste0(Sys.time(),": Producing FRAME-EUROPE input file for ",species," ",year,"..."))
  
  ##############################################################################
  ####    Create emissions surfaces for NH3/NOx/SOx for EMEP grid for 2017  ####
  ####   Using new (2019) EMEP netcdfs with distributions from 1990 - 2017  ####
  ####     These EMEP netcdfs are 0.1 degree cell totals (all countries)    ####
  ####   FRAME still uses the old FRAME 50km grid so requires conversion    ####
  ##############################################################################
  
    ##~~ OCT 2020: ISSUE ~~## 
    ## EMEP 50km2 projection no longer works under PROJ6 - must use manual reprojection (!)
    
  # domain
  euro.frame.ext <- suppressWarnings(raster(xmn = -20, xmx = 180, ymn = -40, ymx = 200, resolution=1))
  
  ######## Process EMEP 0.1 degree euro data into FRAME inputs ########
  
  ## emissions sector table with injection heights - GNFR categories
  sectors <- data.table(sector = LETTERS[1:13], GNFR = c("publicpower","industry","otherstationarycomb","fugitive","solvents","roadtransport","shipping","aviation","offroad","waste","agrilivestock","agriother","other"), top.height = c(12,10,16,6,6,2,5,5,2,2,1,1,1), bot.height = c(10,8,1,1,1,1,1,1,1,1,1,1,1))
  
  ## conversion factor of NOx/SOx/NH3 to N or S
  if(species=="nox")         { AV <- (1/46)*14 
  } else if(species=="sox")  { AV <- (1/64)*32 
  } else if(species=="nh3")  { AV <- (1/17)*14 }
  
  ## IF THERE IS A NETCDF PRESENT FOR THE LATEST EMISSIONS YEAR: USE THIS
  ## IF THERE IS NOT, USE THE TEXT FILES TO MAKE A SERIES OF RASTERS
  #  netcdf downloaded from EMEP (https://www.ceip.at/the-emep-grid/gridded-emissions)
  ## OCT 2020: EMEP 50km2 proj4string is defunct in PROJ6, using manual reprojection below
  
  ## name of latest netcdf downloaded from EMEP
  ## this is the latest version of gridden emissions and ahsa all years in it
  ncname <- paste0("./../EMEP_diffuse_data/",species,"/",toupper(substr(species,1,2)),substr(species,3,3),"_",year+2,"_GRID_1990_to_",year,".nc")
  
  ## for every sector in the table above, loop through, rasterise the data for year, and stack
  ## plus put the totals into a list
  spec.st <- stack()
  emep.tots <- list()
  frame.tots <- list()
  
  if(file.exists(ncname)){
    
    ## THIS IF LOOP IS IF NETCDF EXISTS
    
  # open netcdf connection
  nc <- nc_open(ncname)
  
  for(S in sectors[,GNFR]){
    
    # time band is chosen by converting integer in netcdf to date and grepping chosen year.
    # create raster of 0.1 degrees in LL
    # data is in Tg cell-1 yr-1, all countries summed (if >1 in a cell)
    r <- raster(ncname, varname = S, band = grep(latest.emep.netcdf,as.Date(nc$dim$time$vals, origin = '1850-01-01')))
    crs(r) <- LL
    
    # EMEP emissions surface total, convert to kt from Mt, and list
    r.tot <- cellStats(r, sum) * 1000
    emep.tots[[paste0(species,"_",S)]] <- data.frame(paste0(species,"_",S),r.tot)
    
    # convert to points and reproject to EMEP 50km, for FRAME-Europe
    r.dot <- rasterToPoints(r, spatial=TRUE)
    
    r.dt <- as.data.table(as.data.frame(r, centroids=TRUE,xy=TRUE,na.rm=FALSE))
    
    radconv <- pi/180
    M <- (6370/50)*(1+sin(60*radconv))
    
    r.dt[, x_EMEP := 8 + ( M * tan((pi/4)-((y*radconv)/2)) * sin((x*radconv) - (-32*radconv)) ) ]
    r.dt[, y_EMEP := 110 - ( M * tan((pi/4)-((y*radconv)/2)) * cos((x*radconv) - (-32*radconv)) ) ]
    
    #r.dot.emep <- spTransform(r.dot, CRS(EMEP50))
    r.dt[,c("x","y") := NULL]
    r.dot.emep <- copy(r.dt)
    coordinates(r.dot.emep) <- ~x_EMEP + y_EMEP
    
    # rasterize points back to FRAME EU grid and remove ID attribute
    #r.50 <- rasterize(r.dot.emep, euro.frame.ext, fun=sum)
    r.50 <- rasterize(r.dot.emep, euro.frame.ext, fun=sum)
    r.50 <- dropLayer(r.50, c(1))
    
    # FRAME input emissions surface total, convert to kt from Mt, and list
    r.f.tot <- cellStats(r.50, sum) * 1000
    frame.tots[[paste0(species,"_",S)]] <- data.frame(paste0(species,"_",S),r.f.tot)
    
    # name and add to raster stack 
    names(r.50) <- paste0(species,"_",S)
    spec.st <- stack(spec.st, r.50)
    
  }
  
  # close netcdf connection
  nc_close(nc)
  
  
  }else{
    
    ## THIS IF LOOP IS IF TEXT FILES ARE USED (NO NETCDF)
    
    emis.files <- list.files(paste0("./../EMEP_diffuse_data/",species,"/",toupper(substr(species,1,2)),substr(species,3,3),"_",year+2,"_GRID_",year), full.names = T)
    
    for(S in sectors[,GNFR]){
      
      ## pick out required file from text file list
      if(S =="other"){
        file.needed <- paste0("./../EMEP_diffuse_data/",species,"/",toupper(substr(species,1,2)),substr(species,3,3),"_",year+2,"_GRID_",year,"/",toupper(substr(species,1,2)),substr(species,3,3),"_M_Other_",year+2,"_GRID_",year,".txt")
      }else{
        file.needed <- emis.files[grep(S, emis.files, ignore.case = T)]
      }
      
      dat <- fread(file.needed, skip = 4)
      
      # convert to Mt to be in line with netcdf
      dat[, EMISSION := EMISSION / 1000000]
      
      # EMEP emissions surface total, convert to kt from t, and list
      dat.tot <- sum(dat[,EMISSION])/1000
      emep.tots[[paste0(species,"_",S)]] <- data.frame(paste0(species,"_",S),dat.tot)
      
      ## immediately produce EMEP 50km2 coordinates
      setnames(dat, c("LONGITUDE","LATITUDE"), c("x","y"))
      
      radconv <- pi/180
      M <- (6370/50)*(1+sin(60*radconv))
      
      dat[, x_EMEP := 8 + ( M * tan((pi/4)-((y*radconv)/2)) * sin((x*radconv) - (-32*radconv)) ) ]
      dat[, y_EMEP := 110 - ( M * tan((pi/4)-((y*radconv)/2)) * cos((x*radconv) - (-32*radconv)) ) ]
      
      
      #r.dot.emep <- spTransform(r.dot, CRS(EMEP50))
      dat <- dat[,c("x_EMEP","y_EMEP","EMISSION")]
      r.dot.emep <- copy(dat)
      coordinates(r.dot.emep) <- ~x_EMEP + y_EMEP
      
      # rasterize points back to FRAME EU grid and remove ID attribute
      #r.50 <- rasterize(r.dot.emep, euro.frame.ext, fun=sum)
      r.50 <- rasterize(r.dot.emep, euro.frame.ext, fun=sum)
      r.50 <- dropLayer(r.50, c(1))
      
      # FRAME input emissions surface total, convert to kt from Mt, and list
      r.f.tot <- cellStats(r.50, sum)
      frame.tots[[paste0(species,"_",S)]] <- data.frame(paste0(species,"_",S),r.f.tot)
      
      # name and add to raster stack 
      names(r.50) <- paste0(species,"_",S)
      spec.st <- stack(spec.st, r.50)
      
    }
    
    
  }
  
  
  ### convert to data table and format to FRAME inupt file ###
  dtr <- as.data.table(as.data.frame(spec.st, centroids=TRUE,xy=TRUE,na.rm=FALSE))
  
  # make sums of all sectors following data table conversion (in kt)
  dtr.sums <- as.data.frame(t(dtr[, lapply(.SD, sum, na.rm=TRUE) ][,-c(1:2)] * 1000))
  
  # read in the table that gives the exact sizes (in km2) of the EMEP grid cells
  EMEPgrid <- fread("./../EMEP_diffuse_data/Exact_EMEP_areas.csv")
  EMEPgrid[,ha := `AREA OF GRID-SQUARE (km2)` * 100]
  EMEPgrid <- EMEPgrid[,c(1,2,9)]
  names(EMEPgrid) <- c("x", "y","area_ha")
  EMEPgrid[ , c("x", "y") := list(x-0.5,y-0.5)]
  
  # join EMEP grid size (in ha) and scale sectors by real-world cell size (to get ha-1)
  dtr.ha <- EMEPgrid[dtr, on=.(x,y)]
  # remove NA cells for hectarage (and therefore some emissions data?)
  dtr.ha <- dtr.ha[!is.na(area_ha)]
  # resummarise the emissions following the hectarage addition (in kt)
  dtr.ha.sums <- as.data.frame(t(dtr.ha[, lapply(.SD, sum, na.rm=TRUE) ][,-c(1:3)] * 1000))
  
  # change all coluns to per hectare
  dtr.ha[, setdiff(colnames(dtr.ha), c("x","y","area_ha")) := lapply(.SD, function(x) x / area_ha), .SDcols = setdiff(colnames(dtr.ha), c("x","y","area_ha"))]
  
  # convert all sectors into kgs (from kt) and to N or S
  cols <- setdiff(colnames(dtr.ha), c("x","y","area_ha"))
  dtr.ha[, (cols) := lapply(.SD, function(x) x * AV), .SDcols = cols]
  dtr.ha[, (cols) := lapply(.SD, function(x) x * 1000000000), .SDcols = cols]
  
  # rework coordinates to FRAME (from EMEP) and drop the hectares column
  dtr.ha[, c("x","y") := list(x * 50000, y * 50000) ]
  dtr.ha[is.na(dtr.ha)] <- 0
  
  ## Create header lines
  dtr.ha[ , area_ha := NULL ]
  
  print(paste0(Sys.time(),":       Writing input file and summary file..."))
  
  line1 <- paste0("EMEP 50km model emission data for the FRAME-Europe model [kg ",toupper(substr(species,1,1))," / ha / y]")
  line2 <- paste0("Version....................................: ",date)
  line3 <- paste0("Author_and_source..........................: ",author)
  line4 <- paste0("Year.......................................: ",year)
  line5 <- "Coordinate_grid_reference..................: EMEP 50km (i,j)"
  line6 <- paste0("Number_of_sectors_for_emissions............: ",nrow(sectors))
  line7 <- paste0("Top_height_of_emissions_for_each_sector....: ",paste(sectors[,top.height],collapse=","))
  line8 <- paste0("Bottom_height_of_emissions_for_each_sector.: ",paste(sectors[,bot.height],collapse=","))
  line9 <- paste0("0.1 degree netcdf data for ",year,", published by EMEP in ", year+2)
  header <- c(line1,line2,line3,line4,line5,line6,line7,line8,line9)
  
  ## CREATE FOLDERS
  dir.create(file.path(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/inputs")), showWarnings = T, recursive = T)
  dir.create(file.path(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/raw_output")), showWarnings = T, recursive = T)
  
  
  ## WRITE
  file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/inputs/NFC_",species,"_",year,"_FRAMEeuro.csv"))
  newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/inputs/NFC_",species,"_",year,"_FRAMEeuro.csv")
  zz <- file(newfile,"w+")
  writeLines(header, newfile)
  write.table(dtr.ha,newfile,append=TRUE,row.names=FALSE,sep=",") 
  close(zz)
  
  
  ## SUMMARY DATA
  
  # collapse and join together all summary tables 
  # totals from EMEP
  emep.tots.df <- as.data.table(rbindlist(emep.tots))
  names(emep.tots.df) <- c("Sector","EMEP.data.kt")
  setkey(emep.tots.df, Sector)
  
  # totals converted to FRAME grid
  frame.tots.df <- as.data.table(rbindlist(frame.tots))
  names(frame.tots.df) <- c("Sector","regrid.FRAME.kt")
  setkey(frame.tots.df, Sector)
  
  # totals converted to data table
  dtr.sums.df <- data.table(Sector = row.names(dtr.sums), kt = dtr.sums, row.names = NULL)
  names(dtr.sums.df) <- c("Sector","tabled.FRAME.kt")
  setkey(dtr.sums.df, Sector)
  
  # totals following conversion to ha and loss of rows
  dtr.ha.sums.df <- data.table(Sector = row.names(dtr.ha.sums), kt = dtr.ha.sums, row.names = NULL)
  names(dtr.ha.sums.df) <- c("Sector","hectared.FRAME.kt")
  setkey(dtr.ha.sums.df, Sector)
  
  # totals of input file following backward conversion from kg ha-1 yr-1
  ## THIS DOES NOT USE EXACT HECTARAGE ##
  data.file <- as.data.frame(t(((dtr.ha[, lapply(.SD, sum, na.rm=TRUE) ][,-c(1:2)]/1000000*250000)/AV)))
  data.file$Sector <- row.names(data.file)
  data.file <- as.data.table(data.file)
  names(data.file) <- c("input.file.FRAME.kt","Sector")
  setkey(data.file, Sector)
  
  # JOIN
  summary <- emep.tots.df[frame.tots.df,][dtr.sums.df,][dtr.ha.sums.df,][data.file,]
  summary[,input.EMEP.ratio := input.file.FRAME.kt / EMEP.data.kt]
  
  line1 <- "EMEP INPUT DATA by GNFR THROUGH PROCESSING CHAIN in kt"
  line2 <- paste0("This is ",year," EMEP netcdf data, produced in ",year+2,", for the ",year," NFC run")
  line3 <- "Col2: Original EMEP surface"
  line4 <- "Col3: EMEP surface regridded to FRAME 50km"
  line5 <- "Col4: FRAME 50km converted to tables"
  line6 <- "Col5: FRAME 50km table converted to hectares"
  line7 <- "Col6: FRAME-Europe input file back converted from kg N/S ha-1 yr-1"
  line8 <- "Data is lost when converting to hectares as many cells dont have hectare count - these cells are Russia etc"
  header.summary <- c(line1, line2, line3, line4, line5, line6, line7, line8)
  
  ## WRITE
  file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/inputs/NFC_",species,"_",year,"_FRAMEeuro_SUMMARY.csv"))
  newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/inputs/NFC_",species,"_",year,"_FRAMEeuro_SUMMARY.csv")
  zz <- file(newfile,"w+")
  writeLines(header.summary, newfile)
  suppressWarnings(write.table(summary,newfile,append=TRUE,row.names=FALSE,sep=","))
  close(zz)
  
  print(paste0(Sys.time(),":       COMPLETE."))
  
  # Europe inputs done
  }
}
