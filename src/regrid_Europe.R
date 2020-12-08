
regrid.frame.europe <- function(){
  
  print(paste0(Sys.time(),": Re-gridding FRAME-EUROPE raw output files for ",year,"..."))
  
  ##############################################################################
  ####   Take FRAME-EUROPE concnetration outputs @50km and re-grid them     ####
  ####   for use in the FRAME-UK 1km model                                  ####
  ##############################################################################
  
  
  ##~~ NOV 2019: ISSUE ~~## 
  ##   The re-gridded data needs to be on 5km grid, as FRAME 1km still takes 5km boundary conditions
  
  ##~~ OCT 2020: ISSUE ~~## 
  ## EMEP 50km2 projection no longer works under PROJ6 - must use manual reprojection (!)
  ## manually reproject EMEP output to Lat Lon, rasterize and reproject to BNG
  
  # re-grid resolution (m)
  regrid.res <- 5000
  
  # extended FRAME domain includes both UK & Eire
  NFC.uk.frame <- raster(xmn=-230000,xmx=700000,ymn=0,ymx=1250000, res=regrid.res, crs=BNG, vals=NA)
  tempLL.ext <- raster(xmn=-122,xmx=60,ymn=20,ymx=90, res=0.1, crs=LL, vals=NA)
  # reproject the extent to LL
  r.frame.LL <- projectRaster(NFC.uk.frame,crs=LL,res=0.01)
  NFC.uk.ext.LL <- extent(r.frame.LL)
  
  
  ######################################
  ##  regrid FRAME Europe output to   ## 
  ##   UK NFC boundary conditions     ##
  ######################################

# loop through all 8 directional European concentration files;
      # Read in directional concentration file, adjust coordinates
      # for every pollutant, clip out domain for UK FRAME and stack
      # make new FRAME-UK boundary condition file
  
  
  for(d in 1:8){
    print(paste0(Sys.time(),":       processing and writing direction ",d,"..."))
    
    # read in data
    f.eur.out <- suppressWarnings(fread(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/raw_output/frame_EU_",run.name,"_9-15-0_",year,"_000",d,".dat")))
    
    f.eur.out[,V63 := NULL]
    
    # adjust to FRAME-EUROPE 50km grid
    f.eur.out[,c("ie_OS","in_OS") := list(ie_OS/50000,in_OS/50000)]
    
    # blank stack to populate with rasters of pollutants
    spec.stack <- stack()
    
    # loop through all pollutant columns, rasterize, reproject to BNG and add to stack
    for(spec in names(f.eur.out)[grep("air",names(f.eur.out))]){
      
      cols <- c("ie_OS","in_OS",spec)
      f.eursub <- f.eur.out[, ..cols]
      
      ## dont have a conversion from EMEP anymore in proj6
      ## converting the EMEP coords to LL will create very coarse data, simply rasterizing is wrong
      ## can convert EMEP coords to LL and subset these to fit inside LL FRAME extent
      ## THEN interpolate (akima package) these onto a finer resolution
      ## then rasterize and reproject to BNG
      
      ## MANUAL lat lon conversion
      radconv <- pi/180
      M <- (6370/50)*(1+sin(60*radconv))
      
      f.eursub[, Lon := -32 + (180/pi * atan((ie_OS-8)/(110-in_OS)) ) ]
      f.eursub[, Lat := 90 - (360/pi * atan( (sqrt((ie_OS-8)^2 + (in_OS-110)^2)) / M) ) ]
      
      # cut down the new points to be inside LatLon extent - minimises data for resampling
      f.eursub <- f.eursub[Lon >= NFC.uk.ext.LL[1] & Lon <= NFC.uk.ext.LL[2]]
      f.eursub <- f.eursub[Lat >= NFC.uk.ext.LL[3] & Lat <= NFC.uk.ext.LL[4]]
      f.eursub[,c("ie_OS","in_OS") := NULL]
      
      # interpolate to finer res
      f.eursub.interp <- interp(x = f.eursub[, Lon], y = f.eursub[, Lat], z = f.eursub[[spec]], 
                                xo = seq(min(f.eursub[, Lon]), max(f.eursub[, Lon]), length = ncol(r.frame.LL)), 
                                yo=seq(min(f.eursub[, Lat]), max(f.eursub[, Lat]), length = nrow(r.frame.LL)), 
                                linear=T)
      
      # rasterize
      r.interp <- rasterFromXYZ(cbind(expand.grid(f.eursub.interp$x, f.eursub.interp$y),as.vector(f.eursub.interp$z)))
      crs(r.interp) <- LL
      
      # reproject to NFC UK - the target raster (NFC domain) contains the resolution as well as extent etc.
      rc.uk.res <- projectRaster(r.interp, NFC.uk.frame, method = "bilinear")
      names(rc.uk.res) <- spec
      
      # stack
      spec.stack <- stack(spec.stack,rc.uk.res)
    }
    
    ## turn into a data table, getting ready for FRAME input file
    dtr <- as.data.table(as.data.frame(spec.stack, centroids=TRUE,xy=TRUE,na.rm=FALSE))
    setnames(dtr,c("x","y"),c("ie","in"))
    
    dtr[is.na(dtr)] <- 0
    
    ## header text
    line1 <- paste0("FRAME ",regrid.res/1000,"km UK & Eire directional concentrations generated from FRAME-Europe for NFC ",year)
    line2 <- paste0("FRAME simulation time: ",Sys.time())
    line3 <- paste0("Emissions Year: ",year)
    line4 <- paste("Direction: ",d)
    line5 <- "units: Air concentration ug/m3"
    line6 <- ""
    line7 <- ""
    line8 <- ""
    header <- c(line1,line2,line3,line4,line5,line6,line7,line8)
    
    ## CREATE FOLDERS
    suppressWarnings(dir.create(file.path(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/re-gridded")), showWarnings = T, recursive = T))
    
    # create the output file, write header first, append the dataframe
    
    file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/re-gridded/NFC",substr(year,3,4),"_",regrid.res/1000,"k_EU_000",d,".csv"))
    newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_Europe/re-gridded/NFC",substr(year,3,4),"_",regrid.res/1000,"k_EU_000",d,".csv")
    zz <- file(newfile,"w+")
    writeLines(header, newfile)
    suppressWarnings(write.table(dtr,newfile,append=TRUE,row.names=FALSE,sep=","))
    close(zz) 
  }

  print(paste0(Sys.time(),":       COMPLETE."))
  
}
