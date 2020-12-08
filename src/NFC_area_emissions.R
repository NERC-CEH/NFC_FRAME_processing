
NFC.area.source.inputs <- function(){
  
  uk.frame <<- raster(xmn = -230000, xmx=700000, ymn = 0, ymx = 1250000, res = 1000, crs=BNG)
  
  ##############################################################################
  ####   To create diffuse source emission data for FRAME-UK input file     ####
  ##############################################################################
  
  ## THIS DATA IS now processed under the UK Emissions model
  ## \\nercbuctdb.ad.nerc.ac.uk\projects1\NEC03642_Mapping_Ag_Emissions_AC0112\NAEI_data_and_SNAPS\UK_emissions_model
  ## that model creates flat csv emission files and rasters, to be formatted here to FRAME inputs
  ## the spatial subtraction for Ireland has already taken place
  
  # go through each nominated pollutant to process
  
  for(species in pollutants){
    
    print(paste0(Sys.time(),": Creating new FRAME diffuse file for ",species," ",year,"..."))
    
    
    # emissions injection heights
    secs <- ifelse(species=="nh3",9,11)
    sec.top <- ifelse(species=="nh3",paste0("5,2,6,2,2,2,4,2,4"),paste0("7,7,7,7,7,6,2,4,1,1,1"))
    sec.bot <- ifelse(species=="nh3",paste0("1,1,3,1,1,1,1,1,1"), paste0("6,6,6,6,6,1,1,1,1,1,1"))
  
    # multiplication factor N and S:
    AV <- ifelse(species=="nh3",(1/17*14), ifelse(species=="nox",(1/46*14),(1/64*32)))
   
    if(species!="nh3"){
    
    ## read in the pre-made data csv
    premade.dat <- fread(paste0("./../Emissions_grids_plain/BNG/",species,"/diffuse/",year,"/",species,"_diff_",year,"_UKEIRE_SNAP_t_1km_BNG_",naei.map.year,"NAEImap.csv"))
    premade.dat[, V1 := NULL]
    
    ## crop down the pre-processed data to fit the FRAME model
    premade.dat <- premade.dat[x > -230000 & x < 700000]
    premade.dat <- premade.dat[y > 0 & y < 1250000]
    
    }else{
    
    ## if it's NH3, it has to be made differently for FRAME, agric must be in its own columns
    # read in pre made data for UK, remove agric. read in data for Ireland. 
    
      if(year!=2018){stop("ERROR YEAR IS NOT 2018 - INVESTIGATE")}
      
      print(paste0(Sys.time(),": NH3 requires specific agric sectors, processing..."))
      
      
      ## read in the pre-made data csv
      premade.uk <- fread(paste0("./../Emissions_grids_plain/BNG/",species,"/diffuse/",year,"/",species,"_diff_",year,"_UK_SNAP_t_1km_BNG_",naei.map.year,"NAEImap.csv"))
      premade.uk[, V1 := NULL]
      # remove agric
      premade.uk[, S10 := NULL]
      
      premade.eire <- fread(paste0("./../Emissions_grids_plain/BNG/",species,"/diffuse/",year,"/",species,"_diff_",year,"_EIRE_SNAP_t_1km_BNG_",naei.map.year,"NAEImap.csv"))
      premade.eire[, V1 := NULL]
      premade.eire.ag <- copy(premade.eire)
      premade.eire[, S10 := NULL]
      
      # process UK to one non-agric surface
      col.names <- names(premade.uk)[-c(1:2)]
      premade.uk[, total := rowSums(.SD, na.rm=T), .SDcols = col.names]
      
      uk.non.ag <- premade.uk[,c("x","y","total")]
      nh3.non.ag <- crop(extend(rasterFromXYZ(uk.non.ag), uk.frame), uk.frame)
      names(nh3.non.ag) <- "non_agric"
      
      ### NH3 Agricultural Data - kg NH3 cell-1 ###
      print(paste0(Sys.time(),": Preparing Agricultural emissions from AENEID for ", species, "..."))
      
      dat18.temp <- fread("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/AC0112-Agric_Inventory/2018_Agric_Inventory/2020-03-27_scaled_2018_maps/FRAME_input.csv")
      
      
      # cattle
      cattle <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Cattle")]), uk.frame), uk.frame)
      
      # sheep
      sheep <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Sheep")]), uk.frame), uk.frame)
      
      # pigs
      pigs <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Pigs")]), uk.frame), uk.frame)
      
      # poultry
      poultry <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Poultry")]), uk.frame), uk.frame)
      
      # horses, goats and deer
      horses.goats.deer <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Horses, Goats and Deer")]), uk.frame), uk.frame)
      names(horses.goats.deer) <- "horses.goats.deer"
      
      # Fertilisers
      fert <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Fertiliser")]), uk.frame), uk.frame)
      sludge <- crop(extend(rasterFromXYZ(dat18.temp[,c("EASTING","NORTHING","Sewage_Sludge")]), uk.frame), uk.frame)
      fertiliser.st <- stack(fert, sludge)
      fertiliser <- calc(fertiliser.st, fun=sum, na.rm=T)
      
      names(fertiliser) <- "fertiliser"
      
      # stack up and divide to tonnes
      all.agric <- stack(cattle, sheep, pigs, poultry, horses.goats.deer, fertiliser)
      all.agric <- all.agric/1000
      
      
      ## IRISH surfaces
      
      non.ag.names <- names(premade.eire)[-c(1:2)]
      premade.eire[, total := rowSums(.SD, na.rm=T), .SDcols = non.ag.names]
      
      eire.non.ag <- premade.eire[,c("x","y","total")]
      eire.nh3.non.ag <- crop(extend(rasterFromXYZ(eire.non.ag), uk.frame), uk.frame)
      names(eire.nh3.non.ag) <- "Irish.non_agric"
      
      eire.ag <- premade.eire.ag[,c("x","y","S10")]
      eire.nh3.ag <- crop(extend(rasterFromXYZ(eire.ag), uk.frame), uk.frame)
      names(eire.nh3.ag) <- "Irish.Agric"
      
      ## combine to stack
      ## stack together Ag and Non-Ag for NH3
      uk.diffuse <- stack(all.agric, nh3.non.ag, eire.nh3.ag, eire.nh3.non.ag)
      crs(uk.diffuse) <- BNG
      
      premade.dat <- as.data.table(as.data.frame(uk.diffuse, xy=T, centroid=T))
      premade.dat[is.na(premade.dat)] <- 0
    }
    
    
    ## convert to kgs per hectare
    ## convert to N or S
    col.names <- names(premade.dat)[-c(1:2)]
    premade.dat[, (col.names) := lapply(.SD, function(x) (x * AV) * 1000 / 100), .SDcols = col.names]
    
    dtr <- copy(premade.dat)
    
    print(paste0(Sys.time(),":       Writing ",species," 1km FRAME input file for UK & EIRE..."))
      # create sums for header info
      snapsum <- round(((colSums(dtr[,-(1:2)]))*100/1000000),digits=2)
      snapsum2 <- paste(round(snapsum,digits=2), collapse=",",sep="")
      
      # create header info
      line1 <- paste0(species,"-",substr(species,1,1)," emission data for the FRAME-UK model [kg /ha /y]")
      line2 <- "1km resolution over United Kingdom (NAEI) and Republic of Ireland (MapEire/EMEP) - BNG"
      line3 <- paste("Version....................................:",Sys.Date())
      line4 <- paste("Author_and_source..........................:","Sam Tomlinson")
      line5 <- paste("Year.......................................:",year)
      line6 <- paste("Project...................................:",run.name)
      line7 <- paste0("Number_of_sectors_for_emissions............: ",secs)
      line8 <- paste0("Top_height_of_emissions_for_each_sector....: ",sec.top)
      line9 <- paste0("Bottom_height_of_emissions_for_each_sector.: ",sec.bot)
      header.1 <- c(line1,line2,line3,line4,line5,line6,line7,line8,line9)
      line10 <- paste0("Total emissions (Gg ",species,"-",toupper(substr(species,1,1)),") = ",sum(snapsum),",",snapsum2)
      
      header.2 <- c(header.1,line10)
      
      # create the output file, write header first, append the dataframe
      # Input file max 20 chars with ".csv"
      input.file <- paste0(species,"-",substr(species,1,1),"_UKI_1_d_",substr(year,3,4))
      
      file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv"),overwrite=T)
      newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv")
      zz <- file(newfile,"w+")
      writeLines(header.2, newfile)
      suppressWarnings(write.table(dtr,newfile,append=TRUE,row.names=FALSE,sep=","))
      close(zz)
      
    
  
  print(paste0(Sys.time(),":       COMPLETE."))
  
 
  } # end of pollutant loop

} # end of function




