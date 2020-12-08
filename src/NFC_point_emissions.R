
NFC.point.source.inputs <- function(uk.write, eire.write, uk.eire.write, run.name){
  
  ##############################################################################
  ####   Use NAEI point source emission data to create point input file     ####
  ####   for FRAME-UK; important for plume emissions                        ####
  ##############################################################################
  
  # go through each nominated pollutant to process
  for(species in pollutants){
    
    print(paste0(Sys.time(),": Creating new FRAME point file for ",species," ",year,"..."))
  
  ######### Directories #############
  NAEI.pt.wd <- "//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/NAEI_data/point"
  EIRE.wd <- "//nercbuctdb.ad.nerc.ac.uk/projects1/NEC03642_Mapping_Ag_Emissions_AC0112/NAEI_data_and_SNAPS/EIRE_data"
  
  # scaling factor for species to N or S
  AV <- ifelse(species=="nh3",(1/17*14), ifelse(species=="nox",(1/46*14),(1/64*32)))
  
  # read the latest point sector to SNAP coversion table, subset to pollutant
  pts.sec.to.snap <- fread("./../lookups/points_sectors_to_SNAP.csv")
  pts.sec.to.snap <- pts.sec.to.snap[Pollutant == species]
  pts.sec.to.snap[,c("Pollutant","Sector","Group","NFR") := NULL]
  
  # create a FRAME input file for given year & species, for UK and Ireland
  
  #############################################################
  #### UK POINTS from NAEI data - input to FRAME is Gg N/S ####
  #############################################################
  
    print(paste0(Sys.time(),":       Reading, subsetting and formatting UK data..."))
  
    # NAEI point data - in tonnes
    NAEI.pts <- as.data.table(read_excel(paste0(NAEI.pt.wd,"/raw_data/NAEIPointsSources_",year,".xlsx"), sheet = "Data"))
  
    # subset NAEI data
    NAEI.species <- ifelse(species == "nox","NOx",ifelse(species == "sox", "SO2", "NH3"))
    NAEI.sub <- NAEI.pts[Pollutant == NAEI.species]
    
    # sum all NOx point emissions in NAEI
    naei.tot <- data.table(Area = "UK", Sector = "Total", download.kt = sum(NAEI.sub[,Emission/1000]), row.names = NULL)    
    setkey(naei.tot, Area,Sector)
    
    # read in stack info for some sites, for emission plumes
    ##~~ ISSUE: STACK HEIGHTS NEED UPDATING ~~##
    stack.info <- fread(paste0(NAEI.pt.wd,"/raw_data/all_point_stack_data.csv"))[,1:5]
    
    # join any stack height info to NAEI data
    NAEI.stack <- stack.info[NAEI.sub, on=c("PlantID")]
    
    # join SNAP data based on SectorID
    NAEI.stack.SNAP <- NAEI.stack[pts.sec.to.snap, on = "SectorID"]
    
    # sum all NOx point emissions by SNAP following joins
    naei.snap.tots <- cbind(rep("UK",length(unique(NAEI.stack.SNAP[,SNAP]))),NAEI.stack.SNAP[,sum(Emission)/1000,by=SNAP])
    names(naei.snap.tots) <- c("Area","Sector","sectored.kt")
    naei.snap.tots <- rbind(naei.snap.tots, data.frame(Area = "UK", Sector="Total",t(colSums(naei.snap.tots[,3]))))
    setkey(naei.snap.tots, Area,Sector)
    
    # convert units from Mg to Gg (tonne to kt) and to N/S
    NAEI.stack.SNAP[, paste0("Gg ",toupper(substr(species,1,1))," yr-1") := (Emission / 1000) * AV ]
    
    # set Area and set NAs
    NAEI.stack.SNAP[,AREA:="UK"]
    for(col in c("Height","Velocity","Temp","Diameter")) set(NAEI.stack.SNAP, which(is.na(NAEI.stack.SNAP[[col]])), col, 9999)
    
    # format to final shape
    col.keep <- c("Easting","Northing",paste0("Gg ",toupper(substr(species,1,1))," yr-1"),"Height","Velocity","Temp","Diameter","PlantID","SNAP","AREA")
    NAEI.final <- NAEI.stack.SNAP[,..col.keep]
    
    
  ###################################################
  #### IRISH POINTS - same input file as UK data ####
  ###################################################
    
    #print(paste0(Sys.time(),":       Reading, subsetting and formatting Eire data..."))
    
    # taking from the EPRTR database - the European point source register
    # extract of EPRTR database must have already taken place via;
    # ...\NAEI_data_and_SNAPS\E-PRTR_point_data\EPRTR_extract.R
    
    # read E-PRTR data - data is in kg
    #EPRTR.pts <- fread(paste0("./../E-PRTR_point_data/",species,"/",species,"_pts_",year,"_E-PRTR.csv"))
    
    # subset for Ireland
    #EPRTR.Ire <- EPRTR.pts[CountryName == "Ireland"]
    
    # sum all NOx point emissions in E-PRTR
    #eprtr.tot <- data.table(Area = "Eire", Sector = "Total", download.kt = sum(EPRTR.Ire[,TotalQuantity/1000000]), row.names = NULL)    
    #setkey(eprtr.tot, Area,Sector)
    
    # add columns
    #EPRTR.Ire[, c("Height","Velocity","Temp","Diameter") := 9999]
    
    # sum all NOx point emissions in E-PRTR by GNFR grouping (EMEP version of SNAP)
    #eprtr.snap.tots <- cbind(rep("Eire",length(unique(EPRTR.Ire[,GNFR14]))),EPRTR.Ire[,sum(TotalQuantity)/1000000,by=GNFR14])
    #names(eprtr.snap.tots) <- c("Area","Sector","sectored.kt")
    #eprtr.snap.tots <- rbind(eprtr.snap.tots, data.frame(Area = "Eire", Sector="Total",t(colSums(eprtr.snap.tots[,3]))))
    #setkey(eprtr.snap.tots, Area,Sector)
    
    # convert units
    #EPRTR.Ire[, paste0("Gg ",toupper(substr(species,1,1))," yr-1") := (TotalQuantity / 1000000) * AV ]
    
    # set Area
    #EPRTR.Ire[,AREA:="Eire"]
    
    # change some column names
    #setnames(EPRTR.Ire,c("Xbng","Ybng","GNFR14","FacilityID"), c("Easting","Northing","SNAP","PlantID"))
    
    # format to final shape
    #col.keep <- c("Easting","Northing",paste0("Gg ",toupper(substr(species,1,1))," yr-1"),"Height","Velocity","Temp","Diameter","PlantID","SNAP","AREA")
    #EIRE.final <- EPRTR.Ire[,..col.keep]
  
  ###################################################
  ### Combine UK & EIRE and write for FRAME input ###
  ###################################################
    
    if(uk.write == T){
      
      
      print(paste0(Sys.time(),":       writing UK point data..."))
      
      # create header info
      line1 <- paste0("UK ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1))," point source emissions for the FRAME-UK model [Gg ",toupper(substr(species,1,1))," yr-1]")
      line2 <- run.name
      line3 <- paste("Version....................................:",Sys.Date())
      line4 <- paste("Author_and_source..........................: Sam Tomlinson (CEH)")
      line5 <- paste("Year.......................................:",year)
      line6 <- "Coordinate_grid_reference..................: BNG"
      line7 <- "Number_of_sectors_for_emissions............: 11"
      line8 <- "Comment...................................: 9999 no data"
      line9 <- paste("Comment...................................: UK emissions from NAEI")
      header.1 <- c(line1,line2,line3,line4,line5,line6,line7,line8,line9)
      
      # create final header line10 (based on totals of each SNAP in kt of N or S)
      
      line10 <- paste0("Total UK emissions (Gg ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1)),"),",round(sum(NAEI.final[,get(paste0("Gg ",toupper(substr(species,1,1))," yr-1"))]),1))
      header.all <- c(header.1,line10)
      
      ## CREATE FOLDERS
      dir.create(file.path(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs")), showWarnings = T, recursive = T)
      dir.create(file.path(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output")), showWarnings = T, recursive = T)
      dir.create(file.path(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raw_output")), showWarnings = T, recursive = T)
      
      # create the output file, write header first, append the dataframe
      # FRAME input files are MAX 20 characters (WITH ".csv")
      
      input.file <- paste0(species,"-",substr(species,1,1),"_UK_pts_",substr(year,3,4))
      
      file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv"))
      newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv")
      zz <- file(newfile,"w+")
      writeLines(header.all, newfile)
      suppressWarnings(write.table(NAEI.final,newfile,append=TRUE,row.names=FALSE,sep=","))
      close(zz) 
      
      
    }else{}
    
    if(eire.write == T){
      
      
      print(paste0(Sys.time(),":       writing EIRE point data..."))
      
      # create header info
      line1 <- paste0("EIRE ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1))," point source emissions for the FRAME-UK model [Gg ",toupper(substr(species,1,1))," yr-1]")
      line2 <- run.name
      line3 <- paste("Version....................................:",Sys.Date())
      line4 <- paste("Author_and_source..........................: Sam Tomlinson (CEH)")
      line5 <- paste("Year.......................................:",year)
      line6 <- "Coordinate_grid_reference..................: BNG"
      line7 <- "Number_of_sectors_for_emissions............: 11"
      line8 <- "Comment...................................: 9999 no data"
      line9 <- paste("Comment...................................: Eire from E-PRTR database")
      header.1 <- c(line1,line2,line3,line4,line5,line6,line7,line8,line9)
      
      # create final header line10 (based on totals of each SNAP in kt of N or S)
      
      line10 <- paste0("Total Eire emissions (Gg ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1)),"),",round(sum(EIRE.final[,get(paste0("Gg ",toupper(substr(species,1,1))," yr-1"))]),1))
      header.all <- c(header.1,line10)
      
      # create the output file, write header first, append the dataframe
      # FRAME input files are MAX 20 characters (WITH ".csv")
      
      input.file <- paste0(species,"-",substr(species,1,1),"_I_pts_",substr(year,3,4))
      
      file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv"))
      newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv")
      zz <- file(newfile,"w+")
      writeLines(header.all, newfile)
      suppressWarnings(write.table(EIRE.final,newfile,append=TRUE,row.names=FALSE,sep=","))
      close(zz) 
      
      
    }else{}
    
    if(uk.eire.write == T){
    
    print(paste0(Sys.time(),":       Combining UK & Eire data and writing..."))
    
    # combine datasets
    #all.pts <- rbindlist(list(NAEI.final,EIRE.final))
      all.pts <- NAEI.final
    
    # create header info
    line1 <- paste0("UK & Eire ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1))," point source emissions for the FRAME-UK model [Gg ",toupper(substr(species,1,1))," yr-1]")
    line2 <- run.name
    line3 <- paste("Version....................................:",Sys.Date())
    line4 <- paste("Author_and_source..........................: Sam Tomlinson (CEH)")
    line5 <- paste("Year.......................................:",year)
    line6 <- "Coordinate_grid_reference..................: BNG"
    line7 <- "Number_of_sectors_for_emissions............: 11"
    line8 <- "Comment...................................: 9999 no data"
    line9 <- paste("Comment...................................: UK emissions from NAEI & Eire from E-PRTR database")
    header.1 <- c(line1,line2,line3,line4,line5,line6,line7,line8,line9)
    
    # create final header line10 (based on totals of each SNAP in kt of N or S)
    
    line10 <- paste0("Total UK emissions (Gg ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1)),"),",round(sum(all.pts[AREA=="UK",get(paste0("Gg ",toupper(substr(species,1,1))," yr-1"))]),1),",Total Eire emissions (Gg ",toupper(substr(species,1,2)),substr(species,3,3),"-",toupper(substr(species,1,1)),"),",round(sum(all.pts[AREA=="Eire",get(paste0("Gg ",toupper(substr(species,1,1))," yr-1"))]),1))
    
    
    header.all <- c(header.1,line10)
    
    # create the output file, write header first, append the dataframe
    # FRAME input files are MAX 20 characters (WITH ".csv")
    
    input.file <- paste0(species,"-",substr(species,1,1),"_UKI_pts_",substr(year,3,4))
    
    file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv"))
    newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",input.file,".csv")
    zz <- file(newfile,"w+")
    writeLines(header.all, newfile)
    suppressWarnings(write.table(all.pts,newfile,append=TRUE,row.names=FALSE,sep=","))
    close(zz) 
    
    }else{}
    
    
    
    
    
    ## SUMMARY DATA
    
    all.snap.tots <- all.pts[,sum(get(paste0("Gg ",toupper(substr(species,1,1))," yr-1"))) / AV, by=.(AREA,SNAP)]
    names(all.snap.tots) <- c("Area","Sector","processed.kt")
    
    all.snap.tots <- rbind(all.snap.tots, data.frame(Area = "UK", Sector = "Total", all.snap.tots[Area=="UK",list(processed.kt = sum(processed.kt))]))
    #all.snap.tots <- rbind(all.snap.tots, data.frame(Area = "Eire", Sector = "Total", all.snap.tots[Area=="Eire",list(processed.kt = sum(processed.kt))]))
    all.snap.tots <- rbind(all.snap.tots, data.frame(Area = "Eire", Sector = "Total", processed.kt=0))
    setkey(all.snap.tots, Area,Sector)
    
    # join summary tables
    summary.uk <- naei.tot[naei.snap.tots]
    #summary.eire <- eprtr.tot[eprtr.snap.tots]
    #summary <- rbindlist(list(summary.uk, summary.eire))
    summary <- summary.uk
    summary <- summary[all.snap.tots, on = c("Area","Sector")]
    
    line1 <- "EMISSIONS POINT DATA by SECTOR THROUGH PROCESSING CHAIN in kt"
    line2 <- "Col2: Sector - SNAP or GNFR or TOTAL"
    line3 <- "Col3: NAEI or EPRTR values"
    line4 <- "Col4: NAEI or EPRTR values after sector grouping"
    line5 <- "Col5: NAEI or EPRTR values after processing to input file"
    header.summary <- c(line1,line2,line3,line4,line5)
    
    ## write summary
    file.create(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",species,"-",substr(species,1,1),"_UKI_pts_",substr(year,3,4),"_SUMMARY.csv"))
    newfile <- paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/inputs/",species,"-",substr(species,1,1),"_UKI_pts_",substr(year,3,4),"_SUMMARY.csv")
    zz <- file(newfile,"w+")
    writeLines(header.summary, newfile)
    suppressWarnings(write.table(summary,newfile,append=TRUE,row.names=FALSE,sep=","))
    close(zz)
  
    print(paste0(Sys.time(),":       COMPLETE."))
  
  
  } # end of pollutant loop

}