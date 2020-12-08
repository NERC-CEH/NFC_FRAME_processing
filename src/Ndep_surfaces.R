produce.ndep.surfaces <- function(){
  
  ##########################################################################
  ####  Function to create N dep surfaces from FRAME output - total N   ####
  ####      plus wet & dry, oxidised and reduced (kg ha-1 yr-1)         ####
  ##########################################################################
  
  print(paste0(Sys.time(),": Creating single year N-Dep surfaces from FRAME for NFC ",year,"..."))
  print(paste0(Sys.time(),":      Reading ",year," deposition data..."))
  
  for(habitat in c("grid","forest","moorland")){
    
    print(paste0("             ",habitat,"..."))
    
    FRAME.hab <- ifelse(habitat == "grid", "grd", ifelse(habitat == "forest", "for", "mor"))
    
    ## read in the data files
    NFCdat <- fread(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raw_output/dep_",run.name,"_9-15-0_",year,"unc.dat"), skip=7)
    
    #### loop through the 4 N dep components; create a raster for each one, plus a total
    
    for(nd in c("NHx_dry","NHx_wet","NOy_dry","NOy_wet","totalN")){
      
      ## make rasters, convert to N
      if(nd == "totalN"){
        
        cols.keep <- c("easting(m)","northing(m)",paste0(FRAME.hab,"_",c("NHx_dry","NHx_wet","NOy_dry","NOy_wet")))
        NFCdat.allN <- NFCdat[,..cols.keep]
        NFCdat.allN[,sumN := rowSums(.SD,na.rm=T), .SDcols = (paste0(FRAME.hab,"_",c("NHx_dry","NHx_wet","NOy_dry","NOy_wet")))]
        FRAME.r <- rasterFromXYZ(NFCdat.allN[,c(1,2,7)])
        FRAME.r.N <- FRAME.r * 14  
        
      }else{
        
        cols.keep <- c("easting(m)","northing(m)",paste0(FRAME.hab,"_",nd))
        FRAME.r <- rasterFromXYZ(NFCdat[,..cols.keep])
        FRAME.r.N <- FRAME.r * 14
      
      }
      
      
      writeRaster(FRAME.r.N, paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output/",toupper(FRAME.hab),"_NDep_",nd,"_NFC_",year,"_KgNha_1km_",run.name,".tif"), overwrite=T)
      
      
    } # end of N dep component loop
    
    
  } # end of habitat loop
  
  print(paste0(Sys.time(),": Writing complete"))
  
}