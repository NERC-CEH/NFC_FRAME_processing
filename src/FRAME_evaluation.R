
evaluate.frame.output <- function(){

################################################################################
####   Evaluate the FRAME output data using measurement data from UKEAP     ####
####   https://uk-air.defra.gov.uk/networks/network-info?view=ukeap         ####
################################################################################

  print(paste0(Sys.time(),": Evaluating the FRAME-UK NFC outputs for ",year,"..."))
  uk <- readOGR("./../country_shps","UK")
  # list of the compoumds to evaluate
  compounds <- c("NH3","SO2","NO2","HNO3","SO4","NO3","NH4","SO4","NO3","NH4") # air/gas  = NH3, SO2, NO2, HNO3. wetdep = SO4, NO3, NH4. Aerosol = SO4, NO3, NH4
  result.types <- c("gas","gas","gas","gas","precipitation","precipitation","precipitation","aerosol","aerosol","aerosol") # 'precipitation' = wet dep or precip conc, 'gas' = air concentration, 'aerosol' = particulate aerosol
  
  conc.dat <- fread(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raw_output/conc_",run.name,"_9-15-0_",year,".dat"))
  prconc.dat <- fread(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raw_output/prcon_",run.name,"_9-15-0_",year,"unc.dat"))
  
  
  stat.list <- list()
  
  # loop through all compounds
  for(i in 1:length(compounds)){

  compound <- compounds[i]
  result.type <- result.types[i]
  
  print(paste0(Sys.time(),":      ",compound,"..."))
  
  
  #########################################
  #### FRAME data: 
  # choose the required data
  #file.prefix <- ifelse(result.type == "precipitation", "prcon", "conc")
  #frame.file <- list.files(paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raw_output"), pattern = paste0("^",file.prefix), full.names = T)
  
  # read the file and subset
  
  if(result.type == "precipitation"){
    frame.dat <- copy(prconc.dat)
    cols.keep <- c("easting(m)","northing(m)",names(frame.dat)[grep(glob2rx(paste0("grd_",substr(compound,1,2),"*_wet")),names(frame.dat))])
  }else{
    frame.dat <- copy(conc.dat)
    cols.keep <- c("ie_OS","in_OS",paste0("air_",compound))
  }
  
  frame.sub <- frame.dat[,..cols.keep]
  
  # rasterize the FRAME data
  r <- rasterFromXYZ(frame.sub)
  
  #########################################
  #### Measurement data for evaluation:
  # read and subset to correct measurements
  measure <- fread(paste0("./../UKEAP_measurements/measurements_UKEAP_",year,".csv"))
  measure.com <- measure[Compound == compound & Form == result.type]
  
  # convert coordinates to metres
  measure.com[, c("E", "N") := list(Easting.km * 1000, Northing.km * 1000)]
  
  # promote to points
  coordinates(measure.com) <- ~E+N
  crs(measure.com) <- BNG
  
  #########################################
  #### Extract and compare model and measurements
  # extract model values by point source measurements
  extract.data <- extract(r, measure.com, sp=T)
  #writeOGR(extract.data, ".", "NH4_points", driver = "ESRI Shapefile")
  evaluation <- data.table(as.data.frame(extract.data))
  
  # some stats
  evaluation[, mod.min.mes := get(cols.keep[3]) - Value]
  evaluation[, mes.min.mod := Value - get(cols.keep[3])]
  evaluation[, mod.over.min := get(cols.keep[3]) / Value]
  
  model.bias <- mean(evaluation[,mod.min.mes])
  
  # normalised mean bias is just about identical to (1 - 0.5FB) / (1 + 0.5FB)
  fractional.bias <- ((mean(evaluation[,Value]) - mean(evaluation[,get(cols.keep[3])]))/(mean(evaluation[,Value]) + mean(evaluation[,get(cols.keep[3])]))) * 2
  
  norm.mean.bias <- (1 - 0.5*fractional.bias) / (1 + 0.5*fractional.bias)
  
  rmse <- sqrt( mean(evaluation[,mod.min.mes]^2) )
  
  nmse <- 1/nrow(evaluation) * (sum ( (evaluation[, mes.min.mod]^2) / ( mean(evaluation[, Value]) * mean(evaluation[, get(cols.keep[3])]) ) ))
  
  fac2 <- sum(evaluation[, mod.over.min] >= 0.5 & evaluation[, mod.over.min] <= 2) / nrow(evaluation)
  
  # calculate median bias for correction value (NH3 gas only)
  if(compound == "NH3" & result.type == "gas"){
  
  evaluation[, bias := (get(cols.keep[3]) - Value) / Value ]
  med.bias <- median(evaluation[, bias])
  calib.fac <- 1 / (1 + med.bias)
  evaluation[, paste0(cols.keep[3],"_adj") := get(cols.keep[3]) * calib.fac]
  
  #write.csv(calib.fac,paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/NH3_gas_calib_",year,"22.csv"),row.names = F)
  
  
  }else{print("Median bias correction not applicable")}
  
  m <- lm(Value ~ get(cols.keep[3]), evaluation)
  r2 <- format(summary(m)$r.squared, digits = 2)
  
  # stats data frame
  df <- data.frame(Statistic = c("Points (n)","R2","Factor of 2","Normalised Mean Bias","Normalised Mean Square Error"), Value = c(nrow(evaluation), r2, round(fac2,2),round(norm.mean.bias,2),round(nmse,2)))
  df.all <- data.frame(Variable = rep(paste0(compound,".",result.type),5),Statistic = c("Points (n)","R2","Factor of 2","Normalised Mean Bias","Normalised Mean Square Error"), Value = c(nrow(evaluation), r2, round(fac2,2),round(norm.mean.bias,2),round(nmse,2)))
  
  stat.list[[i]] <- df.all
  
  #########################################
  ### write out a raster surface of the output data
  
  units <- ifelse(result.type == "precipitation","ueql","ugm3")
  
  if(compound == "NH3" & result.type == "gas"){
    rc <- r * calib.fac
  }else{
    rc <- r
  }
  
  if(compound == "NH3" & result.type == "gas"){
    writeRaster(r, paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output/NH3_gas_conc_non.calib_NFC_",year,"_",units,"_1km_",run.name,".tif"), overwrite=T)
    writeRaster(rc, paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output/NH3_gas_conc_calib_NFC_",year,"_",units,"_1km_",run.name,".tif"), overwrite=T)
  }else{
    writeRaster(rc, paste0("./NAEI_NFC_",year,"_forFRAME/FRAME_UK/raster_output/",compound,"_",result.type,"_conc_NFC_",year,"_",units,"_1km_",run.name,".tif"), overwrite=T)
  }
  
  
  #########################################
  ### plot files
  save.image.filename <- paste0(compound,"_",result.type,"_eval_NFC_",year,"_",run.name,".png")
  
  # determine if model or measure has bigger range (and therefore 600 pixels dimension)
  img.height <- 600 / ((max(evaluation[,get(cols.keep[3])]) - min(evaluation[,get(cols.keep[3])])) / (max(evaluation[,Value]) - min(evaluation[,Value])))
  img.height <- ifelse(img.height < 250, 300, img.height)
  
  png(file = paste0("./NAEI_NFC_",year,"_forFRAME/evaluation/",save.image.filename),width=600,height=img.height + 400)
  
  # table grob theme
  mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 1.4)),
  colhead = list(fg_params=list(cex = 1.6)),
  rowhead = list(fg_params=list(cex = 1.6)))
  
  # write an expression dependent on compound and type
  # air/gas  = NH3, SO2, NO2, HNO3. wetdep = SO4, NO3, NH4. Aerosol = SO4, NO3, NH4
  if(result.type == "precipitation"){
  x.axis.name <- bquote("Modelled"~.(substr(compound,1,(nchar(compound) - 1)))[.(substr(compound,nchar(compound),nchar(compound)))]~"("*"in"~.(result.type)*"):"~micro~eq~l^-1)
  y.axis.name <- bquote("Measured"~.(substr(compound,1,(nchar(compound) - 1)))[.(substr(compound,nchar(compound),nchar(compound)))]~"("*"in"~.(result.type)*"):"~micro~eq~l^-1)
  }else{
  x.axis.name <- bquote("Modelled"~.(substr(compound,1,(nchar(compound) - 1)))[.(substr(compound,nchar(compound),nchar(compound)))]~"("*"in"~.(result.type)*"):"~mu*g~m^-3)
  y.axis.name <- bquote("Measured"~.(substr(compound,1,(nchar(compound) - 1)))[.(substr(compound,nchar(compound),nchar(compound)))]~"("*"in"~.(result.type)*"):"~mu*g~m^-3)
  }

  x.limits <- ifelse(max(evaluation[,get(cols.keep[3])]) < 10 ,ceiling(max(evaluation[,get(cols.keep[3])])),ceiling(max(evaluation[,get(cols.keep[3])])/10)*10 )
  y.limits <- ifelse(max(evaluation[,Value]) < 10, ceiling(max(evaluation[,Value])), ceiling(max(evaluation[,Value])/10)*10 )
  
  # plot
  g1 <- ggplot()+
  #annotation_custom(tableGrob(df, theme = mytheme, rows = NULL), xmin=max(evaluation[,get(cols.keep[3])])*0.1, xmax=max(evaluation[,get(cols.keep[3])])*0.3, ymin=max(evaluation[,Value])*0.75, ymax=max(evaluation[,Value]))+
  geom_point(data=evaluation, aes_string(x = paste0(cols.keep[3]), y = "Value"))+
  #geom_point(data=evaluation, aes_string(x = paste0("air_",compound,"_adj"), y = "Value"),colour="red")+
  geom_abline(slope=1, intercept=0)+
  geom_abline(slope=2, intercept=0,linetype=2)+
  geom_abline(slope=0.5, intercept=0,linetype=2)+
  coord_fixed()+
  expand_limits(x = 0, y = 0)+
  scale_x_continuous(expand = c(0,0),limits=c(0,x.limits),name=x.axis.name)+
  scale_y_continuous(expand = c(0,0),limits=c(0,y.limits),name=y.axis.name)+
  theme(panel.background = element_blank(),
        axis.text = element_text(size=14,colour="black"),
        axis.title = element_text(size=16,colour="black"),
        axis.line = element_line(size=1,colour = "black"),
        axis.ticks = element_line(size=1,colour = "black"),
        plot.margin=unit(c(-1,1,-2,1), "cm"))

  grid.arrange(
  g1,
  tableGrob(df, theme = mytheme, rows = NULL),
  nrow = 2,
  heights = c(4/5, 1/5),
  clip = FALSE
  )
  
  dev.off()
  
  #########################################
  ### fairly quick, non-classified plot of results ###
  
  
  
  ## maps
  rc.c <- mask(rc, uk, inverse=F)
  
  rc.c[rc.c==0] <- NA

  # set alpha, pstrip and axes font sizes. 
  alphs <- c(0.9)
  
  
  if(result.type == "precipitation"){
    plot.units <- parse(text = "mu~eq.~l^-1")[[1]]
  }else{
    plot.units <- parse(text = "mu*g~m^-3")[[1]]
  }

  plot.title <- list(bquote(.(substr(compound,1,nchar(compound)-1))[.(substr(compound,nchar(compound),nchar(compound)))]~.(result.type)~conc~.(year)~(.(plot.units))), cex=3)
  
  # filename
  save.image.filename <- paste0(compound,"_",result.type,"_conc_NFC_",year,"_",run.name,".png")
  
  # png plot
  png(file = paste0("./NAEI_NFC_",year,"_forFRAME/plots/",save.image.filename),width=800,height=1000)
  
  pj <- print(levelplot( rc.c, margin=F, alpha.regions=alphs, main = plot.title, scales=list(draw=FALSE),colorkey=list(labels=list(cex=2.5),height=0.7, width=2))  + layer(sp.polygons(uk,lwd=0.1)))
  
  grid.text(label= paste0("Max = ",round(cellStats(rc.c,max),2)), x = 0.25, y = 0.9, just = "left",gp = gpar(cex=2.2))
  
  dev.off()
  
  # end of map
  
  } # end of compound loop

## write out all of the stats summary
nfut.stats <- do.call("rbind", stat.list)
nfut.stats.w <- dcast(nfut.stats, formula = Statistic ~ Variable)

write.csv(nfut.stats.w,paste0("//nercbuctdb.ad.nerc.ac.uk/projects1/NEC07358_NitrogenFutures/WP1/FRAME_modelling/FRAME_outputs/UK/evaluation/all_stats_NF_",scenario,substr(year,3,4),"_",frame.version,".csv"))

#

print(paste0(Sys.time(),":      COMPLETE."))


} # end of function