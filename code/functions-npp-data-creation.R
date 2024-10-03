
em.phase_min_max <- function(glx) {
  glx@other$ind.metrics[, c("trip", "period")] %>% 
    mutate(loco = ymd(trip)) %>% 
    group_by(period) %>% 
    # find min and max
    summarise(min = min(loco)-months(2),
              max = max(loco)) 
}


em.reproject_layer <- function(layer){
  newproj <- '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
  r <- layer
  layer1 <- projectRaster(r, crs=newproj)
  
  return(layer1)
}


em.no_negatives <- function(layer){ 
  values(layer)[values(layer)<0] <- NA  
  return(layer)
}



em.landscape_dates <- function(gl){
  
  rain <- em.rain_caps_phase() %>% 
    mutate(trip = factor(trip))
  gl@other$ind.metrics <- gl@other$ind.metric %>% 
    left_join(rain[, c('trip', 'phase', 'period')])
  pop(gl) <- gl@other$ind.metrics$period
  glphase <- seppop(gl) 
  
  prange <-lapply(glphase, em.phase_min_max)
  
  n <-  lapply(glphase, function(x) table(factor(x@other$ind.metrics$trip)))
  sample_dates  <- reduce(lapply(n, function(x) data.frame(x)),rbind) %>% 
    mutate(sampled = rep(names(n), sapply(n, length)),
           trip = ymd(Var1)) %>% 
    dplyr::select(-Var1)
  
  sample_range <- lapply(prange,
                         function(x) data.frame(trip = seq(x$min,x$max, 
                                                           by = "months"),
                                                name = paste0("X",
                                                              gsub("-", "\\.",
                                                                   seq(x$min,
                                                                       x$max,
                                                              by = "months"))),
                                                period = x$period)) %>% 
    reduce(rbind) %>% 
    left_join(sample_dates)
  
  return(sample_range) 
}

em.landscape_load <- function(phaseSeq){
  tifFiles <- list.files("../data/monthly_npp_simpson/", "tif$", 
                         full.names = TRUE)
  
  ifirst <- grep(min(phaseSeq$trip), tifFiles)
  ilast <- grep(max(phaseSeq$trip), tifFiles)
  
  
  
  stackTifs <- tifFiles[ifirst:ilast]
  monNPP <- stack(stackTifs)
  
  stackdates <- stringr::str_extract_all(stackTifs, "\\d+")
  layerNames <- lapply(stackdates, function(x) paste(x, collapse = "-"))
  names(monNPP) <- do.call("c", layerNames)
  
  sepIndex <- lapply(split(phaseSeq, factor(phaseSeq$period)), function(x) x$name)
  
  NPPphase <- lapply(sepIndex, function(x) monNPP[[x]])
  
  return(NPPphase)
}


em.npp_layers <- function(ph){
  
  layer_sequence <- em.landscape_dates(ph)
  npplayers <- em.landscape_load(layer_sequence)

  npp <- lapply(npplayers, function(x) calc(x, fun = mean, na.rm = T))
  
  
  layersagg3 <- lapply(npp, function(x) aggregate(x, 3))
  layersR <- lapply(layersagg3, em.reproject_layer)
  layersAnalyse <- lapply(layersR, em.no_negatives)
  
  
  b <- brick(layersAnalyse)
  
  bf <- writeRaster(b, './output/npp_layers.tif',
              options="INTERLEAVE=BAND",
              overwrite=TRUE)
  
  nppMean <- data.frame(period = names(layersAnalyse),
                        npp=cellStats(stack(layersAnalyse), "mean"),
                        row.names = NULL)
  write.csv(nppMean, './output/npp_means.csv', row.names = FALSE)    
  
  return(nppMean)
  
}
