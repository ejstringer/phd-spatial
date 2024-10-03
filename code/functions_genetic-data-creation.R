# functions
em.filter.mac <- function(gl, threshold = 4){
  
  t <- threshold + 1
  gen <- as.matrix(gl)
  thresh <- colSums(!is.na(gen))*2 - t
  nall <- colSums(gen, na.rm = T)
  lowmac <- nall[nall < t | nall > thresh]

  gl2 <- gl.drop.loc(gl, loc.list = names(lowmac), verbose = 0)
  return(gl2)
}

# filtering
em.remove.recaps <- function(pherm, threshold = 0.97){
  
  pherm@other$ind.metrics$id <- paste0('x', pherm@other$ind.metrics$id)
  pherm@ind.names <- paste0('x', pherm@ind.names)
  pherm@other$ind.metrics$id == pherm@ind.names
  
  pop(pherm) <- pherm@other$ind.metrics$species
  phpop <- seppop(pherm)
  
  
  # find related inidividuals
  # setup 
  removeList <- list()
  triprelated <- NA
  
  
  r <- 1 #remove list
  
  for(y in 1:length(phpop)){
    # choose trip
    tripx <- phpop[[y]]
    related <- FALSE
    # relativeness
    if(nInd(tripx)>1){ 
      ibdecent <- gl.grm(tripx,plotheatmap = F, verbose = 0) 
      related <- sum(as.dist(ibdecent) >= threshold) != 0
      
    }
    # if related add to remove list
    if(related){
      df <- data.frame(ibdecent)
      for(k in 1:nrow(df)) df[k,k] <- NA
      df[(df) < threshold] <- NA
      df[1:6, 1:5]
      
      
      indRel <- apply(df, MARGIN = 1, function(x) sum(x, na.rm = T))
      
      sibs <- names(indRel[indRel>0])
      
      removeInd <- NA
      x = 1
      df2 <- df[sibs, sibs]
      
      while(sum(df2 >=threshold, na.rm = T) >0){
        
        nsibs <- rowSums(!is.na(as.matrix(df2)))
        topsib <-  which.max(nsibs)
        
        #print(x)
        removeInd[x] <- names(topsib)
        x = x+1
        df2 <- df2[-topsib, -topsib]
        
      }
      
      removeList[[r]] <- removeInd
      triprelated[r] <- names(phpop)[y]
      
      r = r+1
    }
    
  }
  pherm2 <- pherm
  nRemove <- 0
  if(length(removeList)>0){
    
    names(removeList) <- triprelated
    removeALL <- do.call('c', removeList)
    pherm2 <- gl.drop.ind(pherm, ind.list = removeALL, verbose = 0)
    nRemove <- length(removeALL)
  }
  
  pherm2@other$ind.metrics$id <- sub('x','', pherm2@other$ind.metrics$id)
  pherm2@ind.names<- sub('x','', pherm2@ind.names)
  
  cat("number of recaptures removed:", nRemove, "\n", 
      nInd(pherm2), "individuals remaining", '\n')
  return(pherm2)
  
}

em.filtering <- function(species = 'ph') {
  if(species == 'ph') glO <- readRDS("../data/pherm_genotypes.rds")
  if(species == 'sy') glO <- readRDS("../data/syoung_genotypes.rds")
  
  print(ifelse(species == 'ph',
               "loading pherm_genotypes from data folder...",
               "loading syoung_genotypes from data folder..."))
  
  glx <- glO
  glx <- gl.filter.callrate(glx, method = "ind", threshold = 0.5, verbose = 0)
  glx <- gl.filter.rdepth(glx, verbose =0, lower = 5) 
  glx <- gl.filter.reproducibility(glx,threshold = 0.95, verbose = 0)
  glx <- gl.filter.monomorphs(glx, verbose = 0)
  
  glx <- gl.filter.callrate(glx, method = "loc", threshold = 0.95, verbose = 0)
  glx <- gl.filter.callrate(glx, method = "ind", threshold = 0.95, verbose = 0)
  
  
  glx <- gl.filter.secondaries(glx, verbose = 0, method = 'best')
  
  glx <- em.filter.mac(glx, threshold = 4)
  
  glx <- em.remove.recaps(glx)
  
  cat("\n", "pre filtering loc:", nLoc(glO), "and ind:", 
      nInd(glO), "\n",
      "post filtering loc:", nLoc(glx), "and ind:", 
      nInd(glx), "\n")
  
  filterfile <- paste0("./output/", glx@other$ind.metrics$species[1], 
                       ' filtered genotypes.rds')
  saveRDS(glx, gsub(' ', '_', filterfile))
  
  return(glx)
}


em.rain_caps_phase <- function(){
  
  # phase dates
  impactful_rain <- ymd(c("2007-01-01", "2010-02-01", "2015-12-01"))
  
  increase_dates <- impactful_rain + months(2)
  
  decrease_dates <- ymd(c("2008-03-01", "2011-09-01", "2016-11-01")) 
  
  low_dates <- c(ymd("2000-12-01"), decrease_dates[-3] + c(months(13), months(9)))
  
  # session info
  session <- read.csv("../data/em_sessiontable.csv") %>% 
    group_by(trip) %>%
    summarise(rain = mean(rain),
              captures = mean(phAbundance, na.rm = T),
              capturesSy = mean(syAbundance, na.rm = T)) %>%
    ungroup() %>% 
    mutate(year = lubridate::year(lubridate::ymd(trip)),
           trip = ymd(trip)) %>% 
    bind_rows(data.frame(trip = max(.$trip)+months(1:12)))
  
  ## phase data
  phase <- data.frame(low =  low_dates, 
                      decrease = ymd(decrease_dates), 
                      increase = increase_dates) %>% 
    pivot_longer(cols = c(low, decrease, increase)) %>% 
    filter(complete.cases(.)) %>% arrange(value) %>% 
    mutate(no = rep(1:3, each = 3),
           period = paste('period', 1:9),
           phaseNo = paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3)))
  
  session$phase <- NA
  session$period <- NA
  session$phaseNo <- NA
  for(i in 1:nrow(phase)){
    index <- which(session$trip >= phase$value[i])
    session$phase[index] <- phase$name[i]
    session$period[index] <- phase$period[i]
    session$phaseNo[index] <- phase$phaseNo[i]
    
  }
  
  return(session)
}

em.distance_grids <- function (dfcoor, between = "gridId", proj = "+proj=longlat", 
                               transform = "+init=epsg:32754") {
  dfcoor <- filter(dfcoor, coor == 'grid')
  cord.dec = sp::SpatialPoints(cbind(dfcoor$lon, dfcoor$lat), 
                               proj4string = sp::CRS(proj))
  cord.UTM <- sp::spTransform(cord.dec, sp::CRS(transform))
  mat <- as.matrix(dist(cord.UTM@coords))
  matNames <- dfcoor[, between]
  rownames(mat) <- matNames
  colnames(mat) <- matNames
  distmat <- as.dist(mat)
  dfdist <- otuSummary::matrixConvert(distmat, colname = c("pop1", 
                                                           "pop2", "metres"))
  dfdist_pairs <- dfdist %>%
    rowwise() %>% 
    mutate(pairs = paste(sort(c(pop1, pop2)), collapse = "-"),
           km = metres/1000)
  
  return(dfdist_pairs)
}


em.individual_relatedness <- function(ph){
  kinship <- popkin(t(as.matrix(ph)))
  grm <- gl.grm(ph, plotheatmap = FALSE, verbose = 0)
  euc <- dist(as.matrix(ph)-1)
  
  
  phkin <- matrixConvert(kinship, colname = c("id.x", "id.y", "kinship")) %>% 
    mutate(id.x = as.character(id.x),
           id.y = as.character(id.y))
  
  
  phibd <- matrixConvert(grm, colname = c("id.x", "id.y", "ibdecent")) %>% 
    mutate(id.x = as.character(id.x),
           id.y = as.character(id.y))
  
  pheuc <- matrixConvert(euc, colname = c("id.x", "id.y", "euclidean")) %>% 
    mutate(id.x = as.character(id.x),
           id.y = as.character(id.y))
  
  rain <- em.rain_caps_phase() %>% 
    mutate(trip = factor(trip))
  meta <- ph@other$ind.metric[,c('id','pop','lat', 'lon', 'gridId', 'trip',
                                 'species', 'sex', 'mass', 'age')] %>% 
    left_join(rain[, c('trip', 'phase', 'period')])
  
  grids <- read.csv("../data/em_gridcoortable.csv") %>%
    filter(coor == 'grid')
  eucGrid <- em.distance_grids(grids)
  
  individual_metrics <- phkin  %>% 
    left_join(phibd) %>% 
    left_join(pheuc) %>% 
    left_join(meta, by = c('id.x' = 'id')) %>% 
    left_join(meta, by = c('id.y' = 'id'))  
  
  individual_df <- individual_metrics %>% 
    rowwise() %>% 
    mutate(pairs = paste(sort(c(as.character(gridId.x),
                                as.character(gridId.y))), collapse = "-"),
           idpairs = paste(sort(c(as.character(id.x),
                                  as.character(id.y))), collapse = "-"),
           mc = ifelse(grepl('MC', gridId.x) & grepl('MC', gridId.y),
                       TRUE, FALSE),
           within = (period.x == period.y)) %>%
    filter(!duplicated(idpairs)) %>% 
    left_join(eucGrid[,c('pairs', 'metres')], by = 'pairs') %>% 
    mutate(species = species.x,
           metres = ifelse(gridId.x == gridId.y & gridId.x %in% grids$gridId, 0, metres),
           km = metres/1000,
           metres.log = log(metres + 1)) %>%
    dplyr::select(-c(species.x, species.y))
  
  relatednessfile <- paste0("./output/", individual_df$species[1], 
                            ' individual relatedness.csv')
  write.csv(individual_df, gsub(' ', '_', relatednessfile), row.names = F)
  
  return(individual_df)
}

em.gl_assign_pop <- function (glx, define.pop = "gridId") 
{
  if (length(define.pop) == 1) {
    if (!define.pop %in% names(glx@other$ind.metrics)) {
      stop(paste(define.pop, "not in ind.metrics; names(glx@other$ind.metrics)"))
    }
    pop(glx) <- glx@other$ind.metrics[, define.pop]
  }
  else {
    pop(glx) <- define.pop
  }
  return(glx)
}


em.gl_min <- function (glx, define.pop = "gridId", s.size = 4, equal = FALSE) 
{
  my.fun <- function(x, equal) {
    if (length(x@ind.names) >= s.size) {
      aa <- x@ind.names
      if (equal) 
        aa <- sample(x$ind.names, s.size, replace = FALSE)
    }
    else {
      aa <- NULL
    }
    return(aa)
  }
  glx <- em.gl_assign_pop(glx, define.pop)
  glsep <- seppop(glx)
  indexNames <- unlist(lapply(glsep, function(x) my.fun(x, 
                                                        equal)))
  glSample <- NULL
  if (length(indexNames) > 0) 
    glSample <- gl.keep.ind(glx, ind.list = indexNames)
  return(glSample)
}



em.phase_fst <- function(ph){
  
    rain <- em.rain_caps_phase() %>% 
      mutate(trip = factor(trip))
    ph@other$ind.metrics <- ph@other$ind.metric %>% 
      left_join(rain[, c('trip', 'phase', 'period')])
    
    ph <- em.remove.recaps(ph, threshold = 0.5)
    pop(ph) <- ph@other$ind.metrics$period
    phasePh <- seppop(ph)
    phasePh <- lapply(phasePh, em.gl_assign_pop)
    phasePh <- lapply(phasePh, em.gl_min)
    phasePh <- phasePh[!sapply(phasePh, is.null)]
    
    fst <- pbapply::pblapply(phasePh, function(x) gl.fst.pop(x, nboots = 1))
    
    fun.x <- function(x, gdist = fst) {
      ff <- otuSummary::matrixConvert(as.dist(gdist[[x]]), c("pop1", "pop2", "fst"))
      ff$i <- x
     # ff$phase <- sub("[[:digit:]]+", "", names(fst)[x])
      ff$period <- factor(names(fst), levels = names(fst))[x]
      return(ff)
    }
    
    fstdf <- lapply(1:length(phasePh), function(x) fun.x(x, gdist = fst))
    fstdf2 <- do.call("rbind", fstdf) %>% 
      left_join(unique(ph@other$ind.metrics[,c('phase', 'period')]))
    
    fst2 <-fstdf2 %>%
      rowwise() %>% 
      mutate(pairs = paste(sort(c(as.character(pop1),
                                  as.character(pop2))), collapse = "-"))
    grids <- read.csv("../data/em_gridcoortable.csv")
    eucGrid <- em.distance_grids(grids)
    
    fst2$pop1 %in% grids$gridId
    fst2$pairs %in% eucGrid$pairs
    fsteuc <- left_join(fst2, eucGrid[c("pairs","metres", "km")], by = "pairs")
    fsteuc$species <- phasePh[[1]]@other$ind.metrics$species[1]
    
    fstfile <- paste0("./output/", fsteuc$species[1], ' fst.csv')
    write.csv(fsteuc, gsub(' ', '_', fstfile),
              row.names = F)
    
    

  return(fsteuc)
}


em.ne_prep <- function(gl){
  
  # glsib <- em.remove.perc.relatives(gl, threshold = 0.25,
  #                                   removePercent = 0.75)
  
  gl
  rain <- em.rain_caps_phase() %>% 
    mutate(trip = factor(trip))
  
  gl@other$ind.metrics <- gl@other$ind.metrics %>% left_join(rain)
  
  #gl <- em.remove.recaps(gl, threshold = 0.5)
  pop(gl) <- gl@other$ind.metrics$period
  glpop <- seppop(gl)
  sapply(glpop, nInd)
  popMaf <- list() 
  for(i in 1:length(glpop)){
    glx <- glpop[[i]]
    
    glxx <- em.filter.mac(glx, threshold = 1) #remove singletons
    #glxx <- gl.filter.callrate(glxx, threshold = 0.99, verbose = 0)
    popMaf[[i]] <- glxx
    
  }
  names(popMaf) <- gsub(' ', '_', names(glpop))
  
  sapply(1:9, function(x) table(colSums(as.matrix(popMaf[[x]])))[1:5])
  # remove snps that are only heterozygotes
  gg <- popMaf[[1]] 
  
  
  em.he.removal <- function(gg){
    gmat <- as.matrix(gg)
    allhet <- colSums(gmat == 1, na.rm = T)/(nrow(gmat)-is.na(gmat)) == 1
    nhet <- sum(allhet)
    cat(nhet)
  } # none found so did not expand function to removing loci

  sapply(popMaf, em.he.removal)
  sapply(popMaf, nInd)
  sapply(popMaf, nLoc)
  
  neprepfile <- paste0('./output/', popMaf[[1]]@other$ind.metrics$species[1],
                      ' ne prep FOR dungog.rds')
  saveRDS(popMaf, gsub(' ', '_', neprepfile))
  return(popMaf)
}

em.ne_summarise <- function(){
  
  
  fx <- function(x, spp = 'ph') {
    n <- names(x)
    n <- sub('_P_herm', '', n)
    d <- x[[1]]
    d$phaseNo <- n
    d$phase <- gsub('[[:digit:]]+', '', n)
    d$species <- spp
    return(d)
  }
  
  
  neEstim
  c(lapply(neEstimate, fx),
    lapply(neEstimateSy, fx, spp = 'sy')) %>% 
    do.call('rbind', .) %>% 
    dplyr::select(-`Frequency 1`) %>% 
    pivot_wider(names_from = Statistic, 
                values_from = `Frequency 2`) -> ne_df_NOsingletons
  
  c(lapply(neEstimate, fx),
    lapply(neEstimateSy, fx, spp = 'sy')) %>% 
    do.call('rbind', .) %>% 
    dplyr::select(-`Frequency 2`) %>% 
    pivot_wider(names_from = Statistic, 
                values_from = `Frequency 1`) -> ne_df_singletons
  
  
  
  syn <- sampleSize$sy %>% sapply(nInd) %>% 
    data.frame(n = ., phaseNo = names(.), species = 'Pseudomys hermannsburgensis')
  phn <- sampleSize$ph %>% sapply(nInd) %>%
    data.frame(n = ., phaseNo = names(.), species = 'Sminthopsis youngsoni')
  
  nesummarised <- ne %>% 
    mutate(species = ifelse(species == 'ph',
                            'Pseudomys hermannsburgensis',
                            'Sminthopsis youngsoni')) %>% 
    dplyr::select(phaseNo, phase, species,
                  Estimated.Ne., CI.low.Parametric, CI.high.Parametric) %>% 
    left_join(rbind(syn,phn)) %>% 
    rename(ne = Estimated.Ne., 
           lower = CI.low.Parametric,
           upper = CI.high.Parametric)
  
  
  
}