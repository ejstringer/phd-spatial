

em.to_mat <- function(lcf, matname = 'fst'){
  rowsAcols <- sort(unique(c(lcf$pop1,lcf$pop2)))
  
  mat <-matrix(data = NA, nrow = length(rowsAcols),
               ncol = length(rowsAcols))
  rownames(mat) <- rowsAcols
  colnames(mat) <- rowsAcols
  
  take2 <- lcf[,c('pop1','pop2', matname)]
  names(take2) <- c('pop2', 'pop1', matname)
  lcd <- bind_rows(lcf, take2)
  for(i in 1:nrow(mat)) mat[i,i] <- 0
  
  for(i in 1:nrow(lcf)){
    a <- which(rownames(mat)== lcf$pop1[i])
    b <- which(rownames(mat)== lcf$pop2[i])
    mat[a,b] <- unlist(lcf[i,matname])
    mat[b,a] <- unlist(lcf[i,matname])
    
    
  }
  
  return(mat)  
}


em.ibd.period.correlog <- function(kin2, test = c("spearman", 'pearson')){ 
  
  kin2$gendist <- kin2$euclidean
  kin2$pop1 <- kin2$id.x
  kin2$pop2 <- kin2$id.y
  
  
  newMan <- list()
  
  
  for(indx in levels(kin2$phaseNo)){
    
    kin_i <- filter(kin2, phaseNo == indx)
    g <- em.to_mat(kin_i, matname = 'gendist')
    d <- em.to_mat(kin_i, matname = 'km.log')
    bks <- c(0.1,1,2.1,4.5,9.5, 20.1,42.5, 85) %>% log
    
    man <- mantel(as.dist(d), as.dist(g))  
    
    
    mancorrelog <- mantel.correlog(D.eco = as.dist(g), D.geo = as.dist(d),
                                   cutoff = FALSE, r.type = test,
                                   nperm = 10000,
                                   break.pts = bks) 
    mancorrelog
    
    # summary(mancorrelog)
    # plot(mancorrelog)
    # title(main = indx)
    matx <- data.frame(mancorrelog$mantel.res)
    matx$dist <- exp(matx$class.index);matx
    matx$breaksto <- bks[-1]
    matx$phaseNo <- indx
    matx$sig <- matx$Pr.corrected. < 0.05 & matx$Mantel.cor > 0 | matx$class.index < 0
    i <- grep(indx, levels(kin2$phaseNo))
    
    
    newMan[[i]] <- matx
    
  }
  allmans <- do.call('rbind', newMan) %>% 
    mutate(phase = str_sub(phaseNo, 1,1)) %>% # %>% filter(sig) %>% 
    mutate(phase2 = ifelse(phase=='L', 'bust', 'boom'),
           phase = factor(phase, levels = c('L',
                                            'I',
                                            'D')),
           species = kin2$species[1],
           sex_pairs = ifelse('unknown' %in% kin2$sex_pairs,'all',
                              kin2$sex_pairs[1])) %>% 
    left_join(kin2[,c('phaseNo', 'npp')] %>%  unique())
  
  return(allmans)
}