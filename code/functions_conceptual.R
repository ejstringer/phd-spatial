
em.visualisation_subset <- function(rain_caps, ph){
  
  max_y <- max(rain_caps$captures, na.rm = T)
  
  subset_dates <- rain_caps %>% 
    filter(trip >= min(ymd(ph@other$ind.metrics$trip))) %>%
    mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                     rep(1:3, each = 3)))) %>% 
    group_by(phase, phaseNo) %>% 
    summarise(trip = min(trip)) %>% 
    filter(complete.cases(phase)) %>% 
    arrange(trip) %>% 
    rename(subset = phaseNo)
  
  
  vis_dates <- c(subset_dates$trip, 
                 max(ymd(ph@other$ind.metrics$trip))+months(2))
  
  vis_comb <- combn(1:length(vis_dates),2) 
  vis_comb_seq <-t(vis_comb[,vis_comb %>% diff == 1])
  
  subset_data_visualise <- data.frame(subset = subset_dates$subset,
                                      trip1 = ymd(vis_dates[vis_comb_seq[,1]]),
                                      trip2 = ymd(vis_dates[vis_comb_seq[,2]])) %>% 
    mutate(trip3 = ymd(trip1) - days(1),
           trip4 = ymd(trip2) - days(1)) %>% 
    pivot_longer(cols = trip1:trip4, values_to = "trip") %>% 
    arrange(subset, trip) %>% 
    mutate(rain = rep(c(0, rep(max_y, 2), 0), nrow(subset_dates))) %>% 
    dplyr::select(-name) %>% 
    mutate(phase = sub("[[:digit:]]+", "", subset),
           phase = factor(ifelse(phase == 'L', 'low', ifelse(phase == 'I',
                                                             'increase', 
                                                             'decrease')),
                          levels = c('low','increase','decrease'))) %>% 
    relocate(phase)
  
  return(list(subset = subset_data_visualise,
              dates = subset_dates))
  
}


em.conceptual_phase <- function(startdate = "2004-01-01",
                                rep123 = c(0.07, 0.402,0.615)){
  
  ph <- readRDS('../data/pherm_genotypes.rds')
  rain_caps <- em.rain_caps_phase()
  rain_caps$phase <- factor(rain_caps$phase, 
                            levels = c('low','increase', 'decrease'))
  subphase <- em.visualisation_subset(rain_caps, ph)
  
  subset_data_visualise <- subphase$subset 
  start_dates <- subphase$dates 
  
  myfillsubtle <- c(low = "#FFF3CD", increase = "#D2E7FA",
                    decrease = "#F7D1DF")
  
  myfill <- rev(c(increase = '#1E88E5',
                  decrease = '#D81B60', 
                  low = '#FFC107'))
  
  myfill <- c(decrease = terrain.colors(10)[4],
                increase = terrain.colors(10)[1],
                low = terrain.colors(10)[7],
              stable = 'grey40')
  
  

  grob_repeats <- grobTree(textGrob('Periods:',
                                    x=0.002,  
                                    y= 0.93,
                                    hjust=0, rot = 0,
                                    gp=gpar(col="black", 
                                            fontsize=11)))
  
  # grob_repeats <- grobTree(textGrob(paste("period"),
  #                                   x=rep123[1] +0.005,  
  #                                   y= 0.98,
  #                                   hjust=0, rot = 0,
  #                                   gp=gpar(col="black", 
  #                                           fontsize=9, 
  #                                           fontface="italic")))
  
  dd <- filter(rain_caps, 
         trip > ymd(startdate) & complete.cases(captures)) %>% 
    pivot_longer(cols = captures:capturesSy, names_to = 'species', values_to = 'captures') %>% 
    mutate(trip = ymd(trip),
           phaseNo = factor(phaseNo, levels = unique(rain_caps$phaseNo)))
  
   sdv <- subset_data_visualise %>% 
     filter(rain > 0) %>% 
     mutate(phaseNo = factor(subset, levels = unique(rain_caps$phaseNo)),
            period = phaseNo,
            species = 'captures') %>%
     group_by(phaseNo, phase, period, species) %>% 
     summarise(rain = mean(rain),
               trip2 = as.Date(mean(trip), origin = "1970-01-01") - days(125),
               trip3 = as.Date(mean(trip), origin = "1970-01-01")-days(220)) %>% 
     arrange(trip2) %>% 
     filter(!duplicated(phaseNo)) %>% 
     ungroup()  %>% 
     mutate(trip = ifelse(phase == 'increase', 
                          as.character(trip2), 
                          as.character(trip3)),
            trip = ymd(trip))
   subset_data_visualise$species <- 'captures'
   subset_data_visualise2 <- subset_data_visualise
   subset_data_visualise2$species <- 'capturesSy'
   subset_data_visualise2$rain <- ifelse(subset_data_visualise2$rain > 0,
                                         10, 0)
   
   phaseVis <- rbind(subset_data_visualise,
                     subset_data_visualise2) %>% 
     mutate(phase = as.character(phase),
            phase = ifelse(species == 'capturesSy', 'stable', phase),
            phase = factor(phase, levels = names(myfill)))
 #  sdv$period <- (sdv$phaseNo)
    #sdv2 <- rbind(sdv[1,], sdv)
   # sdv2$period[1] <- 'Period: 1'
   # sdv2$trip[2] <- ymd('2005-11-15')
   # sdv2$trip[1] <- ymd('2003-12-01')
   # sdv2$species <- dd$species[1]
   # sdv2 <- sdv2[-2,]
  conceptual <- ggplot(dd,
                       aes(ymd(trip), captures))+
    geom_area(data = phaseVis, 
              aes(x = trip, y = rain*1.07, group = subset, fill = phase), 
              alpha = 0.60) +
    geom_area(alpha = 1, fill = "grey90", col = "grey80")+
    geom_bar(stat = "identity", alpha = 0.8, col="grey20", fill = "black")+
    scale_fill_manual(values = myfill, 
                      breaks = c('low', 'increase', 'decrease', 'stable'),
                      name = 'population phase') +
    #facet_wrap(~species, nrow = 2, )+
    facet_grid(rows = vars(species),  scale = 'free_y', space = 'free_y')+
    # geom_text(data = sdv, aes(x = trip, y = rain + 5, 
    #                            label = period), hjust = -0.5, size = 4.5)+
    theme_classic()+
    xlab("Year") +
    ylab("Captures (per 100 trap nights)")+
    geom_vline(xintercept = start_dates$trip,#[seq(1,9,3)], 
               lty = 2, col = "black", alpha = 0.4, lwd = 0.5)+
   # annotation_custom(grob_repeats)+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
         # legend.position = c(0.11, 0.7),
          legend.position = 'bottom',
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
    scale_x_date(date_breaks = "years" , date_labels = "%Y")
  return(conceptual)
}


