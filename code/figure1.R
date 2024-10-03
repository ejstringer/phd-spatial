source('./code/libraries.R')
source('./code/functions_conceptual.R')
source('./code/functions_genetic-data-creation.R')

# pictures ------------------
speciesName <- grobTree(textGrob('P. hermannsburgensis',
                                 x=0.02,  
                                 y= 0.12,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=9, fontface="italic")))
mammalName <- grobTree(textGrob('rodent',
                                x=0.05+0.25,  
                                y= 0.220,
                                hjust=0, rot = 0,     #was fontsize = 12
                                gp=gpar(col="black", fontsize=9, fontface = 'bold')))


speciesName2 <- grobTree(textGrob('S. youngsoni',
                                  x=0.2,  
                                  y= 0.06,
                                  hjust=0, rot = 0,     #was fontsize = 12
                                  gp=gpar(col="black", fontsize=9, fontface="italic")))


mammalName2 <- grobTree(textGrob('marsupial',
                                 x=0.2+0.05,  
                                 y= 0.16,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=9, fontface = 'bold')))

imgph <- readPNG("../data/pics/pherm_close.png")

imgsy <- readPNG("../data/pics/syoung_close.png")



pic <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(imgph, interpolate=TRUE))+
  annotation_custom(speciesName)+
  annotation_custom(mammalName)
pic2 <- ggplot() + theme_void() +
  annotation_custom(rasterGrob(imgsy, interpolate=TRUE)) +
  annotation_custom(speciesName2)+
  annotation_custom(mammalName2)

# conceptual -------
myfill <- c(decrease = terrain.colors(10)[4],
            increase = terrain.colors(10)[1],
            low = terrain.colors(10)[7],
            stable = 'grey90')



con <- em.conceptual_phase('2003-11-01', c(0.24, 0.507,0.68)) 
#  theme(legend.position = 'none')

phaseCol <- c(decrease = terrain.colors(10)[4],
              increase = terrain.colors(10)[1],
              low = terrain.colors(10)[7])


# grids ------
grids <- read.csv('../data/em_gridcoortable.csv')
phxy <- rgdal::project(as.matrix(grids[,4:3]), 
                       proj = '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
gridxy <- data.frame(x = phxy[,1], y = phxy[,2],
                     gridId = grids$gridId)


gridsUsed <- unique(c(fstdata$pop1, fstdata$pop2))
gridxy 
#notph <- gridsUsed[!gridsUsed %in% ph@other$ind.metrics$gridId]

nrow(gridxy)


gridsph <- split(fstdata[grepl('herm', fstdata$species),], 
             fstdata$phaseNo[grepl('herm', fstdata$species)]) %>% 
  lapply(function(x) unique(c(x$pop1, x$pop2)))

gridssy <- split(fstdata[grepl('you', fstdata$species),], 
                 fstdata$phaseNo[grepl('you', fstdata$species)]) %>% 
  lapply(function(x) unique(c(x$pop1, x$pop2)))

gridfst <- data.frame(gridId = do.call('c', gridsph), 
           variable = rep(names(gridsph),
                           sapply(gridsph, length)),
           species = 'ph') %>% 
  bind_rows(data.frame(gridId = do.call('c', gridssy), 
                       variable = rep(names(gridssy),
                                       sapply(gridssy, length)),
                       species = 'sy')) %>% 
  filter(complete.cases(gridId)) %>% 
  left_join(gridxy) %>% 
  mutate(variable = factor(variable, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                     rep(1:3, each = 3))),
         grid = 'grid')




# npp -----------
npp <- stack('./output/npp_layers.tif') #readRDS('./output/npp_layers.rds')

names(npp) <- c("low1", "increase1", "decrease1", 
                "low2", "increase2", "decrease2",     
                "low3",  "increase3", "decrease3")

phaseCol <- c(low = terrain.colors(10)[7],
              increase = terrain.colors(10)[1],
              decrease = terrain.colors(10)[4])
phaseColmute <- scales::muted(phaseCol, l = 80, c = 70)
scales::show_col(c(phaseCol, phaseColmute))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


npp2 <- npp
names(npp2) <- paste0(rep(c('L', 'I', 'D'), 3),
                      rep(1:3, each = 3)) # paste0('period ', 1:9)

my.labs = names(npp2)#paste('period', 1:9)
names(my.labs) <- names(npp2)
gplot(npp2) + 
  geom_raster(aes(fill = value)) +
  facet_grid(~ variable, #switch = 'x',
             labeller = labeller(variable = my.labs)) +
  geom_point(data = gridfst, aes(x = x, y = y, colour = ' grid'), 
             alpha = 0.6, size = 1)+
  scale_fill_gradientn(colours = rev(terrain.colors(225)[c(seq(1,80,2), seq(81,225,10))]),
                       na.value = 'white',
                       name = "NPP") +
  scale_colour_manual(values = c('black'),name = NULL)+
  coord_equal()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0, "cm", data = NULL),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(colour = 'black'),
        legend.text = element_text(size = 8),
      #  legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, vjust = -1),
        strip.background = element_rect(fill = 'white', colour = NA),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        legend.key.width = unit(0.25, 'cm'),
        legend.spacing = unit(0.0, 'cm'),
        legend.position = 'right',
        plot.margin = margin(0, 0, 0, 0, "cm"))-> p;p

# together ----------
g_legend<-function(a.gplot){
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
l1 <- g_legend(con)
l2 <- g_legend(p)

lay <- rbind(
             c(rep(5,11),1),
             c(rep(5,11),1),
             c(rep(5,11),1),
             c(rep(5,11),1),
             c(rep(5,11),1),
             c(rep(4,10),NA,NA),
             c(rep(4,10),2,2),
             c(rep(4,10),2,2),
             c(rep(4,10),2,2),
             c(rep(4,10),2,2),
             c(rep(4,10),NA,NA),
             c(rep(4,10),3,3),
             c(rep(4,10),3,3),
             c(rep(4,10),3,3),
             c(rep(c(4,4),each = 5),3,3),
             c(rep(c(4,4),each = 5),NA,NA))
fig1 <- grid.arrange(l2,
                     pic, pic2, 
                     
                     con + #labs(tag = 'A)') +
                       theme(legend.position = 'none',
                                 plot.margin = margin(0, 0, 0, 0, "cm")),
                     p + #labs(tag = 'B)') + 
                       theme(legend.position = 'none'),
                     nrow = 2, layout_matrix = lay)
ggsave('./figures/Richard/fig1_captures_phases_npp_stable.png', plot = fig1,
       height = 12.2, width = 20, units = 'cm', dpi = 600)
