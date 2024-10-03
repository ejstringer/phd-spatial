# analysis 2 - structure and dispersal

# create data

# filter data
source('./code/libraries.R')
source('./code/functions_genetic-data-creation.R')
source('./code/functions-npp-data-creation.R')

ph <- em.filtering('ph')
sy <- em.filtering('sy')

neprep <- lapply(list(pherm = ph, syoung = sy), em.ne_prep)
n <- sapply(neprep, function(x) sapply(x, nInd)) %>% 
  as.data.frame() %>% 
  mutate(period = rownames(.),
         phaseNo = paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3))) %>% 
  pivot_longer(pherm:syoung, names_to = 'species', values_to = 'n') %>% 
  mutate(species = ifelse(species == 'pherm', 'Pseudomys hermannsburgensis',
                          'Sminthopsis youngsoni')) %>% 
  arrange(species)
write.csv(n, './output/sampe_sizes.csv', row.names = FALSE)


system.time(gdist <- lapply(list(pherm = ph, syoung = sy), 
                em.individual_relatedness)) # ~12mins

system.time(fst <- lapply(list(pherm = ph, syoung = sy), 
                          em.phase_fst)) # ~22mins

npp <- em.npp_layers(ph)

nppmean <- read.csv('./output/npp_means.csv') 



periodID <- data.frame(phaseNo = paste0(rep(c('L', 'I', 'D'), 3),
                                        rep(1:3, each = 3)),
                       period = paste0('period ', 1:9))

system.time(individual <- gdist$pherm %>% 
  bind_rows(gdist$syoung) %>% 
    filter(within) %>% 
    left_join(nppmean, by = c('period.x' = 'period')) %>% 
    rename(phase = phase.x, period = period.y) %>% 
    left_join(periodID) %>% 
  dplyr::select(species, period, phaseNo, phase, kinship, ibdecent, euclidean, 
                metres, km, npp, mc, pairs, idpairs, trip.x, trip.y, 
                sex.x, sex.y, id.x, id.y)) 
names(individual)
head(individual)

individual$npp.log = log(individual$npp)
individual$phase <- factor(individual$phase,
                             levels = c('low', 'increase', 'decrease'))

individual$phaseNo <- factor(individual$phaseNo,
                           levels =  periodID$phaseNo)

individual$rep <- as.numeric(str_sub(individual$phaseNo, 2,2))

individual$daydist <- abs(ymd(individual$trip.x) - ymd(individual$trip.y))
individual$sex_pairs <- paste0(individual$sex.x, '-', individual$sex.y)

individual$sex_pairs <- ifelse(individual$sex_pairs %in% c('f-m', 'm-f'), 'f-m',
                               individual$sex_pairs)

individual$sex_pairs <- ifelse(grepl('NA', individual$sex_pairs), 'unknown', 
                               individual$sex_pairs)

individual$sex_pairs %>% table

individual_final <- individual %>% dplyr::select(-c(sex.x, sex.y, trip.x, trip.y)) 
head(individual_final)


ngrids <- c(
  sapply(split(fst$pherm, fst$pherm$period), function(x)
    length(unique(c(x$pop1, x$pop2)))),
  sapply(split(fst$syoung, fst$syoung$period), function(x)
    length(unique(c(filter(x, grepl('you', species))$pop1,
                    filter(x, grepl('you', species))$pop2))))) %>% 
  data.frame(ngrids = ., period = names(.), 
             species = rep(c('Pseudomys hermannsburgensis',
                             'Sminthopsis youngsoni'), each = 9)[-18])

fst_final <- fst$syoung %>% 
  bind_rows(fst$pherm) %>% 
  #right_join(ne) %>%
  left_join(ngrids) %>%
  left_join(nppmean) %>%  
  left_join(periodID)  %>%   
  mutate(npp.log = log(npp),
         fst = ifelse(fst < 0, 0, fst),
         fst.log = log(fst + 0.00001),
         phase = factor(phase, levels = c("low", 'increase','decrease')),
         phaseNo =  factor(phaseNo, levels = periodID$phaseNo),
         rep = as.numeric(str_sub(phaseNo, 2,2)),
         species2 = ifelse(grepl('Pse', species), 'ph', 'sy'),
         mRate = NA, #(1/ifelse(fst > 0, fst, NA)-1)/(4*(ne)),
         mRate.log = NA, # log(mRate),
         nMigrants = (1/ifelse(fst > 0, fst, NA)-1)/4,
         nMigrants.log = log(nMigrants)
         )

write.csv(individual_final, './output/individual_analysis.csv', row.names = F)
write.csv(fst_final, './output/fst_analysis.csv', row.names = FALSE)

# summarise ------
unique(c(fstdata$pop1, fstdata$pop2)) %>% is.na() %>% table
# end -------------
phasedf <- fst_final %>% 
  group_by(species, npp, period, phase)  %>% 
  summarise(fst = mean(fst, na.rm = T),
            n = n(),
            grp = ifelse(grepl('herm', species), as.character(phase), 'na')
            ) 

ggplot(phasedf, aes(period, fst, fill = phase)) +
  geom_boxplot()+
  theme_classic()

mfst <- lmer(fst ~ npp.log * phase + (1|period),
            data = fst_final[fst_final$species2=='ph',])
summary(mfst)
MuMIn::r.squaredGLMM(mfst) # marginal R2 only takes the fixed effects, 
AIC(mfst,mfst2)
sjPlot::tab_model(mfst)

mfst2<- lmer(fst ~ npp.log * phase + (1|period),
            data = fst_final[fst_final$species2=='sy',])
summary(mfst)
MuMIn::r.squaredGLMM(mfst) # marginal R2 only takes the fixed effects, 

sjPlot::tab_model(mfst, mfst2, dv.labels = c("ph", "sy"), digits = 4)



ggplot(fst_final,
       aes(npp, fst, colour = phase,fill = phase, group = phase)) +
  #scale_x_log10()+
  #scale_y_log10()+
 scale_color_manual(values = phaseCol, name = 'population phase')+
  #scale_color_manual(values = c('black', 'grey'), name = 'species')+
  scale_fill_manual(values = phaseCol, name = 'population phase')+
  geom_smooth(method = 'lm', se = F, colour = 'grey', alpha = 0.2, size = 2)+
  geom_jitter(alpha = 0.5, width = 0.02)+
  geom_point(data = phasedf, aes(size = n), show.legend = F,
             alpha = 0.95, pch = 21, colour = 'black')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'right',
        strip.text = element_text(face = 'italic'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )+
  facet_grid(~species)+
  ylab(expression("F"[ST]))+
  xlab(expression("Net primary productivity"[` log-scale`]))


mfst<- lmer(mRate.log ~  npp.log * phase + (1|period),
            data = fst_final[fst_final$species2=='ph',])
summary(mfst)
MuMIn::r.squaredGLMM(mfst) # marginal R2 only takes the fixed effects, 

anova(lm(mRate.log ~ phase,
         data = fst_final[fst_final$species2=='ph',]))

mfst2<- lmer(mRate.log ~ npp.log * phase + (1|period),
            data = fst_final[fst_final$species2=='sy',])
summary(mfst2)
MuMIn::r.squaredGLMM(mfst2) # marginal R2 only takes the fixed effects, 

sjPlot::tab_model(mfst, mfst2, dv.labels = c("ph - migration", "sy - migration"))


anova(lm(mRate.log ~ phase,
         data = fst_final[fst_final$species2=='sy',]))

aov(mRate.log ~ phase,
           data = filter(fst_final,species2 == 'ph',
                         period %in% c('period 7',
                                       'period 8',
                                       'period 9'))) %>% TukeyHSD()

ggplot(fst_final, aes(period, mRate, fill = phase)) +
  geom_boxplot()+
  theme_bw()+
  scale_y_log10()+
  facet_grid(~species)

lm(mRate.log ~ npp.log+ phase,
    data = fst_final[fst_final$species2=='sy',]) %>% summary

ggplot(filter(fst_final, mRate < 1),
       aes(npp, mRate, colour = phase,fill = phase, group = phase)) +
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values = phaseCol, name = 'population phase')+
  #scale_color_manual(values = c('black', 'grey'), name = 'species')+
  scale_fill_manual(values = phaseCol, name = 'population phase')+
  geom_smooth(method = 'lm', se = F, colour = 'grey', alpha = 0.2, size = 2)+
  geom_jitter(alpha = 0.5, width = 0.02)+
  # geom_point(data = phasedf, aes(size = n), show.legend = F,
  #            alpha = 0.95, pch = 21, colour = 'black')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'right',
        strip.text = element_text(face = 'italic'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )+
  facet_grid(~species)+
  ylab(expression("migration rate"))+
  xlab(expression("Net primary productivity"[` log-scale`]))


ggsave('./figures/temp_fst_npp.png')

ggplot(phasedf, aes(period, fst, colour = phase, fill = pairs)) +
  geom_jitter()+
  theme_classic()+
  theme(legend.position = 'none') -> p

plotly::ggplotly(p)

ggplot(filter(gdist$syoung, within, 
                               #pairs %in% c('MC5-MC5','FRS3-FRS3','PP1-PP1')
                              km<5 ), 
                        aes(km, euclidean, colour = pairs, group = species))+
  geom_jitter()+
    theme_classic()+
    facet_wrap(~period.x)+
 # geom_hline(yintercept = 0.5)+
  theme(legend.position = 'none')+
  geom_smooth(se = F)


filter(gdist$syoung,within,  gridId.x == gridId.y) %>% 
  group_by(period.x, pairs) %>% 
  summarise(ibd = median(ibdecent, na.rm = T),
            n = n()) %>%
  mutate(name = 'ibdecent') -> outliers  
  hist(outliers)
plotly::ggplotly(ggplot(filter(outliers, n > 0),
                        aes(y = ibd, x = period.x, 
                            colour = pairs))+
  geom_jitter(aes(size = n))+
  theme_classic()+
  theme(legend.position = 'none'))


# individual
gg <- gdist$pherm %>% filter(ibdecent >= 0.5, within) %>% 
  rename(phase = phase.x,
         period = period.x)

ggplot(gg, aes(period, ibdecent, colour = phase))+
  geom_point()+
  theme_classic()

lm(ibdecent ~ phase, data = gg) %>% anova





