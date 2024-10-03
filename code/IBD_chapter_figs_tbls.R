
# load --------------------------------------------------------------------
em.coef.models <-function(m, spp.ph = T, modelname = 'Interaction'){
  
 df1 <- summary(m)$coefficients %>% 
    as.data.frame() 
 df1$df <- ifelse('df' %in% names(df1), df1$df, m$df.residual)
 df1 %>% 
    mutate(Species = ifelse(spp.ph, 'P. hermanns',
                            'S. youngsoni'),
           t = qt(0.975, df),
           me = t*`Std. Error`,
           lower = round(Estimate - me, 3),
           upper = round(Estimate + me, 3),
           confidence = paste0(round(Estimate,3), ' (',
                               lower, ' - ', upper, ')'),
           Parameter = rownames(.),
           #Species = ifelse(duplicated(Species), NA, Species),
           sig = `Pr(>|t|)` < 0.05,
           R2 = round(MuMIn::r.squaredGLMM(m)[2], 2),
           R2 = ifelse(duplicated(R2), NA, R2),
           model = modelname,
           Model = ifelse(duplicated(model), NA, model)) %>% 
    dplyr::select(Species, Model, Parameter, confidence, R2, sig)%>% 
    rename(`Estimate (95% CI)` = confidence)
}


source('./code/libraries.R')
source('./code/functions_data-tables.R')
source('./code/functions_correlograms.R')
phaseCol <- c(low = terrain.colors(10)[7],
              increase = terrain.colors(10)[1],
              decrease = terrain.colors(10)[4],
              stable = 'grey50')

nppCol <- rev(terrain.colors(225)[c(seq(1,80,2),
                                    seq(81,225,10))])[-1]



ind <- read.csv('./output/individual_analysis.csv') %>% 
  mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                   rep(1:3, each = 3))),
         km.log = ifelse(km == 0, log(km + 0.1), log(km)),
         npp = round(npp, 2),
         phase = ifelse(species == 'Sminthopsis youngsoni', 'stable', phase),
         phase = factor(phase, levels = c('low', 'increase',
                                          'decrease','stable')))
table(ind$sex_pairs)
table(ind$phase)

#________ ------------

# raw ibd 5km ----------------
triangle <- data.frame(x = c(1,1,5),
                       y = c(88,60, 88, 36.5,25,36.5), sex_pairs = 'm-m',
                       npp = 100,species = rep(c('Pseudomys hermannsburgensis',
                                                 'Sminthopsis youngsoni'), each = 3))
trianglec <- data.frame(x = c(1,5),
                        y = c(88, 88,36.5,36.5), sex_pairs = 'm-m',
                        npp = 100,species = rep(c('Pseudomys hermannsburgensis',
                                                  'Sminthopsis youngsoni'), each = 2))
ggplot(ind %>% filter(km < 5.1) %>% 
         arrange(npp),
       aes(km, euclidean,
           fill = npp, group = species))+
  facet_wrap(~species,ncol = 2, scale = 'free')+
  geom_point(aes(pch = sex_pairs), colour = 'grey90',
             size = 3)+
  geom_line(data = triangle, aes(x,y), colour = 'grey55',
            linewidth = 2, lty = 2, alpha = 0.5)+
  scale_shape_manual(values = c(24:21),
                     name = 'Comparison')+ 
  geom_line(data = trianglec, aes(x,y), colour = 'grey55',
            linewidth = 2, lty = 2, alpha = 0.5)+
  guides(shape = guide_legend( 
    override.aes=list(colour = 'black')))+
  scale_fill_gradientn(colours = nppCol,
                       na.value = 'white',
                       name = "NPP")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'italic', size = 10.5))+  
  ylab('Genetic distance')+
  xlab('Geographical distance (km)')+
  theme(panel.grid = element_blank()) -> rawslope; rawslope

tiff(filename = './figures/IBD/fig_ibd_5km.tiff', 
     units = 'cm', height = 8, width = 16, res = 600)
rawslope
dev.off()

# All data ----------
ibdDataAll <- lapply(c('herm', 'young'), 
                  function(x) filter(ind, grepl(x, species),
                                     complete.cases(metres)))

## mantel correlog --------
system.time(ibdcorr <- lapply(ibdDataAll, 
                              em.ibd.period.correlog, 
                              test = 'pearson'))

ibdcorr[[1]] %>%
  bind_rows(ibdcorr[[2]]) %>% 
  mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    d.class = exp(breaksto)
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.)) -> mantel.results2 

xibdcorr <- do.call('rbind', ibdcorr) %>% 
  filter(class.index < 0) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0, 
                             Mantel.cor)) 

# sex bias ----------------

df_index <- data.frame(i = 1:6, x = rep(c('herm', 'young'), each = 3),
                       y = rep(c('f-m', 'm-m', 'f-f'), 2)) %>% 
  split(., .$i)

ibdData <- lapply(df_index, function(x) filter(ind, grepl(x$x, species),
                                               sex_pairs %in% unlist(str_split(x$y, ',')),
                                               complete.cases(metres)))
## correlogram ----------
ibdcorrsexes <- lapply(ibdData[c(2,3,5,6)], em.ibd.period.correlog, test = 'pearson')



do.call('bind_rows', ibdcorrsexes)%>% 
  mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.)) -> mantel.results 

# figure
xibdcorrexes <- do.call('rbind', ibdcorrsexes) %>% 
  filter(class.index < 0) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0, 
                             Mantel.cor),
         sex_pairs = ifelse(sex_pairs == 'f-f', 'Female', 'Male')) 

phaseCol2 <- phaseCol
names(phaseCol2) <- toupper(str_sub(names(phaseCol2), 1,1))

# sex correlograms --------

mantel.results %>%  
  bind_rows(mantel.results2) %>%
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0,
                             Mantel.cor),
         sex_pairs = case_when(
           sex_pairs == 'f-f' ~ 'Females',
           sex_pairs == 'm-m' ~ 'Males',
           sex_pairs == 'all' ~ 'All'
         )
         ) %>%
  ggplot(aes(exp(breaksto), Mantel.cor, 
             fill= sig, colour = npp,
             group = phaseNo))+
  geom_hline(yintercept = 0, #size = 1,
             colour = 'grey', 
             linetype = 'dashed')+
  geom_line(size = 0.75)+
  geom_point( shape = 21, size = 2, colour = 'black')+
  
  scale_x_log10(breaks = unique(exp(mantel.results$breaksto)))+
  #guides(size = 'none')+
  # scale_colour_manual(values = unname(phaseCol),
  #                     labels = c('low', 'increase', 'decrease'))+
  scale_fill_manual(values = c('white', 'pink'),
                    labels = c('No IBD', 'IBD'),
                    name = NULL)+
  scale_colour_gradientn(colours = nppCol,
                         na.value = 'white',
                         name = "NPP")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'italic', size = 10),
        strip.text.y = element_text(size = 12))+
  xlab('Distance classes (km)')+
  ylab(expression(italic('r')))+
  facet_grid(sex_pairs~species, scale = 'free_x') -> mantel.plot.sex

tiff('./figures/IBD/fig_correlogram_sexes.tiff', res = 600,
     units = 'cm', width = 16, height = 18)
mantel.plot.sex
dev.off()

# plot simplified -----------

xibdcorrexes %>% 
  bind_rows(xibdcorr) %>% 
  mutate(sex_pairs = ifelse(grepl('all', sex_pairs), 'All', sex_pairs)) %>% 
ggplot(aes(npp, Mantel.cor, colour = sex_pairs,
                         fill = sex_pairs, group = sex_pairs))+
  scale_shape_manual(values = c(24,22))+
  geom_point(pch = 21, #aes(pch = sex_pairs), 
             size = 3,
            colour = 'black',
             alpha = 0.5
  )+
  theme_bw()+
  guides(shape = 'none')+
  geom_smooth(method = 'lm', se = F, show.legend = F)+
  scale_x_continuous(trans = log_trans(), 
                     breaks = 8 * 2^seq(0,5,1))+
  facet_grid(~species) +
  ylab('mantel correlation (<1km)')+
  scale_fill_manual(values = c('grey', 'red', 'black'), #unname(phaseCol),
                    name = 'Comparison')+
  scale_colour_manual(values = c('grey', 'red', 'black'),
                      name = 'Comparison',
                      #labels = c('females', 'males')
                      )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.89,0.65),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        legend.key.size = unit(0.5, units = 'cm'),
        strip.text.x = element_text(face = 'italic', size = 11),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank())+
  ylab(expression(italic("r")*' (distance class 0 - 1km)'))+
  xlab(expression('Net primary productivity '[log-scale])) -> fig_r_npp_sexes2;fig_r_npp_sexes2


ggsave('./figures/IBD/fig_r_npp_sexes.png', plot = fig_r_npp_sexes2,
       units = "cm", width = 14, height = 9, dpi = 600)

## sex x phase ---------
ibdsexes <- rbind(ibdDataAll[[1]],ibdDataAll[[2]]) %>% 
  mutate(sites = factor(ifelse(km == 0, 'within', 
                               ifelse(km>1, 'beyond', 'between')),
                        levels = c('within', 'between', 'beyond')),
         ibd = factor(ifelse(sites == 'beyond', '> 1km', '< 1km'),
                      levels = c('< 1km', '> 1km')),
         phase = factor(phase, levels = c('low', 'increase', 'decrease'))) %>% 
  filter(sex_pairs %in% c('m-m',
                          # 'f-m', 
                          'f-f'))

ibdsexes<-ibdsexes %>% 
  mutate(phase = factor(phase, levels = c('low', 'increase', 
                                          'decrease', 'stable')))


ibdplotting <- ibdsexes %>% filter(km <= 1) %>% 
  arrange((npp)) %>% 
  mutate(phase = factor(phase,
                        levels = c('low', 'increase', 'decrease', 'stable')),
         hline = ifelse(grepl('herm', species), 92, 39),
         sex_pairs = ifelse(sex_pairs == 'f-f', 'Females', 'Males'),
         species = ifelse(species == 'Pseudomys hermannsburgensis',
                          'P. hermannsburgensis',
                          'S. youngsoni'))  
  (ibdplotting$phase %>% table %>% names)== names(phaseCol)
  ggplot(ibdplotting,aes(km, euclidean, colour = phase,
             fill = npp, group = phase))+
  facet_grid(species~sex_pairs, scale = 'free')+
  geom_point(pch = 21, colour = 'white', alpha = 0.7, fill = 'grey70')+
  scale_shape_manual(values = c(21:24))+ 
  scale_fill_gradientn(colours = nppCol,
                       na.value = 'white',
                       name = "NPP")+
  geom_hline(aes(yintercept = hline), alpha = 0.8,
             lty = 2, colour = 'grey70', lwd = 0.3)+
  geom_smooth(method = 'lm', se = F, lwd = 1)+
  scale_color_manual(values = phaseCol,
                     labels = c('Low','Increase', 'Decrease', 'Stable'),
                     name = 'Population phase')+
  #scale_color_manual(values = phaseCol)+
  theme_bw()+
  scale_x_continuous(limits = c(0,1), labels = c(seq(0,0.8, 0.25),'1'))+
  theme(panel.grid = element_blank(),
        #legend.position = 'right',
        legend.direction="vertical",
        #legend.box = "horizontal",
        legend.position = c(0.85, 0.14),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "grey"),
        legend.key.size = unit(0.35, units = 'cm'),
        legend.title = element_text(size = 10),
        strip.background = element_blank(),
        axis.title = element_text(size = 13),
        strip.text.y = element_text(face = 'italic', size = 11),
        strip.text.x = element_text(face =   'bold',size = 12))+  
  ylab('Genetic distance')+
  xlab('Geographical distance (km)') -> IBD_sexes_phase; IBD_sexes_phase


tiff(filename = './figures/IBD/fig_slope.tiff', 
     units = 'cm', height = 12, width = 14, res = 600)
IBD_sexes_phase
dev.off()


# ems coeffs -----------
lm(Mantel.cor ~ log(npp)*sex_pairs, 
   data = filter(xibdcorrexes, grepl('herm', species))) %>% summary 
lm(Mantel.cor ~ log(npp)*sex_pairs, 
   data = filter(xibdcorrexes, grepl('young', species))) %>% summary

mod_r_ph <- lm(Mantel.cor ~ log(npp)+sex_pairs , #weights = (n.dist),
               data = filter(xibdcorrexes, grepl('herm', species))) 
mod_r_sy <- lm(Mantel.cor ~ log(npp)+sex_pairs, 
               data = filter(xibdcorrexes, grepl('young', species)))

sjPlot::tab_model(mod_r_ph, mod_r_sy, digits = 3,
                  df.method = 'satterthwaite')

modeltbsx <- rbind(em.coef.models(mod_r_ph, T,'P. hermann'),
      em.coef.models(mod_r_sy, F, 'S. youngsoni')) %>% 
  dplyr::select(-Model) %>% 
  mutate(Species = ifelse(duplicated(Species), NA, Species))

sexes_raw_ph <- lmer(euclidean ~ km * sex_pairs*phase + (1|phaseNo),
                     data = filter(ibdsexes, km <= 1,
                                   grepl('herm', species)))
sexes_raw_sy <- lmer(euclidean ~  km * sex_pairs  + (1|phaseNo),
                     data = filter(ibdsexes, km <= 1,
                                   grepl('young', species)))

sjPlot::tab_model(sexes_raw_ph, sexes_raw_sy, digits = 3,
                  df.method = 'satterthwaite')


modeltbsx2 <- rbind(em.coef.models(sexes_raw_ph, T,'P. hermann'),
      em.coef.models(sexes_raw_sy, F, 'S. youngsoni'))%>% 
  dplyr::select(-Model) %>% 
  mutate(Species = ifelse(duplicated(Species), NA, Species))


# flex tb ----------
mtb <- modeltbsx 

mtb$Parameter <- gsub('phase', '', mtb$Parameter)
mtb$Parameter <- gsub('sex_pairs', '', mtb$Parameter)
mtb$Parameter <- gsub(':', ' : ', mtb$Parameter)

mtb$n[complete.cases(mtb$Species)] <- 18

mtb <- mtb %>% relocate(n, .after = 'Species')

mtb[-ncol(mtb)] %>% 
  flextable() %>% 
  bold(part = 'header') %>%
  italic(j = 1) %>% 
  italic(part = 'header', j = 2) %>% 
  border_remove() %>%  
  bold(j = 4, i = which(mtb$sig)) %>% 
  autofit() %>%
  width(j = 3, width = 3, unit = 'cm') %>% 
  hline(i = c(nrow(mtb)),
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = which(mtb$Species == 'S. youngsoni')-1,
        border = fp_border(color = "grey60",
                           width = 2, 
                           style = 'dashed')) %>% 
  compose(
    part = "body", j = 3, i = grep('Intercept', mtb$Parameter),
    value = as_paragraph(('Intercept'), as_sub(' female'))) %>% 
  compose(
    part = "body", j = 3, i = grep('Male', mtb$Parameter),
    value = as_paragraph(('Male'), as_sub(' difference'))) %>% 
  compose(
    part = "body", j = 3, i = grep('log', mtb$Parameter),
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) -> fxtb_r_npp;fxtb_r_npp

flextable::save_as_docx(fxtb_r_npp, path = './figures/IBD/table1_sexes_r.docx')





mtb <- modeltbsx2

mtb$Parameter <- gsub('phase', '', mtb$Parameter)
mtb$Parameter <- gsub('sex_pairs', '', mtb$Parameter)
mtb$Parameter <- gsub(':', ' : ', mtb$Parameter)
mtb$n[complete.cases(mtb$Species)] <- c(4901, 561)

mtb <- mtb %>% relocate(n, .after = 'Species')

mtb[,-ncol(mtb)] %>% 
  flextable()  %>% 
  bold(i = (1:nrow(mtb))[mtb$sig], j = 4) %>% 
  bold(part = 'header') %>% 
  #italic(part = 'header') %>% 
  autofit() %>% 
  #width(j = 1, unit = 'cm', width = 4.2) %>% 
  italic(j = c(1), i = which(complete.cases(mtb$Species))) %>% 
  border_remove() %>% 
  hline(i = nrow(mtb), 
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(i = which(mtb$Species == 'S. youngsoni')-1, 
        border = fp_border(color = "grey40", width = 2)) %>%
  hline(part = 'header', 
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(part = 'body', 
        i = which(mtb$Parameter %in% c('decrease','m-m : decrease')), 
        j = 3:(ncol(mtb)-1),
        border = fp_border(color = "grey60", 
                           width = 1,style = 'dashed')) %>% 
  compose(
    part = 'body', j = 3, i = which(mtb$Parameter == 'Intercept'), 
    value = as_paragraph(('Intercept'), as_sub('  low f-f'))
  ) %>% 
  compose(
    part = 'body', j = 3, i = which(mtb$Parameter == 'km'), 
    value = as_paragraph(('km'), as_sub('  low f-f'))
  ) %>% 
  compose(
    part = 'header', j = ncol(mtb)-1, 
    value = as_paragraph(('R'), as_sup('2'))
  ) %>%
  compose(
    part = 'header', j = 5, 
    value = as_paragraph(('R'), as_sup('2'))
  ) %>% 
  italic(part = 'header', j = 2)-> flextable_sex_models;flextable_sex_models

flextable::save_as_docx(flextable_sex_models, path = './figures/IBD/table2_sexes_ibd.docx')

# table S1 ---------------
# ems coeffs -----------
lm(Mantel.cor ~ log(npp)*sex_pairs, 
   data = filter(xibdcorrexes, grepl('herm', species))) %>% summary 
lm(Mantel.cor ~ log(npp)*sex_pairs, 
   data = filter(xibdcorrexes, grepl('young', species))) %>% summary

mod_r_ph <- lm(Mantel.cor ~ log(npp)*sex_pairs , #weights = (n.dist),
               data = filter(xibdcorrexes, grepl('herm', species))) 
mod_r_sy <- lm(Mantel.cor ~ log(npp)*sex_pairs, 
               data = filter(xibdcorrexes, grepl('young', species)))

sjPlot::tab_model(mod_r_ph, mod_r_sy, digits = 3,
                  df.method = 'satterthwaite')

modeltbsxAppendix <- rbind(em.coef.models(mod_r_ph, T,'P. hermann'),
                   em.coef.models(mod_r_sy, F, 'S. youngsoni')) %>% 
  dplyr::select(-Model) %>% 
  mutate(Species = ifelse(duplicated(Species), NA, Species))



# flex tb ----------
mtb <- modeltbsxAppendix 

mtb$Parameter <- gsub('phase', '', mtb$Parameter)
mtb$Parameter <- gsub('sex_pairs', '', mtb$Parameter)
mtb$Parameter <- gsub(':', ' : ', mtb$Parameter)

mtb$n[complete.cases(mtb$Species)] <- 18

mtb <- mtb %>% relocate(n, .after = 'Species')

mtb[-ncol(mtb)] %>% 
  flextable() %>% 
  bold(part = 'header') %>%
  italic(j = 1) %>% 
  italic(part = 'header', j = 2) %>% 
  border_remove() %>%  
  bold(j = 4, i = which(mtb$sig)) %>% 
  autofit() %>%
  width(j = 3, width = 4.25, unit = 'cm') %>% 
  hline(i = c(nrow(mtb)),
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = which(mtb$Species == 'S. youngsoni')-1,
        border = fp_border(color = "grey60",
                           width = 2, 
                           style = 'dashed')) %>% 
  compose(
    part = "body", j = 3, i = grep('Intercept', mtb$Parameter),
    value = as_paragraph(('Intercept'), as_sub(' female'))) %>% 
  compose(
    part = "body", j = 3, i = grep('Male', mtb$Parameter),
    value = as_paragraph(('Male'), as_sub(' difference'))) %>% 
  compose(
    part = 'header', j = ncol(mtb)-1, 
    value = as_paragraph(('R'), as_sup('2'))
  ) %>% 
  compose(
    part = "body", j = 3, i = grep('log', mtb$Parameter),
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) %>% 
  compose(
    part = "body", j = 3, i = grep(':', mtb$Parameter),
    value = as_paragraph(('log'), as_sub('e'), ('NPP : Male'), as_sub(' interaction'))) -> fxtb_r_npp_Sup;fxtb_r_npp_Sup

flextable::save_as_docx(fxtb_r_npp_Sup, 
                        path = './figures/IBD/tableS1_sexes_r_interaction.docx')


