
# load -------
source('./code/libraries.R')
source('./code/functions_data-tables.R')
source('./code/functions_genetic-data-creation.R')
#source('./code/functions_correlograms.R')

phaseCol <- c(
              increase = terrain.colors(10)[1],
              decrease = terrain.colors(10)[4],
              low = terrain.colors(10)[7],
              stable = 'grey70')

nppCol <- rev(terrain.colors(225)[c(seq(1,80,2),
                                    seq(81,225,10))])[-1]

fstdata <- read.csv('./output/fst_analysis.csv') %>% 
  mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                   rep(1:3, each = 3))),
         phase = ifelse(grepl('young', species), 'stable', phase),
         phase = factor(phase, levels = c('increase', 'decrease','low', 'stable')))

# ind <- read.csv('./output/individual_analysis.csv') %>% 
#   mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
#                                                    rep(1:3, each = 3))),
#          km.log = ifelse(km == 0, log(km + 0.1), log(km)),
#          npp = round(npp, 2),
#          phase = factor(phase, levels = c('low', 'increase', 'decrease')))


fstdata$phase %>% table(., useNA = 'always')
# data setup -----------------

caps <- em.rain_caps_phase() %>% 
  mutate(phaseNo = factor(phaseNo, levels =paste0(rep(c('L', 'I', 'D'), 3),
                                                  rep(1:3, each = 3))))

caps %>% 
  # group_by(phaseNo, trip) %>% 
  # summarise(captures = mean(captures, na.rm = T),
  #           capturesSy = mean(capturesSy, na.rm = T)) %>% 
  # ungroup() %>% 
  group_by(phaseNo, phase) %>% 
  summarise(captures = mean(captures, na.rm = T),
            capturesSy = mean(capturesSy, na.rm = T)) %>% 
  filter(complete.cases(phaseNo)) %>% 
  arrange(phase) %>% 
  left_join(unique(fstdata[,c('phaseNo', 'npp')]))

fstdata %>% dplyr::select(phaseNo, phase, npp) %>% 
  unique() %>% mutate(npp = round(npp)) %>% 
  arrange(phase)

captures <- cbind(tapply(caps$captures, caps$phaseNo, mean, na.rm = T),
                  tapply(caps$capturesSy, caps$phaseNo, mean, na.rm = T))



colnames(captures) <- c('ph', 'sy')
captures2 <- as.data.frame(captures) %>% 
  mutate(period = paste('period', 1:9)) %>% 
  pivot_longer(cols = ph:sy, names_to = 'species2', values_to = 'captures')

captures3 <- as.data.frame(captures) %>% 
  mutate(period = paste('period', 1:9),
         phaseNo = rownames(.)) %>% 
  pivot_longer(cols = ph:sy, names_to = 'species2', values_to = 'captures')



fstdata2 <- left_join(fstdata, captures2)%>% 
  mutate(popsize = case_when(
    phase == 'low' ~ 'small', 
    phase ==  'decrease' ~ 'medium',
    phase ==  'increase' ~ 'large'
  ),
  popsize = factor(popsize, levels = c('small', 'medium', 'large')))

## fst data  --------

dat <- fstdata2 %>% 
  mutate(ph = (grepl('herm', species))) %>% 
  #mutate(fst = log(fst+0.001)) %>%
  dplyr::select(fst, species,ph, phaseNo, npp, npp.log, phase, captures)
dat %>% head
plot(dat$npp.log, scale(dat$npp.log))
plot(dat$captures,scale(dat$captures))

datScaled <- dat %>% 
  mutate(captures = scale(captures),
         npp.log = scale(npp.log)) 

# captures vs npp ----------------
cn <- dat |>
  dplyr::select(species,ph, phaseNo, npp.log,npp, phase, captures) |>
  unique() |>
  glimpse() %>% 
  arrange(species)

cn %>% 
  group_by(species) %>% 
  summarise(cor = cor(captures, npp.log)) 


figcor <- ggplot(cn, aes(y = captures, x = npp, group = species)) + 
  stat_smooth(data = filter(cn, ph),
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(fill = phase), size = 5, pch = 21) +
  stat_smooth(method = lm, se = F, colour = 'grey30', lty = 1) +
  facet_wrap(~species, ncol = 1, scale = "free_x") +
  theme_bw()+
  scale_x_log10()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.78,0.3),
        legend.key.size = unit(0.5, units = 'cm'),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  xlab(expression('Net primary productivity '[log-scale]))+
  scale_fill_manual(values = phaseCol,
                    labels = str_to_title(names(phaseCol)),
                    name = 'Population phase')+
  theme(legend.position = 'right',
        legend.background = element_blank())+
  ylab('Mean captures'); figcor
figcor 
ggsave('./figures/FST/caps_npp_correlation_legend.png', plot = figcor,
       units = "cm", width = 14, height = 14, dpi = 600)

m <- (lm(captures ~ npp.log, data = filter(cn, species == "Pseudomys hermannsburgensis")))
m2 <-(lm(captures ~ npp.log, data = filter(cn, species == "Sminthopsis youngsoni")))
sjPlot::tab_model(m, m2, dv.labels = c('P. hermannsburgensis', 'S. youngsoni'),
                  title = 'Captures')


#-------------------------------------------------------------------------------
# fst by captures
# mean fst
mfc <- dat |>
  group_by(species, phaseNo, phase, captures) |>
  summarise(mean.fst = mean(fst)) |>
  glimpse()

figCaps <- ggplot(dat, aes(y = fst, x = captures, fill = phase, group = species)) +
  stat_smooth(data = filter(dat, ph),
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(colour = phase, fill = phase), size = 3, alpha = 0.2) +
  geom_point(data = mfc, aes(y = mean.fst), colour = 'black', 
             shape = 23, size = 4) +
  stat_smooth(method = lm, se = F, colour = 'grey30') +
  facet_wrap(~ species, ncol = 1, scale = "free_x") +
  theme_bw()+
  xlab('Mean captures')+
  ylab(expression(italic('F')[ST]))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
         legend.position = 'none',#c(0.85,0.73),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase',
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase');figCaps

datScaled$x <- datScaled$captures
m3 <- (lmer(fst ~ x  * phase + (1|phaseNo),
                   data = filter(datScaled,
                                 species == "Pseudomys hermannsburgensis")))

summary(m3)
sjPlot::tab_model(m3,dv.labels = c('Captures'),
                  title = 'P. hermannsburgensis: Fst', collapse.ci = F,
                  df.method = "satterthwaite", digits = 3)
# fst by npp
# mean fst
mfn <- dat |>
  group_by(species, phaseNo, phase, npp.log, npp) |>
  summarise(mean.fst = mean(fst),
            x = npp.log) |>
  glimpse()

figNPP <- ggplot(dat, aes(y = fst, x = npp, fill = phase)) +
  stat_smooth(data = filter(dat, ph), show.legend = F,
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(colour = phase, fill = phase), size = 3, show.legend = F,
             alpha = 0.2) +
  geom_point(data = mfn, aes(y = mean.fst), 
             colour = 'black', shape = 23, size = 4) +
  stat_smooth(method = lm, se = F, colour = 'grey30', show.legend = F,
              aes(group = species)) +
  facet_wrap(~ species, ncol = 1, scale = "free_x") +
  scale_colour_discrete(na.translate = F)+
  theme_bw()+
  scale_x_log10()+
  xlab(expression('Net primary productivity '[log-scale]))+
  ylab(expression(italic('F')[ST]))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
         legend.position = c(0.83,0.28),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase');figNPP

datScaled$x <- datScaled$npp.log
m5 <- (lmer(fst ~ x  * phase + (1|phaseNo), data = filter(datScaled, species == "Pseudomys hermannsburgensis")))

sjPlot::tab_model(m5,dv.labels = c('NPP'),
                  title = 'P. hermannsburgensis: Fst', collapse.ci = F,
                  df.method = "satterthwaite", digits = 3)

round(summary(m5)$coefficients,3) %>% 
  cbind(round(summary(m3)$coefficients,3)) %>% 
  as.data.frame

# figures linear --------

figCapsNPP <- grid.arrange(figCaps, 
             figNPP + theme(legend.position = c(0.7,0.3),
                            legend.key.size = unit(0.5, units = 'cm')),
             ncol = 2)

ggsave(filename = './figures/FST/caps_npp_fst_spp.png', figCapsNPP,
       units = "cm", width = 16, height = 14, dpi = 600)

# table simple ---------

coefcaps <- summary(m3)$coefficients[,1]
coefnpp <- summary(m5)$coefficients[,1]

modeltbl <- data.frame(Species = rep(c('P. hermannsburgensis',
                           'S. youngsoni'), c(3,1)),
           Phase = c('Increase', 'Decrease', 'Low', 'Stable'),
           Captures = c(coefcaps[2], 
                        sum(coefcaps[c(2,5)]),
                        sum(coefcaps[c(2,6)]),
                        summary(m3)$coefficients[2]),
           NPP = c(coefnpp[2], 
                   sum(coefnpp[c(2,5)]),
                   sum(coefnpp[c(2,6)]),
                   summary(m5)$coefficients[2])) %>% 
  mutate(Captures = round(Captures, 3),
         NPP = round(NPP, 3),
         Species = ifelse(duplicated(Species), NA, Species)) %>% 
  rename(#`Net Primary Productivity` = NPP,
         `Model Slopes` = Phase) %>% 
  filter(`Model Slopes` != 'Stable')
 
fxPhasemodelslopes <- flextable(modeltbl) %>% 
  italic(j =1) %>% 
  bold(part = 'header') %>%
  border_remove() %>%  
   italic(i = 2:3, j = 3) %>% 
  autofit() %>% 
 #  width(j = 1:4, unit = 'cm', width = c(4.25,2.75,2.25,5)) %>% 
  hline(i = nrow(modeltbl), 
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = grepl('young', modeltbl$Species), 
        border = fp_border(color = "grey40",
                           width = 1, 
                           style = 'dashed'))  %>% 
   compose(
     part = "body", j =2, i = 1,
     value = as_paragraph(('Increase'), as_sub(' boom'))
   )%>% 
   compose(
     part = "body", j =2, i = 2,
     value = as_paragraph(('Decrease'), as_sub(' bust'))
   )%>% 
   compose(
     part = "body", j =2, i = 3,
     value = as_paragraph(('Low'), as_sub(' bust'))
   ) 
  
fxPhasemodelslopes



# multiple regression ----
m7 <- (lmer(fst ~ captures * npp.log + (1|phaseNo), data = filter(dat, species == "Pseudomys hermannsburgensis")))
m8 <- (lmer(fst ~ captures * npp.log + (1|phaseNo), data = filter(dat, species == "Sminthopsis youngsoni")))

sjPlot::tab_model(m7, m8, dv.labels = c('P. hermannsburgensis', 'S. youngsoni'),
                  title = 'Fst', df.method = "satterthwaite",
                  digits = 3,collapse.ci = T)

# multiple regression
m7 <- (lmer(fst ~ captures + (1|phaseNo), data = filter(dat, species == "Pseudomys hermannsburgensis")))
m9 <- (lmer(fst ~ captures * npp.log + (1|phaseNo), data = filter(dat, species == "Pseudomys hermannsburgensis")))
m8 <- (lmer(fst ~  npp.log + (1|phaseNo), data = filter(dat, species == "Pseudomys hermannsburgensis")))


sjPlot::tab_model(m7, m8,m9, # dv.labels = c('P. hermannsburgensis', 'S. youngsoni'),
                  title = 'Fst', df.method = "satterthwaite",
                  digits = 3,collapse.ci = F)

m10 <- (lmer(fst ~ captures + (1|phaseNo), data = filter(dat, species == "Sminthopsis youngsoni")))
m11 <- (lmer(fst ~ npp.log + (1|phaseNo), data = filter(dat, species == "Sminthopsis youngsoni")))
m12 <- (lmer(fst ~ captures * npp.log + (1|phaseNo), data = filter(dat, species == "Sminthopsis youngsoni")))

sjPlot::tab_model(m10, m11, m12, # dv.labels = c('P. hermannsburgensis', 'S. youngsoni'),
                  title = 'Fst', df.method = "satterthwaite",
                  digits = 3,collapse.ci = F)



em.coef.models <-function(m, spp.ph = T, modelname = 'Interaction'){
  
  summary(m)$coefficients %>% 
  as.data.frame() %>% 
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

em.coef.models(m8)



modeltbl2 <- rbind(em.coef.models(m7, T,'Captures'),
                   em.coef.models(m8, T, 'NPP'),
                   em.coef.models(m9, T, 'Interaction'),
                   em.coef.models(m10, F,'Captures'),
                   em.coef.models(m11, F, 'NPP'),
                   em.coef.models(m12, F, 'Interaction')) %>%
  mutate(
         AIC = rep(c(AIC(m7, m8, m9)[,2],AIC(m10, m11, m12)[,2]),
                   c(2,2,4,2,2,4)),
         AIC = round(ifelse(duplicated(AIC), NA, AIC)),
         ) %>% 
  group_by(Species) %>% 
  mutate(deltaAIC = AIC-min(AIC, na.rm = T),
         Species = ifelse(duplicated(Species), NA, Species)) %>% 
  relocate(sig, .after = deltaAIC) %>% 
  dplyr::select(-AIC)

AIC(m7, m8, m9) %>% 
  mutate(delta = c(2.9, 0, 4),#AIC-min(AIC, na.rm = T),
         likelihood = exp(-0.5*delta),
         weight = likelihood/sum(likelihood))
modeltbl2

fxInteractionmodel <- modeltbl2[,-c(ncol(modeltbl2))] %>% 
  flextable() %>% 
  italic(j =1) %>% 
  bold(part = 'header') %>%
  border_remove() %>%  
  bold(j = 4, i = which(modeltbl2$sig)) %>% 
  autofit() %>% 
  width(j = 5:6, unit = 'cm', width = 2) %>% 
  hline(i = nrow(modeltbl2), 
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = which(complete.cases(modeltbl2$Model))[-1]-1,
        border = fp_border(color = "grey40",
                           width = 1, 
                           style = 'dashed')) %>% 
  hline(i = grep('young', modeltbl2$Species)-1, 
        border = fp_border(color = "grey40",
                           width = 2, 
                           style = 'solid'))  %>% 
  compose(
    part = "body", j =3, i = which(modeltbl2$Parameter== 'npp.log'),
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))
  ) %>% 
  compose(
    part = "body", j = 3, i = grep('Intercept', modeltbl2$Parameter),
    value = as_paragraph(('intercept'))
  ) %>% 
  compose(
    part = "body", j =3, i = which(modeltbl2$Parameter == 'captures:npp.log'),
    value = as_paragraph(('captures:log'), as_sub('e'), ('NPP'))
  ) %>% 
  compose(
    part = "header", j =5,
    value = as_paragraph(('R'), as_sup('2'))
  ) %>%
  compose(
    part = "header", j =6,
    value = as_paragraph(('AIC'), as_sub('delta'))
  ) %>%
  colformat_num( j = 6, big.mark = '')
fxInteractionmodel

flextable::save_as_docx(fxInteractionmodel, path = './figures/FST/table2_fst.docx')


# flex table 1 --------

modeltbl1 <- em.coef.models(m3) %>% 
  rbind(c(NA, NA, 'R2', 
          round(MuMIn::r.squaredGLMM(m3), 2)[2], 
          NA, T)) %>% 
  rename(`Captures (95% CI)` = `Estimate (95% CI)`,
         R2captures = R2, sigcap = sig) %>% 
  left_join(em.coef.models(m5)) %>% 
  mutate(Species = 'P. hermannsburgensis',
         Species = ifelse(duplicated(Species), NA, Species)) %>% 
  dplyr::select(-Model, -R2captures, -R2) %>% 
  rename(`NPP (95% CI)` = `Estimate (95% CI)`,
         signpp = sig) %>% 
  relocate(sigcap, .before = signpp) %>% 
  mutate(sigcap = ifelse(sigcap == 'TRUE', T, F),
         signpp = ifelse(signpp == 'TRUE', T, F),
         signpp = ifelse(is.na(signpp), T, signpp))

modeltbl1$`NPP (95% CI)`[nrow(modeltbl1)] <-  round(MuMIn::r.squaredGLMM(m5), 2)[2]

modeltbljoin <- modeltbl
names(modeltbljoin) <- names(modeltbl1[1:4])

modeltbl1join <- rbind(modeltbl1[1:4],c(NA, 'Model Slopes',NA,NA)) %>% 
  rbind(modeltbljoin)

fxPhasemodels <- modeltbl1[,-c(ncol(modeltbl1), ncol(modeltbl1)-1)] %>% 
  flextable() %>% 
  bold(part = 'header') %>%
  italic(j = 1) %>% 
  border_remove() %>%  
  bold(j = 3, i = which(modeltbl1$sigcap)) %>% 
  bold(j = 4, i = which(modeltbl1$signpp)) %>% 
  autofit() %>% 
  hline(i = c(nrow(modeltbl1)),
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = nrow(modeltbl1)-1,
        border = fp_border(color = "grey40",
                           width = 1, 
                           style = 'solid')) %>% 
  compose(
    part = "header", j =4,
    value = as_paragraph(('NPP (95% CI)'), as_sub(' [x]'))
  ) %>%
  compose(
    part = "header", j =3, 
    value = as_paragraph(('Captures (95% CI)'), as_sub(' [x]'))
  ) %>%
  compose(
    part = "body", j = 2, i = grep('Intercept', modeltbl1$Parameter),
    value = as_paragraph(('Intercept'))
  ) %>% 
  compose(
    part = "body", j =2, i = which(modeltbl1$Parameter== 'phasedecrease'),
    value = as_paragraph(('Phase'), as_sub(' Decrease'))
  ) %>% 
  compose(
    part = "body", j =2, i = which(modeltbl1$Parameter== 'phaselow'),
    value = as_paragraph(('Phase'), as_sub(' Low'))
  ) %>% 
  compose(
    part = "body", j =2, i = which(modeltbl1$Parameter== 'x:phasedecrease'),
    value = as_paragraph(('x:Phase'), as_sub(' decrease'))
  ) %>% 
  compose(
    part = "body", j =2, i = which(modeltbl1$Parameter== 'x:phaselow'),
    value = as_paragraph(('x:Phase'), as_sub(' Low'))
  ) %>% 
  compose(
    part = "body", j =2, i = which(modeltbl1$Parameter== 'R2'),
    value = as_paragraph(('Conditional R'), as_sup('2'))
  );fxPhasemodels

fxPhasemodels


flextable::save_as_docx(fxPhasemodels, 
                        fxPhasemodelslopes,
                        path = './figures/FST/table1_fst.docx')


# interaction -----


em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
em.fstequ <- function(ne,m) 1/(4*m*ne+1)




fst <- dat$fst[1:200]+0.0001
fst <- dat$fst[sample(1:length(dat$fst), 200)]+0.0001
N <- rep(seq(25,250, 25), 20)
m <- em.mRate(fst, N)

N <- rep(seq(10,2000, 100), 10)
m <- rep(c(0.00001, 0.00001,0.00002, 0.00003), 50)/0.001
fst<- log(em.fstequ(N, m))

plot(log(m), fst)
plot(N, fst)



mm  <- lm(fst ~ m);mN  <- lm(fst ~ N);mNm <-lm(fst ~ N*m)   


mm  %>% summary
mN  %>% summary
mNm %>% summary

AIC(mN,mm, mNm) %>% 
  mutate(delta = AIC -min(AIC))

