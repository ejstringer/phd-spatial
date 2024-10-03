
# run analysis file
n.sample <- read.csv('./output/sampe_sizes.csv')
# summary table ----------------------------------------------
table(fstdata$species, fstdata$phaseNo)

k <- 1:40

sumCaps <- caps %>% 
  group_by(phaseNo) %>%
  summarise(captures = mean(captures, na.rm = T),
            capturesSy = mean(capturesSy, na.rm = T)) %>%
  filter(complete.cases(phaseNo)) %>%
  pivot_longer(cols = captures:capturesSy, names_to = 'species',
               values_to = 'captures') %>%
  mutate(species = ifelse(species == 'captures',
                          'P. hermanns', 'S. youngsoni')) %>%
  arrange(species, phaseNo) %>%
  ungroup %>%
  as.data.frame
fstdata %>% names

fstdata_n <- right_join(fstdata, n.sample, by = c('phaseNo', 
                                                 'species'))
summarytb <- em.period_table(fstdata_n) %>%
  rename_all(str_to_title) %>%
  right_join(sumCaps, by = c('Period' ='phaseNo', 'Species' = 'species')) %>% 
   mutate(Phase = ifelse(is.na(Phase), 'stable', as.character(Phase)),
          Ngrid = ifelse(is.na(Ngrid), 0, Ngrid),
          Npp = ifelse(is.na(Npp), 80.5, Npp)) %>% 
  rename(NPP = Npp, n = N, gPair = Ngrid, Nm = Migrants) %>% 
  #rowwise()  %>% 
  mutate(gPair = ifelse(is.na(gPair), 0, gPair),
         #ngrids = which(k*(k-1)/2==gPair),
         indPairs = n*(n-1)/2,) %>% 
  relocate(Nm, n,indPairs, .before = gPair) %>% 
  dplyr::select(-Nm) %>% 
  
  mutate(Period = factor(Period, 
                         levels = paste0(rep(c('L', 'I', 'D'), 3), 
                                         rep(1:3, each = 3)))) %>% 
  arrange(Species, Period) %>% 
  mutate(Species = ifelse(duplicated(Species),
                          NA, Species),
         captures = round(captures,1)) %>% 
  dplyr::select(-indPairs) %>% 
    relocate(Captures = captures, .after = NPP)

flex_sumtb <- summarytb %>% 
  flextable() %>% autofit %>% 
  italic(j = 1) %>% 
  border_remove() %>% 
  hline(i = c(9), 
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = c(18), 
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  bold(part = 'header') %>% 
  width(j = 7, width = 2.5, unit = 'cm') %>% 
  #bg(i = seq(1,18,2), bg = 'grey97') %>% 
  hline(i = c(3,6,12,15), j = 2:8,
        border = fp_border(color = "grey40", 
                           width = 1.5, style = "dashed")) %>% 
  compose(
    part = "header", j = 8,
    value = as_paragraph(as_b(as_i("F")), as_sub("ST"))
  ) %>% 
  compose(
    part = "header", j = 7,
    value = as_paragraph(as_b("Grid pairs"))
  );flex_sumtb  
  save_as_docx(flex_sumtb, path = './figures/FST/table_S1_summary.docx')


# fst temporal ----------------

fstdata_sum <- fstdata %>% 
  group_by(species, npp, phaseNo, phase)  %>% 
  summarise(fstmedian = median(fst, na.rm = T),
            meanFst = mean(fst, na.rm = T),
            n = n(),
            se = sd(fst, na.rm = T)/sqrt(n),
            meanlower = meanFst - se,#*2,
            meanupper = meanFst + se#*2
  ) %>% 
  filter(complete.cases(fstmedian)) %>% 
  rename(fst = fstmedian)

fst_diff <- ggplot(fstdata_sum, aes(phaseNo, meanFst, 
                                    fill = phase,colour = species,
                                    group = species))+
  theme_classic()+
  geom_line(size = 1)+

  geom_errorbar(aes(ymin = meanlower, ymax = meanupper),
                width = 0)+
  geom_point(aes(size = n), pch = 21, colour = 'black')+
  theme(
    panel.spacing = unit(0.1, "cm", data = NULL),
    panel.border = element_rect(fill = NA),
    strip.background = element_blank()
  ) +
  guides(size = 'none',
         linetype = 'none',
         fill = guide_legend(override.aes = list(size = 4)))+
  #facet_wrap(~species, ncol = 2) +
  ylab(expression('Mean '* italic('F')[ST]))+
  xlab('Sampling period')+
  scale_fill_manual(values = phaseCol,
                    labels = c('Increase', 'Decrease', 'Low', 'Stable'),
                    name = 'Population phase')+
  scale_colour_manual(values = c('pink3', 'turquoise4'),
                      name  = 'Species',
                      labels = c(expression(italic('P. hermanns')),
                                 expression(italic('S. youngsoni'))));fst_diff

ggsave('./figures/FST/figS1_fst_temp_colourful.png', plot = fst_diff,
       units = "cm", width = 16, height = 9, dpi = 600)

