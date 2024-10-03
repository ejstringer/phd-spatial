
em.period_table <- function(fstdata){ 
period_summaries <- fstdata %>%
  group_by(species, phaseNo, phase, npp, n, ngrids) %>% 
  summarise(nPairs = n(),
            sd = sd(fst),
            fst = mean(fst),
            se = sd/sqrt(nPairs),
            lowerf = round(fst - se*1.96,3),
            upperf = round(fst + se*1.96,3),
            fstsum = paste0(round(fst, 3), " (", lowerf, ' - ' ,upperf, ")" ),
            fstsum = ifelse(is.na(se), as.character(round(fst,3)), fstsum)) %>%
  dplyr::select(-lowerf, -upperf) %>% 
  ungroup() %>% 
  mutate(npp.log = log(npp)) %>% 
  mutate(#mRate = (1/ifelse(fst > 0, fst, NA)-1)/(4*(ne)),
         nMigrants = (1/ifelse(fst > 0, fst, NA)-1)/4,
         #mRate = round(mRate, 3),
         nMigrants = round(nMigrants),
         npp = round(npp,1),
         nPairs = ifelse(is.na(ngrids), NA, nPairs)) %>% 
  arrange(species, phaseNo)  %>% 
  mutate(fst = round(fst, 3)) %>% 
  dplyr::select(species,phaseNo,  phase, npp, nPairs,
                #fst,
                fstsum,
                n,
                nMigrants) ; period_summaries


periodft <- period_summaries %>% 
  mutate(species = ifelse(grepl('hermann', species),
                          'P. hermanns',
                          'S. youngsoni'))

colnames(periodft) <- c('species', 'period','phase',  'NPP', 
                        'ngrid', 'FST', 'n', 'migrants')

return(periodft)
}
