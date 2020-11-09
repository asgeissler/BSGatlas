source('analysis/00_load.R')


load('data/01_refseq.rda')
load('data/01_bsubcyc.rda')
load('data/01_rfam.rda')
load('data/01_nicolas.rda')
load('data/01_dar-riboswitches.rda')

# Script collects the individual gene annotation resources
# that should be merged and make a summary table

# Hierarchy
# refseq coding

# bsubcyc coding

# rfam conservative
# refseq noncoding
# 
# bsubcyc noncoding
# dar term seq
# 
# rfam medium
# nicolas all
# 
# rfam sensetive

data <- bind_rows(
  transmute(refseq$coding,
            priority = 0, src = 'RefSeq coding',
            name,
            id = locus, start, end, strand, type) %>%
    mutate(type = ifelse(type == 'putative', 'putative-coding', type)),
  
  transmute(bsubcyc$coding,
            priority = 1, src = 'BsubCyc coding',
            id = locus, start, end, strand,
            name,
            type = ifelse(
              str_detect(description, 'putative|predicted'),
              'putative-coding',
              'CDS'
            )),
  
  transmute(rfam$conservative,
            name,
            priority = 2, src = 'Rfam conservative',
            id, start, end, strand, type),
  transmute(refseq$noncoding,
            name,
            priority = 2, src = 'RefSeq non-coding',
            id = locus, start, end, strand, type) %>%
    mutate(type = ifelse(type %in% c('unclear', 'putative'),
                         'putative-non-coding',
                         type)),
  
  transmute(bsubcyc$noncoding,
            name,
            priority = 3, src = 'BsubCyc non-coding',
            id = locus, start, end, strand, type),
  transmute(dar_riboswitches,
            name,
            priority = 3, src = 'Dar et al. riboswitches',
            id = paste0('TERMseq-predicted-riboswitch_', 1:n()),
            start, end, strand, type) %>%
    mutate(name = ifelse(name == 'novel', 'TERMseq-predicted-riboswitch', name)),
  
  transmute(rfam$medium,
            name,
            priority = 4, src = 'Rfam medium',
            id , start, end, strand, type),
  nicolas$all.features %>%
    filter(startsWith(type, 'indep')) %>%
    transmute(name, priority = 4, src = 'Nicolas et al. predicted',
              id = locus, start, end, strand),
  
  transmute(rfam$relaxed,
            priority = 5, src = 'Rfam relaxed',
            name,
            id, start, end, strand, type)
) %>%
  mutate(type = case_when(
    type == 'small' ~ 'sRNA',
    type == 'unclear' ~ 'putative-non-coding',
    type == 'snoRNA'  ~ 'putative-non-coding',
    is.na(type)       ~ 'putative-non-coding',
    TRUE ~ type
  )) %>%
  mutate_at(c('start', 'end', 'priority'), as.integer) %>%
  mutate(seqnames = 'dummychr')

save(data, file = 'analysis/01_tomerge.rda')

