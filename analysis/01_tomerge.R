source('analysis/00_load.R')


load('analysis/01_refseq.rda')
load('analysis/01_bsubcyc.rda')
load('analysis/01_rfam.rda')
load('analysis/01_nicolas.rda')
load('analysis/01_dar-riboswitches.rda')


# Hierarchy
# refseq coding

# bsubcyc coding

# nicolas lit review
# nicolas trusted
# rfam conservative
# refseq noncoding
# 
# bsubcyc noncoding
# dar term seq
# 
# rfam medium
# nicolas lower
# 
# rfam sensetive

data <- bind_rows(
  transmute(refseq$coding,
            priority = 0, src = 'refseq coding',
            name,
            id = locus, start, end, strand, type) %>%
    mutate(type = ifelse(type == 'putative', 'putative-coding', type)),
  
  transmute(bsubcyc$coding,
            priority = 1, src = 'bsubcyc coding',
            id = locus, start, end, strand,
            name,
            type = ifelse(
              str_detect(description, 'putative|predicted'),
              'putative-coding',
              'CDS'
            )),
  
  transmute(nicolas$high$their.lit_review,
            priority = 2, src = 'nicolas lit review',
            name,
            id = name, start, end, strand),
  transmute(nicolas$high$their.trusted,
            priority = 2, src = 'nicolas trusted',
            name,
            id = name, start, end, strand),
  transmute(rfam$conservative,
            name,
            priority = 2, src = 'rfam conservative',
            id, start, end, strand, type),
  transmute(refseq$noncoding,
            name,
            priority = 2, src = 'refseq noncoding',
            id = locus, start, end, strand, type) %>%
    mutate(type = ifelse(type %in% c('unclear', 'putative'),
                         'putative-non-coding',
                         type)),
  
  transmute(bsubcyc$noncoding,
            name,
            priority = 3, src = 'bsubcyc noncoding',
            id = locus, start, end, strand, type),
  transmute(dar_riboswitches,
            name,
            priority = 3, src = 'dar riboswitches',
            id, start, end, strand, type),
  
  transmute(rfam$medium,
            name,
            priority = 4, src = 'rfam medium',
            id = paste('row', 1:n()), start, end, strand, type),
  nicolas$lower %>%
    filter(startsWith(type, 'indep')) %>%
    transmute(name, priority = 4, src = 'nicolas lower',
              id = locus, start, end, strand),
  
  transmute(rfam$sensetive,
            priority = 5, src = 'rfam sensetive',
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

tab <- data %>%
  mutate(src = case_when(
    src == "bsubcyc coding"     ~ 'BsubCyc',
    src == "bsubcyc noncoding"  ~ 'BsubCyc',
    src == "dar riboswitches"   ~ 'Dar et al. term-seq',
    src == "nicolas lit review" ~ 'Nicolas et al.\'s literature review',
    src == "nicolas lower"      ~ 'Nicolas et al. predictions',
    src == "nicolas trusted"    ~ 'Nicolas et al. trusted predictions',
    src == "refseq coding"      ~ 'RefSeq',
    src == "refseq noncoding"   ~ 'RefSeq',
    src == "rfam conservative"  ~ 'rfam (conservative)',
    src == "rfam medium"        ~ 'rfam (medium)',
    src == "rfam sensetive"     ~ 'rfam (sensetive)'
  )) %>%
  count(src, type) %>%
  spread(src, n) %>%
  mutate_all(replace_na, '-')

total <- tab %>%
  mutate_all(as.numeric) %>%
  map(sum, na.rm = TRUE)  %>%
  as.tibble %>%
  mutate(type = 'Total') %>%
  mutate_all(as.character)

tab %<>% bind_rows(total)


to.merge <- list(
  todo = data,
  stat = tab
)

save(to.merge, file = 'analysis/01_tomerge.rda')
