
source('analysis/00_load.R')



load('analysis/01_refseq.rda')
load('analysis/01_bsubcyc.rda')
load('analysis/01_rfam.rda')
load('analysis/01_nicolas.rda')


# Hierarchy
# refseq coding

# bsubcyc coding

# nicolas lit review
# nicolas trusted
# rfam conservative
# 
# refseq noncoding
# bsubcyc noncoding
# 
# rfam medium
# nicolas lower
# 
# rfam sensetive

data <- bind_rows(
  transmute(refseq$coding,
            priority = 0, src = 'refseq coding',
            id = locus, start, end, strand, type),
  
  transmute(bsubcyc$coding,
            priority = 1, src = 'bsubcyc coding',
            id = locus, start, end, strand,
            type = ifelse(
              str_detect(description, 'putative|predicted'),
              'putative',
              'CDS'
            )),
  
  transmute(nicolas$high$their.lit_review,
            priority = 2, src = 'nicolas lit review',
            id = name, start, end, strand),
  transmute(nicolas$high$their.trusted,
            priority = 2, src = 'nicolas trusted',
            id = name, start, end, strand),
  transmute(rfam$conservative,
            priority = 2, src = 'rfam conservative',
            id = paste('row', 1:n()), start, end, strand, type),
  
  transmute(refseq$noncoding,
            priority = 3, src = 'refseq noncoding',
            id = locus, start, end, strand, type),
  transmute(bsubcyc$noncoding,
            priority = 3, src = 'bsubcyc noncoding',
            id = locus, start, end, strand, type),
  
  transmute(rfam$medium,
            priority = 4, src = 'rfam medium',
            id = paste('row', 1:n()), start, end, strand, type),
  transmute(nicolas$lower,
            priority = 4, src = 'nicolas lower',
            id = name, start, end, strand),
  
  transmute(rfam$sensetive,
            priority = 4, src = 'rfam sensetive',
            id = paste('row', 1:n()), start, end, strand, type)
) %>%
  mutate(type = case_when(
    type == 'small' ~ 'sRNA',
    type == 'putative' ~ 'putative-coding',
    type == 'unclear' ~ 'putative-non-coding',
    type == 'snoRNA'  ~ 'putative-non-coding',
    is.na(type)       ~ 'putative-non-coding',
    TRUE ~ type
  )) %>%
  mutate_at(c('start', 'end', 'priority'), as.integer) %>%
  mutate(seqnames = 'dummychr')

data %>%
  count(src, priority, type) %>%
  transmute(foo = sprintf('%s (%s)', src, priority), type, n) %>%
  spread(foo, n) %>%
  View

# Level-wise find the union
data %>%
  filter(priority == 2) %>%
  # head(n=100) %>%
  select(seqnames, start, end, strand, id) %>%
  plyranges::as_granges() %>%
  plyranges::find_overlaps_directed(., .) %>%
  as.tibble %>% arrange(strand, start, desc(end)) %>% View
  
