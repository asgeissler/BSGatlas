# 05_tss-term_map.R until tss.near computation
source('scripts/biostrings.R')

load('data/01_nicolas.rda')
load('data/01_bsubcyc.rda')
load('data/03_dbtbs_xml.rda')
load('data/01_refseq.rda')
genome <- refseq$seq[[1]]

dbtbs_xml$promoter %>%
  filter(sigma_factor != '') %>%
  pull(binding_sequence) %>%
  str_remove_all('[^ACTG]') %>%
  unique %>%
  find_pattern(., genome, 0) 

DBTBS$promoter %<>%
  mutate(clean.seq = str_remove_all(binding_sequence, '[^ACTG]'))

dat <- bind_rows(
  dbtbs$tss %>%
    select(TSS, strand) %>%
    unique %>%
    transmute(
      id = paste('DBTBS TSS', 1:n(), sep = '_'),
      strand,
      start = TSS,
      end = TSS
    )
  dbtbs_xml$promoter %>%
    filter(sigma_factor != '')
    select(TSS, strand) %>%
    unique %>%
    transmute(
      id = paste('DBTBS TSS', 1:n(), sep = '_'),
      strand,
      start = TSS,
      end = TSS
    )
)

list(c('DBTBS'), c('BsubCyc'), c('DBTBS', 'BsubCyc')) %>%
  set_names(., .) %>%
  map(function(i) {
    tss.dat %>%
      select(src, TSS, strand) %>%
      filter(src %in% c(i, interest)) %>%
      mutate(src = ifelse(src == interest, 'up', 'other')) %>%
      unique %>%
      mutate(
        id = 1:n(),
        src_id = sprintf('%s_%s', src, id)
      ) %>%
      mutate(start = TSS, end = TSS) %>%
      select(id = src_id, start, end, strand) -> tss.dat2
    
    tss.dat2 %>%
      distance_matching(., .) %>%
      filter(!antisense) %>%
      filter(str_detect(y, 'up_')) %>%
      filter(!str_detect(x, 'up_')) -> foo
      
    
    seq(10, 50, 5) %>%
    # c(15, 25, 50) %>%
      map(function(j) { 
        list(
          x = j,
          near = foo %>%
            filter(abs(distance) <= j) %>%
            select(x) %>%
            unique %>%
            nrow,
          have.over = foo %>%
            select(x) %>%
            unique %>%
            nrow,
          total = tss.dat2 %>%
            filter(str_detect(id, 'other_')) %>%
            nrow
        )
      }) %>%
      bind_rows %>%
      mutate(
        over.r = near / have.over * 100,
        total.r = near  / total * 100
      )
  }) %>%
  map(select, x, near, total, total.r) -> foo
    
foo %>%
  map2(c('DBTBS', 'BSubCyc', 'Both'),
       ~ mutate(.x, lab = .y)) %>%
  bind_rows %>%
  ggplot(aes(x, total.r, color = lab)) +
  geom_line() +
  geom_point(size = 2)
