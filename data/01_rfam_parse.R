source('analysis/00_load.R')


files <- c('conservative', 'medium', 'relaxed')

rfam <- sprintf('data-raw/rfam14.1/%s.bed', files) %>%
  set_names(files) %>%
  map(read_delim, delim = '\t',
      col_names = c('chr', 'start', 'end', 'desc', 'score', 'strand')) %>%
  map(separate, 
      col = desc,
      into = c('tools', 'type', 'family', 'name', 'evalue', 'score.dbl',
              'str.chr', 'str.start', 'str.end', 'str.strand'),
      sep = '[|!]'
  ) %>%
  map(select, - starts_with('str.')) %>%
  map(mutate_at, c('evalue', 'score.dbl'),
      function (i) {i %>% strsplit('=') %>% map(2) %>% unlist}) %>%
  map(select, - score) %>%
  map(rename, score = score.dbl) %>%
  map(mutate_at, c('evalue', 'score'), as.double) %>%
  map(mutate_at, c('start', 'end'), as.integer)


rfam %<>%
  map(function(i) {
    mutate(i,
           type = case_when(
             startsWith(type, 'Cis-reg') ~ 'riboswitch',
             name == 'Bacteria_large_SRP' ~ 'SRP',
             name == 'tmRNA' ~ 'tmRNA',
             name %in% c('6S', 'FsrA') ~ 'sRNA',
             type == 'Gene_ribozyme' ~ 'ribozyme',
             # I choose to classify that one miRNA homolog from highest 
             # sensetivity as a asRNA aswell
             type %in% c('Gene_antisense', 'Gene_miRNA') ~ 'asRNA',
             type == 'Gene_rRNA' ~ 'rRNA',
             type == 'Gene_tRNA' ~ 'tRNA',
             type == 'Gene_snRNA_snoRNA_CD-box' ~ 'snoRNA',
             type == 'Gene_sRNA' ~ 'sRNA',
             type == 'Intron' ~ 'intron'
           ))
  })

rfam %<>%
  map(select, family, name, type, start, end, strand, evalue, score)

rfam %>%
  map(count, family) -> foo
rfam %>%
  map2(names(.), ~ left_join(.x, foo[[.y]], 'family')) %>%
  map(group_by, family) %>%
  map(arrange, `start`) %>%
  map(mutate, id = ifelse(
    n == 1,
    family,
    paste0(family, '_match_', 1:n(), '_of_', n())
  )) %>%
  map(select, - n) -> rfam

save(file = 'data/01_rfam.rda', rfam)
