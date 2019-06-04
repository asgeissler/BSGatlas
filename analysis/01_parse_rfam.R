source('analysis/00_load.R')


files <- list(
  conservative = 'Basu_168_ASM904v1.1e-06.out.filtered.ee.bed.gz',
  medium = 'Basu_168_ASM904v1.noscore.1e-06.out.filtered.ee.bed.gz',
  sensetive = 'Basu_168_ASM904v1.noscore.0.001.out.filtered.ee.bed.gz'
)

rfam <- file.path('data-raw/rfam14.1/', files) %>%
  set_names(names(files)) %>%
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

save(file = 'analysis/01_rfam.rda', rfam)
