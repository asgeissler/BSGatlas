
source('analysis/00_load.R')

source('scripts/overlap_matching.R')



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
# refseq noncoding
# 
# bsubcyc noncoding
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
            id = paste('row', 1:n()), start, end, strand, type),
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
  
  transmute(rfam$medium,
            name,
            priority = 4, src = 'rfam medium',
            id = paste('row', 1:n()), start, end, strand, type),
  transmute(nicolas$lower,
            name,
            priority = 4, src = 'nicolas lower',
            id = name, start, end, strand),
  
  transmute(rfam$sensetive,
            priority = 5, src = 'rfam sensetive',
            name,
            id = paste('row', 1:n()), start, end, strand, type)
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

data %>%
  count(src, priority, type) %>%
  transmute(foo = sprintf('%s (%s)', src, priority), type, n) %>%
  spread(foo, n) %>%
  View

# Merging

# Quarantee that the helper separator has no conflicts
assertthat::assert_that(!any(str_detect(data$id, '%')))
assertthat::assert_that(!any(str_detect(data$src, '%')))

# ids are unique
assertthat::are_equal(data %>% nrow,
                      data %>% unique %>% nrow)

search <- data %>%
  transmute(
    id = paste(src, id, sep = '%'),
    start, end, strand, type, name, priority
  )

over <- overlap_matching(search, search) %>%
  filter(
    # only same sense comparisons
    ! antisense,
    # don't compare identity
    x != y,
    # remove symetry
    #x < y
  )


over %>%
  mutate(
    `jaccard similarity` = cut(jaccard, seq(0, 1, 0.1), include.lowest = TRUE)
  ) %>%
  mutate_at(c('x.priority', 'y.priority'), function(i) {
    i %>% 
      as.character() %>%
      as.factor %>%
      fct_recode(
        'refseq coding' = '0',
        'bsubcyc coding' = '1',
        'conservative ncRNA\n(incl. refseq)' = '2',
        'bsubcyc ncRNA' = '3',
        'mediuam ncRNA' = '4',
        'sensetive rfam' = '5'
    )
  }) %>%
  ggplot(aes(x = `jaccard similarity`, fill = `jaccard similarity`)) +
  # ggsci::scale_fill_ucscgb() +
  scale_fill_brewer(palette = 'RdYlBu', direction = -1) +
  geom_bar() +
  xlab(NULL) +
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  # scale_y_log10() +
  # xlim(0, 1) +
  # geom_vline(xintercept = c(0.7, 0.8, 0.9)) +
  facet_grid(x.priority ~ y.priority, scales = 'free') +
  ggtitle('Comparison confidence levels', 
          'Similarities between each overlapping gene pair (coding and non-coding), identity is ignored')
          

ggsave(filename = 'analysis/02_level-overlaps.pdf',
       width = 25, height = 25, units = 'cm')


over %>%
  filter(jaccard > 0.8) %>%
  mutate(
    `jaccard similarity` = cut(jaccard, seq(0.8, 1, length.out = 10), include.lowest = TRUE)
  ) %>%
  mutate_at(c('x.priority', 'y.priority'), function(i) {
    i %>% 
      as.character() %>%
      as.factor %>%
      fct_recode(
        'refseq coding' = '0',
        'bsubcyc coding' = '1',
        'conservative ncRNA\n(incl. refseq)' = '2',
        'bsubcyc ncRNA' = '3',
        'mediuam ncRNA' = '4',
        'sensetive rfam' = '5'
    )
  }) %>%
  ggplot(aes(x = `jaccard similarity`, fill = `jaccard similarity`)) +
  scale_fill_brewer(palette = 'RdYlBu', direction = -1) +
  geom_bar() +
  xlab(NULL) +
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  facet_grid(x.priority ~ y.priority, scales = 'free') +
  ggtitle('Comparison confidence levels', 
          'Similarities between each overlapping gene pair (coding and non-coding), identity is ignored')

ggsave(filename = 'analysis/02_level-overlaps_high.pdf',
       width = 25, height = 25, units = 'cm')


# Decision, merge if Jaccard is 0.9
over %>%
  filter(jaccard > 0.9) %>%
  mutate(merges = sprintf('priority %s - priority %s', x.priority, y.priority)) %>%
  count(merges) %>%
  separate(merges, c('a', 'b'), sep = ' - ') %>%
  spread(b, n) 

library(tidygraph)

# 1. Identify group of who will be merged
nodes <- transmute(search, n = 1:n(), name = id, priority)
edges <- over %>%
  filter(jaccard > 0.9) %>%
  select(from = x, to = y) %>%
  mutate(row = 1:n()) %>%
  gather('key', 'node', from, to) %>%
  left_join(nodes %>% 
              select(node = name, n),
            'node') %>%
  select(- node) %>%
  spread(key, n) %>%
  select(- row)

grph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(group = group_components())

grph %>%
  activate(nodes) %>%
  filter(group == 1) %>%
  plot

grph %>%
  activate(nodes) %>%
  as_tibble %>%
  count(priority, group) %>%
  spread(priority, n, fill = 0) %>%
  count(`0`, `1`, `2`, `3`, `4`, `5`) %>%
  arrange(desc(n)) %>%
  rename(
    'refseq coding' = `0`,
    'bsubcyc coding' = `1`,
    'conservative ncRNA\n(incl. refseq)' = `2`,
    'bsubcyc ncRNA' = `3`,
    'mediuam ncRNA' = `4`,
    'sensetive rfam' = `5`
  ) %>%
  View

merging_groups <- grph %>%
  activate(nodes) %>%
  as_tibble
  

# 2. Identify highest priority gene(s)
top_prio <- merging_groups %>%
  group_by(group) %>%
  # !top_n provides more if there is a tie
  # negative to get numerically smalles
  top_n(-1, priority) %>%
  ungroup
  
  
# 3. Union if there are multiples
todo <- top_prio %>%
  left_join(search, c('name' = 'id', 'priority'))

# guarantee there is no strand confusion
assertthat::are_equal(
  0,
  todo %>%
    select(group, strand) %>%
    unique %>%
    count(group) %>%
    filter(n > 1) %>%
    nrow
)

merged_coordinates <- todo %>%
  group_by(group) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = first(strand),
    name = invoke(str_c, name.y, sep = ';')
  ) %>%
  arrange(start, desc(end)) %>%
  mutate(merged_name = paste0('BSGatlas-gene-', 1:n()))


# 4. clear map who participated to a merged gene
merge_map <- merged_coordinates %>%
  select(merged_name, group) %>%
  left_join(merging_groups, 'group') %>%
  select(merged_name, src_locus = name) %>%
  left_join(search, c('src_locus' = 'id')) %>%
  separate(src_locus, c('src', 'locus'), sep = '%', remove = FALSE)

# get rid of graph group id
merged_coordinates %<>% select(- group)

# 5. clarify type of merged one
merge_map %>%
  select(merged_name, type) %>%
  unique %>%
  count(merged_name) %>%
  count(n) %>%
  rename(`different types per merged entry` = n, `count` = nn)
# up to 3 (twice) mostly only 1 (coding) and 2 types (~1200)

merge_map %>%
  select(merged_name, type) %>%
  arrange(merged_name, type) %>%
  unique %>%
  group_by(merged_name) %>%
  summarize(type = invoke(str_c, type, sep = ';')) %>%
  count(type) %>%
  View

# anything containing 'putative' should be discarded
# only problamatic cases are

# investigate <- c(
#   'asRNA;putative-non-coding;sRNA',
#   'asRNA;sRNA',
# resolved:
  ## 'putative-coding;putative-non-coding;riboswitch'
# )
# merge_map %>%
#   select(merged_name, type) %>%
#   arrange(merged_name, type) %>%
#   unique %>%
#   group_by(merged_name) %>%
#   summarize(type = invoke(str_c, type, sep = ';')) %>%
#   filter(type %in% investigate) %>%
#   left_join(merge_map, 'merged_name') %>%
#   View

# Decision: prefer asRNA over sRNA

# merge_map %>%
#   select(merged_name, type) %>%
#   filter(!str_detect(type, 'putative')) %>%
#   unique -> foo
# foo %>%
#   count(merged_name) %>%
#   filter(n > 1) %>%
#   select( - n) %>%
#   left_join(foo) 

type.helper <- function(i) {
  if (identical(sort(i), c('asRNA', 'sRNA'))) {
    'sRNA'
  } else {
    assertthat::are_equal(1, length(i))
    i
  }
}

refined_type <- merge_map %>%
  select(merged_name, type) %>%
  filter(!str_detect(type, 'putative')) %>%
  unique %>%
  group_by(merged_name) %>%
  summarize_at('type', type.helper) %>%
  ungroup

assertthat::are_equal(
  0,
  count(refined_type, merged_name) %>%
    filter(n > 1) %>%
    nrow
)

gene_types <- merge_map %>%
  select(merged_name, type.old = type) %>%
  left_join(refined_type) %>%
  mutate(type = ifelse(is.na(type), type.old, type)) %>%
  select(- type.old) %>%
  unique

# 6. get one gene name, here from the highest priority

merged_genes <- merged_coordinates %>%
  mutate(name = str_remove(name, 'tRNA;')) %>%
  # filter(str_detect(name, ';')) %>% View
  mutate(name = name %>%
           str_replace_all('_', '-') %>%
           str_replace_all(';', '_')
    ) %>%
  # and unify with types from last step
  left_join(gene_types) %>%
  select(merged_id = merged_name, merged_name = name, type, start, end, strand)

merged_src <- merge_map %>%
  select(merged_id = merged_name, src, priority, locus, name, type, start, end,
         strand, src_locus)
  
merging <- list(
  merged_genes = merged_genes,
  merged_src = merged_src
)

save(merging, file = 'analysis/02_merging.rda')

# 7. Make statistics
prots <- c('CDS', 'putative-coding')
stat <- merged_src %>%
  left_join(merged_genes, 'merged_id', suffix = c('', '.merged')) %>%
  mutate(type.equal  = type == type.merged,
         # indicate of change from putative-coding to RNA classification
         type.reclassRNA = (type %in% prots) & !(type.merged %in% prots),
         # and the other way
         type.reclassProt = !(type %in% prots) & (type.merged %in% prots),
         coord.equal = (start == start.merged) & (end == end.merged)) %>%
  transmute(merged_id, src, priority, type.equal, coord.equal,
            type.reclassRNA, type.reclassProt,
            is.coding = type.merged %in% prots,
           `gene type` = ifelse(is.coding, 'coding', 'non-coding')) %>%
  mutate(src = case_when(
    startsWith(src, 'bsubcyc') ~ 'BsubCyc',
    startsWith(src, 'refseq') ~ 'RefSeq',
    startsWith(src, 'rfam') ~ 'rfam',
    src %in% c("nicolas trusted", "nicolas lit review") ~ 'literature review',
    src == 'nicolas lower' ~ 'Nicolas et al predictions'
  ))

p.src_count <- stat %>%
  count(merged_id, `gene type`) %>%
  count(n, `gene type`) %>%
  mutate_at(c('gene type', 'n'), as.factor) %>%
  right_join(expand(., `gene type`, `n`)) %>% 
  replace_na(list(nn = 0)) %>%
  ggplot(aes(x = n, y = nn,
             fill =  `gene type`,
             group = `gene type`)) +
  ggsci::scale_fill_jama(name = 'Gene Type') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  xlab('Number of sources contributing to a merged gene') +
  ylab('Count of occurences')

p.type_count <- merged_genes %>%
  mutate(type = case_when(
    type %in% c('riboswitch', 'CDS', 'putative-coding',
              'putative-non-coding') ~ type,
    TRUE ~ 'non-coding with known type (tRNA, sRNA, ...)'
  ) %>%
    fct_infreq) %>%
  pull(type) %>%
  fct_count(prop = TRUE) %>% mutate(cs = rev(cumsum(rev(p)))) %>%
  rename(type = f) %>%
  mutate(
    lab = round(p * 100),
    lab = ifelse(lab >= 0,
                 sprintf('%s %%', as.character(lab)),
                 ''),
    at = cs * sum(n) - n / 2
  ) %>%
  ggplot(aes(x = '', y = n, fill = type)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama(name = NULL) +
  coord_polar(theta = 'y') +
  geom_label(aes(y = at,# + c(0, cumsum(n)[-length(n)]), 
                x = c(1, 1, 1.3, 1.05, 0.8),
                label = lab), size=3,
            color = 'white') +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text.x=element_blank()
  ) +
  guides(fill = guide_legend(ncol = 2))
p.type_count


cowplot::plot_grid(
  p.type_count, p.src_count,
  labels = 'auto'
) %>% ggsave(filename = 'analysis/02_merge_grid.pdf',
             width = 30, height = 10, units = 'cm')

# how often is a gene specific in a source
stat.src_uniq <- stat %>%
  count(merged_id, `gene type`) %>%
  filter(n == 1) %>%
  left_join(stat) %>%
  count(`gene type`, src) %>%
  rename(uniq = n)

# how often has source type the same as merged
stat.src_type <- stat %>%
  count(src, `gene type`, type.equal) %>%
  mutate(type.equal = ifelse(type.equal, 'same type', 'diff type')) %>%
  spread(type.equal, n, fill = 0)

# total number of genes per resource
stat.src_total <- stat %>%
  count(src, `gene type`) %>%
  rename(Total = n)
  

# assert: both reclass can't happen at the same time
assertthat::are_equal(
  0, stat %>% filter(type.reclassProt, type.reclassRNA) %>% nrow
)
stat.src_reclass <- stat %>%
  count(src, type.reclassRNA, type.reclassProt) %>%
  filter(type.reclassProt | type.reclassRNA) %>%
  gather(... = c(type.reclassProt, type.reclassRNA)) %>%
  filter(value) %>%
  select(- value) %>%
  mutate(key = ifelse(key == 'type.reclassProt',
                      'non-coding%putative ncRNA reclassified as coding',
                      'coding%putative coding reclassified as non-coding')) %>%
  separate(key, c('gene type', 'desc'), sep = '%') %>%
  spread(desc, n)


# how often has source coord the same as merged
stat.src_coord <- stat %>%
  count(src, `gene type`, coord.equal) %>%
  mutate(coord.equal = ifelse(coord.equal, 'same coord', 'diff coord')) %>%
  spread(coord.equal, n, fill = 0)

# general providence stat
stat.src_genes <- merged_src %>%
  mutate(src = case_when(
    startsWith(src, 'bsubcyc') ~ 'BsubCyc',
    startsWith(src, 'refseq') ~ 'RefSeq',
    startsWith(src, 'rfam') ~ 'rfam',
    src %in% c("nicolas trusted", "nicolas lit review") ~ 'literature review',
    src == 'nicolas lower' ~ 'Nicolas et al predictions'
  )) %>%
  select(merged_id, src, type)
  
# as here is a multi source aggregation, I need to work on the type
foo <- stat.src_genes %>%
  select(merged_id, src, type) %>%
  filter(!str_detect(type, 'putative')) %>%
  unique %>%
  group_by(merged_id, src) %>%
  summarize_at('type', type.helper) %>%
  ungroup
  
# now the actual stat
stat.src_genes %<>%
  left_join(foo, c('merged_id', 'src')) %>%
  mutate(type = ifelse(is.na(type.y), type.x, type.y)) %>%
  count(type, src) %>%
  rename(count = n, `gene type` = type)


# summary <-
tmp <- stat.src_coord %>%
  left_join(stat.src_type,  c('src', 'gene type')) %>%
  left_join(stat.src_total,  c('src', 'gene type')) %>%
  left_join(stat.src_uniq,  c('src', 'gene type')) %>%
  left_join(stat.src_reclass,  c('src', 'gene type'))


# make sure numbers add up before continuing
assertthat::assert_that(with(tmp, `diff coord` + `same coord` == Total) %>% all)
assertthat::assert_that(with(tmp, `diff type` + `same type` == Total) %>% all)

# The coding part
coding.part <- tmp %>%
  filter(`gene type` == 'coding') %>%
  left_join(
    # stat on who is putative
    stat.src_genes %>%
      transmute(src, count, specific = `gene type`,
                `gene type` = ifelse(`gene type` %in% prots, 'coding', 'non-coding')) %>%
      filter(`gene type` == 'coding') %>%
      spread(specific, count)
  )
assertthat::assert_that(with(coding.part, `CDS` + `putative-coding` == Total) %>%
                          all(na.rm = TRUE))
  
pct.helper <-function(i, total) {
  sprintf('%s (%s %%)',
          i,
          round(i / total * 100) %>% as.character
  )
}
coding.part2 <- coding.part %>%
  replace_na(list(
    'putative-coding' = 0
  )) %>%
  # special case in which a putative ncRNA has been re.classified as coding
  filter(src != 'Nicolas et al predictions') %>%
  transmute(
    src,
    'Protein Coding Genes' = Total,
    'of these hypothetical' = pct.helper(`putative-coding`, Total),
    'Merging refined coordinates of' = pct.helper(`diff coord`, Total),
    'Merging removed hypothetical status for' = pct.helper(`diff type`, Total),
    'Genes uniq to this resource' = ifelse(is.na(uniq),
                                           '-',
                                           pct.helper(`uniq`, Total))
  ) %>%
  # transpose
  gather(... = - src) %>%
  spread(src, value) %>%
  arrange(desc(1:n()))
  
  
  
  
  

# The non-coding part
nc.part <- tmp %>%
  filter(`gene type` == 'non-coding') %>%
  left_join(
    # stat on who is putative
    stat.src_genes %>%
      transmute(src, count, specific = `gene type`,
                `gene type` = ifelse(`gene type` %in% prots, 'coding', 'non-coding')) %>%
      filter(`gene type` == 'non-coding') %>%
      spread(specific, count)
  )
assertthat::assert_that(with(
  nc.part, pmap(list(
    `putative-non-coding`, `putative ncRNA reclassified as coding`,
    `asRNA`, `intron`, `putative-non-coding`, `riboswitch`, `ribozyme`,
    `rRNA`, `sRNA`, `SRP`, `tmRNA`, `tRNA` 
  ),  sum
  ) == Total) %>%
    all(na.rm = TRUE)
)
  
nc.part2 <- nc.part %>%
  mutate_at(c(
    'putative-non-coding', 'asRNA', 'intron', 'riboswitch', 'ribozyme', 'rRNA',
    'sRNA', 'SRP', 'tmRNA', 'tRNA'
    ),
    replace_na, 0
  ) %>%
  transmute(
    src,
    'Non-Coding Structures and Genes' = Total,
    'of these predecited/hypothetical' = pct.helper(`putative-non-coding`, Total),
    'specific ncRNA type known' = pmap(list(
        `asRNA`, `intron`, `riboswitch`, `ribozyme`,
        `rRNA`, `sRNA`, `SRP`, `tmRNA`, `tRNA` 
      ),  sum) %>%
      unlist %>%
      pct.helper(Total),
    'ribosomal RNA (rRNA)' = pct.helper(`rRNA`, Total),
    'transfer RNA (tRNA)' = pct.helper(`tRNA`, Total),
    'small regulatory RNA (sRNA)' = pct.helper(`sRNA`, Total),
    'regulatory antisense RNA (asRNA)' = pct.helper(`asRNA`, Total),
    'riboswitch' = pct.helper(`riboswitch`, Total),
    'self-splicing intron' = pct.helper(intron, Total),
    'other (ribozyme, SRP, tmRNA)' = pct.helper(ribozyme + SRP + tmRNA,
                                                Total),
    
    
    'Merging refined coordinates of' = pct.helper(`diff coord`, Total),
    'Merging removed hypothetical status for' = pct.helper(`diff type`, Total),
    'Merging reclassified putative ncRNA as coding' = ifelse(
     is.na(`putative ncRNA reclassified as coding`),
     '-',
     `putative ncRNA reclassified as coding`
    ),
    'Genes uniq to this resource' = ifelse(is.na(uniq),
                                           '-',
                                           pct.helper(`uniq`, Total))
  )

# transpose and sort columns correctly
nc.part3 <- nc.part2 %>%
  gather(... = - src) %>%
  spread(src, value) %>%
  mutate_at('key', fct_relevel,
            names(nc.part2) %>%
              discard(`==`, 'src') %>%
              unlist) %>%
  arrange(key) %>%
  mutate_at('key', as.character)
  

# stat for the final product
foo <- merged_genes %>%
  select(type) %>%
  mutate(gene.type = ifelse(type %in% prots, 'total.coding', 'total.noncoding'))

stat.bsgatlas <- foo %>%
  count(gene.type) %>%
  set_names(c('type', 'n')) %>%
  bind_rows(count(foo, type)) %>%
  spread(type, n) %>%
  transmute(
    src = 'BSGatlas',
    'Protein Coding Genes' = total.coding,
    'of these hypothetical' = pct.helper(`putative-non-coding`, total.coding),
    
    'Non-Coding Structures and Genes' = total.noncoding,
    'of these predecited/hypothetical' = pct.helper(`putative-non-coding`, total.noncoding),
    'ribosomal RNA (rRNA)' = pct.helper(`rRNA`, total.noncoding),
    'transfer RNA (tRNA)' = pct.helper(`tRNA`, total.noncoding),
    'small regulatory RNA (sRNA)' = pct.helper(`sRNA`, total.noncoding),
    'regulatory antisense RNA (asRNA)' = pct.helper(`asRNA`, total.noncoding),
    'riboswitch' = pct.helper(`riboswitch`, total.noncoding),
    'self-splicing intron' = pct.helper(intron, total.noncoding),
    'other (ribozyme, SRP, tmRNA)' = pct.helper(ribozyme + SRP + tmRNA,
                                                total.noncoding)
  ) %>%
  gather(... = -src) %>%
  spread(src, value) %>%
  rename(description = key)

# Show full table
bind_rows(coding.part2, nc.part3) %>%
  mutate_all(replace_na, '-') %>%
  # names
  select(
    description = key,
    RefSeq, BsubCyc, 
    `Literature Review` = `literature review`,
    `rfam (various sensetivity levels)` = `rfam`,
    `Nicolas et al predictions`
  ) %>%
  left_join(stat.bsgatlas) %>%
  replace_na(list(BSGatlas = '-')) %>%
  kable('latex', booktabs = TRUE) %>%
  add_indent(c(2, 7, 9:15)) %>%
  row_spec(5, hline_after = TRUE) %>%
  kable_styling() -> tbl_tex

as_image(tbl_tex)

tbl_tex %>%
  str_split('\\n') %>%
  unlist %>%
  discard(str_detect, pattern = '\\{table\\}') %>%
  write_lines(path = 'analysis/02_stat.tex')

