
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

helper <- function(i) {
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
  summarize_at('type', helper) %>%
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
  merged_strc = merged_src
)

# 7. Make statistics
stat <- merged_src %>%
  left_join(merged_genes, 'merged_id', suffix = c('', '.merged')) %>%
  mutate(type.equal  = type == type.merged,
         coord.equal = (start == start.merged) & (end == end.merged)) %>%
  transmute(merged_id, src, priority, type.equal,
           is.coding = type.merged %in% c('CDS', 'putative-non-coding'))

p.src_count <- stat %>%
  count(merged_id, is.coding) %>%
  count(n, is.coding) %>%
  mutate(`gene type` = ifelse(is.coding, 'coding', 'non-coding')) %>%
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
    str_detect(type, 'putative') ~ type,
    type == c('riboswitch', 'CDS') ~ type,
    TRUE ~ 'non-coding with known type (tRNA, sRNA, ...)'
  ) %>%
    fct_infreq) %>%
  count(type) %>%
  mutate(lab = round(n / nrow(merged_genes) * 100),
         lab = ifelse(lab > 10,
                      sprintf('%s %%', as.character(lab)),
                      '')) %>%
  ggplot(aes(x = '', y = n, fill = type)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama(name = NULL) +
  coord_polar('y') +
  geom_text(aes(y = n/2.5 + c(0, cumsum(n)[-length(n)]), 
                label = lab), size=5,
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
