
source('analysis/00_load.R')

source('scripts/overlap_matching.R')

load('analysis/01_tomerge.rda')

data <- to.merge$todo

# Quarantee that the helper separator has no conflicts
assertthat::assert_that(!any(str_detect(data$id, '%')))
assertthat::assert_that(!any(str_detect(data$src, '%')))

# ids are unique
assertthat::are_equal(data %>% nrow,
                      data %>% unique %>% nrow)

search <- data %>%
  # ignore senseitve rfam
  filter(priority <= 4) %>%
  transmute(
    id = paste(src, id, sep = '%'),
    start, end, strand, type, name, priority
  )

over <- overlap_matching(search, search) %>%
  filter(
    # only same sense comparisons
    ! antisense,
    # don't compare identity
    x != y
  )

# Jaccard similarities when the loci match
# load('analysis/01_refseq.rda')
# bind_rows(
#   refseq$coding %>%
#     transmute(
#       general.type = 'coding',
#       refseq.locus = locus,
#       bsubcyc.locus = old_locus
#     ),
#   refseq$noncoding %>%
#     transmute(
#       general.type = 'noncoding',
#       refseq.locus = locus,
#       bsubcyc.locus = old_locus
#     )
# ) %>%
#   mutate(bsubcyc.locus = ifelse(is.na(bsubcyc.locus),
#                                 refseq.locus,
#                                 bsubcyc.locus),
#          refseq.locus = sprintf('refseq %s%%%s', general.type, refseq.locus),
#          bsubcyc.locus = sprintf('bsubcyc %s%%%s', general.type, bsubcyc.locus)
#          ) %>%
#   inner_join(over,
#              c('refseq.locus' = 'x', 'bsubcyc.locus' = 'y')) %>%
#   # pull(jaccard) %>%
#   # summary
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.4536  1.0000  1.0000  0.9989  1.0000  1.0000
# # On average really high similarity
#   # filter(jaccard < 0.9) %>%
#   # nrow
# # only 11 cases below 0.9
#   filter(jaccard < 0.8) %>%
#   nrow
# # only 5 cases below 0.8


over %>%
  mutate(
    `jaccard similarity` = cut(jaccard, seq(0, 1, 0.1), include.lowest = TRUE)
  ) %>%
  mutate_at(c('x.priority', 'y.priority'), function(i) {
    i %>% 
      as.character() %>%
      as.factor %>%
      fct_recode(
        'RefSeq coding' = '0',
        'BsubCyc coding' = '1',
        'conservative ncRNA set\n(incl. refseq)' = '2',
        'BsubCyc ncRNA\nterm-seq riboswitches' = '3',
        'medium rfam\nRemaining Nicolas predictions' = '4',
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
        'RefSeq coding' = '0',
        'BsubCyc coding' = '1',
        'conservative ncRNA set\n(incl. refseq)' = '2',
        'BsubCyc ncRNA\nterm-seq riboswitches' = '3',
        'medium rfam\nRemaining Nicolas predictions' = '4',
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


# Decision, merge if Jaccard is 0.8
over %>%
  filter(jaccard > 0.8) %>%
  mutate(merges = sprintf('priority %s - priority %s', x.priority, y.priority)) %>%
  count(merges) %>%
  separate(merges, c('a', 'b'), sep = ' - ') %>%
  spread(b, n) 

library(tidygraph)

# 1. Identify group of who will be merged
nodes <- transmute(search, n = 1:n(), name = id, priority)
edges <- over %>%
  filter(jaccard > 0.8) %>%
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
    'asRNA'
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
