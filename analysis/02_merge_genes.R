# this script merged the gene annotation, as collected in `01_tomerge.rda`

source('analysis/00_load.R')

source('scripts/overlap_matching.R')

load('analysis/01_tomerge.rda')

# Guarantee that the helper separator has no conflicts
assertthat::assert_that(!any(str_detect(data$id, '%')))
assertthat::assert_that(!any(str_detect(data$src, '%')))

# ids are unique
assertthat::are_equal(data %>% nrow,
                      data %>% unique %>% nrow)

# find overlapping pairs
search <- data %>%
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
    # # special case do not merge known CDS and riboswitch
    # NO don't filter that here, only later when extracting the graph
    # ! (paste0(x.type, y.type) %in% c('CDSriboswitch', 'riboswitchCDS'))
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

#############################################################################

foo.help <- function(i) {
  # Replace last space with line break
  str_replace(i, ' ([^[:space:]]*)$', '\n\\1')
}

data %>%
  filter(priority <= 4) %>%
  select(src, priority) %>%
  unique %>%
  arrange(priority) %>%
  pull(src) %>%
  foo.help -> src.ord

over %>%
  mutate(
    `jaccard similarity` = cut(jaccard, seq(0, 1, 0.1), include.lowest = TRUE)
  ) %>%
  separate(x, c('src.x', 'id.x'), sep = '%') %>%
  separate(y, c('src.y', 'id.y'), sep = '%') %>%
  mutate_at('src.x', foo.help) %>%
  mutate_at('src.y', foo.help) %>%
  mutate_at('src.x', fct_relevel, src.ord) %>%
  mutate_at('src.y', fct_relevel, src.ord) %>%
  filter(as.integer(src.x) >= as.integer(src.y)) -> baz

foo <- c("RefSeq\ncoding", "BsubCyc\ncoding")

baz %>%
  filter(src.x %in% foo) %>%
  filter(src.y %in% foo) %>%
  ggplot(aes(x = `jaccard similarity`)) +
  geom_bar() +
  xlab('Jaccard Index') +
  theme_bw(18) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_grid(src.x ~ src.y,
             switch = 'both',
             scales = 'free_y', drop = TRUE) -> A

baz %>%
  filter(!(src.x %in% foo)) %>%
  filter(src.y %in% foo) %>%
  ggplot(aes(x = `jaccard similarity`)) +
  geom_bar() +
  xlab('Jaccard Index') +
  theme_bw(18) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_grid(src.y ~ src.x,
             switch = 'both',
             scales = 'free_y', drop = TRUE) -> B
baz %>%
  filter(!(src.x %in% foo)) %>%
  filter(!(src.y %in% foo)) %>%
  ggplot(aes(x = `jaccard similarity`)) +
  geom_bar() +
  xlab('Jaccard Index') +
  theme_bw(18) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_grid(src.x ~ src.y,
             switch = 'both',
             scales = 'free_y', drop = TRUE) -> C

cowplot::plot_grid(
  cowplot::plot_grid(A, B,
                     nrow = 1,
                     rel_widths = c(1, 2.5),
                     labels = c('(a)', '(b)'), label_size = 24),
  C,
  nrow = 2,
  rel_heights = c(1, 1.5),
  labels = c(NA, '(c)'), label_size = 24
)

ggsave(filename = 'analysis/02_overlaps.pdf',
       width = 4 * 4, height = 5 * 4,
       units = 'in')

#############################################################################


library(tidygraph)

# 1. Identify groups of genes that will be merged
nodes <- transmute(search, n = 1:n(), name = id, priority)
prots <- c('CDS', 'putative-coding')
edges <- over %>%
  # ignore protein ribo overlaps
  filter(!((x.type %in% prots) & (y.type == 'riboswitch'))) %>%
  filter(!((y.type %in% prots) & (x.type == 'riboswitch'))) %>%
  # relaxed contion, RNA-RNA overlap or fully contained
  mutate(
    # ! (paste0(x.type, y.type) %in% c('CDSriboswitch', 'riboswitchCDS'))
    rnarna = ! ((x.type %in% prots) | (y.type %in% prots)),
    relaxed = ((jaccard > 0.5) | str_detect(mode, 'contain')) & rnarna
  ) %>%
  filter((jaccard > 0.8) | relaxed) %>%
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

# grph %>%
#   activate(nodes) %>%
#   filter(group == 1) %>%
#   plot

# grph %>%
#   activate(nodes) %>%
#   as_tibble %>%
#   count(priority, group) %>%
#   spread(priority, n, fill = 0) %>%
#   count(`0`, `1`, `2`, `3`, `4`) %>%
#   arrange(desc(n)) %>%
#   View

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
  mutate(id = paste0('tmp-', 1:n()))

###############################################################################
# Carry over IDs from last version
'data-gff/BSGatlas_v1.0.gff' %>%
  rtracklayer::import.gff3() %>%
  as_tibble() %>%
  filter(type == 'gene') %>%
  select(start, end, strand, id = ID) -> last.ids

merged_coordinates %>%
  select(id, start, end, strand) %>%
  overlap_matching(last.ids) %>%
  filter(!antisense) -> cmp

# Step 1: equal positions
cmp %>%
  filter(mode == 'equal') %>%
  select(x, y) -> eq

# Step 2: Max on rest, excluding form the pool the equal ones
cmp %>%
  anti_join(eq, 'x') %>%
  anti_join(eq, 'y') %>%
  group_by(x) %>%
  top_n(1, jaccard) %>%
  ungroup %>%
  # count(x) %>%
  # count(n)
  # => clear match
  bind_rows(eq) %>%
  select(tmp = x, id = y) -> carry

cmp %>%
  filter(is.na(x)) %>%
  pull(y) -> obsolete
  
# two new ids, start counting at last pos
last.ids$id[nrow(last.ids)] %>%
  strsplit('-') %>%
  map(3) %>%
  unlist %>%
  as.integer -> foo
cmp %>%
  filter(is.na(y)) %>%
  transmute(
    tmp = x,
    id = paste0('BSGatlas-gene-', (foo + 1):(foo + n()))
  ) -> new.entries

new.entries %>%
  bind_rows(carry) -> look

# Guarantee 1:1 lookup
assertthat::are_equal(
  look %>%
    select(tmp) %>%
    unique %>%
    nrow,
  nrow(look)
)
assertthat::are_equal(
  look %>%
    select(id) %>%
    unique %>%
    nrow,
  nrow(look)
)

merged_coordinates %<>%
  left_join(look, c('id' = 'tmp')) %>%
  select(-id) %>%
  rename(merged_name = id.y)
###############################################################################


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
# merge_map %>%
#   select(merged_name, type) %>%
#   unique %>%
#   count(merged_name) %>%
#   count(n) %>%
#   rename(`different types per merged entry` = n, `count` = nn)

merge_map %>%
  select(merged_name, type) %>%
  arrange(merged_name, type) %>%
  unique %>%
  group_by(merged_name) %>%
  summarize(type = invoke(str_c, type, sep = ';')) %>%
  count(type) %>%
  View

# anything containing 'putative' should be discarded
# only 2 problamatic cases are

investigate <- c(
  'asRNA;putative-non-coding;sRNA',
  'asRNA;sRNA',
  'riboswitch;sRNA',
  'putative-non-coding;riboswitch;sRNA'
# resolved:
# 'putative-coding;putative-non-coding;riboswitch'
)
merge_map %>%
  select(merged_name, type) %>%
  arrange(merged_name, type) %>%
  unique %>%
  group_by(merged_name) %>%
  summarize(type = invoke(str_c, type, sep = ';')) %>%
  filter(type %in% investigate) %>%
  left_join(merge_map, 'merged_name') %>%
  View

# Decision: prefer asRNA over sRNA, and prefer sRNA over riboswitch

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
  } else if (identical(sort(i), c('riboswitch', 'sRNA'))) {
    'sRNA'
  } else {
    stopifnot(
      # an ugly way to make the message appear
      assertthat::are_equal(1, length(i)) |
        {print(sprintf('unchanged type: %s', i)) ; FALSE}
    )
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
  select(merged_name, priority, type.old = type) %>%
  left_join(refined_type) %>%
  mutate(type = ifelse(is.na(type), type.old, type)) %>%
  select(- type.old) %>%
  group_by(merged_name) %>%
  top_n(-1, priority) %>%
  ungroup %>%
  select(- priority) %>%
  unique

assertthat::are_equal(
  0,
  gene_types %>%
    count(merged_name) %>%
    filter(n > 1) %>%
    nrow
)

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

###############################################################################
# 6b Additional query to prefer the names from SubtiWiki
load('data/03_subtiwiki.rda')


###############################################################################

merged_src <- merge_map %>%
  select(merged_id = merged_name, src, priority, locus, name, type, start, end,
         strand, src_locus)
  
merging <- list(
  merged_genes = merged_genes,
  merged_src = merged_src
)

save(merging, file = 'analysis/02_merging.rda')

###############################################################################

