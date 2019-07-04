# Infer additional isoforms based on the just computed UTRs

# Idea: compute a directed graph 
# Nodes are genes and UTRs
# (+ special case of 'empty' UTRs to model very short UTRs that were filtered out)
# Edges indicate transcriptional direction
# Isoforms including UTR implied new ones are given by maximal paths
# between nodes with only outgoing and only ingoing edges

source('scripts/graph_paths.R')
source('scripts/frame_helpers.R')
source('analysis/00_load.R')

load('analysis/02_merging.rda')
load('analysis/03_tus.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/06_raw.utrs.rda')


# Collect transcriptional units, including genes and internal UTRs
tus %>%
  select(id, genes) %>%
  separate_rows(genes, sep = ';') %>%
  left_join(merging$merged_genes, c('genes' = 'merged_id')) %>%
  select(- merged_name, - type) %>%
  rename(tu = id, node = genes) %>%
  bind_rows(
    UTRs$internal_UTR %>%
      select(tu = boundary, node = id, start, end, strand) %>%
      separate_rows(tu, sep = ';')
  ) %>%
  # sort in chromosome direction
  arrange(tu, start, end) %>%
  # add link to next node 
  group_by(tu) %>%
  mutate(next.node = lead(node)) %>%
  ungroup %>%
  drop_na %>%
  # infer edge direction based on strand
  transmute(
    from = ifelse(strand == '+', node, next.node),
    to = ifelse(strand == '+', next.node, node)
  ) -> tu.edges

# edges of genes to 3' UTR and 5' UTR to genes
UTRs[c("5'UTR", "3'UTR")] %>%
  bind_rows() %>%
  mutate(len = end - start + 1) %>%
  filter(len > 15) -> us

us %>%
  transmute(
    from = ifelse(type == "5'UTR", id, gene),
    to = ifelse(type == "5'UTR", gene, id)
  ) -> utr.edges
# edges of TSS and terminators to UTRs
us %>%
  transmute(
    from = ifelse(type == "5'UTR", boundary, id),
    to = ifelse(type == "5'UTR", id, boundary)
  ) -> tssterm.edges

# add edges for too short UTRs, or those that would have negativel lengths
# (eg overlaps)

raw.utrs %>%
  mutate(utr.len = end.utr - start.utr + 1) %>%
  filter(utr.len <= 15) %>%
  # similiar to UTRs, edge direction depends on TSS or Terminator type
  transmute(
    from = ifelse(bound_type == 'TSS', bound, gene),
    to = ifelse(bound_type == 'TSS', gene, bound)
  ) -> empty.edges


# Challange: some Tus have neither TSS/term not UTR to them
# -> find these and add dummies, such that the paths can find all TUs

tus %>%
  select(id, genes) %>%
  separate_rows(genes, sep = ';') %>%
  left_join(merging$merged_genes, c('genes' = 'merged_id')) %>%
  select(- merged_name, - type) %>%
  group_by(id) %>%
  arrange(start, end) %>%
  # first/last gene in tu by chromosome order
  summarize(
    start = first(genes),
    end = last(genes),
    strand = first(strand)
  ) %>%
  # first/last one by transcription order
  mutate(
    swap = start,
    start = ifelse(strand == '+', start, end),
    end = ifelse(strand == '+', end, swap)
  ) %>%
  select(-swap) -> tu.firstlast.gene

tu.firstlast.gene %>%
  anti_join(
    raw.utrs %>%
      filter(bound_type == 'TSS'),
    c('start' = 'gene')
  ) %>%
  transmute(
    from = paste('dummy-TSS', 1:n()),
    to = start
  ) -> dummy.tss
tu.firstlast.gene %>%
  anti_join(
    raw.utrs %>%
      filter(bound_type == 'terminator'),
    c('end' = 'gene')
  ) %>%
  transmute(
    from = end,
    to = paste('dummy-terminator', 1:n())
  ) -> dummy.term

# !!!!!!!!!!!!
# The paths must be filterd to only keep dummy terminators corresponding
# to the original TU they were designed for
# !!!!!!!!!!!!

# Collect all nodes
bind_rows(
  UTRs %>%
    map(select, id, type) %>%
    bind_rows,
  merging$merged_genes %>%
    transmute(id = merged_id, type = 'a gene'),
  bsg.boundaries %>%
    map2(names(.), ~ transmute(.x, id, type = .y)) %>%
   bind_rows,
  dummy.tss %>%
    transmute(id = from, type = 'dummy-TSS'),
  dummy.term %>%
    transmute(id = to, type = 'dummy-Terminator')
) %>%
  mutate(number = 1:n()) %>%
  select(number, name = id, type) -> nodes


edges.full <- bind_rows(
  tu.edges,
  utr.edges,
  tssterm.edges,
  empty.edges,
  dummy.term,
  dummy.tss
) %>%
  unique 
  # filter(!complete.cases(.))
edges.full %>%
  mutate(row = 1:n()) %>%
  gather('key', 'name', from, to) %>%
  left_join(nodes, 'name') %>%
  select(row, key, number) %>%
  spread(key, number) %>%
  select(from, to) %>%
  drop_na -> edges

iso.graph <- tbl_graph(nodes, edges, directed = TRUE)

paths <- find_paths(iso.graph)

paths %<>%
  group_by(path) %>%
  mutate(transcription.order = 1:n()) %>%
  ungroup

save(paths, file = 'analysis/07_paths.rda')

###############################################################################

paths %>%
  select(path, transcription.order, node = value) %>%
  left_join(select(nodes, node = name, type), 'node') %>%
  count(path, type) %>%
  spread(type, n, fill = 0) %>%
  group_by_at(vars(- path)) %>%
  count %>%
  select(
    TSS, `5'UTR`, `a gene`, internal_UTR,
    `3'UTR`, `terminator`, count = n
  ) %>%
  arrange(desc(count)) -> paths.stat

assertthat::assert_that(
  with(paths.stat, all(TSS <= 1)),
  with(paths.stat, all(terminator <= 1)),
  with(paths.stat, all(`5'UTR` <= 1)),
  with(paths.stat, all(`3'UTR` <= 1)),
  with(paths.stat, all(internal_UTR < `a gene`)),
  with(paths.stat, all(`a gene` >= 1))
)

View(paths.stat)

###############################################################################
# Which paths have which TU?
# How are the operons?
# How many genes have still no transcript
# What are the spans?
# How many are due to the UTRs
# How many paths have the same genes, but just differen 5/3 UTRs


# > nrow(paths)
# [1] 28133

paths %>%
# Which paths have which TU?
  select(path, merged_id = value) %>%
  semi_join(merging$merged_genes, 'merged_id') %>%
  group_by(path) %>%
  summarize(genes = merged_id %>%
              sort %>%
              clean_paste) -> genes.sorted

genes.sorted %>%
  left_join(tus, 'genes') %>%
  select(path, tu = id, src) -> path.tu

# check each path has only one TU
assertthat::are_equal(
  path.tu %>%
    count(path) %>%
    count(n) %>%
    filter(n > 1) %>%
    nrow,
  0
)

path.tu %>%
  drop_na %>%
  count(tu) %>%
  count(n) %>%
  rename('paths per tu' = n, count = nn)
# `paths per tu` count
# 1  1385
# 2   647
# 3    85
# 4   105
# 5     1
# 6    27
# 8     2
# 9     1

path.tu %>%
  drop_na %>%
  select(tu) %>%
  unique %>%
  nrow
# 2253 number of TUs used? is less?

# Why are these tus not used?
tus %>%
  anti_join(path.tu, c('id' = 'tu')) %>%
  select(id, src, genes) %>%
  # nrow
  # 221
  head

# interest <- c('BSGatlas-gene-21', 'BSGatlas-gene-23')
# edges.full %>%
#   filter((from %in% interest) | (to %in% interest)) %>%
#   View
# okay, should have been expected, because some TUs don't have both explicit 
# TSS/terminator
# Revist computation with placeholder for those ->

# How are the operons?
# How many genes have still no transcript
# What are the spans?
# How many are due to the UTRs
# How many paths have the same genes, but just differen 5/3 UTRs