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
    to = start,
    tu = id
  ) -> dummy.tss
tu.firstlast.gene %>%
  anti_join(
    raw.utrs %>%
      filter(bound_type == 'terminator'),
    c('end' = 'gene')
  ) %>%
  transmute(
    from = end,
    to = paste('dummy-terminator', 1:n()),
    tu = id
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
  dummy.term %>% select(- tu),
  dummy.tss %>% select(- tu)
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
# Which paths have which TU?

# !!!! Also filter erraneous dummy matching !!!!

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
  select(path, tu = id, src) -> paths.tu

# Make blacklist, these are TUs that contain a dummy not ment for them
paths %>%
  filter(str_detect(value, 'dummy')) %>%
  left_join(paths.tu, 'path') %>%
  # accept correct TSS map
  anti_join(dummy.tss, c('value' = 'from', 'tu')) %>%
  # accept correct Term map
  anti_join(dummy.term, c('value' = 'to', 'tu')) %>%
  select(path) %>%
  unique -> blacklist

# > nrow(blacklist)
# [1] 205

filtered.paths <- paths %>% anti_join(blacklist)
filtered.paths <- paths %>% anti_join(blacklist)

filtered.tu <- paths.tu %>% anti_join(blacklist)

###############################################################################
# How are the operons?
# How many genes have still no transcript
# What are the spans?
# How many are due to the UTRs
# How many paths have the same genes, but just differen 5/3 UTRs

# > nrow(filtered.paths)
# [1] 29597

filtered.paths %>%
  select(path, transcription.order, node = value) %>%
  left_join(select(nodes, node = name, type), 'node') %>%
  count(path, type) %>%
  spread(type, n, fill = 0) %>%
  group_by_at(vars(- path)) %>%
  count %>%
  select(
    `dummy-TSS`,
    TSS, `5'UTR`, `a gene`, internal_UTR,
    `3'UTR`, `terminator`, `dummy-Terminator`,
    count = n
  ) %>%
  arrange(desc(count)) -> paths.stat

# path constraints
assertthat::assert_that(
  with(paths.stat, all(`dummy-TSS` <= 1)),
  with(paths.stat, all(TSS <= 1)),
  with(paths.stat, !(TSS & `dummy-TSS`)) %>% all,
  with(paths.stat, all(terminator <= 1)),
  with(paths.stat, all(`dummy-Terminator` <= 1)),
  with(paths.stat, !(terminator & `dummy-Terminator`)) %>% all,
  with(paths.stat, all(`5'UTR` <= 1)),
  with(paths.stat, all(`3'UTR` <= 1)),
  with(paths.stat, all(internal_UTR < `a gene`)),
  with(paths.stat, all(`a gene` >= 1))
)

View(paths.stat)

# check each path has only one TU
assertthat::are_equal(
  filtered.tu %>%
    count(path) %>%
    count(n) %>%
    filter(n > 1) %>%
    nrow,
  0
)

filtered.tu %>%
  drop_na %>%
  count(tu) %>%
  count(n) %>%
  rename('paths per tu' = n, count = nn)
# `paths per tu` count
# 1  1560
# 2   690
# 3    87
# 4   106
# 5     1
# 6    27
# 8     2
# 9     1

assertthat::are_equal(
  nrow(tus),
  filtered.tu %>%
    drop_na %>%
    select(tu) %>%
    unique %>%
    nrow
)
# found all, dummy stuff worked

filtered.paths %>%
  filter(!str_detect(value, 'dummy')) %>%
  left_join(filtered.tu, 'path') %>%
  select(path, transcription.order, value, path, tu, src) %>%
  arrange(path, transcription.order) -> trans

###############################################################################
# How are the operons?

operon.nodes <- trans %>%
  select(value, path) %>%
  gather('key', 'node', path, value) %>%
  unique %>%
  transmute(i = 1:n(), node)

trans %>%
  select(path, value) %>%
  mutate(row = 1:n()) %>%
  gather('key', 'node', path, value) %>%
  left_join(operon.nodes, 'node') %>%
  select(- node) %>%
  spread(key, i) %>%
  select(from = path, to = value) -> operon.edges

operon.grph <- tbl_graph(nodes = operon.nodes, edges = operon.edges, directed = FALSE)

operon.grph %>%
  activate(nodes) %>%
  mutate(operon = group_components()) %>%
  as_tibble %>%
  mutate(path = as.integer(node)) %>%
  select(operon, path) %>%
  drop_na -> operons


###############################################################################
# Make new TUs, compute spans, and provide nicer names

trans %>%
  inner_join(merging$merged_genes, c('value' = 'merged_id')) %>%
  mutate_at('src', replace_na, 'BSGatlas') %>%
  group_by(path) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand),
    src = clean_paste(src),
    genes =  value %>%
      fct_reorder(transcription.order) %>%
      clean_paste()
  ) %>%
  group_by(genes) %>%
  summarize_all(clean_paste) %>%
  ungroup() -> trans.span

assertthat::assert_that(with(trans.span, all(strand %in% c('+', '-'))))

trans.span %>%
  arrange(start, end) %>%
  mutate(id = sprintf('BSGatlas-TU-%s', 1:n())) %>%
  mutate_at(c('start', 'end'), as.integer) -> new.tus
  

operons %>%
  left_join(
    new.tus %>%
      separate_rows(path, sep = ';') %>%
      mutate_at('path', as.integer),
  'path') %>%
  group_by(operon) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand),
    TUs = id %>% clean_paste
  ) -> operons.span

assertthat::assert_that(with(operons.span, all(strand %in% c('+', '-'))))

operons.span %>%
  arrange(start, end) %>%
  mutate(id = sprintf('BSGatlas-operon-%s', 1:n())) %>%
  left_join(
    operons %>% 
      group_by(operon) %>%
      summarize_all(clean_paste),
    'operon'
  ) -> new.operons

# Get full transcript lengths, including UTRs
# but without TSSs or Terminators (negative lengths)
# These are added later on

bind_rows(
  trans %>%
    select(path, value) %>%
    inner_join(merging$merged_genes, c('value' = 'merged_id')),
  trans %>%
    select(path, value) %>%
    inner_join(
      bind_rows(UTRs),
      c('value' = 'id')
    )
) %>%
  select(path, start, end, strand) %>%
  group_by(path) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand)
  ) -> full.span
  
assertthat::assert_that(with(full.span, all(strand %in% c('+', '-'))))


full.span %>%
  arrange(start, end) %>%
  mutate(id = sprintf('BSGatlas-transcript-%s', 1:n())) -> new.transcript

new.transcript %>%
  left_join(
    trans %>%
      filter(str_detect(value, 'TSS')) %>%
      select(path, TSS = value),
    'path'
  ) %>%
  left_join(
    trans %>%
      filter(str_detect(value, 'terminator')) %>%
      select(path, Terminator = value),
    'path'
  ) %>%
  left_join(
    trans %>%
      filter(str_detect(value, 'gene') | str_detect(value, 'UTR')) %>%
      group_by(path) %>%
      summarize(
        features =  value %>%
          fct_reorder(transcription.order) %>%
          clean_paste()
      ),
    'path'
  ) -> new.transcript.assoc

new.operons %>%
  select(id, path) %>%
  separate_rows(path, sep = ';') %>%
  mutate_at('path', as.integer) %>%
  left_join(
    new.transcript.assoc %>%
      select(path, transcripts = id),
    'path'
  ) %>%
  group_by(id) %>%
  summarize(transcripts = clean_paste(transcripts)) %>%
  left_join(new.operons, 'id') -> new.operons.assoc

###############################################################################
# just some final column rearrangment, then save
# (statistics are in separate file)


isoforms <- list(
  operons = new.operons.assoc %>%
    select(id, start, end, strand, transcripts, TUs),
  
  transcripts =  new.transcript.assoc %>%
    select(id, path) %>%
    left_join(
      select(new.tus, path, TU = id) %>%
        separate_rows(path, sep = ';') %>%
        mutate_at('path', as.integer),
      'path'
    ) %>%
    group_by(id) %>%
    summarize(TUs = clean_paste(TU)) %>%
    left_join(new.transcript.assoc, 'id') %>%
    select(
      id, start, end, strand, TSS, Terminator, features, TUs
    ),
  tus = new.tus %>%
    select(id, start, end, strand, src, genes)
)


save(isoforms, file = 'analysis/07_isoforms.rda')

