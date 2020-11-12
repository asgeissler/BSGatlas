# Infer additional isoforms based on the just computed UTRs

# Idea: compute a directed graph 
# Nodes are genes and UTRs
# (+ special case of 'empty' UTRs to model very short UTRs that were filtered out)
# Edges indicate transcriptional direction
# Isoforms including UTR implied new ones are given by maximal paths
# between nodes with only outgoing and only ingoing edges

source('scripts/graph_paths.R')
source('scripts/frame_helpers.R')
source('scripts/overlap_matching.R')
source('analysis/00_load.R')

load('analysis/02_merging.rda')
load('analysis/03_tus.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/06_raw.utrs.rda')


# Collect transcriptional units, including genes and internal UTRs
tus %>%
  # !!!!!  DO NOT use the sigK special case
  filter(id != 'tmp-tu-1544') %>%
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
  bind_rows() -> us


us %>%
  left_join(
    raw.utrs %>%
      select(start = start.utr, end = end.utr, strand = strand.utr, gene, type),
    c('start', 'end', 'strand', 'type')
  ) %>%
  separate_rows(gene, sep = ';') %>%
  transmute(
    from = ifelse(type == "5'UTR", id, gene),
    to = ifelse(type == "5'UTR", gene, id)
  ) -> utr.edges
# edges of TSS and terminators to UTRs
us %>%
  left_join(
    raw.utrs %>%
      select(start = start.utr, end = end.utr, strand = strand.utr, type,
             boundary = bound),
    c('start', 'end', 'strand', 'type')
  ) %>%
  transmute(
    from = ifelse(type == "5'UTR", boundary, id),
    to = ifelse(type == "5'UTR", id, boundary)
  ) -> tssterm.edges

# add edges for too short UTRs, or those that would have negative lengths
# (eg overlaps)

raw.utrs %>%
  mutate(utr.len = end.utr - start.utr + 1) %>%
  filter(utr.len <= 15) %>%
  # similiar to UTRs, edge direction depends on TSS or Terminator type
  separate_rows(gene, sep = ';') %>%
  transmute(
    from = ifelse(bound_type == 'TSS', bound, gene),
    to = ifelse(bound_type == 'TSS', gene, bound)
  ) -> empty.edges


# Challange: some Tus have neither TSS/term not UTR to them
# -> find these and add dummies, such that the paths can find all TUs

tus %>%
  filter(id != 'tmp-tu-1544') %>%
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
      filter(bound_type == 'TSS') %>%
      separate_rows(gene, sep = ';'),
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
      filter(bound_type == 'terminator') %>%
      separate_rows(gene, sep = ';'),
    c('end' = 'gene')
  ) %>%
  transmute(
    from = end,
    to = paste('dummy-terminator', 1:n()),
    tu = id
  ) -> dummy.term

# !!!!!!!!!!!!
# The paths must be filterd to only keep dummy terminators corresponding
# to the original TU they were designed for (done after computation)
# !!!!!!!!!!!!

# Collect all nodes
bind_rows(
  UTRs %>%
    map(select, id, type) %>%
    bind_rows,
  merging$merged_genes %>%
    transmute(id = merged_id, type = 'a gene'),
  bsg.boundaries %>%
    `[`(c('TSS', 'terminator')) %>%
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

x <- Sys.time()
paths <- find_paths(iso.graph)
Sys.time() - x

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
  left_join(tus ,
            'genes') %>%
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
# [1] 206

filtered.paths <- paths %>% anti_join(blacklist)

filtered.tu <- paths.tu %>% anti_join(blacklist)

###############################################################################
# How are the operons?
# How many genes have still no transcript
# What are the spans?
# How many are due to the UTRs
# How many paths have the same genes, but just different 5/3 UTRs

# > nrow(filtered.paths)
# [1] 24801

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
  # max one TSS
  with(paths.stat, all(`dummy-TSS` <= 1)),
  with(paths.stat, all(TSS <= 1)),
  with(paths.stat, !(TSS & `dummy-TSS`)) %>% all,
  # max one Term
  with(paths.stat, all(terminator <= 1)),
  with(paths.stat, all(`dummy-Terminator` <= 1)),
  with(paths.stat, !(terminator & `dummy-Terminator`)) %>% all,
  # at most one 5'/3' TUR
  with(paths.stat, all(`5'UTR` <= 1)),
  with(paths.stat, all(`3'UTR` <= 1)),
  # note more internal UTRs then genes
  with(paths.stat, all(internal_UTR < `a gene`)),
  # at least one gene
  with(paths.stat, all(`a gene` >= 1))
)


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
# 1   1692 (-1 for sigK)
# 2    614
# 3     69
# 4     81
# 5      3
# 6     21
# 7      3

assertthat::are_equal(
  nrow(tus) - 1,
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
  transmute(i = 1:n(), node, node.type = key)

trans %>%
  select(path, value) %>%
  mutate(row = 1:n()) %>%
  gather('key', 'node', path, value) %>%
  left_join(operon.nodes, 'node') %>%
  select(row, key, i) %>%
  spread(key, i) %>%
  select(from = path, to = value) -> operon.edges

operon.grph <- tbl_graph(nodes = operon.nodes, edges = operon.edges, directed = FALSE)

operon.grph %>%
  activate(nodes) %>%
  mutate(operon = group_components()) %>%
  as_tibble %>%
  filter(node.type == 'path') %>%
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
  mutate_at(c('start', 'end'), as.integer) %>%
  arrange(start, end) %>%
  mutate(id = sprintf('tmp-TU-%s', 1:n())) -> new.tus
  
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
  mutate(id = sprintf('tmp-transcript-%s', 1:n())) -> new.transcript

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


# Get Operons spanS
# once for TUs (gene part) and once via transcripts for UTR lengths

# 1. TU
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

# 2. with UTRs
operons %>%
  left_join(
    new.transcript.assoc %>%
      separate_rows(path, sep = ';') %>%
      mutate_at('path', as.integer),
  'path') %>%
  group_by(operon) %>%
  summarize(
    utr.start = min(start),
    utr.end = max(end)
  ) -> operons.utrs


assertthat::assert_that(with(operons.span, all(strand %in% c('+', '-'))))

operons.utrs %>%
  left_join(operons.span, 'operon') %>%
  arrange(utr.start, utr.end) %>%
  mutate(id = sprintf('tmp-operon-%s', 1:n())) %>%
  left_join(
    operons %>% 
      group_by(operon) %>%
      summarize_all(clean_paste),
    'operon'
  ) -> new.operons


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
# just some final column rearrangment
# (statistics are in separate file)


isoforms <- list(
  operons = new.operons.assoc %>%
    select(id, utr.start, utr.end, strand, transcripts,
           tu.start = start, tu.end = end, TUs),
  
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
###############################################################################
# The last step before saving -> carry over legacy names

legacy <- new.env(parent = emptyenv())
load('data-raw/BSGatlas_v1_isoforms.rda', legacy)

legacy$isoforms  %>%
  map(pull, id) %>%
  map(strsplit, '-') %>%
  map(map, 3) %>%
  map(unlist) %>%
  map(as.integer) %>%
  map(max) %>%
  cbind %>%
  unlist %>%
  set_names(c('operon', 'transcript', 'TU')) -> counter.status


# Step 1: Candidates from overlaps
helper <- function(bar) {
  bind_rows(
    bar$operons %>%
      transmute(id, start = tu.start, end = tu.end, strand, type = 'operon'),
    bar$tus %>%
      transmute(id, start, end, strand, type = 'TU'),
    bar$transcripts %>%
      transmute(id, start, end, strand, type = 'transcript')
  )
}
overlap_matching(helper(isoforms), helper(legacy$isoforms)) %>%
  filter(!antisense, x.type == y.type) -> cmp

# Separate out he clear cut cases by content
helper2 <- function(bar) {
  bind_rows(
    bar$operons %>%
      transmute(id, tu = TUs, type = 'operon') %>%
      separate_rows(tu, sep = ';') %>%
      left_join(
        bar$tus %>%
          select(tu = id, content = genes) %>%
          separate_rows(content, sep = ';'),
        'tu'
      ) %>%
      group_by(id, type) %>%
      summarize_at('content', clean_paste),
    bar$tus %>%
      transmute(id, content = genes, type = 'TU'),
    bar$transcripts %>%
      transmute(id, content = paste(TSS, features, Terminator),
                type = 'transcript')
  )
}

inner_join(
  helper2(isoforms),
  helper2(legacy$isoforms),
  c('type', 'content')
) %>%
  select(x = id.x, type, y = id.y) -> clear

# Figure 
cmp %>%
  inner_join(clear)  %>%
  ggplot(aes(jaccard, color = x.type)) +
  stat_ecdf() +
  ylim(c(0, 0.1)) +
  xlim(c(0.8, 1))

#.9 is an ok cut-off
cmp %>%
  anti_join(clear, 'x') %>%
  anti_join(clear, 'y') %>%
  ggplot(aes(jaccard, color = x.type)) +
  stat_ecdf() +
  geom_vline(xintercept = .9)

cmp %>%
  anti_join(clear, 'x') %>%
  anti_join(clear, 'y') %>%
  filter(jaccard >= .9) %>%
  group_by(x) %>%
  top_n(1, jaccard) %>%
  slice(1) %>%
  group_by(y) %>%
  top_n(1, jaccard) %>%
  arrange(x) %>%
  slice(1) %>%
  ungroup %>%
  # count(x) %>% count(n)
  # count(y) %>% count(n)
  select(x, type = x.type, y) %>%
  bind_rows(clear) -> look

assertthat::are_equal(
  look %>%
    unique %>%
    count(x) %>%
    filter(n > 1) %>%
    nrow,
  0
)
assertthat::are_equal(
  look %>%
    unique %>%
    count(y) %>%
    filter(n > 1) %>%
    nrow,
  0
)
  
# update ids
map2(c('operon', 'transcript', 'TU'),
     isoforms,
     ~ mutate(.y, type = .x)) %>%
  map(left_join, look, c('id' = 'x', 'type')) %>%
  map(mutate, n = cumsum(is.na(y)) + counter.status[type]) %>%
  map(mutate, n = sprintf('BSGatlas-%s-%s', type, n)) %>%
  map(rename, tmp.id = id) %>%
  map(mutate, id = ifelse(is.na(y), n, y)) %>%
  map(select, - c(y, n, type)) %>%
  set_names(names(isoforms)) -> ongoing


ongoing %>%
  map(select, x = tmp.id, y = id) %>%
  bind_rows -> look.full


# make sure that links are coherent

ongoing$transcripts %>%
  rename(old.TUs = TUs) %>%
  left_join(select(ongoing$tus, old.TUs = tmp.id, TUs = id), 'old.TUs') %>%
  select(- old.TUs) -> ongoing$transcripts

ongoing$operons %>%
  select(id, transcripts, TUs) %>%
  gather('k', 'x', transcripts, TUs) %>%
  separate_rows('x', sep = ';') %>%
  left_join(look.full, 'x') %>%
  group_by(id, k) %>%
  summarize_at('y', str_c, collapse = ';') %>%
  spread(k, y) %>%
  ungroup %>%
  left_join(
    ongoing$operons %>%
      select(- transcripts, - TUs),
    'id'
  ) -> ongoing$operons

ongoing %>%
  map(select, - tmp.id) -> isoforms


###############################################################################


save(isoforms, file = 'analysis/07_isoforms.rda')

