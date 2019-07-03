# Infer additional isoforms based on the just computed UTRs

# Idea: compute a directed graph 
# Nodes are genes and UTRs
# (+ special case of 'empty' UTRs to model very short UTRs that were filtered out)
# Edges indicate transcriptional direction
# Isoforms including UTR implied new ones are given by maximal paths
# between nodes with only outgoing and only ingoing edges

source('scripts/graph_paths.R')
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
    to = ifelse(type == "5'UTR", gene, gene)
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
    to = ifelse(bound_type == 'TSS', gene, bound),
  ) -> empty.edges

# Collect all nodes
bind_rows(
  UTRs %>%
    map(select, id, type) %>%
    bind_rows,
  merging$merged_genes %>%
    transmute(id = merged_id, type = 'a gene'),
  bsg.boundaries %>%
    map2(names(.), ~ transmute(.x, id, type = .y)) %>%
    bind_rows
) %>%
  mutate(number = 1:n()) %>%
  select(number, id, type) -> nodes


edges <- bind_rows(
  tu.edges,
  utr.edges,
  tssterm.edges,
  empty.edges
) %>%
  # filter(!complete.cases(.))
  mutate(row = 1:n()) %>%
  gather('key', 'id', from, to) %>%
  left_join(nodes, 'id') %>%
  select(row, key, number) %>%
  spread(key, number) %>%
  select(from, to) %>%
  drop_na

iso.graph <- tbl_graph(nodes, edges, directed = TRUE)
