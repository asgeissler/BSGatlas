source('analysis/00_load.R')
load('analysis/03_subtiwiki.rda')
load('analysis/03_dbtbs.rda')
load('analysis/01_bsubcyc.rda')

load('analysis/02_merging.rda')

source('scripts/frame_helpers.R')
source('scripts/overlap_matching.R')

#####################################
# 1. collect transcripts /operons

bind_rows(
  subtiwiki$transcripts %>%
    transmute(id = transcript, src = 'subti', gene, incomplete.flag),
  
  dbtbs$operon %>%
    transmute(id = operon, src = 'dbtbs', gene = merged_id, incomplete.flag) %>%
    separate_rows(gene, sep = ';'),
  
  bsubcyc$transunit %>%
    transmute(id = as.character(1:n()), src = 'bsubcyc',
              locus = genes, incomplete.flag = FALSE) %>%
    separate_rows(locus, sep = ';') %>%
    left_join(merging$merged_src %>%
                filter(str_detect(src, 'bsubcyc')) %>%
                select(gene = merged_id, locus),
              'locus') %>%
    select(- locus)
) -> dat

# dat %>%
#   select(id, src, incomplete.flag) %>%
#   unique %>%
#   count(src, incomplete.flag)
# # A tibble: 5 x 3
# src     incomplete.flag     n
# 1 bsubcyc FALSE            1602
# 2 dbtbs   FALSE            1107
# 3 dbtbs   TRUE               16
# 4 subti   FALSE            2249
# 5 subti   TRUE               10

dat %<>% mutate(idsrc = paste0(id, '%', src))

#####################################
# 2. complement by overlap

genes <- merging$merged_genes %>%
  rename(id = merged_id)
span <- dat %>%
  left_join(genes, c('gene' = 'id')) %>%
  group_by(idsrc) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand)
  )

# no confusion about strand (hints at possible id mapping problem)
assertthat::are_equal(
  span %>% filter(!strand %in% c('+', '-')) %>% nrow,
  0
)

cmp <- overlap_matching(rename(span, id = idsrc), genes) %>%
  filter(!antisense)

cmp %>%
  count(mode)
# mode                n
# 1 3-5_overlap        63
# 2 5-3_overlap        79
# 3 contained_by       10
# 4 contains         9486
# 5 equal            3004
# 6 without_overlap   432 (genes not contained in any operon)

cmp %>%
  transmute(
    x, y, mode,
    overlap.x = overlap / x.length * 100,
    overlap.y = overlap / y.length * 100,
    jaccard.cut = cut(jaccard * 100, seq(0, 100, length.out = 5),
                      include.lowest = TRUE)
  ) %>%
  mutate(
     mode = case_when(
       mode == '3-5_overlap' ~ "3' of transcript overlaps 5' of gene",
       mode == '5-3_overlap' ~ "5' of transcript overlaps 3' of gene",
       TRUE ~ NA_character_
     )
  ) %>%
  drop_na %>%
  ggplot(aes(x = overlap.x, y = overlap.y, col = jaccard.cut)) +
  ggsci::scale_color_jama(name = 'Jaccard Similarity') +
  geom_point() + 
  xlab('overlap rel. gene length [%]') +
  ylab('overlap rel. transcript length [%]') +
  facet_wrap(~ mode, scales = 'free')

cmp %>%
  filter(str_detect(mode, 'contained')) %>%
  View

# Rules
# contains or contained by or equal
# spans at least 70 % of transcript

complement <- cmp %>%
  filter((mode %in% c("equal", "contains", "contained_by")) |
           (overlap /x.length >= 0.7)) %>%
  select(idsrc = x, gene = y) %>%
  unique

# complement %>%
#   left_join(dat, c('idsrc', 'gene')) %>%
#   mutate(new = is.na(src)) %>%
#   select(idsrc, gene, new) %>%
#   separate(idsrc, c('id', 'src'), '%', remove = FALSE) %>%
#   count(src, new)
# src     new       n
# 1 bsubcyc FALSE  3007
# 2 bsubcyc TRUE     46
# 3 dbtbs   FALSE  2176
# 4 dbtbs   TRUE    753
# 5 subti   FALSE  4547
# 6 subti   TRUE   1968

# complement %>%
#   left_join(dat, c('idsrc', 'gene')) %>%
#   mutate(new = is.na(src)) %>%
#   select(idsrc, gene, new) %>%
#   separate(idsrc, c('id', 'src'), '%', remove = FALSE) %>%
#   group_by(idsrc, src) %>%
#   summarize(any.new = any(new)) %>%
#   ungroup %>%
#   count(src, any.new)
# src     any.new     n
# 1 bsubcyc FALSE    1585
# 2 bsubcyc TRUE       17
# 3 dbtbs   FALSE    1102
# 4 dbtbs   TRUE       21
# 5 subti   FALSE    2209
# 6 subti   TRUE       50

# complement %>%
#   left_join(dat, c('idsrc', 'gene')) %>%
#   mutate(new = is.na(src)) %>%
#   filter(new) %>%
#   select(idsrc) %>%
#   left_join(dat) %>%
#   select(idsrc, incomplete.flag) %>%
#   unique %>%
#   separate(idsrc, c('id', 'src'), '%', remove = FALSE) %>%
#   count(src, incomplete.flag)
# # src     incomplete.flag     n
# # bsubcyc FALSE              17
# # dbtbs   FALSE              17 of 1107
# # dbtbs   TRUE                4 of 16
# # subti   FALSE              49 of 2249
# # subti   TRUE                1 of 10
# 
# dat %>%
#   select(idsrc, incomplete.flag) %>%
#   unique %>%
#   separate(idsrc, c('id', 'src'), '%', remove = FALSE) %>%
#   count(src, incomplete.flag)
# # src     incomplete.flag     n
# # bsubcyc FALSE            1602
# # dbtbs   FALSE            1107
# # dbtbs   TRUE               16
# # subti   FALSE            2249
# # subti   TRUE               10


# There was no gene supposedly part of transcript, but it did not overlap
assertthat::are_equal(dat %>% anti_join(complement) %>% nrow,
                      0)

#####################################
# 3. Unify

complement %>%
  unique %>%
  left_join(dat %>%
              select(src, id, incomplete.flag, idsrc) %>%
              unique,
            'idsrc') %>%
  group_by(idsrc, src, incomplete.flag) %>%
  summarize(genes = gene %>%
              sort %>%
              clean_paste) -> short

short %>%
  group_by(genes) %>%
  summarize(n.flags = sum(incomplete.flag),
            all.flags = all(incomplete.flag),
            ratio.flags = n.flags / n() * 100,
            src = src %>% sort %>% clean_paste) -> un

un %>%
  count(n.flags, all.flags)
# n.flags all.flags     n
# 0       FALSE      2450
# 1       FALSE        15
# 1       TRUE         11

un %>%
  transmute(genes, row = 1:n()) %>%
  separate_rows(genes, sep = ';') %>%
  left_join(genes, c('genes' = 'id')) %>%
  group_by(row) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand)
  ) -> short.span

short.span %>%
  left_join(mutate(un, row = 1:n()), 'row') %>%
  arrange(start, desc(end)) %>%
  mutate(
    id = sprintf('BSGatlas-transcript-%s', 1:n()),
    possibly.incomplete = n.flags > 0,
    src = src %>%
      str_replace('subti', 'SubtiWiki') %>%
      str_replace('bsubcyc', 'BsubCyc') %>%
      str_replace('dbtbs', 'DBTBS')
  ) %>%
  select(id, start, end, strand, src, possibly.incomplete, genes) %>%
  ungroup -> transcripts


#####################################
# 4. Infer operons/isoforms

library(tidygraph)

edges <- transcripts %>%
  select(from = id, to = genes) %>%
  separate_rows(to, sep = ';')
labels <- c(edges$from, edges$to) %>% unique
nodes <- tibble(id = 1:length(labels), label = labels)
# transform edges to ids
edges %<>%
  mutate(row = 1:n()) %>%
  gather('col', 'label', from, to) %>%
  left_join(nodes, 'label') %>%
  arrange(row) %>%
  select( - label) %>%
  spread(col, id) %>%
  select( - row)

grph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

grph %>%
  activate(nodes) %>%
  mutate(operon = group_components()) %>%
  as_tibble -> group

group %>%
  filter(str_detect(label, 'transcript')) %>%
  left_join(transcripts, c('label' = 'id')) %>%
  select(operon, start, end, strand, label) %>%
  group_by(operon) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand),
    isoforms = clean_paste(label)
  ) %>%
  ungroup %>%
  arrange(start, desc(end)) %>%
  mutate(id = sprintf('BSGatlas-operon-%s', 1:n())) %>%
  select(id, isoforms, start, end, strand) -> operons

# count(operons, strand)
# no awkwards strands, check

operons <- list(operon = operons, transcrip = transcripts)
save(operons, file = 'analysis/03_operons.rda')

#####################################
# 5. Make statistics and nice plots


# ratio of genes not described by any operon