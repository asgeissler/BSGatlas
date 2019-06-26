source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

source('scripts/frame_helpers.R')

load('analysis/05_bsg_boundaries.rda')
load('analysis/03_operons.rda')
load('analysis/02_merging.rda')

# The data to work with

foo <- bsg.boundaries
foo$TSS %<>% mutate(start = TSS, end = TSS)
foo %>%
  map(select, id, start, end, strand) %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows -> bounds

merging$merged_genes %>%
  rename(id = merged_id) -> genes

###########################################################################
# Questions to guide the analysis:
# - Distance of TSSs relative to genes
# - How many genes are part of a TU
# - How often is it the first one?
# - What are new isoforms
# - Which resource causes the creation of the new isoform?
# - Should I have a cut-off for limiting it
# - How many TSSs/terminators do not have an association to a gene?
# - Is there a nicolas UTR to it?

# 1. Compute distances to genes
distance_matching(bounds, genes) %>%
  filter(!antisense) %>%
  left_join(genes %>%
              select(y = id, gene_type = type), 
            'y') %>%
  left_join(
    bounds %>% select(x = id, bound_type = type),
    'x'
  ) %>%
  select(
    x, y,
    bound_type, gene_type,
    distance, distance_mode = mode
  ) -> cmp_dist

overlap_matching(bounds, genes) %>%
  filter(!antisense) %>%
  select(
    x, y, 
    bound_type = x.type, gene_type = y.type,
    overlap, overlap_mode = mode,
    x5.dist, x3.dist, x.length
  ) -> cmp_over

full_join(cmp_dist, cmp_over,
          c('x', 'y', 'bound_type', 'gene_type')) %>%
  drop_na(x, y) %>%
  mutate_if(is.numeric, as.integer) %>%
  mutate(
    interest_mode = case_when(
      bound_type == 'terminator' & overlap_mode == 'contained_by' ~
        'term_contained',
      bound_type == 'terminator' & overlap_mode == '5-3_overlap' ~
        'term_5_overlap',
      bound_type == 'terminator' & distance_mode == 'x.after.y' ~
        'term_after',
      bound_type == 'TSS' & distance_mode == '5prime.equal' ~
        'tss_exact',
      bound_type == 'TSS' & overlap_mode == 'contained_by' ~
        'tss_contained',
      bound_type == 'TSS' & distance_mode == 'x.before.y' ~
        'tss_before',
      TRUE ~ 'other'
    ),
    interest_dist = case_when(
      interest_mode == 'term_contained' ~ x3.dist + x.length,
      interest_mode == 'term_5_overlap' ~ overlap,
      interest_mode == 'term_after' ~ distance,
      interest_mode == 'tss_contained' ~ x5.dist,
      interest_mode == 'tss_exact' ~ distance,
      interest_mode == 'tss_before' ~ distance
    )
  ) -> cmp_full

cmp_full %>%
  # group_by(x) %>%
  # top_n(-1, interest_dist) %>%
  # count(interest_mode) 
  select(x, interest_mode) %>%
  mutate(foo = TRUE) %>%
  unique %>%
  spread(interest_mode, foo, fill = '-') %>%
  group_by_at(vars(- x)) %>%
  count %>%
  View

cmp_full %>%
  select(x, y, bound_type, gene_type, interest_mode, interest_dist) %>%
  mutate(
    rel_dist = ifelse(
      interest_mode %in% c('term_contained', 'term_5_overlap', 'tss_before'),
      -interest_dist,
      interest_dist
    )
  ) %>%
  drop_na(rel_dist) -> cmp_rel

cmp_rel %>%
  group_by(x) %>%
  top_n(-1, interest_dist) %>%
  ungroup %>%
  filter(interest_dist < 1e3) %>%
  # filter(between(rel_dist, -500, 300)) %>%
  ggplot(aes(x = rel_dist)) +
  geom_histogram() +
  geom_vline(xintercept = 0, color = 'red') +
  facet_wrap(~ bound_type, scales = 'free_y')

# Cut-off 1e3, similar to nicolas

cmp_rel %>%
  filter(interest_dist < 1e3) %>%
  select(x, interest_mode) %>%
  unique %>%
  mutate(foo = 'yes') %>%
  spread(interest_mode, foo, fill = '-') %>%
  group_by_at(vars(-x)) %>%
  count %>% View

# Prefer
# term_after < term_5_overlap < term_contained
# TSS exact < before < contained
# in each category the min
cmp_rel %>%
  filter(interest_dist < 1e3) %>%
  mutate(prio = case_when(
    interest_mode == 'term_after' ~     0,
    interest_mode == 'term_5_overlap' ~ 1,
    interest_mode == 'term_contained' ~ 2,
    interest_mode == 'tss_exact' ~     0,
    interest_mode == 'tss_before' ~    1,
    interest_mode == 'tss_contained' ~ 2
  )) %>%
  group_by(x) %>%
  top_n(-1, prio) %>%
  ungroup -> bound_map
# count(interest_mode)
# interest_mode      n
# term_5_overlap   627
# term_after      1456
# term_contained    42
# tss_before      2986
# tss_contained     74
# tss_exact         82


bound_map %>%
  ggplot(aes(x = rel_dist)) +
  geom_histogram() +
  geom_vline(xintercept = 0, color = 'red') +
  facet_wrap(~ bound_type, scales = 'free_y')

bound_map %>%
  count(x) %>%
  count(n)

bounds %>%
  select(id, type) %>%
  left_join(bound_map, c('id' = 'x')) %>%
  select(id, type, y) %>%
  mutate(has.partner = !is.na(y)) %>%
  select(id, type, has.partner) %>%
  unique %>%
  count(type, has.partner) %>%
  mutate(has.partner = ifelse(has.partner, 'partner', 'none')) %>%
  spread(has.partner, n)
# type        none partner
# terminator   150    2117
# TSS          248    3142

bound_map %>%
  # count(bound_type, interest_mode, gene_type) %>%
  count(bound_type, gene_type) %>%
  spread(bound_type, n, fill = 0) %>%
  arrange(desc(TSS))
# gene_type           terminator   TSS
# CDS                       1996  2795
# putative-non-coding         43   160
# riboswitch                  30    98
# putative-coding             28    30
# sRNA                         8    24
# tRNA                        14    15
# rRNA                         3    13
# asRNA                        3     3
# tmRNA                        0     2
# ribozyme                     0     1
# SRP                          0     1

###########################################################################
# An initial good map -> minor adjustment: Chain riboswitches

distance_matching(genes, genes) %>%
  filter(!antisense) %>%
  left_join(genes, c('x' = 'id')) %>%
  select(x, type, y, mode, distance) %>%
  filter(type == 'riboswitch', x != y) %>%
  filter(mode == 'x.before.y') %>%
  group_by(x) %>%
  top_n(-1, distance) %>%
  ungroup %>%
  select(y = x, gene = y) %>%
  right_join(bound_map, 'y') %>%
  transmute(bound = x, bound_type,
            # chain riboswitches but only for TSSS
            gene = ifelse(is.na(gene) | (bound_type != 'TSS'),
                          y, gene)) %>%
  left_join(genes %>% select(gene = id, gene_type = type),
            'gene') -> bound_chain


bound_chain %>%
  count(bound_type, gene_type) %>%
  spread(bound_type, n, fill = 0) %>%
  arrange(desc(TSS))
# gene_type           terminator   TSS
# CDS                       1996  2894
# putative-non-coding         43   160
# putative-coding             28    30
# sRNA                         8    24
# tRNA                        14    15
# rRNA                         3    13
# asRNA                        3     3
# tmRNA                        0     2
# ribozyme                     0     1
# SRP                          0     1
# riboswitch                  31     0

###########################################################################
# Helper track for adhoc viz in browser
load('analysis/01_refseq.rda')
genome.size <- refseq$seq$`168` %>% length

bound_chain %>%
  left_join(bounds, c('bound' = 'id')) %>%
  left_join(genes, c('gene' = 'id'),
            suffix = c('.bound', '.gene')) %>%
  mutate_if(is.numeric, as.integer) %>%
  mutate(
    strand.utr = strand.bound,
    start.utr = case_when(
      (bound_type == 'terminator') & (strand.utr == '+') ~
        end.gene,
      (bound_type == 'terminator') & (strand.utr == '-') ~
        start.bound,
      (bound_type == 'TSS') & (strand.utr == '+') ~
        end.bound,
      (bound_type == 'TSS') & (strand.utr == '-') ~
        end.gene
    ),
    end.utr = case_when(
      (bound_type == 'terminator') & (strand.utr == '+') ~
        end.bound,
      (bound_type == 'terminator') & (strand.utr == '-') ~
        start.gene,
      (bound_type == 'TSS') & (strand.utr == '+') ~
        start.gene,
      (bound_type == 'TSS') & (strand.utr == '-') ~
        end.bound
    )
  ) %>%
  select(start = start.utr, end = end.utr, strand = strand.utr) %>%
  arrange(start) %>%
  transmute(
    chr = 'basu168',
    start1 = start - 1, end1 = end,
    primary.name = paste0('test-', 1:n()),
    score = 0, strand,
    start2 = start1, end2 = end1,
    rgb = ifelse(strand == '+', '0,255,0', '255,0,0')
  ) -> foo


special <- foo %>%
  filter(end1 > genome.size)
foo %>%
  filter(end1 <= genome.size) %>%
  bind_rows(
    special %>%
      mutate(end1 = genome.size,
             end2 = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    special %>%
      mutate(start1 = 0,
             end1 = end1 %% genome.size,
             start2 = 0,
             end2 = end1,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start1) %>%
  write_tsv('~/Downloads/test.bed', col_names = FALSE)

###########################################################################
# Implications of this map for isoforms

overlap_matching(genes, operons$transcript) %>%
  filter(!antisense) %>%
  filter(mode %in% c('equal', 'contained_by')) %>%
  mutate(
    start = x5.dist == 0,
    ends = x3.dist == 0,
    tu_mode = case_when(
      start & ends ~ 'mono-cistronic',
      start ~ 'gene starts TU',
      ends ~ 'gene ends TU',
      TRUE ~ 'gene in middle (novel isoform)'
    )
  ) %>%
  select(gene = x, tu = y, part, start, ends, tu_mode) -> gene.tu

bound_chain %>%
  left_join(gene.tu, 'gene') %>%
  mutate_at('tu_mode', replace_na, 'without TU') %>%
  count(bound_type, tu_mode) %>%
  spread(bound_type, n)
