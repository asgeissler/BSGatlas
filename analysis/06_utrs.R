source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

source('scripts/frame_helpers.R')

load('analysis/02_merging.rda')
load('analysis/03_tus.rda')
load('analysis/05_bsg_boundaries.rda')

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
  filter(interest_dist < 3e3) %>%
  # filter(between(rel_dist, -500, 300)) %>%
  ggplot(aes(x = rel_dist)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-3e3, +3e3, 1e3)) +
  geom_vline(xintercept = 0, color = 'red') +
  facet_wrap(~ bound_type, scales = 'free_y')

# Cut-off 2e3, similar to nicolas
BOUND <- 2e3

cmp_rel %>%
  filter(interest_dist < BOUND) %>%
  select(x, interest_mode) %>%
  unique %>%
  mutate(foo = 'yes') %>%
  spread(interest_mode, foo, fill = '-') %>%
  group_by_at(vars(-x)) %>%
  count %>% View

# Find closesed partner within cut-off
cmp_rel %>%
  filter(interest_dist < BOUND) %>%
  group_by(x) %>%
  top_n(-1, interest_dist) %>%
  ungroup %>%
  unique -> bound_map

# count(bound_map, interest_mode)
# # interest_mode      n
# # term_5_overlap   979
# # term_after      1116
# # term_contained   454
# # tss_before      2924
# # tss_contained    208
# # tss_exact         82

bound_map %>%
  mutate(bound_type = bound_type %>%
           fct_recode(Terminator = 'terminator') %>%
           fct_relevel('TSS', 'Terminator')) %>%
  ggplot(aes(x = rel_dist)) +
  geom_histogram() +
  geom_vline(xintercept = 0, color = 'red') +
  xlab("Distance in [nt] relative to 5' start of a gene (3' end for terminators)") +
  facet_wrap(~ bound_type, scales = 'free_y')

ggsave('analysis/06_distances.pdf')

# bound_map %>%
#   count(x) %>%
#   count(n)

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
# terminator   121    2545
# TSS          184    3213

bound_map %>%
  # count(bound_type, interest_mode, gene_type) %>%
  count(bound_type, gene_type) %>%
  spread(bound_type, n, fill = 0) %>%
  arrange(desc(TSS))
# gene_type           terminator   TSS
# CDS                       2311  2837
# putative-non-coding         83   176
# riboswitch                  71   104
# putative-coding             27    34
# sRNA                        22    26
# tRNA                        24    16
# rRNA                         2    12
# asRNA                        6     3
# tmRNA                        1     3
# intron                       0     1
# ribozyme                     1     1
# SRP                          1     1

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
# CDS                       2311  2941
# putative-non-coding         83   176
# putative-coding             27    34
# sRNA                        22    26
# tRNA                        24    16
# rRNA                         2    12
# asRNA                        6     3
# tmRNA                        1     3
# intron                       0     1
# riboswitch                  71     1
# ribozyme                     1     1
# SRP                          1     1

###########################################################################
# 2. Compute UTR elements

# a) 5'/3' UTRs based on the bound_chain from step 1.
bound_chain %>%
  left_join(bounds, c('bound' = 'id')) %>%
  left_join(genes, c('gene' = 'id'),
            suffix = c('.bound', '.gene')) %>%
  mutate_if(is.numeric, as.integer) %>%
  # make relative to 5' end for TSS and 3' end for Terminators
  mutate(
    strand.utr = strand.bound,
    start.utr = case_when(
      (bound_type == 'terminator') & (strand.utr == '+') ~
        end.gene,
      (bound_type == 'terminator') & (strand.utr == '-') ~
        start.bound,
      (bound_type == 'TSS') & (strand.utr == '+') ~
        start.bound,
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
    ),
    type = ifelse(bound_type == 'TSS', "5'UTR", "3'UTR")
  ) -> raw.utrs

# needed for isoforms
save(raw.utrs, file = 'analysis/06_raw.utrs.rda')

raw.utrs %>%
  # remove empyt entries (from exact or overlapping cases)
  filter(start.utr < end.utr) %>%
  group_by(type) %>%
  arrange(start.utr, end.utr) %>%
  mutate(id = sprintf('BSGatlas-%s-%s', type, 1:n())) %>%
  ungroup %>%
  select(id, type, start = start.utr, end = end.utr, strand = strand.utr,
         gene, boundary = bound) %>%
  group_by(type) %>%
  do(foo = set_names(list(.), first(.$type))) %>%
  pull(foo) %>%
  invoke(.f = c) -> UTRs


# b) internal UTRs from the TU list
tus %>%
  select(id, tu.start = start, tu.end = end, strand, genes) %>%
  separate_rows(genes, sep = ';') %>%
  left_join(
    select(genes, genes = id, start, end),
    'genes'
  ) %>%
  group_by(id) %>%
  do(foo = list(.)) %>%
  ungroup %>%
  pull(foo) %>%
  invoke(.f = c) -> jobs

helper <- function(job) {
  job %>%
    select(start = tu.start, end = tu.end, strand, id) %>%
    unique %>%
    mutate(seqnames = 'foo') %>%
    plyranges::as_granges() -> tu
  
  job %>%
    select(start, end, strand) %>%
    mutate(seqnames = 'foo') %>%
    plyranges::as_granges() %>%
    GenomicRanges::gaps() %>%
    GenomicRanges::intersect(tu) %>%
    as_tibble() %>%
    select(- seqnames) %>%
    mutate(tu = tu$id)
}

library(parallel)
cluster <- makeCluster(detectCores() - 1, type = 'FORK')

internals <- parLapply(jobs, cl = cluster, fun = safely(helper))

internals %>%
  map('error') %>%
  map(is.null) %>%
  unlist %>%
  table
# all done, yeah
internals %>%
  map('result') %>%
  bind_rows -> internal.utrs

# The sigK related UTRs should be removed, because of its unique 
# regulation mechanism (these computed up to 10k+ nt UTRs are not transcribed)

internal.utrs %>%
  filter(tu != 'BSGatlas-tu-1538') %>%
  ggplot(aes(x = width)) + geom_histogram() +
  scale_x_log10(breaks = c(10, 20, 30, 50, 100, 500, 1e3)) +
  xlab("Length predicted internal UTRs")

ggsave('analysis/06_length_internal.pdf',
       width =  10, height = 7, units = 'cm')

# we chose min length 15
internal.utrs %>%
  filter(tu != 'BSGatlas-tu-1538') %>%
  filter(width > 15) %>%
  select(- width) %>%
  group_by(start, end, strand) %>%
  summarize(tus = clean_paste(tu)) %>%
  ungroup %>%
  arrange(start, end) %>%
  mutate(id = sprintf('BSGatlas-internal_UTR-%s', 1:n())) -> int.utr

int.utr %>%
  rename(boundary = tus) %>%
  mutate(type = 'internal_UTR') -> UTRs$internal_UTR

save(UTRs, file = 'analysis/06_utrs.rda')

###########################################################################
###########################################################################
# Comparison with the Nicolas et al. UTRs
# and general assessment of quality

# a) comparison of lengths

load('data/01_nicolas.rda')

nicolas$all.features %>%
  drop_na(type) %>%
  count(type)
  filter(!str_detect(type, 'indep')) %>%
  # count(type)
  transmute(
    id = locus, start, end, strand, 
    type = case_when(
      type == 'intra' ~ 'internal',
      type == 'inter' ~ 'intergenic',
      startsWith(type, '3') ~ "3' UTR",
      TRUE ~ "5' UTR"
    )
  ) -> nic.utrs


# UTRs %>% map(nrow) %>% unlist
# 3'UTR        5'UTR internal_UTR 
# 2095         2943         1126 
# 
# nic.utrs %>% count(type)
# type           n
# 3' UTR       249
# 5' UTR       676
# intergenic   319
# internal     186

bind_rows(
  UTRs %>%
    bind_rows %>%
    mutate(len = end - start + 1) %>%
    filter(len > 15) %>%
    mutate(src = 'BSGatlas'),
  nic.utrs %>%
    filter(type != 'intergenic') %>%
    mutate(len = end - start + 1) %>%
    mutate(src = 'Nicolas et al.')
) %>%
  # filter(len > 15) %>%
  ggplot(aes(x = len)) +
  geom_histogram() +
  scale_x_log10(breaks = c(15, 30, 50, 100, 500, 1e3)) +
  facet_wrap(src ~ type, scales = 'free_y')

ggsave('analysis/06_length_comparison.pdf', 
       width = 25, height = 15, units = 'cm')

###########################################################################
# b) comparison of overlaps
UTRs %>%
  bind_rows %>%
  mutate(len = end - start + 1) %>%
  filter(len > 15) %>%
  select(id, type, start, end, strand) %>%
  overlap_matching(nic.utrs) -> cmp


# cmp %>% 
#   filter(!antisense) %>%
#   drop_na(overlap) %>%
#   ggplot(aes(x = overlap / y.length)) + geom_histogram()


cmp %>%
  filter(!antisense) %>%
  mutate(
    # ratio = overlap / y.length > 0.5,
    s = overlap > 25,
    mode2 = case_when(
      (mode == 'without_overlap') & is.na(x) ~ 'missed Nicolas et al. annotation',
      (mode == 'without_overlap') & is.na(y) ~ 'new in BSGatlas',
      # ratio ~ 'match',
      s ~ 'match',
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(mode2) %>%
  select(x, y, mode2, bsgatlas = x.type, nicolas = y.type) -> possible.matches


possible.matches %>%
  mutate(y = ifelse(is.na(y), x, y)) %>%
  select(-x ) %>%
  unique %>%
  drop_na(mode2) %>%
  mutate(foo = TRUE) %>%
  spread(bsgatlas, foo, fill = FALSE) %>%
  select(- `<NA>`) %>%
  mutate_at('nicolas', replace_na, 'N/A') %>%
  group_by_at(vars(- y)) %>%
  count %>% View
# first good over-view, needs better structuring
# separate newness and missed ones, with indication of Nicolas S6.3 table
# only list those with problamatic matching

possible.matches %>%
  mutate_at('nicolas', str_remove, ' ') %>%
  mutate(nicolas = ifelse(nicolas %in% c('intergenic', 'internal'),
                          'internal_UTR',
                          nicolas)) %>%
  filter(bsgatlas == nicolas) -> good.matches


# The good matches
nic.utrs %>%
  semi_join(good.matches, c('id' = 'y')) %>%
  select(nicolas = type,
         id) %>%
  # count(id) %>% count(n)
  count(nicolas) %>%
  mutate(bsgatlas = ifelse(nicolas != 'intergenic',
                           'same',
                           'internal_UTR')) -> mat.match

possible.matches %>%
  drop_na(y) %>%
  # anti_join(good.matches, 'x') %>%
  anti_join(good.matches, 'y') %>%
  arrange(y) %>%
  # View
  # select(y, bsgatlas) %>% unique %>% count(y) %>% count(n)
  # select(y, bsgatlas) %>% unique %>% count(y) %>% filter(n>1)
  # might need grouping for it
  select(x, y, bsgatlas, nicolas) %>%
  drop_na -> assoc.reclass

assoc.reclass %>%
  select(-x) %>%
  unique %>%
  group_by(y, nicolas) %>%
  summarize_at('bsgatlas', clean_paste) %>%
  ungroup %>%
  count(bsgatlas, nicolas) -> mat.reclass

nic.utrs %>%
  anti_join(good.matches, c('id' = 'y')) %>%
  anti_join(assoc.reclass, c('id' = 'y')) %>%
  count(type) %>%
  rename(nicolas = type) %>%
  mutate(bsgatlas = 'missed') -> mat.missed

UTRs %>%
  bind_rows() %>%
  select(id, type) %>%
  unique %>%
  left_join(
    bind_rows(
      select(good.matches, id = x),
      select(assoc.reclass, id = x)
    ) %>%
      unique %>%
      mutate(novel = FALSE),
    'id'
  ) %>%
  mutate_at('novel', replace_na, TRUE) %>%
  count(type, novel) %>%
  spread(novel, n) %>%
  mutate(Total = `TRUE` + `FALSE`) %>%
  mutate_at('type', fct_recode, 'internal UTR' = 'internal_UTR') %>%
  select(new = `TRUE`, Total, type) -> bsg.mat

bind_rows(
  nic.utrs %>%
    count(type) %>%
    rename(nicolas = type) %>%
    mutate(bsgatlas = 'Total Annotated'),
  mat.match,
  mat.missed,
  mat.reclass
) %>%
  arrange(nicolas, desc(n)) %>%
  mutate(nicolas = fct_recode(
    nicolas,
    "5'UTR" = "5' UTR",
    "3'UTR" = "3' UTR",
    "internal UTR" = 'internal',
    'intergenic' = 'intergenic'
  ) %>%
    fct_relevel("5'UTR", "3'UTR", "internal UTR", "intergenic")) %>%
  arrange(nicolas, desc(n)) %>%
  left_join(bsg.mat, c('nicolas' = 'type')) %>%
  mutate_all(replace_na, 'N/A') %>%
  mutate_all(str_replace, '(\\d)(\\d{3})', '\\1,\\2') %>%
  mutate_at(c('new', 'Total'),
            ~ ifelse(.x != lag(.x) | is.na(lag(.x)), .x, '')) %>%
  select(type = nicolas, Description = bsgatlas,
         'Nicolas et al' = n, 'BSGatlas' = Total,
         'new UTRs' = new) %>%
  mutate(Description = case_when(
    Description == 'Total Annotated' ~ 'Annotated',
    Description == 'same' ~ 'Overlap with same type',
    Description == 'missed' ~ 'Without/Low overlap',
    TRUE ~ paste('Overlaps', Description)
  )) -> utr.stat

utr.stat %>%
  select(-type) %>%
  kable('latex', booktabs = TRUE) %>%
  kable_styling(latex_options = 'scale_down') -> k

utr.stat %>%
  mutate(row = 1:n()) %>%
  group_by(type) %>%
  summarize(from = min(row), to = max(row)) %>%
  rowwise %>%
  do({
    k <<- kableExtra::group_rows(k, .$type, .$from, .$to)
    NULL
  })
k %>%
  strsplit('\n') %>%
  unlist %>%
  `[`(2:(length(.) - 1)) %>%
  write_lines('analysis/06_utr_comparison.tex')


###########################################################################
# why missed
# nic.utrs %>%
#   anti_join(good.matches, c('id' = 'y')) %>%
#   anti_join(assoc.reclass, c('id' = 'y')) %>%
#   select(y = id) %>%
#   left_join(cmp) %>%
#   filter(!antisense) %>%
#   mutate(ratio = overlap / y.length) %>%
#   arrange(desc(ratio)) %>%
#   View
#   with(overlap / x.length * 100) %>%
#   round %>%
#   replace_na('X') %>%
#   table
#   summary
#   View
#   head
###########################################################################