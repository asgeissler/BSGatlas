source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

source('scripts/frame_helpers.R')

load('analysis/02_merging.rda')
load('analysis/03_tus.rda')
load('analysis/05_bsg_boundaries.rda')

# The data to work with

foo <- bsg.boundaries
foo$obsolete <- NULL
foo %>%
  map(select, id, start, end, strand, without.tu = wihtout.tu) %>%
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
    bounds %>% select(x = id, bound_type = type, without.tu),
    'x'
  ) %>%
  select(
    x, y,
    bound_type, gene_type,
    without.tu,
    distance, distance_mode = mode
  ) -> cmp_dist

# distance_match does not overlaps as 0,
# accquiore more specific values with 5/3 distances from overlap_matching
overlap_matching(bounds, genes) %>%
  filter(!antisense) %>%
  select(
    x, y, 
    bound_type = x.type, gene_type = y.type,
    without.tu = x.without.tu,
    overlap, overlap_mode = mode,
    x5.dist, x3.dist, x.length
  ) -> cmp_over

full_join(cmp_dist, cmp_over,
          c('x', 'y', 'bound_type', 'without.tu', 'gene_type')) %>%
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
  ) %>%
  select(x, y, bound_type, without.tu,
         interest_mode, interest_dist, gene_type) -> cmp_full


# assign negative numbers for improved interpretability
cmp_full %>%
  mutate(
    rel_dist = ifelse(
      interest_mode %in% c('term_contained', 'term_5_overlap', 'tss_before'),
      -interest_dist,
      interest_dist
    )
  ) %>%
  drop_na(rel_dist) -> cmp_rel

# Cut-off 2e3, similar to nicolas
BOUND <- 2e3


cmp_rel %>%
  group_by(x) %>%
  top_n(-1, abs(interest_dist)) %>%
  ungroup %>%
  mutate(bound_type = bound_type %>%
           fct_recode(Terminator = 'terminator') %>%
           fct_relevel('TSS', 'Terminator')) %>%
  mutate_at('without.tu', as.character) %>%
  mutate_at('without.tu', fct_recode, yes = 'FALSE', no = 'TRUE')  %>%
  ggplot(aes(x = rel_dist, color = without.tu)) +
  stat_ecdf(size = 1.5) +
  ggsci::scale_color_jama(name = 'Is DBTBS/BsubCyc annotation or associated\nwith a transcribed region (Nicolas et al.)') +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(-BOUND, BOUND, 500),
                     minor_breaks = seq(-BOUND, BOUND, 100),
                     limits = c(-BOUND, BOUND)) +
  geom_vline(xintercept = 0, color = 'red') +
  ylab('Empirical density') +
  xlab("Distance in [nt] relative to 5' start of a gene (3' end for terminators)") +
  facet_grid(~ bound_type, scales = 'free_y') +
  theme_bw(18) +
  theme(legend.position = 'bottom' ) -> p1

# Alternative cut-offs
# !without.tu ? 2k : 200

# Find closesed partner within cut-off
cmp_rel %>%
  filter(abs(interest_dist) < ifelse(without.tu, 200, BOUND)) %>%
  group_by(x) %>%
  top_n(-1, interest_dist) %>%
  ungroup %>%
  unique -> bound_map


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
# 1 terminator   209    2357
# 2 TSS          269    3121

bound_map %>%
  count(bound_type, gene_type) %>%
  spread(bound_type, n, fill = 0) %>%
  arrange(desc(TSS))
#   gene_type           terminator   TSS
# 1 CDS                       2165  2779
# 2 putative-non-coding         57   140
# 3 riboswitch                  67   103
# 4 putative-coding             23    36
# 5 sRNA                        17    26
# 6 tRNA                        22    16
# 7 rRNA                         2    13
# 8 asRNA                        5     3
# 9 tmRNA                        1     3
# 10 intron                       0     1
# 11 ribozyme                     1     1
# 12 SRP                          1     1

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
#   gene_type           terminator   TSS
# 1 CDS                       2165  2882
# 2 putative-non-coding         57   140
# 3 putative-coding             23    36
# 4 sRNA                        17    26
# 5 tRNA                        22    16
# 6 rRNA                         2    13
# 7 asRNA                        5     3
# 8 tmRNA                        1     3
# 9 intron                       0     1
# 10 riboswitch                  67     1
# 11 ribozyme                     1     1
# 12 SRP                          1     1


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
  # also, set min length to 15
  filter(end.utr - start.utr + 1 > 15) %>%
  group_by(type) %>%
  arrange(start.utr, end.utr) %>%
  # count(start.utr, end.utr) %>% filter(n>1)
  # # only affects 3 x 3' UTR
  group_by(type, strand.utr, start.utr, end.utr) %>%
  summarize(gene = str_c(gene, collapse = ';')) %>%
  ungroup %>%
  mutate(id = sprintf('tmp-%s-%s', type, 1:n())) %>%
  select(id, type, start = start.utr, end = end.utr, strand = strand.utr) %>%
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
  filter(tu != 'tmp-tu-1544') %>%
  # arrange(desc(width))
  ggplot(aes(x = width)) + geom_histogram() +
  scale_x_log10(breaks = c(10, 20, 30, 50, 100, 500, 1e3)) +
  xlab("Length predicted internal UTRs")

# we chose min length 15
internal.utrs %>%
  filter(tu != 'tmp-tu-1544') %>%
  filter(width > 15) %>%
  select(- width) %>%
  group_by(start, end, strand) %>%
  summarize(tus = clean_paste(tu)) %>%
  ungroup %>%
  arrange(start, end) %>%
  mutate(id = sprintf('tmp-internal_UTR-%s', 1:n())) -> int.utr

int.utr %>%
  rename(boundary = tus) %>%
  mutate(type = 'internal_UTR') -> UTRs$internal_UTR


#########
# Match with old ids again

'data-gff/BSGatlas_v1.0.gff' %>%
  rtracklayer::import.gff3() %>%
  as_tibble %>%
  filter(str_detect(type, 'UTR')) %>%
  select(id = ID, start, end, strand, type) %>%
  mutate(
    type = id %>%
      strsplit('-') %>%
      map(2) %>%
      unlist) -> legacy

legacy %>%
  mutate(i = id %>%
           strsplit('-') %>%
           map(3) %>%
           unlist() %>%
           as.integer) %>%
  group_by(type) %>%
  summarize_at('i', max) %>%
  with(set_names(i, type)) -> count.status

UTRs %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows %>%
  select(id, start, end, strand, type) %>%
  overlap_matching(legacy) %>%
  filter(!antisense, x.type == y.type) -> cmp

cmp %>%
  filter(mode == 'equal') %>%
  # count(x.type)
  select(x, y) -> clear

cmp %>%
  anti_join(clear, 'x') %>%
  anti_join(clear, 'y') %>%
  filter(jaccard > .9) %>%
  group_by(y) %>%
  top_n(1, jaccard) %>%
  ungroup %>%
  group_by(x) %>%
  top_n(1, jaccard) %>%
  ungroup %>%
  # count(y) %>% count(n)
  # good, none
  # count(x) %>% count(n)
   # count(x) %>% filter(n > 1) %>% left_join(cmp) %>% View
  # the special cases in which the two cotG genes had two UTRs
  group_by(x) %>%
  top_n(1, y) %>%
  ungroup %>%
  select(x, y) %>%
  bind_rows(clear) -> look
  # count(y) %>% count(n)

# The ones that map
map2(
  UTRs %>%
    map(inner_join, look, c('id' = 'x')) %>%
    map(select, - id) %>%
    map(rename, id = y),
  # The ones with new numbers
  UTRs %>%
    map(anti_join, look, c('id' = 'x')) %>%
    map2(names(.), ~ mutate(.x, id = sprintf('BSGatlas-%s-%s',
                                             .y,
                                             count.status[.y] + 1:n()))),
  bind_rows
) %>%
  map(arrange, start)  -> UTRs

#########

save(UTRs, file = 'analysis/06_utrs.rda')

###########################################################################
###########################################################################
# Comparison with the Nicolas et al. UTRs
# and general assessment of quality

# a) comparison of lengths

load('data/01_nicolas.rda')

nicolas$all.features %>%
  drop_na(type) %>%
  # count(type)
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
# 1449         2760         1125
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
  geom_vline(xintercept = 47, color = 'red') +
  theme_bw(18) +
  facet_wrap(src ~ type, scales = 'free_y') -> p2


# a2) KS test for UTR lengths
bind_rows(
  UTRs %>%
    bind_rows %>%
    mutate(len = end - start + 1) %>%
    mutate(src = 'BSGatlas'),
  nic.utrs %>%
    filter(type != 'intergenic') %>%
    mutate(len = end - start + 1) %>%
    mutate(src = 'Nicolas et al.')
) %>%
  filter(len >= 47) %>%
  mutate_at('type', str_remove, ' ') %>%
  mutate_at('type', str_remove, '_UTR') %>%
  select(type, src, len) %>%
  group_by(type, src) %>%
  do(i = list(.$len)) %>%
  group_by(type) %>%
  do(j = set_names(.$i, .$src)) %>%
  with(set_names(.$j, .$type)) %>%
  map(map, 1) -> lens

lens %>%
  map(function (x) {
    ks.test(x$`Nicolas et al.`,
            x$BSGatlas)
  })

bind_rows(
  UTRs %>%
    bind_rows %>%
    mutate(len = end - start + 1) %>%
    mutate(src = 'BSGatlas'),
  nic.utrs %>%
    filter(type != 'intergenic') %>%
    mutate(len = end - start + 1) %>%
    mutate(src = 'Nicolas et al.')
) %>%
  filter(len < 47) %>%
  # pull(len) %>% summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 16.00   26.00   34.00   32.94   41.00   46.00
  # count(type, src)
# 1 3'UTR        BSGatlas   707
# 2 5'UTR        BSGatlas  1161
# 3 internal_UTR BSGatlas   475
  nrow
# 2333


###########################################################################
# b) comparison of overlaps
nicolas$all.features %>%
  select(id = locus, start, end, strand, type) %>%
  mutate(
    type = case_when(
      type == "5'" ~ "5'UTR",
      type  == "3'UTR" ~ "3'UTR",
      str_detect(type, '3') ~ "3'UTR (unclear termination)",
      # type %in% c("inter", "intra") ~ "internal_UTR + intergenic"
      type %in% c("inter") ~ "internal_UTR ",
      type %in% c("intra") ~ "intergenic"
    )
  ) %>%
  drop_na -> over.ref

UTRs %>%
  bind_rows %>%
  mutate(len = end - start + 1) %>%
  # pull(len) %>% summary
  select(id, type, start, end, strand) %>%
  overlap_matching(over.ref) %>%
  filter(!antisense) %>%
  drop_na(x, y) -> cmp

cmp %>%
  select(id = y, x.type) %>%
  unique %>%
  group_by(id) %>%
  arrange(x.type) %>%
  summarize(over.type = str_c(x.type, collapse = ',')) %>%
  ungroup %>%
  right_join(over.ref, 'id') %>%
  mutate_at('over.type', replace_na, 'without overlap') %>%
  count(type, over.type) %>%
  spread(type, n, fill = 0) -> foo

over.ref %>%
  count(type) %>%
  with(set_names(n, type)) %>%
  as.list -> bar
bar$over.type <- 'Total'

bind_rows(bar, foo) %>%
  select(over.type, everything()) %>%
  rename('Overlapping BSGatlas UTRs' = over.type) -> utr.stat

write_tsv(utr.stat, 'analysis/06_utrstat.tsv')
View(utr.stat)

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
#   View
#   head
##########################################################################
# combined cowplot for publication

cowplot::plot_grid(
  p1, p2,
  ncol = 1,
  labels = 'AUTO',
  label_size = 24
)

ggsave('analysis/06_dist_len.pdf',
       width = 35, height = 35, units = 'cm')


###########################################################################