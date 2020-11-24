path <- '/home/projects/nextprod/results/CRISPR/gw_crispr_cas9_grnas/20190321_AEB284x_assemblies/runs/ASM904v1/'

library(tidyverse)
library(plyranges)

library(conflicted)

guides <- read_tsv('data-raw/CRISPRoff/00_guides.tsv.gz')

# assign integer gid to save space
guides %>%
  select(guide) %>%
  unique %>%
  arrange(guide) %>%
  mutate(
    gid = 1:dplyr::n(),
    suffix = str_extract(guide, '^....')
  ) -> gids
write_tsv(gids, 'data-raw/CRISPRoff/01_gids.tsv.gz')

# lookup table to identify on-targets
guides %>%
  select(start, end, strand) %>%
  mutate(is.on.target = TRUE) -> on.targets

# Load the ref-annotation
rtracklayer::import.gff3('data-gff/BSGatlas_v1.1.gff') %>%
  dplyr::filter(str_detect(ID, 'gene') | str_detect(ID, 'operon')) %>%
  dplyr::filter(!str_detect(ID, '_transcribed$')) %>%
  select(ID, type, Name, Parent) %>%
  mutate(rid = 1:plyranges::n()) -> ref

# The most compact version
ref %>%
  select(rid) -> ref.range

# Prettyfied names
ref %>%
  as_tibble() %>%
  select(rid, ID, type, Parent, Name) %>%
  rowwise() %>%
  mutate_all(str_c, collapse = ';') %>%
  ungroup  -> foo

foo %>%
  left_join(foo, c('ID' = 'Parent')) %>%
  select(rid = rid.x, ID, type = type.x, Name.x, Name.y) %>%
  mutate(name = ifelse(is.na(Name.y), Name.x, Name.y)) %>%
  select(-Name.x, - Name.y) %>%
  group_by(rid, ID, type) %>%
  summarize_at('name', str_c, collapse = ', ') %>%
  ungroup %>%
  mutate_at('rid', as.integer) %>%
  arrange(rid) -> ref.names
  
assertthat::are_equal(length(ref), nrow(ref.names))

# A helper to more clearly highlighting mismatches
viz.helper <- function(ref, xs) {
  rs <- strsplit(ref, '')[[1]]
  xs %>%
    strsplit('') %>%
    map(map2, rs, ~ ifelse(.x == .y, '.', .x)) %>%
    map(str_c, collapse = '') %>%
    unlist()
}
assertthat::are_equal(
  viz.helper('xxx', c('xxx', 'TCx')),
  c('...', 'TC.')
)


# Look-up for each guide the off-targets
# and match these against the reference annotation
# This is the computationally expensive and is thus parallilzed
worker <- function(guide, gid, suffix) {
  # guide <- 'AAAAAAAAAAAAAAAAACAAGGG' ; gid <- 1L ; suffix <- 'AAAA'
  
  sprintf('%s/%s_GCF_000009045.1_ASM904v1_genomic.suff/%s.CRISPRoff.tsv',
          path, suffix, guide) %>%
    read_tsv(comment = '#',
             col_names = c('chr', 'start', 'end', 'seq', 'CRISPRspec', 'strand')) %>%
    transmute(
      gid = gid,
      start, end, seq, CRISPRspec, strand,
      mismatches = stringdist::stringdist(guide, seq, 'hamming'),
      cmp.pam = viz.helper(guide, seq) %>%
        str_replace('(?=...$)', ' ')
    ) %>%
    select( - seq) %>%
    mutate(cut.pos = ifelse(strand == '+', end - 6 - 1, start + 6)) %>%
    left_join(on.targets, c('start' , 'end', 'strand')) %>%
    mutate_at('is.on.target', replace_na, FALSE) -> bindings
  
  bindings %>%
    mutate(start = cut.pos, end = cut.pos) %>%
    select(start, end, strand, is.on.target) %>%
    mutate(seqnames = 'basu168') %>%
    as_granges() %>%
    plyranges::join_overlap_left(ref.range) %>%
    as_tibble %>%
    select(is.on.target, rid) %>%
    unique %>%
    left_join(ref.names, 'rid') %>%
    unique %>%
    select(-rid) -> targets
  
  return(list(bindings = bindings, targets = targets))
}

worker <- safely(worker)
worker('AAAAAAAAAAAAAAAAACAAGGG', 1L, 'AAAA')

# microbenchmark::microbenchmark({
#   gids %>%
#     head(n = 10) %>%
#     pmap(worker)
# }, times = 10L)
# Unit: seconds
#                                        expr      min       lq     mean   median       uq      max neval
# {     gids %>% head(n = 10) %>% pmap(worker) } 2.706565 2.856593 3.186819 3.140312 3.420395 3.892114    10

# Worst-case estimate
# > nrow(guides) / 10 * 4 / 60 / 60
# [1] 42.89167 hours, single core