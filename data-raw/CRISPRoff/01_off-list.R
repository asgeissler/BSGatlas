path <- '/home/projects/nextprod/results/CRISPR/gw_crispr_cas9_grnas/20190321_AEB284x_assemblies/runs/ASM904v1/'

library(rslurm)
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
  mutate_all(function(i) {
    res <- str_c(i, collapse = ';')
    if (identical(res, character(0))) {
      ''
    } else {
      res
    }
  }) %>%
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

worker.safe <- safely(worker)
worker('AAAAAAAAAAAAAAAAACAAGGG', 1L, 'AAAA')

microbenchmark::microbenchmark({
  gids %>%
    head(n = 10) %>%
    pmap(worker)
}, times = 10L)
# on server rstudio worst-case 8sec

# Worst-case estimate
# nrow(guides) / 10 * 8 / 60 / 60
# ~ 86 h single core

slurm_apply(
  worker, gids,
  jobname = 'offtargets',
  nodes = 50,
  cpus_per_node = 1,
  add_objects = ls(),
  slurm_options = list(
    partition = 'panic',
    time = '1-00:00:00',
    mem = '8GB'
  ), submit = FALSE
)
#> _rslurm_offtargets
#> sbatch submit.sh

#
'_rslurm_offtargets/results_*.RDS' %>%
  Sys.glob() %>%
  head() %>%
  map(readRDS) -> rds

rds[[1]]$AAAAAAAAAAAAAAAAACAAGGG$bindings %>%
  count(is.on.target, mismatches)


guides %>%
  dplyr::filter(guide == 'AAAAAAAAAAAAAAAAACAAGGG')


rds[[1]]$AAAAAAAAAAAAAAAAACAAGGG$bindings %>%
  dplyr::filter(is.on.target)  %>%
  arrange(mismatches)

guides %>%
  dplyr::filter(start ==  15149)

rds[[1]]$AAAAAAAAAAAAAAAAACAAGGG$bindings %>%
  dplyr::filter(is.on.target)  %>%
  select(start, end, strand, CRISPRspec, mismatches, cmp.pam) %>%
  left_join(on.targets)

rds[[1]]$AAAAAAAAAAAAAAAAACAAGGG$bindings %>%
  arrange(desc(CRISPRspec))

guides %>%
  dplyr::filter(start ==  3255452)

rds[[1]] %>%
  map('bindings') %>%
  map(~ count(.x, gid, pos = CRISPRspec > 0, is.on.target)) %>%
  bind_rows() -> foo


foo %>%
  # select(pos, is.on.target) %>%
  # unique all 4 exists
  mutate(
    m = case_when(
      pos & is.on.target ~ 'clear on.target',
      pos & !is.on.target ~ 'masked on.target',
      !pos & is.on.target ~ 'weired off.target',
      !pos & !is.on.target ~ 'clear off.target'
    )
  )  %>%
  select(gid, m, n) %>%
  spread(m, n, fill = 0) %>%
  left_join(
    foo %>%
      group_by(gid) %>%
      summarize_at('n', sum),
    'gid'
  ) %>%
  dplyr::rename(total = n) %>%
  mutate_at(vars(contains('target')), ~ .x / total * 100) %>%
  select(- total, - gid) %>%
  map(summary)


# > foo$gid %>% unique %>% length
# [1] 7562