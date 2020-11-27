path <- '/home/projects/nextprod/results/CRISPR/gw_crispr_cas9_grnas/20190321_AEB284x_assemblies/runs/ASM904v1/'

library(rslurm)
library(tidyverse)
library(plyranges)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")

source('data-raw/CRISPRoff/01_helper.R')

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

# Flag for multi on-
guides %>%
  count(guide) %>%
  mutate(multi.match = n > 1) %>%
  select(-n) %>%
  left_join(gids) %>%
  select(gid, multi.match) -> multi.flag

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
  guide.short <- str_remove(guide, '...$')
  
  sprintf('%s/%s_GCF_000009045.1_ASM904v1_genomic.suff/%s.CRISPRoff.tsv',
          path, suffix, guide) %>%
    read_tsv(comment = '#',
             col_names = c('chr', 'start', 'end', 'seq', 'CRISPRoff', 'strand')) %>%
    mutate(seq.short = str_remove(seq, '...$')) %>%
    transmute(
      gid = gid,
      start, end, seq, CRISPRoff, strand,
      mismatches = stringdist::stringdist(guide.short, seq.short, 'hamming'),
      cmp.pam = paste(
        viz.helper(guide.short, seq.short),
        str_sub(seq, -3L)
      )
    ) %>%
    select( - seq) %>%
    mutate(cut.pos = ifelse(strand == '+', end - 6 - 1, start + 6)) -> bindings
  
  bindings %>%
    dplyr::filter(mismatches <= 4) %>%
    mutate(start = cut.pos, end = cut.pos) %>%
    select(start, end, strand, cut.pos) %>%
    mutate(seqnames = 'basu168') %>%
    as_granges() %>%
    plyranges::join_overlap_left(ref.range) %>%
    as_tibble %>%
    transmute(gid = gid, cut.pos, rid) %>%
    left_join(ref.names, 'rid') %>%
    drop_na -> targets
  
  return(list(bindings = bindings,
              targets = targets,
              meta = guide.meta(bindings, targets, guide, gid)))
}

# worker('AAAAAAAAAAAAAAAAACAAGGG', 1L, 'AAAA') -> foo
# str(foo)
foo$meta %>% write_lines('foo.html')
# microbenchmark::microbenchmark({
#   gids %>%
#     head(n = 10) %>%
#     pmap(worker)
# }, times = 10L)
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

