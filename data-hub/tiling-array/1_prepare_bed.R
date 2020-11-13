library(tidyverse)

# Tasks
# - export array probe design
# - load raw tiling array signals
# - create bed files
# - liftover to updated assembly


# Probe design
old.genome <- 'NC_000964.2'



#filter to select area between
#!platform_table_begin
#!platform_table_end
custom.filter <- function(xs) {
  i <- 1
  while (xs[i] != '!platform_table_begin') {
    i <- i + 1
  }
  j <- i + 1
  while (xs[j] != '!platform_table_end') {
    j <- j + 1
  }
  i <- i + 1
  j <- j - 1
  return(xs[i:j])
}

'GEO/GPL8486_family.soft.gz' %>%
  read_lines() %>%
  custom.filter %>%
  read_tsv() -> design

design %>%
  transmute(
    chr = old.genome,
    start = RANGE_START - 1,
    end = RANGE_END,
    name = ID %>%
      strsplit(';') %>%
      map(2) %>%
      unlist,
    score =0,
    RANGE_STRAND
  ) %>%
  # filter(start > end)
  # ommitted 6 probes over origin for viz simplisity
  filter(start < end) %>%
  write_tsv('old_genome/probes.bed', col_names = FALSE)
  

# Sample meta info
nicolas <- GEOquery::getGEO(filename = 'GEO/GSE27219_family.soft.gz')
GEOquery::GSMList(nicolas) %>%
  map(GEOquery::Meta) %>%
  map(`[`, c('title', 'description', 'geo_accession')) %>%
  map(as_tibble) %>%
  bind_rows() %>%
  separate(title, c('experiment', 'replicate'), '_', remove = FALSE) %>%
  mutate(
    description = description %>%
      strsplit('_') %>%
      map(1) %>%
      unlist
  ) -> meta


# Load raw data, filter, and write to bedgraph
meta$geo_accession %>%
  # head %>%
  map(function(i) {
    x <- sprintf('GEO/GSE27219_RAW/%s.pair.gz', i)
    x %>%
      read_tsv(comment = '#') %>%
      filter(str_starts(SEQ_ID, 'Bsub_')) %>%
      transmute(
        strand = ifelse(SEQ_ID == 'Bsub_positive_strand', 'fwd', 'rev'),
        pos = POSITION,
        value = log2(PM)
      ) -> foo
    
    foo %>%
      plyr::dlply('strand', list) %>%
      map(1) %>%
      # map(select, - strand) %>%
      map(arrange, pos) %>%
      map(transmute, chr = old.genome,
          # make signal from position up to next probe start
          s = pos - 1,
          # end is open interval yet still minus extra one because of some
          # weird overlap
          e = replace_na(lead(pos), 4214630) - 1,
          nr = 1:n(),
          score = "0",
          strand = ifelse(strand == 'fwd', '+', '-'),
          v = value) %>%
      # map2(names(.), ~ write_tsv(.x, sprintf('old_genome/%s_%s.bedgraph', i, .y),
      map2(names(.), ~ write_tsv(.x, sprintf('old_genome/%s_%s.bed', i, .y),
                                 col_names = FALSE)) -> bar
    x
  })

# Prepare a trackhub addition file

# selection matrix lookups
hash <- list(
  exp = meta %>%
    select(experiment) %>%
    unique %>%
    mutate(rn = sprintf('%02x', 1:n())) %>%
    with(set_names(rn, experiment)),
  rep = meta %>%
    select(replicate) %>%
    unique %>%
    mutate(rn = LETTERS[1:n()]) %>%
    with(set_names(rn, replicate))
)
# hash$cl <-
library(scales)
RColorBrewer::brewer.pal(11, 'Spectral') %>%
  gradient_n_pal() -> foo
foo(seq(0, 1, length.out = length(hash$exp))) %>%
# ggsci::pal_gsea(n = length(hash$exp))(length(hash$exp)) %>%
  colorspace::hex2RGB() %>%
  `@`(coords) %>%
  `*`(255) -> foo
paste(foo[, 1], foo[, 2], foo[, 3], sep = ',') %>%
  set_names(names(hash$exp)) -> hash$cl

# block container
block <- function(strand) {
  sprintf(
"

track tiling_%s
type bigWig
group tiling
container multiWig
shortLabel Tiling array %s
longLabel log2 signals from Nicolas et al. tiling array - %s strand
visibility full
autoScale off
alwaysZero off
viewLimits 2:20
windowingFunction mean
smoothingWindow 10
yLineOnOff on
yLineMark 10
maxHeightPixels 300:250:150
showSubtrackColorOnUi on
aggregate solidOverlay
graphTypeDefault points
priority  50
subGroup1 exp Experiment %s
subGroup2 rep Replicate %s
dimensions dimX=sap dimY=day
sortOrder exp=+ rep=+

", strand, strand, strand,
hash$exp %>%
  map2(names(.), ~ paste(.x, .y, sep = '=')) %>%
  invoke(.f = paste),
hash$rep %>%
  map2(names(.), ~ paste(.x, .y, sep = '=')) %>%
  invoke(.f = paste)
)
}



indiv <- function(strand, title, geo, exp, expL, rep) {
  sprintf(
"

  track tiling_%s_%s
  type bigWig
  parent tiling_%s
  shortLabel %s
  longLabel %s replicate %s (%s-%s)
  color %s
  subGroups rep=%s exp=%s
  bigDataUrl lifted_genome/%s_%s.bw

", strand, geo, strand, title,
expL, rep, title, geo,
hash$cl[exp],
hash$rep[rep], hash$exp[exp],
geo, strand
)
}

indiv('fwd', 'TITLE', "GEO", 'EXZP', 'LONG', "REP")

# build trackDB addition
c(
  '
track tiling_probes
  type bigBed 6 .
  group tiling
  shortLabel Lifted Probes
  longLabel LiftedOver Probe coordinates
  visibility hide
  colorByStrand 0,0,50 50,0,0
  bigDataUrl lifted_genome/probes.bb

  ',
  block('fwd'),
  with(meta, indiv('fwd', title, geo_accession, experiment, description, replicate)),
  block('rev'),
  with(meta, indiv('rev', title, geo_accession, experiment, description, replicate))
) %>%
  write_lines('track_add.txt')



