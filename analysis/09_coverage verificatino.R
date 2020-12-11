# Validate what the coverage of signals on the annotations are
# Challange: Tiling-array signal has gaps (older genome sequence)
# These positions need to be excluded from stat

library(tidyverse)
library(rtracklayer)
library(plyranges)

library(conflicted)

conflict_prefer('filter', 'dplyr')
conflict_prefer('first', 'dplyr')
conflict_prefer('lag', 'dplyr')
conflict_prefer('desc', 'dplyr')
conflict_prefer('rename', 'dplyr')
conflict_prefer("n", "dplyr")

###############################################################################
# Load annotations
list(
  'BSGatlas v1' =  '../BSGatlas/data-gff/BSGatlas_v1.0.gff',
  'BSGatlas v1.1' = '../BSGatlas/data-gff/BSGatlas_v1.1.gff'
) %>%
  map2(names(.), function(path, x) {
    foo <- import.gff3(path)
    y <- foo %>%
      filter(type %in% c('operon', 'transcript', 'gene', 'UTR')) %>%
      select(ID)
    y$type <- y$ID %>%
      strsplit('-') %>%
      map(2) %>%
      set_names(y$ID) %>%
      unlist()
    foo %>%
      filter(type == 'CDS') %>%
      as_tibble %>%
      pull(derives_from) -> cds
    
    y$type <- y %>%
      as_tibble %>%
      with(case_when(
        str_detect(type, 'UTR') ~ type,
        type %in% c('operon', 'transcript') ~ type,
        ID %in% cds ~ 'coding gene',
        TRUE ~ 'non-coding gene'
      ))
    
    y$src <- x
    return(y)
  }) -> annots.bsg
  

## Add nicolas information
load('data/01_nicolas.rda')
nicolas$all.features %>%
  mutate(
    type = case_when(
      startsWith(type, "3'") ~ "3'UTR",
      type == "5'" ~ "5'UTR",
      type %in% c('inter', 'intra') ~ "intra/intergenic"
      # type %in% c('inter', 'intra') ~ "internal_UTR"
    )
  ) %>%
  drop_na(type) %>%
  transmute(
    seqnames = 'basu168',
    start, end, strand,
    ID = paste0('Nicolas-', type, '-', 1:n()),
    type,
    src = 'Nicolas et al.'
  ) %>%
  as_granges() -> nic.annot

# The negative control: Gaps in annotations
annots.bsg$`BSGatlas v1.1` %>%
  c(nic.annot) %>%
  reduce_ranges_directed() %>%
  gaps() -> foo
foo$ID <- as.character(1:length(foo))
foo$type <- 'annotation gap'
foo$src <- 'BSGatlas + Nicolas'
    

annot <- c(annots.bsg, nic.annot, foo) %>%
  map(function(i) {i$ID <- sprintf('%s%%%s', i$src, i$ID) ; i}) %>%
  invoke(.f = c)
annot.lens <- annot %>%
  as_tibble() %>%
  select(ID, type, src, width)
###############################################################################
# The signals
max.c <- import('data-hub/tiling-array/max_coverage.bed') %>%
  select(-name)
###############################################################################
# Gaps in the signal
read_lines('data-hub/tiling-array/genome.sizes') %>%
  strsplit('\t') %>%
  map(2) %>%
  unlist %>%
  as.integer -> genome.length

max.c %>%
  reduce_ranges_directed() %>%
  gaps() -> signal.gaps

signal.gaps %>%
  width %>%
  summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   21.00   23.00   32.06   30.00 1003.00

###############################################################################
# The sum of signal per annotation
join_overlap_intersect_directed(max.c, annot) %>%
  as_tibble %>%
  group_by(ID) %>%
  summarize(signal.sum = sum(score * width),
            signal.len = sum(width)) %>%
  ungroup %>%
  right_join(annot.lens, 'ID') %>%
  # complement missing values
  mutate_at('signal.sum',  replace_na, 0L) -> raw.coverage
###############################################################################
# gaps per nt
join_overlap_intersect_directed(annot, signal.gaps) %>%
  as_tibble %>%
  group_by(ID) %>%
  summarize(gap.len = sum(width)) -> gaps
###############################################################################
# Compute average signal excl. gaps
raw.coverage %>%
  left_join(gaps, 'ID') %>%
  mutate(
    gap.len = replace_na(gap.len, 0L),
    eff.len = width - gap.len,
    avg.coverage = signal.sum / eff.len
  ) %>%
  # ggplot(aes(eff.len, signal.len)) + geom_point()
  # good, they are equal
  # ggplot(aes(eff.len, avg.coverage)) + geom_point() + scale_x_log10() +
  # geom_vline(xintercept = c(15, 25, 50, 66, 100))
  filter(eff.len >= 100) %>%
  mutate(lab = paste(type, '-', src)) %>%
  # filter(str_detect(src, '1.1')) %>%
  filter(type != 'operon') %>%
  filter(type != 'transcript') %>%
  filter(src != 'BSGatlas v1') %>%
  mutate_at('lab', str_remove, ' v1.1$') %>%
  mutate_at('lab', fct_relevel, c(
    'annotation gap - BSGatlas + Nicolas',
    "3'UTR - Nicolas et al.",
    "3'UTR - BSGatlas",
    'intra/intergenic - Nicolas et al.',
    "5'UTR - BSGatlas",
    'non-coding gene - BSGatlas',
    "5'UTR - Nicolas et al.",
    'internal_UTR - BSGatlas',
    'coding gene - BSGatlas'
  )) %>%
  ggplot(aes(avg.coverage, color = lab)) +
  stat_ecdf(size = 1.5) +
  ggsci::scale_color_jco(name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Average coverage by tiling-array log2 signal') +
  ylab('Cumulative Distribution') +
  guides(col = guide_legend(nrow = 3, byrow = FALSE)) +
  theme_minimal(18) +
  theme(legend.position = 'bottom')

ggsave('analysis/09_coverage.pdf',
       width = 11, height = 8)

