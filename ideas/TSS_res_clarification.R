source('analysis/00_load.R')

source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

load('data/01_bsubcyc.rda')
load('data/01_nicolas.rda')
load('data/03_dbtbs.rda')
load('data/03_dbtbs_xml.rda')
load('data/01_refseq.rda')

# Nicolas used binding site, not TSS for benchmark, get them again
genome <- refseq$seq[[1]]
comp <- Biostrings::reverseComplement(genome)
both <- Biostrings::DNAStringSet(list(genome, comp))
names(both) <- c('+', '-')


dbtbs_xml$promoter %>%
  filter(sigma_factor != '') %>%
  pull(binding_sequence) %>%
  str_remove_all('[^ACTG]') %>%
  unique %>%
  discard(~ .x == '') %>%
  discard(~ nchar(.x) < 30) %>%
  length
  map(~ Biostrings::vmatchPattern(.x, both, max.mismatch = 0)) %>%
  map(GRanges) %>%
  invoke(.f = c) %>%
  as_tibble() %>%
  mutate(strand = seqnames) %>%
  select(- seqnames, - width) %>%
  # fix coords of reverse complement
  mutate(
    tmp = start,
    start = ifelse(strand == '+', start,
                   as.integer(length(genome) - end + 1)),
    end   = ifelse(strand == '+', end,
                   as.integer(length(genome) - tmp + 1))
  ) %>%
  select(- tmp) -> dbtbs.full
  
  
  
# set negative strand positions to correctly match the forward

ref <- bind_rows(
  # dbtbs.full %>%
  #   mutate(id = paste0("DBTBS binding site_", 1:n())),
  dbtbs$tss %>%
    select(TSS, strand) %>%
    unique %>%
    transmute(
      id = paste('DBTBS', 1:n(), sep = '_'),
      strand,
      start = TSS,
      end = TSS
    ),
  bsubcyc$TSS %>%
    select(tss, strand) %>%
    unique %>%
    transmute(
      id = paste('BSubCyc', 1:n(), sep = '_'),
      strand,
      start = tss,
      end = tss
    ),
  bind_rows(
    dbtbs$tss %>%
      select(TSS, strand),
    bsubcyc$TSS %>%
      select(TSS = tss, strand)
  ) %>%
    unique %>%
    transmute(
      id = paste('DBTBS+BSubCyc', 1:n(), sep = '_'),
      strand,
      start = TSS,
      end = TSS
    )
)
list(
  nicolas$upshifts %>%
    transmute(id = paste0('Nicolas et al upshift_', id),
              strand,
              start = pos,
              end = pos),
  nicolas$all.features %>%
    filter(startsWith(type, "5'")) %>%
    transmute(id = paste0('Nicolas et al 5\' UTR start_', name),
              strand,
              tmp = ifelse(strand == "+", start, end),
              start = tmp,
              end = tmp) %>%
    select(- tmp)
) -> nics


nics %>%
  map(function(nic) {
    # nic <- nics[[2]]
    q <- nic$id[[1]] %>%
      strsplit('_') %>%
      map(1) %>%
      unlist
    
    qs <- nic %>%
      mutate(
        seqnames = 'placeholder'
      ) %>%
      plyranges::as_granges()
    r <- ref %>%
      mutate(
        seqnames = 'placeholder'
      ) %>%
      plyranges::as_granges()
    
    distanceToNearest(r, qs) %>%
      as_tibble %>%
      # bind_cols(ref) %>%
      # filter(between(abs(distance), 20, 30)) %>%
      # View
      transmute(
        x = ref$id,
        d = distance,
        q = q
      )
  }) %>%
  bind_rows %>%
  separate(x, c('src', 'i'), sep = '_') -> ds
  

ds %>%
  # filter(q == 'Nicolas et al upshift') %>%
  # filter(q == 'Nicolas et al 5\' UTR start') %>%
  ggplot(aes(d, color = src)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95, 1)) +
  scale_x_continuous(breaks = c(0, 15, 25, 50, 100, 150),
  # scale_x_continuous(breaks = c(0, 15, 20, 25, 50, 100, 150),
                     minor_breaks = NULL,
                     limits = c(0, 150)) +
  geom_density(stat = 'ecdf', size = 1.5, alpha = 0.8) +
  xlab("Abs. distance to closest upshift [bp]") +
  ylab("Cumulative distribution") +
  ggsci::scale_color_jama(name = NULL) +
  theme_bw(18) +
  theme(legend.position = 'bottom') +
  facet_wrap(~ q)

ggsave('ideas/TSS_res.pdf', width = 12, height = 7)

## find ±33 example
nics %>%
  map(mutate, seqnames = 'placeholder') %>%
  map(plyranges::as_granges) %>%
  (function(x) distanceToNearest(x[[2]], x[[1]])) %>%
  as_tibble %>%
  mutate(
    utr = nics[[2]]$id[queryHits],
    shift = nics[[1]]$id[subjectHits]
  ) %>%
  arrange(desc(distance)) %>%
  filter(distance > 22) -> utr.no.shift

nic <- nics[[2]]
q <- nic$id[[1]] %>%
  strsplit('_') %>%
  map(1) %>%
  unlist

qs <- nic %>%
  mutate(
    seqnames = 'placeholder'
  ) %>%
  plyranges::as_granges()
r <- ref %>%
  mutate(
    seqnames = 'placeholder'
  ) %>%
  plyranges::as_granges()

distanceToNearest(qs, r) %>%
  as_tibble %>%
  mutate(
    utr = nic$id[queryHits],
    ref = r$id[subjectHits]
  ) %>%
  semi_join(utr.no.shift, 'utr') -> foo

View(foo)
# foo %>%
#   filter(between(distance, 25, 35)) %>%
#   select(utr) %>%
#   left_join(foo) -> bar
# 
# bar %>%
#   anti_join(filter(goo, distance <= 25)) %>%
#   View
