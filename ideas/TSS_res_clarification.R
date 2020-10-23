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
      id = paste('DBTBS TSS', 1:n(), sep = '_'),
      strand,
      start = TSS,
      end = TSS
    ),
  bsubcyc$TSS %>%
    select(tss, strand) %>%
    unique %>%
    transmute(
      id = paste('BSubCyc TSS', 1:n(), sep = '_'),
      strand,
      start = tss,
      end = tss
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
              strand, start, end = start)
) %>%
  map(function(nic) {
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
      transmute(
        x = ref$id,
        d = distance,
        q = q
      )
  }) %>%
  bind_rows %>%
  separate(x, c('src', 'i'), sep = '_') -> ds
  

total <- dat %>%
  separate(id, c('src', 'i'), sep = '_') %>%
  count(src) %>%
  rename(total = n)


ds %>%
  count(src, q, d) %>%
  arrange(src, q, d) %>%
  group_by(src, q) %>%
  mutate(n.cum = cumsum(n)) %>%
  ungroup %>%
  left_join(total, 'src') %>%
  mutate(r = n.cum / total * 100) -> foo
  
      
foo %>%
  # filter(r < 95) %>%
  filter(q == 'Nicolas et al upshift') %>%
  ggplot(aes(d, r, color = src)) +
  geom_line() +
  geom_hline(yintercept = seq(70, 100, 5)) +
  geom_vline(xintercept = c(15, 25, 50, 100)) +
  xlim(c(1, 150)) +
  facet_wrap(~ q, scales = 'free')
    

ds %>%
  # filter(r < 95) %>%
  filter(q == 'Nicolas et al upshift') %>%
  ggplot(aes(d, color = src)) +
  xlim(c(0, 150)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95, 0.99, 1)) +
  geom_vline(xintercept = c(15, 25, 50)) +
  geom_density(stat = 'ecdf')
