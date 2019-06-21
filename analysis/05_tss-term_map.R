source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

load('analysis/01_bsubcyc.rda')
load('analysis/01_nicolas.rda')
load('analysis/03_dbtbs.rda')

tss.dat <- bind_rows(
  nicolas$upshifts %>%
    transmute(src = 'Nicolas et al upshift', id, TSS = pos, strand, sigma) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  dbtbs$tss %>%
    transmute(src = 'DBTBS', id = 1:n(), TSS, strand, sigma) %>%
    mutate_at('id', as.character) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  bsubcyc$TSS %>%
    transmute(src = 'BsubCyc', id = 1:n(), TSS = tss, strand, sigma) %>%
    mutate_at('id', as.character)
) %>%
  mutate(src_id = paste0(src, '%', id))
  

tss.dat %>%
  mutate(start = TSS, end = TSS) %>%
  select(id = src_id, start, end, strand) %>%
  distance_matching(., .) %>%
  filter(x != y, !antisense) %>%
  mutate(
    from = x %>%
      strsplit('%') %>%
      map(1) %>%
      unlist,
    to = y %>%
      strsplit('%') %>%
      map(1) %>%
      unlist
  ) %>%
  mutate(nearest = ifelse(mode == 'x.after.y',
                          - distance,
                          distance),
         abs.near = abs(nearest)) %>%
  group_by(x, from, to) %>%
  top_n(-1, abs.near) %>%
  ungroup -> tss.near

tss.near %>%
  filter(from <= to) %>%
  filter(abs.near < 100) %>%
  ggplot(aes(x = nearest)) +
  geom_histogram() +
  xlab('Distance to closest TSS') +
  facet_wrap(from ~ to, scale = 'free_y')

ggsave(file = 'analysis/05_tss_comparison.pdf',
       width = 7, height = 7, units = 'in')
