source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

load('analysis/01_bsubcyc.rda')
load('analysis/01_nicolas.rda')
load('analysis/03_dbtbs.rda')

##############################################################################
# 1. find TSS resolution

tss.dat <- bind_rows(
  nicolas$upshifts %>%
    transmute(src = 'Nicolas et al upshift', id, TSS = pos, strand, sigma,
              pubmed = '') %>%
    mutate_at('sigma', str_remove, 'Sig'),
  dbtbs$tss %>%
    transmute(src = 'DBTBS', id = 1:n(), TSS, strand, sigma,
              pubmed = reference) %>%
    mutate_at('id', as.character) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  bsubcyc$TSS %>%
    transmute(src = 'BsubCyc', id = 1:n(), TSS = tss, strand, sigma,
              pubmed = cite) %>%
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

##############################################################################
# 2. TSS merging

tss.dat %>%
  transmute(id = src_id, TSS, strand, sigma) %>%
  mutate(res.limit = ifelse(startsWith(id, 'Nicolas'), 22, 0),
         start = TSS - res.limit,
         end = TSS + res.limit) -> tss.win

tss.win %>%
  transmute(i = 1:n(), name = id) -> nodes

overlap_matching(tss.win, tss.win) %>%
  filter(!antisense) %>%
  select(from = x, to = y) %>%
  drop_na %>%
  unique %>%
  mutate(row = 1:n()) %>%
  gather('pair', 'name', from, to) %>%
  left_join(nodes, 'name') %>%
  select(- name) %>%
  spread(pair, i) %>%
  select(- row) -> edges

library(tidygraph)
grp <- tbl_graph(nodes = nodes, edges = edges,
                 directed = FALSE)

grp %>%
  activate(nodes) %>%
  mutate(group = group_components()) %>%
  # filter(group == 1) %>%
  # plot
  as_tibble %>%
  select(- i) %>%
  left_join(tss.win, c('name' = 'id')) %>%
  #tmp substitute YlaC -> C
  mutate(sigma = case_when(
    is.na(sigma) ~ '?',
    sigma == '-' ~ '?',
    sigma == 'YlaC' ~ 'C',
    TRUE ~ sigma
  )) -> tss.merge

tss.merge %>%
  group_by(group) %>%
  filter(res.limit == min(res.limit)) %>%
  ungroup %>%
  arrange(TSS) %>%
  transmute(id = paste0('BSGatlas-TSS-', 1:n()),
            TSS, res.limit, start, end ,strand) -> bsg.tss

cmp <- overlap_matching(bsg.tss, tss.win) %>%
  filter(!antisense) %>%
  select(x, y) %>%
  drop_na %>%
  left_join(tss.dat, c('y' = 'src_id')) %>%
  select(id = x, src, pubmed) %>%
  gather('key', 'value', src, pubmed) %>%
  separate_rows(value, sep = ';') %>%
  unique %>%
  drop_na %>%
  filter(value != '') %>%
  group_by(id, key) %>%
  summarize(value = value %>% sort %>% invoke(.f = paste, sep = ';')) %>%
  spread(key, value, fill = '')
  

bsg.tss %<>% left_join(cmp, 'id')

library(venn)
pdf(file = 'analysis/05_venn_tss.pdf')
bsg.tss %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  mutate_at('src', fct_recode,
            'Nicolas et al.\nupshift' = 'Nicolas et al upshift') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  with(set_names(map(i, 1), src)) %>%
  venn(cexil = 1.3,
       cexsn = 1.5, zcolor = 'style')
dev.off()



##############################################################################
# 2. Terminator merging

dat.term <- bind_rows(
  bsubcyc$terminator %>%
    arrange(start) %>%
    unique %>%
    transmute(id = 1:n(), start, end, strand, energy,
              src = 'BsubCyc', prio = 0) %>%
    mutate_at(c('id', 'energy'), as.character),
  dbtbs$term %>%
    arrange(start) %>%
    transmute(id = 1:n(), start, end, strand, energy = energies,
              src = 'DBTBS', prio = 1,
              pubmed = reference) %>%
    mutate_at('id', as.character),
  nicolas$downshifts %>%
    transmute(id, start = pos, end = pos, strand, energy,
              src = 'Nicolas et al. downshift', prio = 2) %>%
    mutate_at('energy', as.character)
) %>%
  mutate_at('energy', replace_na, '') %>%
  drop_na(start) %>%
  rename(rid = id) %>%
  mutate(id = paste0(src, '%', rid))

distance_matching(dat.term, dat.term) %>%
  filter(!antisense) %>%
  filter(x != y) %>%
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
  ungroup -> term.near
  
term.near %>%
  # filter(from < to) %>%
  filter((from != to) | (x < y)) %>%
  filter(abs.near < 1e2) %>%
  # group_by(from, to) %>% top_n(-5, abs.near)  %>% slice(1:5) %>% View
  ggplot(aes(x = nearest)) +
  geom_histogram() +
  xlab('Distance to closest TSS') +
  facet_wrap(from ~ to, scale = 'free_y')

ggsave(file = 'analysis/05_term_comparison.pdf',
       width = 7, height = 7, units = 'in')
