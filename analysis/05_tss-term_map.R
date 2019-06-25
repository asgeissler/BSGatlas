source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

source('scripts/frame_helpers.R')

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
  select(TSS, res.limit, start, end, strand, sigma) %>%
  group_by(TSS, res.limit, start, end, strand) %>%
  summarize(sigma = sigma %>% unique %>% invoke(.f = paste0)) %>%
  ungroup %>%
  # count(sigma) %>% View
  mutate(id = paste0('BSGatlas-TSS-', 1:n())) %>%
  select(id, everything()) -> bsg.tss

cmp <- overlap_matching(bsg.tss, tss.win) %>%
  filter(!antisense) %>%
  select(x, y) %>%
  drop_na %>%
  left_join(tss.dat, c('y' = 'src_id')) %>%
  select(id = x, src, pubmed) %>%
  unique %>%
  gather('key', 'value', src, pubmed) %>%
  separate_rows(value, sep = ';') %>%
  unique %>%
  drop_na %>%
  filter(value != '') %>%
  group_by(id, key) %>%
  summarize(value = value %>% sort %>% invoke(.f = paste, sep = ';')) %>%
  spread(key, value, fill = '')
  
# cleanup
bsg.tss %<>% left_join(cmp, 'id')
bsg.tss %<>% mutate(sigma = ifelse(sigma == 'C', 'YlaC', sigma))

# sigma distribution

tss.dat %>%
  select(src, sigma) %>%
  bind_rows(
    bsg.tss %>%
      transmute(src = 'BSGatlas', sigma)
  ) %>%
  mutate_at('sigma', replace_na, '?') %>%
  mutate_at('sigma', str_replace, '-', '?') %>%
  mutate_at('sigma', str_replace_all, '([A-Z])(?=[A-Z])', '\\1;') %>%
  separate_rows(sigma, sep = ';') -> sigma.bind

sigma.bind %>%
  count(src, sigma) %>%
  spread(src, n, fill = 0) %>%
  arrange(desc(BSGatlas)) %>%
  mutate(sigma = ifelse(sigma != 'YlaC',
                        paste0('Sig', sigma),
                        sigma),
         group = case_when(
           sigma == 'Sig?' ~ 'unknown',
           1:n() <= 6 ~ sigma,
           TRUE ~ 'other'
        )) %>%
  group_by(group) %>%
  summarize_if(is.numeric, sum) %>%
  gather('src', 'n', - group) %>%
  left_join(
    sigma.bind %>%
      count(src) %>%
      rename(total = n),
    'src'
  ) %>%
  # mutate(prop = n) %>%
  mutate(prop = n / total * 100) %>%
  select(src, group, prop) %>%
              # filter(src == 'BSGatlas') %>%
              # with(reorder(group, - prop)) %>%
              # levels %>% dput
  mutate_at('src', fct_recode, 'Nicolas et al.' = 'Nicolas et al upshift') %>%
  mutate_at('group', fct_relevel,
            "SigA", "other", "SigF", "unknown", "SigE", "SigK", "SigG") %>%
  ggplot(aes(x = src, fill = group, y = prop)) +
  geom_bar(stat = 'identity') +
  xlab(NULL) + ylab('Proprotion [%]') +
  ggsci::scale_fill_jama(name = 'Sigma Factor') +
  theme_minimal(base_size = 16)
ggsave(file = 'analysis/05_sigma_proportion.pdf')
  
bsg.tss %>%
  count(res.limit)
# res.limit     n
# 0           706
#22          2684
bsg.tss %>%
  separate_rows(src, sep = ';') %>%
  count(src) %>%
  bind_rows(tibble(src = 'BSGatlas', n = nrow(bsg.tss))) %>%
# 1 BsubCyc                 556
# 2 DBTBS                   644
# 3 Nicolas et al upshift  3269
# 4 BSGatlas               3390
  mutate(src = fct_reorder(src, - n) %>%
           fct_recode('Nicolas et al.' = 'Nicolas et al upshift'))  %>%
  ggplot(aes(x = src, y = n, fill = src)) + 
  ggsci::scale_fill_jama(name = NULL) +
  geom_bar(stat = 'identity') +
  theme_minimal(base_size = 16) +
  scale_y_continuous(breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500)) +
  ylab('count of sigma factor binding sites') + xlab(NULL) +
  theme(legend.position = 'none')

ggsave(file = 'analysis/05_count_promoter.pdf')

# venn diagram

library(venn)
pdf(file = 'analysis/05_venn_tss.pdf')
bsg.tss %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  mutate_at('src', fct_recode,
            'DBTBS (TSS)' = 'DBTBS',
            'BsubCyc\n(TSS)' = 'BsubCyc',
            'Nicolas et al.\nupshift' = 'Nicolas et al upshift') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  with(set_names(map(i, 1), src)) %>%
  venn(cexil = 1.3,
       cexsn = 1.3,
       zcolor = ggsci::pal_jama()(3))
dev.off()

##############################################################################
# 2. Terminator merging

dat.term <- bind_rows(
  bsubcyc$terminator %>%
    arrange(start) %>%
    unique %>%
    transmute(id = 1:n(), start, end, strand, energy,
              src = 'BsubCyc', prio = 1) %>%
    mutate_at(c('id', 'energy'), as.character),
  dbtbs$term %>%
    arrange(start) %>%
    transmute(id = 1:n(), start, end, strand, energy = energies,
              src = 'DBTBS', prio = 0,
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
  filter(from <= to) %>%
  # filter((from != to) | (x < y)) %>%
  filter(abs.near < 1e2) %>%
  # group_by(from, to) %>% top_n(-5, abs.near)  %>% slice(1:5) %>% View
  ggplot(aes(x = nearest)) +
  geom_histogram() +
  xlab('Distance to closest TSS') +
  facet_wrap(from ~ to, scale = 'free_y')

ggsave(file = 'analysis/05_term_comparison.pdf',
       width = 7, height = 7, units = 'in')


dat.term %>%
  mutate(
    start = start - ifelse(src != 'DBTBS', 50, 0),
    end = end + ifelse(src != 'DBTBS', 50, 0)
  ) -> dat.term2

term.over <- overlap_matching(dat.term2, dat.term2) %>%
  filter(!antisense) %>%
  select(x, y) %>%
  drop_na

nodes <- dat.term2 %>%
  transmute(i = 1:n(), name = id, i)
edges <- term.over %>%
  mutate(row = 1:n()) %>%
  gather('key', 'name', x, y) %>%
  left_join(nodes, 'name') %>%
  select(-name) %>%
  spread(key, i) %>%
  select(-row)

graph <- tbl_graph(nodes, edges, FALSE)
graph %>%
  activate(nodes) %>%
  mutate(group = group_components()) %>%
  as_tibble %>%
  left_join(dat.term2, c('name' = 'id')) %>%
  select(group, id = name, src, energy, prio, pubmed, start, end, strand) -> term.merge


term.merge %>%
  group_by(group) %>%
  top_n(-1, prio) %>%
  summarize(start = min(start), end = max(end),
            strand = clean_paste(strand)) -> term.pos
# count(term.pos, strand)
# ok, no strand problems
term.merge %>%
  select(group, src, energy, pubmed) %>%
  gather('key', 'value', src, energy, pubmed) %>%
  separate_rows(value, sep = ';') %>%
  drop_na %>%
  filter(value != '') %>%
  unique %>%
  group_by(group, key) %>%
  summarize_at('value',
               c(clean_paste,
                 function(i) i %>% as.numeric %>% min)) -> term.meta

term.meta %>%
  mutate(value = ifelse(key == 'energy', as.character(fn2), fn1)) %>%
  select(- starts_with('fn')) %>%
  spread(key, value, fill = '') %>%
  left_join(term.pos, 'group') %>%
  ungroup %>%
  arrange(start) %>%
  unique %>%
  transmute(id = paste0('BSGatlas-terminator-', 1:n()),
            start, end, strand, energy, src) -> bsg.term


pdf(file = 'analysis/05_venn_term.pdf')
bsg.term %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  mutate_at('src', fct_recode,
            'DBTBS (terminators)' = 'DBTBS',
            'BsubCyc\n(term.)' = 'BsubCyc',
            'Nicolas et al.\ndownshift' = 'Nicolas et al. downshift') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  with(set_names(map(i, 1), src)) %>%
  venn(cexil = 1.3,
       cexsn = 1.3,
       zcolor = ggsci::pal_jama()(3))
dev.off()


bsg.term %>%
  separate_rows(src, sep = ';') %>%
  count(src) %>%
  bind_rows(tibble(src = 'BSGatlas', n = nrow(bsg.term))) %>%
  mutate(src = fct_reorder(src, -n ) %>%
           fct_recode('Nicolas et al.' = 'Nicolas et al. downshift')) %>%
  ggplot(aes(x = src, y = n, fill = src)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama() +
  xlab(NULL) + ylab('count of transcription terminators') +
  theme_minimal(base_size = 16) +
  theme(legend.position = 'none')

ggsave(file = 'analysis/05_term_bars.pdf')

  
# save results
bsg.boundaries <- list(TSS = bsg.tss, terminator = bsg.term)
save(bsg.boundaries, file = 'analysis/05_bsg_boundaries.rda')
