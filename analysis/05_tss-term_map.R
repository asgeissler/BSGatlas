source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

source('scripts/frame_helpers.R')

load('data/01_bsubcyc.rda')
load('data/01_nicolas.rda')
load('data/03_dbtbs.rda')

library(tidygraph)

##############################################################################
# 0. Collect input

dat <- bind_rows(
  # Part 1: TSS
  bind_rows(
    nicolas$upshifts %>%
      transmute(src = 'Nicolas et al. upshift', id, TSS = pos, strand, sigma,
                without.tu,
                pubmed = '') %>%
      mutate_at('sigma', str_remove, 'Sig'),
    dbtbs$tss %>%
      transmute(src = 'DBTBS', id = 1:n(), TSS, strand, sigma,
                without.tu = FALSE,
                pubmed = reference) %>%
      mutate_at('id', as.character) %>%
      mutate_at('sigma', str_remove, 'Sig'),
    bsubcyc$TSS %>%
      transmute(src = 'BsubCyc', id = 1:n(), TSS = tss, strand, sigma,
                without.tu = FALSE,
                pubmed = cite) %>%
      mutate_at('id', as.character)
  ) %>%
    mutate(src_id = paste0(src, '%', id)) %>%
    rename(extra = sigma) %>%
    mutate(type = 'TSS', start = TSS, end = TSS) %>%
    select(-TSS),
  # Part 2: Terminators
  bind_rows(
    bsubcyc$terminator %>%
      arrange(start) %>%
      unique %>%
      transmute(id = 1:n(), start, end, strand, energy,
                without.tu = FALSE,
                src = 'BsubCyc') %>%
      mutate_at(c('id', 'energy'), as.character),
    dbtbs$term %>%
      arrange(start) %>%
      transmute(id = 1:n(), start, end, strand, energy = energies,
                without.tu = FALSE,
                src = 'DBTBS',
                pubmed = reference) %>%
      mutate_at('id', as.character),
    nicolas$downshifts %>%
      transmute(id, start = pos, end = pos, strand, energy,
                without.tu,
                src = 'Nicolas et al. downshift') %>%
      mutate_at('energy', as.character)
  ) %>%
    mutate_at('energy', replace_na, '') %>%
    drop_na(start) %>%
    mutate(src_id = paste0(src, '%', id)) %>%
    rename(extra = energy) %>%
    mutate(type = 'Terminator')
)
##############################################################################
# 1. find resolutions

dat  %>%
  select(id = src_id, start, end, strand) %>%
  distance_matching(., .) %>%
  filter(!antisense) %>%
  filter(x != y) %>%
  left_join(select(dat, x = src_id, x.type = type, x.src = src), 'x') %>%
  left_join(select(dat, y = src_id, y.type = type, y.src = src), 'y') %>%
  filter(x.type == y.type) %>%
  rename(type = x.type) %>%
  mutate(mode = ifelse(x.src == y.src, 'within', 'between')) %>%
  select(- y.type) -> dat.cmp


dat.cmp %>%
  filter(((mode == 'between') & (x.src < y.src)) | (x < y)) %>%
  group_by(type, mode, x, x.src, y.src) %>%
  summarize_at('distance', min) %>%
  ungroup -> dist.stat


cowplot::plot_grid(
  dist.stat %>%
    filter(mode == 'within') %>%
    mutate_at('x.src', str_remove, ' [updown]*shift$') %>%
    ggplot(aes(distance, col = x.src)) +
    stat_ecdf(size = 1.5) +
    scale_x_continuous(
      breaks = c(10, 25, 50, 100, 150),
      limits = c(0, 150)
    ) +
    ggsci::scale_color_jama(name = NULL) +
    scale_y_continuous(breaks = c(0, .25, .5, .75, .85, .9, 1)) +
    xlab('Distance to closest neighboring annotation [bp]') +
    ylab('Empirical density') +
    facet_wrap(~ type) +
    theme_bw(18) +
    theme(legend.position = 'bottom'),
  dist.stat %>%
    filter(mode == 'between') %>%
    mutate_at('y.src', str_remove, ' [updown]*shift$') %>%
    mutate(short = sprintf('%s vs\n%s', x.src, y.src)) %>%
    ggplot(aes(distance)) +
    stat_ecdf(size = 1.5) +
    scale_x_continuous(
      breaks = c(10, 25, 50, 100, 150),
      limits = c(0, 150)
    ) +
    ggsci::scale_color_jama(name = 'Resource') +
    scale_y_continuous(breaks = c(0, .25, .5, .75, .85, .9, 1)) +
    geom_hline(yintercept = .9, color = 'red') +
    xlab('Distance to closest annotation in the other resource [bp]') +
    ylab('Empirical density') +
    facet_grid(type ~ short) +
    theme_bw(18) +
    theme(legend.position = 'bottom'),
  labels = 'AUTO',
  ncol = 2,
  rel_widths = c(1.5, 3),
  label_size = 24
)


ggsave('analysis/05_resolution.pdf',
       width = 18, height = 8)
##############################################################################
# 2. Merging

# Make windows for Nicolas
# Then merge by overlaps, similar code as for gene merging
dat %>%
  mutate(
    res.limit = ifelse(str_detect(src, 'Nicolas'), 22, 0),
    start = start - res.limit,
    end = end - res.limit
  ) %>%
  transmute(
    i = 1:n(),
    id = src_id, type, src, 
    start, end, strand,
    extra, without.tu
  ) -> merge.that
  
merge.that %>%
  overlap_matching(., .) %>%
  filter(x != y) %>%
  filter(!antisense) %>%
  # prevent mixture of TSS/Terminators
  filter(x.type == y.type) %>%
  select(from = x.i, to = y.i) %>%
  drop_na -> edges

grp <- tbl_graph(nodes = merge.that,
                 edges = edges,
                 directed = FALSE)

grp %>%
  activate(nodes) %>%
  mutate(group = group_components()) %>%
  filter(group == 1) %>%
  plot(label = id)
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
  mutate_at('src', fct_recode, "Nicolas et al.\nUpshift" = "Nicolas et al upshift") %>%
  mutate_at('src', fct_recode, "Nicolas et al.\n5'UTR start" = "Nicolas et al 5' UTR start") %>%
  mutate_at('group', fct_relevel,
            "SigA", "other", "SigF", "unknown", "SigE", "SigK", "SigG") %>%
  ggplot(aes(x = src, fill = group, y = prop)) +
  geom_bar(stat = 'identity') +
  xlab(NULL) + ylab('Proprotion [%]') +
  ggsci::scale_fill_jama(name = 'Sigma Factor') +
  theme_minimal(base_size = 16)
ggsave(file = 'analysis/05_sigma_proportion.pdf',
       width = 20, height = 12, units = 'cm')
  
bsg.tss %>%
  count(res.limit)
# res.limit     n
#         0   706
#        22  2679
#        33    12
bsg.tss %>%
  separate_rows(src, sep = ';') %>%
  count(src) %>%
  bind_rows(tibble(src = 'BSGatlas', n = nrow(bsg.tss))) %>%
#  BsubCyc                      556
#  DBTBS                        644
#  Nicolas et al 5' UTR start   690
#  Nicolas et al upshift       3264
#  BSGatlas                    3397
  mutate(src = fct_reorder(src, - n)) %>%
  mutate_at('src', fct_recode, "Nicolas et al.\nUpshift" = "Nicolas et al upshift") %>%
  mutate_at('src', fct_recode, "Nicolas et al.\n5'UTR start" = "Nicolas et al 5' UTR start") %>%
  ggplot(aes(x = src, y = n, fill = src)) + 
  ggsci::scale_fill_jama(name = NULL) +
  geom_bar(stat = 'identity') +
  theme_minimal(base_size = 16) +
  scale_y_continuous(breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500)) +
  ylab('count of sigma factor binding sites') + xlab(NULL) +
  theme(legend.position = 'none')

ggsave(file = 'analysis/05_count_promoter.pdf',
       width = 30, height = 12, units = 'cm')

# venn diagram

library(venn)
pdf(file = 'analysis/05_venn_tss.pdf')
bsg.tss %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  mutate_at('src', fct_recode,
            'DBTBS (TSS)' = 'DBTBS',
            'BsubCyc\n(TSS)' = 'BsubCyc',
            "Nicolas et al.\n5'UTR start" = "Nicolas et al 5' UTR start",
            'Nicolas et al.\nupshift' = 'Nicolas et al upshift') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  with(set_names(map(i, 1), src)) %>%
  venn(cexil = 1.3,
       cexsn = 1.3,
       zcolor = ggsci::pal_jama()(4))
dev.off()


# FOO

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
            "Nicolas et al.\n3'UTR end" = "Nicolas et al 3' UTR end",
            'Nicolas et al.\ndownshift' = 'Nicolas et al. downshift') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  with(set_names(map(i, 1), src)) %>%
  venn(cexil = 1.3,
       cexsn = 1.3,
       ilabels = TRUE,
       zcolor = ggsci::pal_jama()(4))
dev.off()


bsg.term %>%
  separate_rows(src, sep = ';') %>%
  count(src) %>%
  bind_rows(tibble(src = 'BSGatlas', n = nrow(bsg.term))) %>%
  mutate(src = fct_reorder(src, -n )) %>%
  mutate_at('src', fct_recode, "Nicolas et al.\nDownshift" = "Nicolas et al. downshift") %>%
  mutate_at('src', fct_recode, "Nicolas et al.\n3'UTR end" = "Nicolas et al 3' UTR end") %>%
  ggplot(aes(x = src, y = n, fill = src)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama() +
  xlab(NULL) + ylab('count of transcription terminators') +
  theme_minimal(base_size = 16) +
  theme(legend.position = 'none')

ggsave(file = 'analysis/05_term_counts.pdf',
       width = 30, height = 12, units = 'cm')

  
# save results
bsg.boundaries <- list(TSS = bsg.tss, terminator = bsg.term)
save(bsg.boundaries, file = 'analysis/05_bsg_boundaries.rda')


bsg.boundaries %>%
  map(separate_rows, src, sep = ';') %>%
  map(count, src) %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() %>%
  mutate(src = case_when(
    str_detect(src, 'UTR') ~ "Nicolas et al. 5'/3' UTR ends",
    str_detect(src, 'Nicolas') ~ 'Nicolas et al',
    TRUE ~ src
  )) %>%
  bind_rows(
    bsg.boundaries %>%
      map(nrow) %>%
      map2(names(.), ~ tibble(type = .y, n = .x)) %>%
      bind_rows %>%
      mutate(src = 'BSGatlas')
  ) %>%
  spread(type, n) %>%
  kable('latex', booktabs = TRUE) %>%
  kable_styling() %>%
  strsplit('\n') %>%
  unlist %>%
  `[`(2:(length(.) - 1)) %>%
  write_lines('analysis/05_overview.tex')
