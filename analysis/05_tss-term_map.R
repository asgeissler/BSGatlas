source('analysis/00_load.R')
source('scripts/distance_matching.R')
source('scripts/overlap_matching.R')

source('scripts/frame_helpers.R')

load('data/01_bsubcyc.rda')
load('data/01_nicolas.rda')
load('data/03_dbtbs.rda')

library(tidygraph)
library(venn)
library(gridGraphics)

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
  transmute(
    i = 1:n(),
    id = paste0(type, '_', src_id),
    type, src, 
    res.limit = ifelse(str_detect(src, 'Nicolas'), 22, 0),
    start = start - res.limit,
    end = end + res.limit,
    strand,
    extra, pubmed, without.tu
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
  morph(to_subgraph, !str_detect(src, 'Nicolas')) %>%
  mutate(sub.group = group_components()) %>%
  unmorph() %>%
  as_tibble %>%
  # fill-up NAs with distinged placeholder to prevent shifts form beeing merged
  mutate(sub.group = ifelse(is.na(sub.group),
                            paste0('placeholder-', 1:n()),
                            sub.group)) -> sub.groups

# Merging depends on type
# Option A - TSS:
#   * Per group, detect and keep those of lowest resolution
#     (also solved via sub-groups)
#   * Per position within group, combine pubmed/sigma info
#   * Add sigma info of the potential overlapping upshift
#
# Option B - Terminator:
#   * Detect sug-groups, when downshift is exlcuded
#   * Unify per sub-group (incl.pubmed)
#   * compute average energy (including downshift info)
# A) unite sigma factor information

f <- partial(str_c, collapse = ';')
sub.groups %>%
  mutate_at('src', str_remove, ' [updown]*shift$') %>%
  group_by(type, group, sub.group) %>%
  summarize(
    src = f(src),
    res.limit = f(res.limit),
    start = min(start),
    end = max(end),
    strand = f(strand),
    wihtout.tu = invoke(all, without.tu),
    extra = f(extra),
    pubmed = f(pubmed)
  ) %>%
  ungroup %>% 
  mutate(id = paste0('tmp-', 1:n())) -> pre.merged

pre.merged %>% count(strand)
# -> no strand confusion, check

# Only keep up/down-shifts if there was no more precice annotation
# -> groups without subgroups

pre.merged %>%
  left_join(
    sub.groups %>%
      transmute(group, x = str_detect(sub.group, 'placeholder')) %>%
      unique %>%
      count(group) %>%
      filter(n == 1) %>%
      transmute(group, allow.shift = TRUE),
    'group'
  ) %>%
  mutate_at('allow.shift', replace_na, FALSE) %>%
  filter(!str_detect(sub.group, 'placeholder') | allow.shift) -> filtered.merged

# Add back info from the upshifts

filtered.merged %>%
  left_join(
    anti_join(pre.merged, filtered.merged) %>%
      select(group, src, nic.extra = extra),
    'group'
  ) %>%
  # Then clean up
  select(- group, - sub.group) %>%
  # okay to do hear, because res.limit and strand are unambiguous
  mutate_at(c('res.limit', 'strand'), str_remove, ';.*$') %>%
  mutate_at('pubmed', replace_na, '') %>%
  gather('k', 'v', src.x, src.y, extra, nic.extra, pubmed) %>%
  mutate_at('k', str_remove, '.x$') %>%
  mutate_at('k', str_remove, '.y$') %>%
  mutate_at('k', str_remove, '^nic.') %>%
  separate_rows(v, sep =';') %>%
  group_by_at(vars(- v)) %>%
  summarise(
    clean = clean_paste(v),
    avg = v %>% as.double() %>% mean(na.rm = TRUE)
  ) %>%
  ungroup %>%
  mutate(v = ifelse(is.na(avg) | (k == 'pubmed'), clean, avg)) %>%
  select(- clean, -avg) %>%
  spread(k, v) -> clean.merge

# Little bit more work for the sigma factors and energies
# But now need to separate TSS and terminators

clean.merge %>%
  group_by(type) %>%
  do(i = list(.)) %>%
  with(set_names(i, type)) %>%
  map(1) %>%
  map(select, - type, - allow.shift) -> merged

merged %>%
  map(count, wihtout.tu)
  
# energies 
merged$Terminator %<>%
  mutate_at('extra', as.double) %>%
  rename(energy = extra) %>%
  select(- res.limit)

# sigma
merged$TSS %>%
  rename(sigma = extra) %>%
  mutate(sigma = ifelse(sigma == '', '?', sigma)) %>%
  mutate_at('sigma', str_replace_all, '-', '?') %>%
  #tmp substitute YlaC -> C
  mutate_at('sigma', str_replace_all, 'YlaC', 'C') %>%
  mutate_at('sigma', str_remove_all, ';') %>%
  mutate_at('sigma', str_replace_all, '([A-Z?])(?=[A-Z?])', '\\1;') %>%
  separate_rows(sigma, sep = ';') %>%
  # and C -> YlaC back transform
  mutate(sigma = ifelse(sigma == 'C', 'YlaC', sigma)) %>%
  unique %>%
  group_by_at(vars(- sigma)) %>%
  summarize_at('sigma', ~ invoke(paste, sort(.x), sep = ';')) %>%
  ungroup -> merged$TSS

# undo res limit
merged$TSS %<>%
  mutate_at('res.limit', as.integer) %>%
  mutate(
    start = start + res.limit,
    end = end - res.limit
  )
##############################################################################
#The naming everything backwards again part

# merged$TSS %>% with(end - start + 1) %>% table()
# sth went wrong?, fixed again

'data-gff/BSGatlas_v1.0.gff' %>%
  rtracklayer::import.gff3() %>%
  as_tibble() %>%
  filter(type %in% c('TSS', 'terminator')) %>%
  select(type, id = ID, start, end, strand) %>%
  mutate_at('type', fct_recode, 'Terminator' = 'terminator') %>%
  mutate_at('type', as.character) -> legacy

legacy %>%
  select(type, i = id) %>%
  mutate(i = i %>%
           strsplit('-') %>%
           map(3) %>%
           unlist %>%
           as.integer) %>%
  group_by(type) %>%
  summarize_all(max) %>%
  with(set_names(i, type)) -> counts.status

merged %>%
  map(select, start, end, strand, id) %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() -> qs

overlap_matching(qs, legacy) %>%
  filter(x.type == y.type, !antisense) %>%
  group_by(y) %>%
  filter(jaccard == max(jaccard)) %>%
  ungroup -> foo

# Only keep unambiguous cases
foo %>%
  semi_join(
    foo %>%
      count(y) %>%
      filter(n == 1),
    'y'
  ) %>%
  semi_join(
    foo %>%
      count(x) %>%
      filter(n == 1),
    'x'
  ) %>%
  # count(mode)
  # also happen to be only the equal cases
  select(id = x, new.id = y) -> look

# continue with incremental counts on rest
merged %>%
  map(left_join, look, 'id') -> bar

map2(
  bar %>%
    map(filter, is.na(new.id)) %>%
    map2(
      names(.),
      ~ mutate(.x, new.id = paste(
        'BSGatlas', .y, counts.status[.y] + 1:n(),
        sep = '-'
      ))) %>%
    map(~ mutate_at(.x, 'new.id', str_replace, 'Terminator', 'terminator')),
  bar %>%
    map(drop_na),
  bind_rows
) %>%
  map(arrange, start) %>%
  map(select, - id) %>%
  map(rename, id = new.id) -> named

legacy %>%
  select(id, type) %>%
  anti_join(look, c('id' = 'new.id')) -> obsolete

##############################################################################
# save results (in a legacy naming of small `t`)
bsg.boundaries <- list(TSS = named$TSS,
                       terminator = named$Terminator,
                       obsolete = obsolete)
save(bsg.boundaries, file = 'analysis/05_bsg_boundaries.rda')

##############################################################################
##############################################################################
##############################################################################
# Some finalizing stat and figure for manuscript

merged$TSS %>%
  count(res.limit)
# res.limit     n
#         0   706
#        22  2684


load('data/01_bsub_raw.rda')

bsg.boundaries[c('TSS', 'terminator')]  %>%
  map2(c('TSS', 'TTS'), ~ mutate(.x, type = .y)) %>%
  map(select, id, type, src) %>%
  bind_rows -> dat

dat %>%
  separate_rows(src, sep = ';') %>%
  unique %>%
  bind_rows(
    dat %>%
      mutate(src = 'BSGatlas')
  ) %>%
  arrange(id) -> dat.full

dat.full %>%
  count(type, src) -> dat.count




list('TSS', 'TTS') %>%
  set_names(., .) %>%
  map(function(i) {
    dat.count %>%
      mutate(src = fct_reorder(src, n, max) %>%
               fct_rev) %>%
      mutate(lab = str_replace(n, '(\\d)(\\d{3})',  '\\1,\\2')) %>%
      filter(type == i) %>%
      ggplot(aes(x = src, y = n, fill = src, group = src)) +
      geom_bar(stat = 'identity',
               position = 'dodge') +
      geom_text(aes(label = lab),
                # size = 6,
                position = position_dodge(width=0.9),
                vjust=-0.25) +
      scale_fill_manual(values = ggsci::pal_jama()(5)[-2]) +
      theme_minimal(14) +
      # theme_bw(14) +
      scale_y_continuous(breaks = NULL) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      ylab(ifelse(i == 'TSS',
                  'Nr. annotated TSSs',
                  'Nr. annotated TTSs')) +
      xlab(NULL) +
      theme(legend.text = element_text(face = "italic"),
            legend.position = 'none',
            legend.justification = c(1, 1))
  }) -> p1

###############################################################################
# comparison via venn diagrams

helper <- function(i) {
  dat.full %>%
    filter(src != 'BSGatlas') %>%
    filter(type == i) %>%
    mutate_at('src', fct_recode, 'Nicolas\net al.' = 'Nicolas et al.') %>%
    group_by(src) %>%
    do(i = list(.$id)) %>%
    with(set_names(map(i, 1), src)) %>%
    `[`(c('Nicolas\net al.', 'BsubCyc', 'DBTBS')) %>%
    venn(ilcs = 1.1,
         sncs = 1.1,
         opacity = 0.6,
         zcolor =  ggsci::pal_jama()(5)[-c(1, 2)])
  grid.echo()
  # grid.ls()
  grid.remove('graphics-plot-1-lines-1')
  size <- 15
  cowplot::ggdraw() + cowplot::draw_grob(
    grid.grab(),
    scale = 0.8
    # width = (size + 1)/2.54,
    # height = (size + 1)/2.54
  )
}

##############################################################################
# comparison of sigma factor distributions

# getting the raw dist from bsubcyc is annoying

bsub_raw$proteins %>%
  select(name = `COMMON-NAME`, other = SYNONYMS,
         box = `RECOGNIZED-PROMOTERS`) %>%
  drop_na(box) %>%
  separate_rows(box, sep = ';') -> bs.box

bsub_raw$promoters %>%
  select(pro.id = `UNIQUE-ID`, box = `PROMOTER-BOXES`) %>%
  separate_rows(box, sep = ';') %>%
  left_join(bs.box, 'box') %>%
  mutate(sigma = str_sub(other, -1)) %>%
  mutate_at('sigma', replace_na, '?') %>%
  transmute(id = pro.id, src = 'BsubCyc', sigma) -> bs

bind_rows(
  nicolas$upshifts %>%
    transmute(src = 'Nicolas et al.', id, sigma) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  dbtbs$tss %>%
    transmute(src = 'DBTBS', id = 1:n(), sigma) %>%
    mutate_at('id', as.character) %>%
    mutate_at('sigma', str_remove, 'Sig'),
  bs,
  bsg.boundaries$TSS %>%
    transmute(id, sigma, src = 'BSGatlas')
) %>%
  mutate_at('sigma', replace_na, '?') %>%
  mutate_at('sigma', str_replace, '-', '?') %>%
  mutate_at('sigma', str_replace_all, '([A-Z])(?=[A-Z])', '\\1;') %>%
  separate_rows(sigma, sep = ';') %>%
  mutate(sigma = ifelse(sigma == 'C', 'YlaC', sigma)) %>%
  mutate(
    group = ifelse(sigma %in% c('?', 'A', 'E','F', 'G', 'K'),
                   sigma,
                   'other')
  ) -> tss.dat


tss.dat %>%
  count(src, group) %>%
  left_join(
    tss.dat %>%
      select(id, src) %>%
      # unique() %>%
      # nope, one position can have multiple binding sites
      count(src) %>%
      rename(total = n),
    'src'
  ) %>%
  mutate(prop = n / total * 100) -> tss.prop

tss.prop %>%
  mutate(group = fct_reorder(group, prop, max) %>% fct_rev) %>%
  mutate(group = fct_recode(group, 'unknown' = '?')) %>%
  mutate_at('src', fct_relevel,
            'BSGatlas', 'Nicolas et al.', 'BsubCyc', 'DBTBS') %>%
  mutate_at('group', fct_relevel,
            'unknown', 'A', 'E', 'F', 'G', 'K', 'other') %>%
  ggplot(aes(x = src, fill = group, y = prop)) +
  geom_bar(stat = 'identity') +
  xlab(NULL) + ylab('Proprotion [%]') +
  # ggsci::scale_fill_jco(name = 'Sigma Factor') +
  scale_fill_brewer(palette = 'RdYlBu', name = 'Sigma Factor') +
  # theme_bw(base_size = 14) -> p2
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) -> p2


##############################################################################
# put it all togehter

cowplot::plot_grid(
  cowplot::plot_grid(
    p1$TSS, p1$TTS,
    labels = c('(a)', '(b)'),
    ncol = 2
  ),
  cowplot::plot_grid(
    helper('TSS'), helper('TTS'), p2,
    # function() helper('TSS'),
    # partial(helper, 'TSS'), partial(helper, 'TTS'), p2,
    labels = c('(c) TSS', '(d) TTS', '(e)'),
    rel_widths = c(1, 1, 2),
    ncol = 3, hjust = 0
  ),
  ncol = 1
)

ggsave('analysis/05_stat.pdf', 
       width = 30, height = 20, units = 'cm')

