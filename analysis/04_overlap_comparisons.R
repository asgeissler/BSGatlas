source('analysis/00_load.R')
source('scripts/overlap_matching.R')
source('scripts/distance_matching.R')
load('analysis/02_merging.rda')
load('analysis/01_nicolas.rda')
load('analysis/03_operons.rda')

# investigation of ammound of refinment due to merging
# and investigation of Nicolas et al UTRs


# 1. Refinment

merged <- merging$merged_genes %>%
  select(id = merged_id, start, end, strand)
indiv <- merging$merged_src %>%
  select(id = src_locus, start, end, strand, merged_id)

cmp <- overlap_matching(merged, indiv) %>%
  # only same direction
  filter(!antisense) %>%
  # keep only pairs of indiv genes and the merged gene
  semi_join(indiv, c('x' = 'merged_id', 'y' = 'id')) %>%
  separate(y, c('src', 'indiv_loci'), sep = '%', remove = FALSE)

# cmp %>%
#   with(x.length + y.length - overlap  ==  (x5.dist + x3.dist + overlap)) %>%
#   summary
# Mode    TRUE 
# logical    9643

cmp %<>%
  mutate(
    total.diff = x5.dist + x3.dist,
    total.diff.cat = cut(total.diff, c(0, 10, 50, 100, 250, 500)) %>%
      fct_explicit_na('No difference') %>%
      fct_relevel('No difference') %>%
      fct_recode('[1,10]' = '(0,10]')
  )

merging$merged_src %>%
  select(src, priority) %>%
  unique() %>%
  mutate(
    nice = case_when(
      src == "nicolas lower"      ~ 'Nicolas et al.\npredictions',
      src == "nicolas trusted"    ~ 'Nicolas et al.\ntrusted predictions',
      src == "dar riboswitches"   ~ 'Dar et al. riboswitches',
      src == "refseq coding"      ~ 'RefSeq Coding',
      src == "bsubcyc coding"     ~ 'BsubCyc Coding',
      src == "rfam conservative"  ~ 'Rfam, conservative',
      src == "rfam medium"        ~ 'Rfam, medium',
      src == "refseq noncoding"   ~ 'RefSeq Non-Coding',
      src == "bsubcyc noncoding"  ~ 'BsubCyc Non-Coding',
      src == "nicolas lit review" ~ 'Nicolas et al.\'s\nliterature review'
    ),
    nice = sprintf('%s (%s)', nice, priority) %>%
      fct_reorder(priority)
  ) -> nice.names

cmp %>%
  left_join(nice.names, 'src') %>%
  ggplot(aes(x = total.diff.cat, y = jaccard)) +
  geom_boxplot() +
  facet_wrap(~ nice) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Difference to corresponding merged gene [nt]') +
  ylab('Jaccard Similarity')

ggsave(file = 'analysis/04_refinment_boxplot.pdf')

cmp %>%
  left_join(nice.names, 'src') %>%
  count(nice, total.diff.cat) %>%
  spread(total.diff.cat, n, fill = 0) -> stat

View(stat)

stat %>%
  mutate_at('nice', str_replace, 'et al.', '\\\\emph{et al.}') %>%
  mutate_at('nice', str_remove, '\\n') %>%
  rename(Resource = nice) %>%
  kable('latex', escape = FALSE, caption = 'foo') %>%
  kable_styling(latex_options = 'scale_down') %>%
  str_split('\\n') %>%
  unlist %>%
  # thousand digit mark
  str_replace_all('(\\d)(\\d{3})', '\\1,\\2') %>%
  # without environment
  `[`(4:(length(.) - 1)) %>%
  write_lines(path = 'analysis/04_refinment_stat.tex')

# 2. UTRs

nicolas$all.features %>%
  drop_na(type) %>%
  filter(!str_detect(type, 'indep')) %>%
  # count(type)
  transmute(
    id = locus, start, end, strand, 
    type = case_when(
      type == 'intra' ~ 'internal',
      type == 'inter' ~ 'intergenic',
      startsWith(type, '3') ~ "3' UTR",
      TRUE ~ "5' UTR"
    )
  ) -> nic.utrs

nic.utrs
cmp.dist <- distance_matching(nic.utrs, operons$operon)
cmp.over <- overlap_matching(nic.utrs, operons$operon)

cmp.over %>%
  filter(!antisense) %>%
  drop_na(x) %>%
  # filter(x.type == 'intergenic') %>%
  group_by(x.type) %>%
  do(foo = cut(.$overlap / .$x.length * 100,
               seq(0, 100, length.out = 5),
               include.lowest = TRUE) %>%
      fct_explicit_na('not overlapping') %>%
      fct_relevel('not overlapping') %>%
      fct_count()
  ) %>%
  unnest %>%
  ggplot(aes(x = f, y = n )) +
  geom_bar(stat = 'identity') +
  xlab('overlap with operon relative to prediction length [%]') +
  ylab('count') +
  facet_wrap(~ x.type, scale = 'free_y')

ggsave(file = 'analysis/04_nic_overlaps.pdf',
       width = 9, height = 6, units = 'in')


# cmp.over %>%
#   filter(!antisense) %>%
#   filter(is.na(y)) %>%
#   select(x, type = x.type) %>%
#   left_join(cmp.dist, 'x') %>%
  # filter(!antisense) %>%
  # group_by(x) %>%
  # top_n(-1) %>%
  # ungroup %>%
#   count(type, mode)
# type       mode           n
# 3' UTR     x.after.y    220 ~95%
# 3' UTR     x.before.y    12
# 5' UTR     x.after.y     43
# 5' UTR     x.before.y   550 ~93%
# intergenic x.after.y    166
# intergenic x.before.y    69
# internal   x.after.y     25
# internal   x.before.y    27
cmp.dist %>%
  filter(!antisense) %>%
  left_join(nic.utrs, c('x' = 'id')) %>%
  group_by(x, type) %>%
  summarize_at('distance', min) %>%
  ungroup %>%
  mutate(
    dist.cut = cut(distance, 
                   c(0, unlist(map(1:4, ~ 10 ** ..1)), Inf)) %>%
      fct_explicit_na('overlaps') %>%
      fct_relevel('overlaps')
  ) %>%
  count(type, dist.cut) %>%
  left_join(count(nic.utrs, type), 'type') %>%
  mutate(nice = sprintf('%s (%s%%)', n.x, round(n.x / n.y * 100, 0))) %>%
  select(type, n = nice, dist.cut) %>%
  spread(type, n, fill = 0) %>%
  mutate(dist.cut = c('Overlapping', '1..10', '10..100', '100..1,000',
                      '1,000..10,000', '10,000+')) %>%
  rename('distance to closest operon' = dist.cut) %>%
  kable('latex', caption = 'foo') %>%
  kable_styling(latex_options = 'scale_down') %>%
  str_split('\\n') %>%
  unlist %>%
  # without environment
  `[`(4:(length(.) - 1)) %>%
  write_lines(path = 'analysis/04_nic_dist_stat.tex')
