source('analysis/00_load.R')
source('scripts/overlap_matching.R')
source('scripts/distance_matching.R')
load('data/01_nicolas.rda')
load('analysis/02_merging.rda')
load('analysis/03_tus.rda')

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

cmp %>%
  with(x.length + y.length - overlap  ==  (x5.dist + x3.dist + overlap)) %>%
  summary
# Mode    TRUE 
# logical    9588

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
    nice = sprintf('%s (%s)', src, priority) %>%
      fct_reorder(priority)
  ) -> nice.names

cmp %>%
  left_join(nice.names, 'src') %>%
  # S708 and S1202 are weired after changing the input to Nicolas et al. table S2
  filter(total.diff <= 500) %>%
  ggplot(aes(x = total.diff.cat, y = jaccard)) +
  geom_boxplot() +
  facet_wrap(~ nice) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Difference in coordinates to corresponding merged gene') +
  ylab('Jaccard Similarity')

ggsave(file = 'analysis/04_refinment_boxplot.pdf',
       width = 6, height = 6)

cmp %>%
  left_join(nice.names, 'src') %>%
  # S708 and S1202 are weired after changing the input to Nicolas et al. table S2
  filter(total.diff <= 500) %>%
  count(nice, total.diff.cat) %>%
  spread(total.diff.cat, n, fill = 0) -> stat

View(stat)

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
cmp.dist <- distance_matching(nic.utrs, merged)
cmp.over <- overlap_matching(nic.utrs, tus)

cmp.dist %>%
  filter(!antisense) %>%
  left_join(nic.utrs, c('x' = 'id')) %>%
  mutate(len = end - start + 1) %>%
  group_by(x, type) %>%
  summarize_at(c('distance', 'len'), min) %>%
  ungroup %>%
  mutate(
    bar = ifelse(distance == 0, 0, len + distance),
    dist.cut = cut(bar, 
                   # c(0, unlist(map(1:3, ~ 10 ** ..1)), Inf)) %>%
                   c(0, 100, 500, 1e3, 2e3, Inf)) %>%
      fct_explicit_na('overlaps') %>%
      fct_relevel('overlaps')
  ) %>%
  count(type, dist.cut) %>%
  left_join(count(nic.utrs, type), 'type') %>%
  mutate(nice = sprintf('%s (%s%%)', n.x, round(n.x / n.y * 100, 0))) %>%
  select(type, n = nice, dist.cut) %>%
  spread(type, n, fill = 0) %>%
  mutate(dist.cut = c('Overlapping', '0..100', '100..500',
                      '500..1,000',  '1,000..2,000',  '2,000+')) %>%
  # mutate(dist.cut = c('Overlapping', '1..10', '10..100', '100..1,000',
  #                     # '1,000..10,000', '10,000+')) %>%
  #                     '1,000+')) %>%
  rename('distance to closest gene' = dist.cut) %>%
  View

