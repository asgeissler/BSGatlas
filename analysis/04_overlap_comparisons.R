source('analysis/00_load.R')
source('scripts/overlap_matching.R')
load('analysis/02_merging.rda')
load('analysis/01_nicolas.rda')

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

cmp %>%
  left_join(nice.names, 'src') %>%
  count(nice, total.diff.cat) %>%
  spread(total.diff.cat, n, fill = 0) %>%
  View

# 2. UTRs
