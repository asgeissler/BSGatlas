
source('analysis/00_load.R')

source('scripts/overlap_matching.R')



load('analysis/01_refseq.rda')
load('analysis/01_bsubcyc.rda')
load('analysis/01_rfam.rda')
load('analysis/01_nicolas.rda')


# Hierarchy
# refseq coding

# bsubcyc coding

# nicolas lit review
# nicolas trusted
# rfam conservative
# refseq noncoding
# 
# bsubcyc noncoding
# 
# rfam medium
# nicolas lower
# 
# rfam sensetive

data <- bind_rows(
  transmute(refseq$coding,
            priority = 0, src = 'refseq coding',
            id = locus, start, end, strand, type),
  
  transmute(bsubcyc$coding,
            priority = 1, src = 'bsubcyc coding',
            id = locus, start, end, strand,
            type = ifelse(
              str_detect(description, 'putative|predicted'),
              'putative',
              'CDS'
            )),
  
  transmute(nicolas$high$their.lit_review,
            priority = 2, src = 'nicolas lit review',
            id = name, start, end, strand),
  transmute(nicolas$high$their.trusted,
            priority = 2, src = 'nicolas trusted',
            id = name, start, end, strand),
  transmute(rfam$conservative,
            priority = 2, src = 'rfam conservative',
            id = paste('row', 1:n()), start, end, strand, type),
  transmute(refseq$noncoding,
            priority = 2, src = 'refseq noncoding',
            id = locus, start, end, strand, type),
  
  transmute(bsubcyc$noncoding,
            priority = 3, src = 'bsubcyc noncoding',
            id = locus, start, end, strand, type),
  
  transmute(rfam$medium,
            priority = 4, src = 'rfam medium',
            id = paste('row', 1:n()), start, end, strand, type),
  transmute(nicolas$lower,
            priority = 4, src = 'nicolas lower',
            id = name, start, end, strand),
  
  transmute(rfam$sensetive,
            priority = 5, src = 'rfam sensetive',
            id = paste('row', 1:n()), start, end, strand, type)
) %>%
  mutate(type = case_when(
    type == 'small' ~ 'sRNA',
    type == 'putative' ~ 'putative-coding',
    type == 'unclear' ~ 'putative-non-coding',
    type == 'snoRNA'  ~ 'putative-non-coding',
    is.na(type)       ~ 'putative-non-coding',
    TRUE ~ type
  )) %>%
  mutate_at(c('start', 'end', 'priority'), as.integer) %>%
  mutate(seqnames = 'dummychr')

data %>%
  count(src, priority, type) %>%
  transmute(foo = sprintf('%s (%s)', src, priority), type, n) %>%
  spread(foo, n) %>%
  View

# Merging

# Quarantee that the helper separator has no conflicts
assertthat::assert_that(!any(str_detect(data$id, '%')))
assertthat::assert_that(!any(str_detect(data$src, '%')))

# ids are unique
assertthat::are_equal(data %>% nrow,
                      data %>% unique %>% nrow)

search <- data %>%
  transmute(
    id = paste(src, id, sep = '%'),
    start, end, strand, type, priority
  )

over <- overlap_matching(search, search) %>%
  filter(
    # only same sense comparisons
    ! antisense,
    # don't compare identity
    x != y,
    # remove symetry
    #x < y
  )


over %>%
  # mutate(
  #   c.sum = map2(str_detect(x.type, 'coding'),
  #                str_detect(y.type, 'coding'),
  #                sum),
  #   comparison = case_when(
  #     c.sum == 2 ~ 'coding-coding overlap',
  #     c.sum == 1 ~ 'coding-RNA overlap',
  #     c.sum == 0 ~ 'RNA-RNA overlap'
  #   )
  # ) %>%
  # drop_na(jaccard) %>%
  # group_by(x) %>%
  # top_n(1, jaccard) %>%
  # ungroup %>%
  mutate(
    `jaccard similarity` = cut(jaccard, seq(0, 1, 0.1), include.lowest = TRUE)
  ) %>%
  # pull(j.cut) %>% levels
  mutate_at(c('x.priority', 'y.priority'), function(i) {
    i %>% 
      as.character() %>%
      as.factor %>%
      fct_recode(
        'refseq coding' = '0',
        'bsubcyc coding' = '1',
        'conservative ncRNA\n(incl. refseq)' = '2',
        'bsubcyc ncRNA' = '3',
        'mediuam ncRNA' = '4',
        'sensetive rfam' = '5'
    )
  }) %>%
  ggplot(aes(x = `jaccard similarity`, fill = `jaccard similarity`)) +
  # ggsci::scale_fill_ucscgb() +
  scale_fill_brewer(palette = 'RdBu', direction = -1) +
  geom_bar() +
  xlab(NULL) +
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  # scale_y_log10() +
  # xlim(0, 1) +
  # geom_vline(xintercept = c(0.7, 0.8, 0.9)) +
  facet_grid(x.priority ~ y.priority, scales = 'free') +
  ggtitle('Comparison confidence levels', 
          'Similarities between each overlapping gene pair (coding and non-coding), identity is ignored')
          

ggsave(filename = 'analysis/02_level-overlaps.pdf',
       width = 25, height = 25, units = 'cm')

