# Create statistics and figures for the computed isoforms.
#
# Other questions to adress are:
# How many genes have still no transcript
# How many are due to the UTRs
# How many paths have the same genes, but just different 5/3 UTRs
###############################################################################

source('scripts/frame_helpers.R')
source('analysis/00_load.R')

load('data/01_bsubcyc.rda')
load('data/03_subtiwiki.rda')
load('analysis/02_merging.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/07_isoforms.rda')

load('analysis/06_utrs.rda')
load('data/01_nicolas.rda')


#####################################


isoforms$tus %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  count(src) %>%
  mutate_at('src', fct_recode, 'Novel TUs' = 'BSGatlas') %>%
  bind_rows(
    tibble(src = 'Combined',
           n = sum(isoforms$tus$src != 'BSGatlas')),
    tibble(src = 'BSGatlas', n = nrow(isoforms$tus))
  ) %>%
  mutate_at('src',
            fct_relevel,
            'DBTBS', 'BsubCyc', 'SubtiWiki',
            'Combined',
            'Novel TUs',
            'BSGatlas') %>%
  arrange(src) %>%
  rename(`Transcriptional Units` = n) -> tu.stat
            
            

isoforms$operons %>%
  select(id, TUs) %>%
  separate_rows(TUs, sep = ';') %>%
  left_join(
    isoforms$tus %>%
      select(TUs = id, src) %>%
      separate_rows(src, sep = ';'),
    'TUs'
  ) %>%
  select(id, src) %>%
  unique -> op.src

op.src %>%
  count(src) %>%
  mutate_at('src', fct_recode, 'Novel TUs' = 'BSGatlas') %>%
  bind_rows(
    tibble(src = 'Combined',
           n = op.src %>%
             filter(src != 'BSGatlas') %>%
             select(id) %>%
             unique %>%
             nrow),
    tibble(src = 'BSGatlas', n = nrow(isoforms$operons))
  ) %>%
  mutate_at('src',
            fct_relevel,
            'DBTBS', 'BsubCyc', 'SubtiWiki',
            'Combined',
            'Novel TUs',
            'BSGatlas') %>%
  arrange(src) %>%
  rename(`Operons` = n) -> op.stat

c(
  BsubCyc = nrow(bsubcyc$transunit),
  # SubtiWiki = subtiwiki$transcripts %>%
  #   select(transcript) %>%
  #   unique %>%
  #   nrow,
  BSGatlas = nrow(isoforms$transcripts)
) %>%
  map2(names(.), ~ tibble(src = .y, Transcripts = .x)) %>%
  bind_rows -> trans.stat

N <- nrow(merging$merged_genes)
isoforms$tus %>%
  select(id, src, genes) %>%
  separate_rows(src, sep = ';') %>%
  separate_rows(genes, sep = ';') %>%
  select(src, genes) %>%
  unique %>%
  count(src) %>%
  mutate_at('src', fct_recode, 'Novel TUs' = 'BSGatlas') %>%
  bind_rows(
    tibble(
      src = 'Combined',
      n = isoforms$tus %>%
        filter(src != 'BSGatlas') %>%
        separate_rows(genes, sep = ';') %>%
        select(genes) %>%
        unique %>%
        nrow
    ),
    tibble(
      src = 'BSGatlas',
      n = isoforms$tus %>%
        separate_rows(genes, sep = ';') %>%
        select(genes) %>%
        unique %>%
        nrow
    )
  ) %>%
  mutate('% genes with transcripts' = 100 * n / N) %>%
  select(-n) -> genome.stat


  
bsg.boundaries %>%
  map(select, id, src) %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows %>%
  separate_rows(src, sep = ';') %>%
  mutate(src = ifelse(str_detect(src, 'Nicolas'), 'Nicolas et al.', src)) %>%
  unique %>%
  count(src, type) %>%
  spread(type, n) %>%
  bind_rows(
    bsg.boundaries %>%
      map(nrow) %>%
      as_tibble %>%
      mutate(src = 'BSGatlas')
  ) -> bound.stat

bind_rows(
  nicolas$all.features %>%
    mutate(type = case_when(
      startsWith(type, "3'") ~ "3'UTR",
      startsWith(type, "5'") ~ "5'UTR",
      type == 'intra' ~ 'internal_UTR'
    )) %>%
    drop_na(type) %>%
    count(type) %>%
    spread(type, n) %>%
    mutate(src = 'Nicolas et al.'),
  UTRs %>%
    map(nrow) %>%
    as_tibble %>%
    mutate(src = 'BSGatlas')
) -> utr.stat

op.stat %>%
  left_join(tu.stat, 'src') %>%
  left_join(trans.stat, 'src') %>%
  left_join(genome.stat, 'src') -> grand.stat

grand.stat %>%
  mutate_at('% genes with transcripts', round, 2) %>%
  mutate_all(as.character) %>%
  mutate_at('src', as_factor) %>%
  mutate_at('Transcripts', replace_na, '-') %>%
  mutate_at(c( "Operons", "Transcriptional Units", "Transcripts"),
            str_replace, '(\\d)(\\d{3})', '\\1,\\2') %>%
  rename(Resource = src) %>%
  knitr::kable('latex', linesep = '', booktabs = TRUE) %>%
  kable_styling() %>%
  strsplit('\\n') %>%
  unlist -> xs

# add midrule line
c(xs[1:10], '\\midrule', xs[11:length(xs)]) %>%
  `[`(3:15) %>%
  write_lines('analysis/08_operon_table.tex')
            
#####################################

bound.stat %>%
  left_join(utr.stat, 'src') %>%
  # mutate_at('src', fct_recode, 
  #           'SubtiWiki' = 'Nicolas et al.') %>%
  gather('key', 'value', - src) %>%
  spread(src, value) %>%
  arrange(desc(BSGatlas)) %>%
  select(key,
         DBTBS, BsubCyc, 
         `Nicolas et al.`,
          BSGatlas) %>%
  mutate_at('key', fct_recode,
            'internal UTR' = 'internal_UTR') %>%
  rename(' ' = key) %>%
  mutate_all(as.character) %>%
  mutate_all(replace_na, '-') %>%
  mutate_all(str_replace, '^(\\d)(\\d{3})$', '\\1,\\2') %>%
  # rename('Nicolas et al. / SubtiWiki' = SubtiWiki) %>%
  knitr::kable('latex', booktabs = TRUE) %>%
  strsplit('\\n') %>%
  unlist %>%
  write_lines('analysis/08_utr_table.tex')

#####################################
#####################################
# Is the operon assumption given?
# Each operon has a TU that contains all genes?
isoforms$operons %>%
  select(id, TUs) %>%
  separate_rows(TUs, sep = ';') %>%
  left_join(isoforms$tus, c('TUs' = 'id')) %>%
  separate_rows(genes, sep = ';') %>%
  select(id, TUs, genes) -> op.genes

op.genes %>%
  unique %>%
  count(id, TUs) %>%
  group_by(id) %>%
  summarise(most.genes = max(n)) %>%
  left_join(
    op.genes %>%
      select(-TUs) %>%
      unique %>%
      count(id),
    'id'
  ) %>%
  mutate(given = most.genes == n) %>%
  count(given) %>%
  mutate(pct = n / nrow(isoforms$operons) * 100)

# given      n    pct
# 1 FALSE    6   0.260
# 2 TRUE   2303  99.7


######################################
# Attempt similar classification to
# Conway T, et al., mBio. 2014


isoforms$operons %>%
  select(id, transcripts) %>%
  separate_rows(transcripts, sep = ';') %>%
  left_join(isoforms$transcripts, c('transcripts' = 'id')) %>%
  mutate(
    has.five = str_detect(features, "5'UTR"),
    has.three = str_detect(features, "3'UTR")
  ) %>%
  group_by(id) %>%
  summarize(
    has.five = any(has.five),
    has.three = any(has.three),
    n.TSS = TSS %>%
      unique %>%
      discard(is.na) %>%
      length,
    n.term = Terminator %>%
      unique %>%
      discard(is.na) %>%
      length,
    n.trans = transcripts %>%
      unique %>%
      length
  ) %>%
  left_join(
    isoforms$operons %>%
      select(id, TUs) %>%
      separate_rows(TUs, sep = ';') %>%
      count(id) %>%
      rename(n.tu = n),
    'id'
  ) %>% 
  left_join(
    isoforms$operons %>%
      select(id, TUs) %>%
      separate_rows(TUs, sep = ';') %>%
      left_join(isoforms$tus, c('TUs' = 'id')) %>%
      separate_rows(genes, sep = ';') %>%
      select(id, genes) %>%
      unique %>%
      count(id) %>%
      rename('# genes' = n),
    'id'
  ) -> indiv.stat


indiv.stat %>%
  mutate(
    class.TSS = case_when(
      n.TSS > 1 ~ 'Multi-TSS',
      n.TSS == 1 ~ 'Single-TSS',
      n.TSS == 0 ~ 'Missing-TSS'
    ),
    class.term = case_when(
      n.term > 1 ~ 'Multi-Term.',
      n.term == 1 ~ 'Single-Term.',
      n.term == 0 ~ 'Missing-Term.'
    ),
    class.tu = case_when(
      n.tu == 1 ~ 'Single TU',
      between(n.tu, 2, 10) ~ 'Multi TUs',
      TRUE ~ '> 10 TUs'
    ),
    class.gene = case_when(
      `# genes` == 1 ~ 'Single gene',
      between(`# genes`, 2, 10) ~ 'Multi genes',
      TRUE ~ '> 10 genes'
    )
  ) %>%
  # select(id, starts_with('class.')) %>%
  # group_by_at(vars(-id)) %>%
  # count %>%
  # View
  mutate(
    class = case_when(
    #   (class.term != 'Multi-Term.') & (class.TSS != 'Multi-TSS') &
    #     (class.gene == 'Single gene') ~ 'Simple',
    #   (class.term != 'Multi-Term.') & (class.TSS != 'Multi-TSS') &
    #     (class.gene != 'Single gene') ~ 'Traditional',
    #   (class.term == 'Multi-Term.') & (class.TSS == 'Multi-TSS') ~ 'Multi-both',
    #   (class.term == 'Multi-Term.') ~ 'Multi-Term.',
    #   (class.TSS == 'Multi-TSS') ~ 'Multi-TSS'
      # (class.term == 'Missing-Term.') | (class.TSS == 'Missing-TSS') ~ 'Misses Promoter and/or Term.',
      (class.term == 'Multi-Term.') & (class.TSS == 'Multi-TSS') ~ 'Multi-TSS\n& Term.',
      (class.term == 'Multi-Term.') ~ 'Multi-Term.',
      (class.TSS == 'Multi-TSS') ~ 'Multi-TSS',
      (class.term != 'Multi-Term.') & (class.TSS != 'Multi-TSS') &
        (class.gene != 'Single gene') ~ 'Traditional',
      (class.term != 'Multi-Term.') & (class.TSS != 'Multi-TSS') &
        (class.gene == 'Single gene') ~ 'Simple',
      TRUE ~ 'Other'
    )
  ) %>%
  # filter(is.na(class)) %>%
  # View
  count(`class`)  %>%
  arrange(desc(n)) %>%
  mutate(ratio = n / sum(n) * 100) -> op.type.stat

op.type.stat %>%
  mutate(class = fct_reorder(class, -n)) %>%
  ggplot(aes(x = class, y = n, fill = class)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama() +
  geom_text(aes(label = paste(round(ratio, 1), '%')),
            size = 6,
            position = position_dodge(width=0.9),
            vjust=-0.25) +
  xlab(NULL) + ylab('Number of Operons') +
  ylim(0, 1100) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'hide')

ggsave('analysis/08_operon_types.pdf',
       width = 20, height = 10, units = 'cm')
  


#####################################
# features distributions similar to existing work

isoforms$operons %>%
  select(id, tus = TUs) %>%
  separate_rows(tus, sep = ';') %>%
  left_join(isoforms$tus, c('tus' = 'id')) %>%
  separate_rows(genes, sep = ';') -> op.tu

op.tu %>%
  select(id, tus, start, end, strand) %>%
  unique %>%
  group_by(id) %>%
  mutate(op.start = min(start), op.end = max(end)) %>%
  ungroup %>%
  mutate(
    first = start == op.start,
    last = end == op.end
  ) %>%
  left_join(
    isoforms$transcripts %>%
      select(tus = TUs, TSS, Terminator),
    'tus'
  ) %>%
  mutate(
    internal.TSS = ifelse(
      ((strand == '+') & !first) | ((strand == '-') & !last),
      TSS,
      NA
    ),
    internal.term = ifelse(
      ((strand == '+') & !last) | ((strand == '-') & !first),
      Terminator,
      NA
    )
  ) %>%
  select(id, internal.TSS, internal.term) %>%
  gather('internal', 'bound', internal.TSS, internal.term) %>%
  drop_na %>%
  unique %>%
  count(id, internal) %>%
  count(internal, n) %>%
  rename(x = n, y = nn) %>%
  mutate(internal = ifelse(internal == 'internal.TSS',
                           'internal TSS',
                           'internal TTS')) %>%
  rename(desc = internal) %>%
  bind_rows(
    op.tu %>%
      select(id, genes) %>%
      unique %>%
      count(id) %>%
      count(n) %>%
      transmute(desc = 'Genes in operons', x = n, y = nn)
  ) -> dat

crossing(desc = unique(dat$desc), x = 1:32) %>%
  left_join(dat) %>%
  mutate_at('y', replace_na, 0) %>%
  mutate(
    lab = ifelse(y > 500, y, ''),
    y = pmin(500, y)
  ) %>%
  ggplot(aes(x = as.integer(x), y = y, fill = desc, group = desc)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label = lab),
            size = 6,
            position = position_dodge(width=0.9),
            vjust=-0.25) +
  scale_x_continuous(breaks = seq(1, 32, 2)) +
  ylim(0, 550) +
  ggsci::scale_fill_jama(name = NULL) +
  xlab('Number of features in operons') +
  ylab('Count of occurences') +
  theme_minimal(18)  +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))
  # scale_y_continuous(breaks = c(seq(0, 140, 20), 500, 1000, 1500))
    
ggsave('analysis/08_feature_dist.pdf',
       width = 16, height = 10, units = 'cm')

######################################
######################################
######################################
######################################

# One single cowplot for the publication


grand.stat %>%
  select(- Transcripts) %>%
  filter(src != 'Novel TUs') %>%
  gather('what', 'value', - src) %>%
  mutate_at('src', fct_relevel,
            "DBTBS", "BsubCyc", "SubtiWiki", "Combined", "BSGatlas") %>%
  mutate_at('what', fct_recode,
            'TUs' = "Transcriptional Units",
            '% genes with TU' = "% genes with transcripts"
            ) %>%
  mutate_at('what', fct_relevel,
            'Operons', 'TUs') %>%
  ggplot(aes(x = src, y = value, fill = fct_rev(src))) +
  ggsci::scale_fill_jama() +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = value %>%
                  round %>%
                  str_replace('(\\d)(\\d{3})', 
                              '\\1,\\2')),
            # size = 6,
            position = position_dodge(width=0.9),
            vjust=-0.25) +
  facet_wrap(~ what, scale = 'free') +
  theme_minimal(14) +
  scale_y_continuous(breaks = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'none') -> p1

tribble(
  ~ `Operon class`, ~ `Streptococcus pneumoniae`,
  ~ `Escherichia coli`, ~ `Bacillus subtilis`,
  'Simple', 47, 45, 43,
  'Traditional', 10, 19, 12,
  'Multi TSS', 1, 21, 21,
  'Multi TTS', 16, 7, 13,
  'Complex', 26, 8, 11
) %>%
  # map_if(is.numeric, sum)
  gather('org', 'pct', - `Operon class`) %>%
  mutate_at('Operon class', fct_relevel,
         'Simple', 'Traditional',
         'Multi TSS', 'Multi TTS', 'Complex') %>%
  ggplot(aes(x = `Operon class`, y = pct, fill = org, group = org)) +
  geom_bar(stat = 'identity',
           position = position_dodge2()) +
  geom_text(aes(label = paste0(pct, '%')),
            # size = 6,
            # position = position_dodge(width=0.9),
            position = position_dodge2(width=0.9),
            vjust=-0.25) +
  theme_minimal(14) +
  scale_y_continuous(breaks = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.text = element_text(face = "italic"),
        legend.justification = c(1, 1), legend.position = c(0.8, 0.9)) +
  ggsci::scale_fill_jco(name = 'Organism') -> p2

cowplot::plot_grid(p1, p2, ncol = 1,
                   labels = c('(a)', '(b)'))

ggsave('analysis/08_operon.pdf',
       width = 30, height = 20, units = 'cm')
