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
            str_replace, '(\\d)(\\d{3})', '\\1,\\2')
            
  

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

# given     n    pct
# 1 FALSE     6  0.261
# 2 TRUE   2293 99.7


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
  mutate(
    class = case_when(
      # (class.term == 'Missing-Term.') | (class.TSS == 'Missing-TSS') ~ 'Misses Promoter and/or Term.',
      (class.term == 'Multi-Term.') & (class.TSS == 'Multi-TSS') ~ 'Multi-Promoter &\nTSS Operon',
      (class.term == 'Multi-Term.') ~ 'Multi-Term. Operon',
      (class.TSS == 'Multi-TSS') ~ 'Multi-Promoter Operon',
      (class.tu == 'Single TU') & (class.gene != 'Single gene') ~ 'Traditional Operon',
      (class.tu == 'Single TU') & (class.gene == 'Single gene') ~ 'Simple Operon',
      TRUE ~ 'Other'
    )
  ) %>%
  count(`class`)  %>%
  arrange(desc(n)) %>%
  mutate(ratio = n / sum(n) * 100) -> op.type.stat

op.type.stat %>%
  mutate(class = fct_reorder(class, -n)) %>%
  ggplot(aes(x = class, y = n, fill = class)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama() +
  geom_text(aes(label = paste(round(ratio, 1), '%')),
            position = position_dodge(width=0.9),
            vjust=-0.25) +
  xlab('') + ylab('Number of Operons') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'hide')

ggsave('analysis/08_operon_types.pdf',
       width = 28, height = 15, units = 'cm')
  


#####################################
# From last template run


load('analysis/03_operons.rda')

library(venn)
genes <- merging$merged_genes %>%
  rename(id = merged_id)

# > nrow(operons$operon)
# [1] 1897
# > nrow(operons$transcript)
# [1] 2474



pdf(file = 'analysis/03_venn_trans.pdf')
operons$transcrip %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  with(set_names(map(i, 1), src)) %>%
  venn(cexil = 1.3,
       cexsn = 1.5, zcolor = 'style')
dev.off()



# ratio of genes not described by any operon
operons$transcript %>%
  select(src, genes) %>%
  separate_rows(genes, sep = ';') %>%
  separate_rows(src, sep = ';') %>%
  mutate(has = TRUE) %>%
  select(id = genes, src, has) %>%
  unique %>%
  spread(src, has) -> has.mat

# STATISTICS 
# ratio of genes with a transcript
select(genes, id) %>%
  left_join(has.mat) %>%
  mutate_all(replace_na, FALSE) %>%
  mutate(BSGatlas = BsubCyc | DBTBS | SubtiWiki) %>%
  select(- id) %>%
  map(sum) %>%
  map(divide_by, nrow(genes)) %>%
  map(multiply_by, 100) %>%
  tibble(`% genes with transcript` = ., src = names(.)) %>%
  unnest %>%
  # Total number of transcripts
  left_join(
    operons$transcript %>%
      select(id, src) %>%
      separate_rows(src, sep = ';') %>%
      unique %>%
      count(src) %>%
      spread(src, n) %>%
      mutate(BSGatlas =  nrow(operons$transcrip) ) %>%
      gather('src', 'number of transcripts'),
    'src'
  ) %>%
  # Total number of operons
  left_join(
    operons$operon %>%
      select(op = id, id = isoforms) %>%
      separate_rows(id, sep = ';') %>%
      left_join(operons$transcrip, 'id') %>%
      select(op, src) %>%
      separate_rows(src, sep = ';') %>%
      unique %>%
      count(src) %>%
      spread(src, n) %>%
      mutate(BSGatlas = nrow(operons$operon)) %>%
      gather('src', 'number of operons')
  ) %>%    
  gather('what', 'value', -src) %>%
  mutate(src = fct_relevel(src, 'DBTBS', 'BsubCyc', 'SubtiWiki', 'BSGatlas')) %>%
  ggplot(aes(x = '', y = value, fill = src)) +
  ggsci::scale_fill_jama(name = 'resource') +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=value %>% round()),
            position=position_dodge(width=0.9), vjust=-0.25) +
  xlab(NULL) +
  ylab(NULL) +
  facet_wrap(~ what, scales = 'free') +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(size = 14))

ggsave(file = 'analysis/03_bars.pdf',
       width = 25, height = 9, units = 'cm')



operons$operon %>%
  select(id, isoforms) %>%
  separate_rows(isoforms, sep = ';') %>%
  left_join(operons$transcript, c('isoforms' = 'id')) %>%
  select(id, gene = genes) %>%
  separate_rows(gene, sep = ';') %>%
  left_join(genes, c('gene' = 'id')) %>%
  mutate(is.prot = type %in% c('putative-coding', 'CDS')) %>%
  select(id, gene, is.prot) %>%
  unique %>%
  group_by(id) %>%
  summarize('#genes' = n(),
            '#coding' = sum(is.prot),
            '%coding' = sum(is.prot) / n() * 100) %>%
  ungroup %>%
  left_join(
    operons$operon %>%
      transmute(id, isoforms, len = end - start + 1) %>%
      separate_rows(isoforms, sep = ';') %>%
      group_by(id) %>%
      summarize('#isoforms' = n(), len = first(len)) %>%
      ungroup,
    'id'
  ) -> operons.stat

save(operons.stat, file = 'analysis/03_operonstat.rda')

library(ggforce)

operons.stat %>%
  pull(`#isoforms`) %>%
  as_factor() %>%
  fct_expand(as.character(1:10)) %>%
  fct_relevel(as.character(1:10)) %>%
  fct_count %>%
  ggplot(aes(x = f, y = n)) +
  geom_bar(stat = 'identity') +
  xlab('isoforms per operon') + 
  ylab('count') +
  facet_zoom(ylim = c(0, 100))

ggsave(file = 'analysis/03_bar_Nisoforms.pdf',
       width = 13, height = 7, units = 'cm')

operons.stat %>%
  ggplot(aes(x = `#genes`)) +
  geom_bar() +
  # geom_freqpoly(bins = 20) +
  scale_x_continuous(breaks = seq(1, 35) %>% keep(. %% 5 == 0 ),
                     labels = as.character) +
  xlab('genes per operon') + 
  ylab('count') +
  facet_zoom(ylim = c(0, 100))

ggsave(file = 'analysis/03_bar_Ngenes.pdf',
       width = 13, height = 7, units = 'cm')

operons.stat %>%
  ggplot(aes(x = `#genes`, y = `%coding`)) +
  geom_point()

