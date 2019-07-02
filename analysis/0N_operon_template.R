source('analysis/00_load.R')


#####################################
# 4. Infer operons/isoforms

library(tidygraph)

edges <- transcripts %>%
  select(from = id, to = genes) %>%
  separate_rows(to, sep = ';')
labels <- c(edges$from, edges$to) %>% unique
nodes <- tibble(id = 1:length(labels), label = labels)
# transform edges to ids
edges %<>%
  mutate(row = 1:n()) %>%
  gather('col', 'label', from, to) %>%
  left_join(nodes, 'label') %>%
  arrange(row) %>%
  select( - label) %>%
  spread(col, id) %>%
  select( - row)

grph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

grph %>%
  activate(nodes) %>%
  mutate(operon = group_components()) %>%
  as_tibble -> group

group %>%
  filter(str_detect(label, 'transcript')) %>%
  left_join(transcripts, c('label' = 'id')) %>%
  select(operon, start, end, strand, label) %>%
  group_by(operon) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand),
    isoforms = clean_paste(label)
  ) %>%
  ungroup %>%
  arrange(start, desc(end)) %>%
  mutate(id = sprintf('BSGatlas-operon-%s', 1:n())) %>%
  select(id, isoforms, start, end, strand) -> operons

# count(operons, strand)
# no awkwards strands, check

operons <- list(operon = operons, transcript = transcripts)
save(operons, file = 'analysis/03_operons.rda')

#####################################
# 5. Make statistics and nice plots


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

# Is the operon assumption given?
# Each operon has a transcript that contains all genes?
operons$operon %>%
  select(id, isoforms) %>%
  separate_rows(isoforms, sep = ';') %>%
  left_join(operons$transcript, c('isoforms' = 'id')) %>%
  separate_rows(genes, sep = ';') %>%
  select(id, isoforms, genes) %>%
  unique %>%
  count(id, isoforms) %>%
  group_by(id) %>%
  summarise(most.genes = max(n)) %>%
  left_join(operons.stat, 'id') %>%
  mutate(given = most.genes == `#genes`) %>%
  count(given) %>%
  mutate(pct = n / nrow(operons$operon) * 100)

# given     n    pct
# FALSE    16  0.843
# TRUE   1881 99.2