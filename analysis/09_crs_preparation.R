# Preparing the annotation in a simplified manner for CRS analysis
# output for all coding sequences max. length of down-/up-stream UTRs
# without distinguishing 5'/3'/internal UTRs

source('analysis/00_load.R')

load('analysis/02_merging.rda')
load('analysis/06_utrs.rda')
load('analysis/07_isoforms.rda')


# list of coding genes
merging$merged_genes %>%
  filter(type %in% c('CDS', 'putative-coding')) -> coding

# for all triplet pairs from isoforms of up/down stream relationships.
# possibly skipping UTRs

isoforms$transcripts %>%
  select(id, strand, features) %>%
  separate_rows(features, sep = ';') %>%
  group_by(id) %>%
  mutate(
    one.next = lead(features),
    two.next = lead(features, n = 2)
  ) %>%
  ungroup() %>%
  select(id, strand, before = features, gene = one.next, after = two.next) %>%
  inner_join(coding, c('gene' = 'merged_id'),
             suffix = c('', '.gene')) %>%
  gather('key', 'utr', before, after) %>%
  drop_na %>%
  filter(str_detect(utr, 'UTR')) %>%
  left_join(bind_rows(UTRs),
            c('utr' = 'id'),
            suffix = c('', '.utr')) %>%
  # make sure the position indication before/after is correct for both
  # strands. Distance below should be either zero or one
  mutate_at(c('start.utr', 'end.utr'), as.integer) %>%
  mutate(
    direction.test = case_when(
      (key == 'before') & (strand == '+') ~ start - end.utr,
      (key == 'before') & (strand == '-') ~ start.utr - end,
      (key == 'after') &  (strand == '+') ~ start.utr - end,
      (key == 'after') &  (strand == '-') ~ start - end.utr
    )
  ) %>%
  # filter(direction.test > 1) %>%
  # View
  # there are due to oberlapping ncRNA some larger values, and one
  # special case that will be ignored here
  # but otherwise the after/before does consider strandness correctly
  filter(direction.test <= 2) %>%
  mutate(utr.len = end.utr - start.utr + 1) %>%
  select(gene, merged_name, type, start, end, strand = strand.gene,
         key, utr.len) %>%
  group_by_at(vars(- utr.len)) %>%
  summarize_all(max) %>%
  spread(key, utr.len, fill = 0) %>%
  select(gene, merged_name, type, start, end, strand,
         "UTR length on 5' end" = before,
         "UTR length on 3' end" = after) %>%
  mutate(
    'UTR length before start' = ifelse(
      strand == '+',
      `UTR length on 5' end`,
      `UTR length on 3' end`
    ),
    'UTR length after end' = ifelse(
      strand == '+',
      `UTR length on 3' end`,
      `UTR length on 5' end`
    )
  ) -> simplified


write_tsv(simplified, 'analysis/09_utr_lengths.tsv')

