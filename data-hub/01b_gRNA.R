library(tidyverse)

'table crispr
"CRISPR guides incl. off-targets"
(
  string chrom;       "Reference sequence chromosome or scaffold"
  uint   chromStart;  "Start position in chromosome"
  uint   chromEnd;    "End position in chromosome"
  string name;        "Name or ID of item, ideally both human readable and unique"
  uint score;         "Score (CRISPRspec rounded)"
  char[1] strand;     "+ or - for strand"
  uint thickStart;    "Start of where display should be thick (start codon)"
  uint thickEnd;      "End of where display should be thick (stop codon)"
  uint reserved;       "color"
  string CRISPRspec;  "CRISPR specificity score (CRISPRspec)"
  string Azimuth;      "CRISPR efficiency score (Azimuth)"
  bigint _offset;        "Offset into tab-sep file for details page"
)' %>%
  write_lines('data-hub/grna.as')
  

read_tsv('data-raw/CRISPRoff/02_offsets.tsv.gz') %>%
  select(guide, offset) %>%
  left_join(
    read_tsv('data-raw/CRISPRoff/00_guides.tsv.gz'),
    'guide'
  ) -> xs

xs %>%
  transmute(
    chr = 'basu168',
    start2 = as.integer(start) - 1,
    end2 = as.integer(end),
    guide,
    s = floor(CRISPRspec),
    strand,
    tS = ifelse(strand == '+', start2, start2 + 3),
    tE = ifelse(strand == '+', end2 - 3, end2),
    cl = '150,150,150',
    CRISPRspec, Azimuth,
    "_offset" = offset
  ) %>%
  write_tsv('data-hub/grna.bed', col_names = FALSE)


