

##############################

parse.off <- function(i) {
  g <- i %>%
    basename() %>%
    str_remove('.CRISPRoff.tsv$')
  i %>%
    read_tsv(comment = '#',
             col_names = c('chr', 'start', 'end', 'seq', 'CRISPRspec', 'strand')) %>%
    transmute(
      guide = g,
      start, end, seq, CRISPRspec, strand,
      mismatches = stringdist::stringdist(g, seq, 'hamming'),
      cmp.pam = viz.helper(g, seq) %>%
        str_replace('(?=...$)', ' ')
    )
  #   mutate(
  #     cut.pos = ifelse(strand == '+', end - 6 - 1, start + 6),
  #     ucsc = sprintf('basu168;%s%s;%.1f', cut.pos, strand, CRISPRspec),
  #     guide = g
  #   ) -> xs
  # 
  # xs %>%
  #   group_by(guide) %>%
  #   summarize(
  #     ucsc = str_c(ucsc, collapse = '|'),
  #     nr.exact.matches = sum(mis == 0),
  #     mm.str = count(., mis) %>%
  #       right_join(tibble(mis = 0:8)) %>%
  #       mutate_at('n', replace_na, 0) %>%
  #       arrange(mis) %>%
  #       pull(n) %>%
  #       str_c(collapse = ',')
  #   ) %>%
  #   ungroup
}


# file.path(path, '*.suff', '*CRISPRoff.tsv') %>%
file.path(path, 'AAAA_GCF_000009045.1_ASM904v1_genomic.suff', '*.CRISPRoff.tsv') %>%
  Sys.glob() %>%
  head %>%
  map(parse.off) %>%
  bind_rows -> offs
##############################
# Write a toy example

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
   uint reserved;       "Doench 2016 / Fusi et al. Score"
   string CRISPRspec;  "CRISPR specificity score"
   string Azimuth;      "Azimuth score"
   bigint _offset;        "Offset into tab-sep file for details page"
   )' %>%
  write_lines('grna.as')
   # bigint _offset;        "Offset into tab-sep file for details page"

guides %>%
  inner_join(offs, 'guide') -> xs
xs %>%
  transmute(
    chr = 'basu168',
    start2 = as.integer(start) - 1,
    end2 = as.integer(end),
    guide,
    s = round(CRISPRspec),
    strand,
    tS = ifelse(strand == '+', start2, start2 + 3),
    tE = ifelse(strand == '+', end2 - 3, end2),
    cl = '150,150,150',
    CRISPRspec, Azimuth,
    mm.str, ucsc
  ) %>%
  write_tsv('grna.full.bed', col_names = FALSE)

# ./cons.py
# module load hg
# bedSort grna.bed grna.bed
# bedToBigBed -tab \
#   -type=bed9+3 -as=grna.as \
#   grna.bed '/home/projects/nextprod/subprojects/BSGatlas_v1.1/data-hub/genome.info' grna.bb

