# assumption: working dir is rproject dir
source('analysis/00_load.R')

load('data/01_refseq.rda')

load('analysis/02_merging.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/07_isoforms.rda')


##############################################################################
# Make a unified color scheme
color.scheme <- read_tsv('data-hub/color_scheme.tsv')

# -1. Make 2bit file of genome
genome.id <- 'basu168'
genome.size <- length(refseq$seq[[1]])


##############################################################################
# 0. track for merged genes

merging$merged_genes %>%
  mutate(rgb.type = case_when(
    type %in% c('putative-coding', 'CDS') ~ 'protein',
    type %in% c('tRNA', 'rRNA', 'riboswitch') ~ type,
    type %in% c('sRNA', 'asRNA') ~ 'shortRNA',
    TRUE ~ 'other'
  )) %>%
  left_join(color.scheme, c('rgb.type' = 'type', 'strand')) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = merged_name, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb, nice.name = merged_id) -> merged.genes
# handle special case of split over origin
special <- merged.genes %>%
  filter(end > genome.size)
merged.genes %>%
  filter(end <= genome.size) %>%
  bind_rows(
    special %>%
      mutate(start2 = start2, end = genome.size,
             thickEnd = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    special %>%
      mutate(start2 = 0, end = end %% genome.size,
             thickStart = 0,
             thickEnd = end %% genome.size,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start2) -> merged.genes

merged.genes %>%
  write_tsv('igv-tracks/genes.bed', col_names = FALSE)

# the autosql part
'table my_genes
"my genes"
    (
    string chrom;          "Reference sequence genome"
    uint   chromStart;     "Start position in genome"
    uint   chromEnd;       "End position in genome"
    string name;           "BSGatlas ID"
    uint score;            "Score, not used"
    char[1] strand;        "+ or - for strand"
    uint thickStart;       "Start of where display should be thick"
    uint thickEnd;         "End of where display should be thick"
    uint reserved;         "User friendly color"
    string ID;             "Gene Name"
    )' %>%
  write_lines('igv-tracks/genes.as')
##############################################################################
# 3. operons
isoforms$operons %>%
  mutate(type = 'operon') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(chr = genome.id,
            start2 = utr.start - 1,
            end2 = utr.end,
            primary.name = id, score = 0,
            strand = strand,
            thickStart = tu.start - 1,
            thickEnd = tu.end,
            rgb) -> oper
oper.special <- oper %>%
  filter(end2 > genome.size)
oper %>%
  filter(end2 <= genome.size) %>%
  bind_rows(
    oper.special %>%
      mutate(start2 = start2,
             end2 = genome.size,
             thickStart = thickStart,
             thickEnd = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    oper.special %>%
      mutate(start2 = 0,
             end2 = end2 %% genome.size,
             thickStart = 0,
             thickEnd = thickEnd %% genome.size,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start2) -> oper

write_tsv(oper, 'igv-tracks/operons.bed', col_names = FALSE)

# the autosql part
'table boundary
"my boundary units"
    (
    string chrom;          "Reference sequence genome"
    uint   chromStart;     "Start position in genome"
    uint   chromEnd;       "End position in genome"
    string name;           "ID"
    uint score;            "Score, no used"
    char[1] strand;        "+ or - for strand"
    uint thickStart;       "Start of where display should be thick"
    uint thickEnd;         "End of where display should be thick"
    uint reserved;        "User friendly color"
    string extra;          "Extra information"
    )' %>%
  write_lines('igv-tracks/extra.as')


##############################################################################
# 5. the unified TSS/terminator map

# the autosql part
'table boundary
"my boundary units"
    (
    string chrom;          "Reference sequence genome"
    uint   chromStart;     "Start position in genome"
    uint   chromEnd;       "End position in genome"
    string name;           "BSGatlas ID"
    uint score;            "Score, no used"
    char[1] strand;        "+ or - for strand"
    uint thickStart;       "Start of where display should be thick"
    uint thickEnd;         "End of where display should be thick"
    uint reserved;       "User friendly color"
    string src;             "Origin of this annotation"
    string extra;          "Extra information"
    )' %>%
  write_lines('igv-tracks/src_extra.as')


bsg.boundaries$TSS %>%
  mutate(type = 'TSS') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(
    chrom = 'basu168',
    start1 = start - res.limit - 1, end1 = end + res.limit,
    name = sprintf('sigma %s', sigma), score = res.limit,
    strand,
    start2 = start - 1, end2 = end,
    rgb, id
  ) %>%
  arrange(start1, desc(end1)) %>%
  mutate(start1 = pmax(start1, 0)) %>%
  mutate_if(is.numeric, ~ pmin(genome.size, .x )) %>%
  write_tsv('igv-tracks/tss.bed', col_names = FALSE)

bsg.boundaries$terminator %>%
  mutate(type = 'terminator') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(
    chrom = 'basu168',
    start3 = start - 1, end,
    name = sprintf('Free energy: %.2f[kcal/mol]', energy),
    score = 0,
    strand,
    start2 = start3, end2 = end,
    rgb, id
  ) %>%
  write_tsv('igv-tracks/terminator.bed', col_names = FALSE)

##############################################################################
# 6. The new UTRs
UTRs %>%
  bind_rows() %>%
  mutate(type = 'UTR') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  arrange(start, end) %>%
  transmute(
    chrom = 'basu168',
    start3 = start - 1,
    end,
    name = id,
    score = 0,
    strand,
    start2 = start3,
    end2 = end,
    rgb, id
  ) %>%
  write_tsv('igv-tracks/utrs.bed', col_names = FALSE)

##############################################################################
# 7. All isoforms and their full lengths


isoforms$transcripts %>%
  left_join(
    isoforms$tus,
    c('TUs' = 'id'),
    suffix = c('.trans', '.tu')
  ) %>%
  arrange(start.trans, end.trans) %>%
  mutate(type = 'transcript') %>%
  left_join(color.scheme, c('type', 'strand.trans' = 'strand')) %>%
  left_join(
    bsg.boundaries$TSS %>%
      select(TSS = id, sigma),
    'TSS'
  ) %>%
  mutate(
    TSS = ifelse(is.na(TSS), 'without TSS',
                 sprintf('%s(sigma: %s)', TSS, sigma)),
    Terminator = ifelse(is.na(Terminator), 'without Terminator', Terminator),
    extra = paste(TSS, Terminator, sep = '<br/>')
  ) %>%
  transmute(
    chrom = 'basu168',
    start2 = start.trans - 1,
    end2 = end.trans,
    name = id,
    score = 0,
    strand.trans,
    start3 = start.trans - 1,
    end3 = end.trans,
    rgb, id
  ) -> trans


special.trans <- trans %>%
  filter(end2 > genome.size)
trans %>%
  filter(end2 <= genome.size) %>%
  bind_rows(
    special.trans %>%
      mutate(start2 = start2,
             end2 = genome.size,
             start3 = start3,
             end3 = genome.size,
             name = paste0(name, '.part1')),
    special.trans %>%
      mutate(start2 = 0,
             end2 = end2 %% genome.size,
             start3 = 0,
             end3 = end3 %% genome.size,
             name = paste0(name, '.part2'))
  ) %>%
  arrange(start2, end2) -> trans

trans %>%
  write_tsv('igv-tracks/transcripts.bed', col_names = FALSE)
