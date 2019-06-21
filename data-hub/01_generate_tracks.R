# assumptions working dir is rproject dir
source('analysis/00_load.R')

load('analysis/01_refseq.rda')
load('analysis/01_bsubcyc.rda')
load('analysis/01_rfam.rda')
load('analysis/01_nicolas.rda')
load('analysis/01_dar-riboswitches.rda')

load('analysis/01_tomerge.rda')


load('analysis/02_merging.rda')
load('analysis/02_mergign_stat.rda')

load('analysis/03_subtiwiki.rda')
load('analysis/03_dbtbs.rda')
load('analysis/03_operons.rda')

load('analysis/05_bsg_boundaries.rda')

# Make a unified color scheme
colors <- tribble(
  ~type, ~normal,
  # blue
  'protein', '#0000FF',
  # bright red
  'TSS', '#FF0000',
  # grey
  'terminator', '#999999',
  # green
  'UTR', '#358000',
  # yellow
  'tRNA', '#FFFF00',
  # pink
  'operon', '#FF00CC',
  # orange
  'transcript', '#FF9900',
  # bright green
  'shortRNA', '#00FF00',
  # purple
  'rRNA', '#CC33FF',
  # light pink
  'TF', '#FFCCCC',
  # cyan
  'riboswitch', '#00FFFF',
  # red
  'other', '#CC0000'
)

colors %>%
  rowwise %>%
  do(foo = colorspace::hex2RGB(.$normal)@coords %>% as_tibble) %>%
  unnest %>%
  multiply_by(255) %>%
  bind_cols(colors) -> colors

colors %>%
  with(rgb2hsv(R, G, B)) %>%
  t %>%
  as_tibble %>%
  bind_cols(colors) %>%
  rowwise %>%
  mutate(darker = hsv(h, 0.8 * s, 0.7 * v )) -> colors


colors %>%
  gather('what', 'rgb', normal, darker) %>%
  mutate(what = fct_recode(what,
                           'forward' = 'normal',
                           'reverse' = 'darker')) %>%
  arrange(type, desc(what)) %>%
  with(set_names(rgb, sprintf('%s (%s)', type, what))) %>%
  map(function(i) {
    function () {scales::show_col(i, labels = TRUE)}
  }) %>%
  invoke(cowplot::plot_grid, ., labels = names(.),
         ncol = 2,
         scale = 0.7, label_size = 10, hjust = 0)

ggsave(file = 'data-hub/color_scheme.pdf',
       width=42, height = 14, units = 'cm')

colors %>%
  select(type, normal, darker) %>%
  gather('strand', 'rgb', normal, darker) %>%
  # convert to ',' separated string for bed export
  rowwise %>%
  mutate(
    html = rgb,
    rgb = colorspace::hex2RGB(html)@coords %>%
           as_tibble %>%
           multiply_by(255) %>%
           with(paste(R, G, B, sep = ',')),
    strand = ifelse(strand == 'normal', '+', '-')
  ) -> color.scheme

write_tsv(color.scheme, path = 'data-hub/color_scheme.tsv')

# -1. Make 2bit file of genome
genome.id <- 'basu168'
genome.size <- length(refseq$seq[[1]])
# refseq$seq %>%
#   set_names(genome.id) %>%
#   Biostrings::writeXStringSet('data-hub/genome.fna')

# excectued manually, because rstudio does not seem to load binary files of
# the conda enviornment?
# > faToTwoBit genome.fna genome.2bit
# > twoBitInfo genome.2bit genome.info
# > rm genome.fna


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
            primary.name = merged_id, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb, nice.name = merged_name) -> merged.genes
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
  write_tsv('data-hub/merged_genes.bed', col_names = FALSE)

# the autosql part
'table my_genes
"my genes"
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
    string ID;             "Gene Name"
    )' %>%
  write_lines('data-hub/genes.as')

# convert to bigbed
# bedToBigBed -type=bed9+1 -tab -as=genes.as \
#   merged_genes.bed                         \
#   genome.info                              \
#   -extraIndex=name,ID                      \
#   merged_genes.bb


# 0b. sensetive rfam genes that would be new
rfam$sensetive %>%
  anti_join(rfam$medium, c('family', 'start', 'end', 'strand')) %>%
  mutate(rgb.type = case_when(
    type %in% c('putative-coding', 'CDS') ~ 'protein',
    type %in% c('tRNA', 'rRNA', 'riboswitch') ~ type,
    type %in% c('sRNA', 'asRNA') ~ 'shortRNA',
    TRUE ~ 'other'
  )) %>%
  left_join(color.scheme, c('rgb.type' = 'type', 'strand')) %>%
  arrange(start) %>%
  transmute(chr = genome.id,
            start2 = start - 1, end,
            primary.name = sprintf('Rfam-sensetive-extra-%s', 1:n()),
            score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb,
            nice.name = paste(family, name, type, sep = '.')) %>%
  write_tsv('data-hub/rfam-extra.bed', col_names = FALSE)

# bedToBigBed -type=bed9+1 -tab -as=genes.as \
#   rfam-extra.bed                         \
#   genome.info                              \
#   -extraIndex=name,ID                      \
#   rfam-extra.bb

# 1. tracks for the individual gene resources
# note: meta page will be linked to main description of the merged set

merging$merged_src %>%
  mutate(rgb.type = case_when(
    type %in% c('putative-coding', 'CDS') ~ 'protein',
    type %in% c('tRNA', 'rRNA', 'riboswitch') ~ type,
    type %in% c('sRNA', 'asRNA') ~ 'shortRNA',
    TRUE ~ 'other'
  )) %>%
  left_join(color.scheme, c('rgb.type' = 'type', 'strand')) %>%
  arrange(start) %>%
  mutate(src = src %>%
           str_replace('^rfam', 'Rfam') %>%
           str_replace('^refseq', 'RefSeq') %>%
           str_replace('^bsubcyc', 'BsubCyc') %>%
           str_replace('nicolas lower', 'nicolas predicitons') %>%
           str_replace('^nicolas', 'Nicolas et al.') %>%
           str_replace('^dar', 'Dar et al.')) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = merged_id, score = priority,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb, nice.name = locus,
            src) -> raw

raw.special <- raw %>%
  filter(end > genome.size)
raw %>%
  filter(end <= genome.size) %>%
  bind_rows(
    raw.special %>%
      mutate(start2 = start2, end = genome.size,
             thickEnd = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    raw.special %>%
      mutate(start2 = 0, end = end %% genome.size,
             thickStart = 0,
             thickEnd = end %% genome.size,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start2) -> raw

raw %>%
  mutate_if(is.character, str_replace_all, ' ', '-') %>%
  group_by(src) %>%
  do(foo = set_names(list(.), first(.$src))) %>%
  pull(foo) %>%
  invoke(.f = c) %>%
  map(select, - src) %>%
  map2(names(.), function(tbl, i) {
    sprintf('data-hub/individual/%s.bed', i) %>%
      str_replace_all(' ', '_') -> to
    write_tsv(tbl, to, col_names = FALSE)
  })

# for i in individual/*.bed ; do
# suff=${i%%[.]bed*}
# echo $suff
# bedToBigBed -type=bed9+1 -tab -as=genes.as \
#   $i                                       \
#   genome.info                              \
#   -extraIndex=name,ID                      \
#   $suff.bb
# done

# 2. transcriptional units
operons$transcript %>%
  mutate(type = 'transcript') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = id, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb, src,
            inc = ifelse(possibly.incomplete, 'yes', 'no')) -> trans
trans.special <- trans %>%
  filter(end > genome.size)
trans %>%
  filter(end <= genome.size) %>%
  bind_rows(
    trans.special %>%
      mutate(start2 = start2, end = genome.size,
             thickEnd = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    trans.special %>%
      mutate(start2 = 0, end = end %% genome.size,
             thickStart = 0,
             thickEnd = end %% genome.size,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start2) -> trans

# the autosql part
'table transunits
"my trans units"
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
    string src;             "Origin of this transcript annotation"
    string incomplete;             "Transcript might be incomplete?"
    )' %>%
  write_lines('data-hub/transunit.as')

write_tsv(trans, 'data-hub/transunit.bed', col_names = FALSE)

# bedToBigBed -type=bed9+2 -tab -as=transunit.as \
#   transunit.bed                            \
#   genome.info                              \
#   -extraIndex=name                         \
#   transunit.bb

# 3. operons
operons$operon %>%
  mutate(type = 'operon') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = id, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb) -> oper
oper.special <- oper %>%
  filter(end > genome.size)
oper %>%
  filter(end <= genome.size) %>%
  bind_rows(
    oper.special %>%
      mutate(start2 = start2, end = genome.size,
             thickEnd = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    oper.special %>%
      mutate(start2 = 0, end = end %% genome.size,
             thickStart = 0,
             thickEnd = end %% genome.size,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start2) -> oper

write_tsv(oper, 'data-hub/operon.bed', col_names = FALSE)

# bedToBigBed -type=bed9   -tab    \
#   operon.bed                     \
#   genome.info                    \
#   -extraIndex=name               \
#   operon.bb


# 3. TSS and terminators

bind_rows(
  dbtbs$tss %>%
    transmute(start = TSS, end = TSS, strand, type = 'TSS',
              desc = ifelse(is.na(name), sigma, paste0(sigma, '_', name)),
              src = 'DBTBS-TSS') %>%
    arrange(start) %>%
    mutate(id = paste0('DBTBS-TSS-', 1:n())),
  dbtbs$factor %>%
    transmute(start, end, strand, type = 'TF',
              desc = paste(factor, mode),
              src = 'DBTBS-TF') %>%
    arrange(start) %>%
    mutate(id = paste0('DBTBS-TF-', 1:n())),
  dbtbs$term %>%
    transmute(start, end, strand, type = 'terminator',
              desc = paste0('energy=', energies),
              src = 'DBTBS-terminators') %>%
    arrange(start) %>%
    mutate(id = paste0('DBTBS-term-', 1:n())),
  #
  nicolas$upshifts %>%
    transmute(start = pos - 22, end = pos + 22, strand,
              desc = sigma, type = 'TSS',
              src ='nicolas-upshift',
              id = paste0('nicolas-TSS-', str_remove(id, '^U'))),
  nicolas$downshifts %>%
    transmute(start = pos - 22, end = pos + 22, strand,
              desc = paste0('ernergy=', energy),
              src ='nicolas-downshift', type = 'terminator',
              id = paste0('nicolas-downshift-', str_remove(id, '^D'))),
  #
  bsubcyc$terminator %>%
    transmute(start, end, strand, type = 'terminator',
              desc = paste0('energy=', energy),
              src = 'BsubCyc-term') %>%
    arrange(start) %>%
    mutate(id = paste0('BsubCyc-term-', 1:n())),
  bsubcyc$TSS %>%
    transmute(start = tss, end = tss, strand, type = 'TSS',
              desc = ifelse(is.na(name),
                            paste0('Sig', sigma),
                            paste0('Sig', sigma, '_', name)),
              src = 'BsubCyc-TSS') %>%
    arrange(start) %>%
    mutate(id = paste0('BsubCyc-TSS', 1:n()))
) %>%
  drop_na(start, end, strand, id) -> dat

dat %>%
  arrange(start) %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = id, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb, desc, src) %>%
  group_by(src) %>%
  do(foo = set_names(list(.), first(.$src))) %>%
  pull(foo) %>%
  invoke(.f = c) %>%
  map(select, - src) %>%
  map2(names(.), function(tbl, i) {
    sprintf('data-hub/transbounds/%s.bed', i) %>%
      str_replace_all(' ', '_') -> to
    write_tsv(tbl, to, col_names = FALSE)
  })

# for i in transbounds/*.bed ; do
# suff=${i%%[.]bed*}
# echo $suff
# bedToBigBed -type=bed9+1 -tab -as=genes.as \
#   $i                                       \
#   genome.info                              \
#   -extraIndex=name,ID                      \
#   $suff.bb
# done


# 4. The seemingly faulty UTRs
nicolas$all.features %>%
  drop_na(type) %>%
  filter(!str_detect(type, 'indep')) %>%
  mutate(type = 'UTR') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  arrange(start) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = locus, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb) %>%
  write_tsv('data-hub/transbounds/nicolas-utrs.bed', col_names = FALSE)

# bedToBigBed -type=bed9         \
#   transbounds/nicolas-utrs.bed \
#   genome.info                  \
#   transbounds/nicolas-utrs.bb

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
  write_lines('data-hub/src_extra.as')

bsg.boundaries$TSS %>%
  mutate(type = 'TSS') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(
    chrom = 'basu168',
    start = TSS - res.limit, end = TSS + res.limit,
    name = id, score = res.limit,
    strand,
    start2 = TSS, end2 = TSS,
    rgb, src,
    extra = sprintf('Resolution limit=%s<br/>PubMed: %s', res.limit, pubmed)
  ) %>%
  arrange(start, desc(end)) %>%
  write_tsv('data-hub/bsgatlas-tss.bed', col_names = FALSE)

bsg.boundaries$terminator %>%
  mutate(type = 'terminator') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(
    chrom = 'basu168',
    start, end,
    name = id, score = 0,
    strand,
    start2 = start, end2 = end,
    rgb, src,
    extra = sprintf('Free energy: %s[kcal/mol]', energy)
  ) %>%
  write_tsv('data-hub/bsgatlas-terminator.bed', col_names = FALSE)
 
# for i in bsgatlas-{tss,terminator} ; do
# bedToBigBed -type=bed9+2 -tab -as=src_extra.as \
#   $i.bed                                       \
#   genome.info                                  \
#   -extraIndex=name                             \
#   $i.bb
# done