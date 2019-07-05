# assumption: working dir is rproject dir
source('analysis/00_load.R')

load('data/01_refseq.rda')
load('data/01_bsubcyc.rda')
load('data/01_rfam.rda')
load('data/01_nicolas.rda')
load('data/01_dar-riboswitches.rda')

load('analysis/01_tomerge.rda')

load('analysis/02_merging.rda')
load('analysis/02_merging_stat.rda')

load('data/03_subtiwiki.rda')
load('data/03_dbtbs.rda')

load('analysis/05_bsg_boundaries.rda')

load('analysis/06_utrs.rda')

load('analysis/07_isoforms.rda')


##############################################################################
# Make a unified color scheme
colors <- tribble(
  ~type, ~normal, ~track,
  # dodger blue
  'protein', '#1E90FF',
  'gene',
  # bright red
  'TSS', '#FF0000',
  'boundaries',
  # grey
  'riboswitch', '#999999',
  'gene',
  'terminator', '#999999',
  'boundaries',
  # green
  'UTR', '#358000',
  'utrs',
  # yellow
  'tRNA', '#FFFF00',
  'gene',
  # bright green
  'shortRNA', '#00FF00',
  'gene',
  # purple
  'rRNA', '#CC33FF',
  'gene',
  # red
  'other', '#CC0000',
  'gene',
  # pink
  'operon', '#FF00CC',
  'operon',
  # orange
  'transcript', '#FF9900',
  'transcripts',
  # light pink
  'TF', '#FFCCCC',
  'boundaries'
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
  arrange(track, type, desc(what)) %>%
  with(set_names(rgb, sprintf('%s (%s track, %s)', type, track, what))) %>%
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
  write_tsv('data-hub/BSGatlas_genes.bed', col_names = FALSE)

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
  write_lines('data-hub/genes.as')


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
            nice.name = sprintf('%s(%s,rfam=%s)', name, type, family)) %>%
  write_tsv('data-hub/rfam-extra.bed', col_names = FALSE)


##############################################################################
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
           str_replace('nicolas lower', 'nicolas predictions') %>%
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
  mutate(nice.name = ifelse(str_detect(nice.name, 'row[- ]\\d+'),
                            'NA', nice.name)) %>%
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


##############################################################################
# 2. transcriptional units
isoforms$tus %>%
  mutate(type = 'transcript') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = id, score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb, src) -> tus
tus.special <- tus %>%
  filter(end > genome.size)
tus %>%
  filter(end <= genome.size) %>%
  bind_rows(
    tus.special %>%
      mutate(start2 = start2, end = genome.size,
             thickEnd = genome.size,
             primary.name = paste0(primary.name, '.part1')),
    tus.special %>%
      mutate(start2 = 0, end = end %% genome.size,
             thickStart = 0,
             thickEnd = end %% genome.size,
             primary.name = paste0(primary.name, '.part2'))
  ) %>%
  arrange(start2) -> tus

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
    uint reserved;         "User friendly color"
    string src;            "Origin of this transcriptional unit (TU) annotation"
    )' %>%
  write_lines('data-hub/transunit.as')

write_tsv(tus, 'data-hub/BSGatlas_tus.bed', col_names = FALSE)


##############################################################################
# 3. operons
isoforms$operons %>%
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

write_tsv(oper, 'data-hub/BSGatlas_operons.bed', col_names = FALSE)



##############################################################################
# 3. TSS and terminators

bind_rows(
  dbtbs$tss %>%
    transmute(start = TSS, end = TSS, strand, type = 'TSS',
              desc = ifelse(is.na(name), sigma, paste0(sigma, '_', name)) %>%
                paste0(ifelse(is.na(reference), '',
                              sprintf(';Pubmed=%s', reference))),
              src = 'DBTBS-TSS') %>%
    arrange(start) %>%
    mutate(id = paste0('DBTBS-TSS-', 1:n())),
  dbtbs$factor %>%
    transmute(start, end, strand, type = 'TF',
              desc = paste(factor, mode) %>%
                paste0(ifelse(is.na(reference), '',
                              sprintf(';Pubmed=%s', reference))),
              src = 'DBTBS-TF') %>%
    arrange(start) %>%
    mutate(id = paste0('DBTBS-TF-', 1:n())),
  dbtbs$term %>%
    transmute(start, end, strand, type = 'terminator',
              desc = paste0('energy=', energies) %>%
                paste0(ifelse(is.na(reference), '',
                              sprintf(';Pubmed=%s', reference))),
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
              desc = paste0('energy=', energy),
              src ='nicolas-downshift', type = 'terminator',
              id = paste0('nicolas-downshift-', str_remove(id, '^D'))),
  #
  bsubcyc$terminator %>%
    transmute(start, end, strand, type = 'terminator',
              desc = paste0('energy=', energy) %>%
                paste0(ifelse(is.na(rho.independent), '',
                              sprintf('rho.idependent=%s', 
                                      ifelse(rho.independent, 'yes', 'no')))),
              src = 'BsubCyc-term') %>%
    arrange(start) %>%
    mutate(id = paste0('BsubCyc-term-', 1:n())),
  bsubcyc$TSS %>%
    transmute(start = tss, end = tss, strand, type = 'TSS',
              desc = ifelse(is.na(name),
                            paste0('Sig', sigma),
                            paste0('Sig', sigma, '_', name)) %>%
                paste0(ifelse(is.na(cite), '',
                              sprintf(';Pubmed=%s', cite))),
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
  write_lines('data-hub/extra.as')


##############################################################################
# 4. The seemingly UTRs by nicolas et al
nicolas$all.features %>%
  drop_na(type) %>%
  filter(!str_detect(type, 'indep')) %>%
  mutate(
    specific.type = case_when(
      type == "5'" ~ "5'UTR",
      type == "inter" ~ 'internal_UTR',
      type == "inter" ~ 'intergenic',
      TRUE ~ "3'UTR"
    ),
    type = 'UTR'
  ) %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  arrange(start) %>%
  transmute(chr = genome.id, start2 = start - 1, end,
            primary.name = sprintf('%s(%s)', specific.type, name),
            score = 0,
            strand = strand,
            thickStart = start2, thickEnd = end,
            rgb) %>%
  write_tsv('data-hub/transbounds/nicolas-utrs.bed', col_names = FALSE)


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
  write_lines('data-hub/src_extra.as')


bsg.boundaries$TSS %>%
  mutate(type = 'TSS') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(
    chrom = 'basu168',
    start = TSS - res.limit - 1, end = TSS + res.limit,
    name = id, score = res.limit,
    strand,
    start2 = TSS - 1, end2 = TSS,
    rgb, src,
    extra = sprintf('Resolution limit=%s<br/>PubMed: %s', res.limit, pubmed)
  ) %>%
  arrange(start, desc(end)) %>%
  mutate(start = pmax(start, 0)) %>%
  write_tsv('data-hub/BSGatlas_tss.bed', col_names = FALSE)

bsg.boundaries$terminator %>%
  mutate(type = 'terminator') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  transmute(
    chrom = 'basu168',
    start3 = start - 1, end,
    name = id, score = 0,
    strand,
    start2 = start3, end2 = end,
    rgb, src,
    extra = sprintf('Free energy: %s[kcal/mol]', energy)
  ) %>%
  write_tsv('data-hub/BSGatlas_terminator.bed', col_names = FALSE)

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
    rgb
  ) %>%
  write_tsv('data-hub/BSGatlas_UTRs.bed', col_names = FALSE)

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
                 sprintf('%s(%s)', TSS, sigma)),
    Terminator = ifelse(is.na(Terminator), 'without Terminator', Terminator),
    extra = paste(TSS, Terminator, sep = '\\n')
  ) %>%
  transmute(
    chrom = 'basu168',
    start2 = start.trans - 1,
    end2 = end.trans,
    name = id,
    score = 0,
    strand.trans,
    start3 = start.tu - 1,
    end3 = end.tu,
    rgb,
    src,
    extra
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
  write_tsv('data-hub/BSGatlas_transcripts.bed', col_names = FALSE)
