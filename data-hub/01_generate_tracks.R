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
  arrange(type, what) %>%
  with(set_names(rgb, sprintf('%s (%s)', type, what))) %>%
  map(function(i) {
    function () {scales::show_col(i, labels = TRUE)}
  }) %>%
  invoke(cowplot::plot_grid, ., labels = names(.),
         ncol = 2,
         scale = 0.7, label_size = 10, hjust = 0)

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


# 0. track for merged genes
genome.size <- length(refseq$seq[[1]])

merging$merged_genes %>%
  mutate(rgb.type = case_when(
    type %in% c('putative-coding', 'CDS') ~ 'protein',
    type %in% c('tRNA', 'rRNA', 'riboswitch') ~ type,
    type %in% c('sRNA', 'asRNA') ~ 'shortRNA',
    TRUE ~ 'other'
  )) %>%
  left_join(color.scheme, c('rgb.type' = 'type', 'strand')) %>%
  transmute(chr = 'ncbi168', start2 = start - 1, end,
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
             primary.name = paste0(primary.name, '.part1')),
    special %>%
      mutate(start2 = 1, end = end %% genome.size,
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
system()

# 0b. sensetive rfam genes that would be new

# 1. tracks for the individual gene resources
# note: meta page will be linked to main dscription of the merged set

# 2. transcripts

# 3. operons

# meta page for operons and transcripts
# meta pages for genes, also link individual genes to it
