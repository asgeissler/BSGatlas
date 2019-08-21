# Convert annotation to a sqlite version for the front-end

source('analysis/00_load.R')

source('scripts/frame_helpers.R')


load('analysis/02_merging.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/07_isoforms.rda')
load('data-gff/01_meta.full.rda')

# Idea structure a file in the style of meta.full
# But with info for all elements

# get Ec nice names
ec.names <- read_tsv('data-gff/ec.tsv') %>%
  mutate_at('ID', str_remove, ' \\(transferred.*$')

#########################################################################

# Split "1 BSGatlas" out of meta.full and tuck togehter with positional
# information
gene.rest <- meta.full %>% filter(Resource != "1 BSGatlas")
gene.part <- meta.full %>% filter(Resource == "1 BSGatlas")

merging$merged_genes %>%
  transmute(merged_id,
            Coordinates = sprintf('%s..%s', start, end),
            Strand = strand,
            `Genomic Size` = sprintf('%s bp', end -start + 1)) %>%
  gather('meta', 'info', - merged_id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  # binding in this order should also preserve the order
  bind_rows(gene.part, gene.rest) %>%
  rename(id = merged_id) -> genes.step

# prettyfy 
genes.step %>% 
  filter(meta == 'Enzyme Classifications') %>%
  separate_rows(info, sep = ';') %>%
  unique %>%
  mutate_at('info', str_remove, '^EC-') %>%
  left_join(ec.names, c('info' = 'ID')) %>%
  # filter(is.na(Name)) %>%
  # View
  drop_na %>%
  mutate(info = sprintf('EC %s: %s', info, Name)) %>%
  select(-Name) %>%
  bind_rows(
    genes.step %>%
    filter(meta != 'Enzyme Classifications')
  ) %>%
  arrange(id, Resource, meta, info) -> genes

#########################################################################

bsg.boundaries$TSS %>%
  select(id, pid = pubmed) %>%
  separate_rows(pid, sep = ';') %>%
  filter(pid != '') %>%
  # break into smaller chunks
  mutate(chunk = 1:n() %/% 100) -> tss.pid

tss.pid %>%
  group_by(chunk) %>%
  do(response = rentrez::entrez_fetch(db = 'pubmed', id = .$pid,
                      rettype = "xml", retmode = "full", parsed = TRUE) %>%
       XML::xmlToList() ) -> queried

queried %>%
  pull(response) %>%
  unlist(recursive = FALSE) %>%
  set_names(NULL) %>%
  # map('PubmedArticle') %>%
  map('MedlineCitation') %>%
  # map('Article')  %>%
  map(safely(function(x) {
    tibble(
      pid     = x$PMID$text,
      lname   = x$Article$AuthorList[[1]]$LastName,
      iname   = x$Article$AuthorList[[1]]$Initials,
      title   = x$Article$ArticleTitle,
      journal = x$Article$Journal$ISOAbbreviation,
      year    = x$Article$Journal$JournalIssue$PubDate$Year
    )
  })) -> convert

# map(convert, 'error') %>% map(is.null) %>% unlist %>% summary
# map(convert, 'error') %>% discard(is.null)
# oh well, just 3

convert %>%
  map('result') %>%
  bind_rows() -> pid.look

tss.pid %>%
  select(-chunk) %>%
  left_join(pid.look) %>%
  mutate(nice = sprintf('%s, %s <i>et al.</i> %s. %s (%s)<br>PUBMED: %s',
                        lname, iname, title, journal, year, pid)) %>%
  select(id, nice) %>%
  mutate(
    Resource = '1 BSGatlas',
    meta = 'Citation'
  ) %>%
  unique %>%
  rename(info = nice) -> tss.cite

bsg.boundaries$TSS %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  rename(info = src) %>%
  mutate(
    Resource = '1 BSGatlas',
    meta = 'Based on'
  ) -> tss.src


isoforms$transcripts %>%
  select(trans = id, id = TSS) %>%
  drop_na %>%
  mutate(trans = sprintf('<a href="details.php?id=%s">%s</a>', trans, trans)) %>%
  group_by(id) %>%
  summarize_at('trans', clean_paste, sep = ', ') %>%
  ungroup %>%
  mutate(
    Resource = '1 BSGatlas',
    meta = 'Associated Transript(s)'
  ) %>%
  rename(info = trans) -> tss.trans

bsg.boundaries$TSS %>%
  transmute(id,
            Coordinates = paste0(TSS, '..', TSS),
            Strand = strand,
            `Sigma Factor` = sigma,
            `Resolution` = sprintf('&#177;%s bp', res.limit)) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(tss.src, tss.cite, tss.trans) %>%
  arrange(id, desc(meta)) %>%
  unique -> tss.meta

#########################################################################

bsg.boundaries$terminator %>%
  select(id, src) %>%
  separate_rows(src, sep = ';') %>%
  rename(info = src) %>%
  mutate(
    Resource = '1 BSGatlas',
    meta = 'Based on'
  ) -> tts.src

isoforms$transcripts %>%
  select(trans = id, id = Terminator) %>%
  drop_na %>%
  mutate(trans = sprintf('<a href="details.php?id=%s">%s</a>', trans, trans)) %>%
  group_by(id) %>%
  summarize_at('trans', clean_paste, sep = ', ') %>%
  ungroup %>%
  mutate(
    Resource = '1 BSGatlas',
    meta = 'Associated Transript(s)'
  ) %>%
  rename(info = trans) -> tts.trans
  
bsg.boundaries$terminator %>%
  transmute(
    id, 
    Coordinates = sprintf('%s..%s', start, end),
    Strand = strand,
    `Free Energy` = sprintf('%s [kcal/mol]', energy)
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(tts.src, tts.trans) %>%
  unique %>%
  arrange(id, desc(meta)) -> tts.meta


#########################################################################

isoforms$transcripts %>%
  select(trans = id, id = features) %>%
  separate_rows(id, sep = ';') %>%
  filter(str_detect(id, 'UTR')) %>%
  drop_na %>%
  mutate(trans = sprintf('<a href="details.php?id=%s">%s</a>', trans, trans)) %>%
  group_by(id) %>%
  summarize_at('trans', clean_paste, sep = ', ') %>%
  ungroup %>%
  mutate(
    Resource = '1 BSGatlas',
    meta = 'Associated Transript(s)'
  ) %>%
  rename(info = trans) -> utr.trans


UTRs %>%
  {.$internal_UTR %<>% select(- boundary) ; .} %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() %>%
  transmute(
    id,
    Coordinates = sprintf('%s..%s', start, end),
    Type = type,
    Strand = strand,
    `Nearest Gene` = gene,
    `Associated TSS/TTS` = boundary
  ) %>%
  gather('meta', 'info', - id) %>%
  drop_na  %>%
  mutate(info = ifelse(meta %in% c('Nearest Gene', 'Associated TSS/TTS'),
                       sprintf('<a href="details.php?id=%s">%s</a>', info, info),
                       info)) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(utr.trans) %>%
  unique %>%
  drop_na %>%
  filter(info != '') %>%
  arrange(id, desc(meta)) -> utr.meta


#########################################################################
isoforms$operons %>% 
  select(id, transcripts) %>%
  separate_rows(transcripts, sep = ';') %>%
  left_join(isoforms$transcripts, c('transcripts' = 'id')) %>%
  select(id, features) %>%
  separate_rows(features, sep = ';') %>%
  inner_join(merging$merged_genes, c('features' = 'merged_id')) %>%
  mutate(info = sprintf('%s (<a href="details.php?id=%s">%s</a>)', merged_name, features, features)) %>%
  group_by(id) %>%
  summarize_at('info', clean_paste, sep = ', ') %>%
  transmute(
    id, info,
    meta = 'Contained Genes',
    Resource = '1 BSGatlas'
  ) -> op.genes

  
isoforms$operons %>%
  transmute(
    id = id,
    Coordinates = sprintf('%s..%s', utr.start, utr.end),
    Strand = strand,
    `List of Isoforms` = transcripts %>%
      str_replace_all(';', ', ') %>%
      str_replace_all("([A-Za-z0-9-_']+)(,?)",
                        '<a href="details.php?id=\\1">\\1</a>\\2')
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(op.genes) %>%
  unique %>%
  arrange(id, desc(meta)) -> op.meta
#########################################################################
isoforms$transcripts %>%
  transmute(
    id = id,
    Coordinates = sprintf('%s..%s', start, end),
    Strand = strand,
    `Transcribed elements` = features %>%
      str_replace_all(';', ', ') %>%
      str_replace_all("([A-Za-z0-9-_']+)(,?)",
                        '<a href="details.php?id=\\1">\\1</a>\\2'),
    `Underlying TU` = TUs,
    `Promoter/TSS` = TSS,
    `Terminator/TTS` = Terminator
  ) %>%
  gather('meta', 'info', - id) %>%
  drop_na %>%
  mutate(info = ifelse(meta %in% c('Underlying TU', 'Terminator/TTS', 'Promoter/TSS'),
                       sprintf('<a href="details.php?id=%s">%s</a>', info, info),
                       info)) %>%
  mutate(Resource = '1 BSGatlas') %>%
  unique %>%
  arrange(id, meta) -> op.trans
#########################################################################
isoforms$tus %>%
  transmute(
    id = id,
    Coordinates = sprintf('%s..%s', start, end),
    Strand = strand,
    `Based on` = src %>%
      str_replace_all(';', ', '),
    `Contained Genes` = genes %>%
      str_replace_all(';', ', ') %>%
      str_replace_all("([A-Za-z0-9-_']+)(,?)",
                        '<a href="details.php?id=\\1">\\1</a>\\2')
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  unique %>%
  arrange(id, desc(meta)) -> op.tu
#########################################################################

all.meta <- bind_rows(genes, op.meta, op.trans, op.tu,
                      tss.meta, tts.meta, utr.meta)

# fix syntax on long strings
all.meta %<>%
  mutate(info = ifelse(str_detect(meta, 'Expression '),
                        str_replace_all(info, ';', ', '),
                        info)) 

save(all.meta, file = 'data-gff/03_meta.rda')

#########################################################################

# Generated cleaned up maps for these gene set types
interest <- c('Category', 'Gene Ontology', 'Functions', 'Enzyme Classifications')

# extra query nice name to make summary table work easier
all.meta %>%
  filter(Resource == '1 BSGatlas') %>%
  filter(str_detect(id, 'gene')) %>%
  filter(meta == 'Name') %>%
  select(gene = id, name = info) -> nice.name
  

# easy part
all.meta %>%
  filter(meta %in% interest) %>%
  select(group = meta, gene = id, set = info) %>%
  left_join(nice.name, 'gene')  %>%
  arrange(group, set, gene) -> p1


# Types but only for the one listed in BSGatlas
all.meta %>%
  filter(Resource == '1 BSGatlas', meta == 'Type') %>%
  mutate_at('info', fct_recode, 'coding-sequence' = 'CDS') %>%
  left_join(nice.name, c('id' = 'gene'))  %>%
  transmute(group = 'Type', gene = id, set = info, name) %>%
  arrange(group, set, gene) -> p2

bind_rows(p1, p2) %>%
  arrange(group, set, tolower(name)) %>%
  unique -> genesets


#########################################################################
# make explicit search index for more performance

searchable <- c(
  "Name", "Alternative Name",
  "Locus Tag", "Alternative Locus Tag"
)

all.meta %>%
  filter(Resource == '4 BsubCyc', meta == 'Description') %>%
  select(id, info) %>%
  unique %>%
  drop_na %>%
  filter(info != '') %>%
  rename(desc = info) -> txt

all.meta %>%
  filter(meta == 'Type', Resource == '1 BSGatlas') %>%
  select(id, type = info) %>%
  unique %>%
  bind_rows(
    bsg.boundaries %>%
      map2(names(.), ~ mutate(.x, type = .y)) %>%
      map(select, id, type),
    isoforms %>%
      set_names(c('Operon', 'Transcript', 'TU')) %>%
      map2(names(.), ~ mutate(.x, type = .y)) %>%
      map(select, id, type)
  ) -> types
  
all.meta %>%
  filter(meta %in% searchable) %>%
  select(- Resource) %>%
  unique %>%
  group_by(id, meta) %>%
  arrange(info) %>%
  summarize(info = info %>% invoke(.f = str_c, sep = ';')) -> look

look %>%
  spread(meta, info, fill = '') %>%
  left_join(txt) %>%
  left_join(types) %>%
  mutate(desc = ifelse(is.na(desc), type, desc)) -> search

#########################################################################
# Make an SQLite version

meta <- all.meta
con <- DBI::dbConnect(RSQLite::SQLite(), 'data-gff/03_meta.sqlite')
copy_to(con, meta, temporary = FALSE)
copy_to(con, genesets, temporary = FALSE)
copy_to(con, search, temporary = FALSE)

DBI::dbListTables(con)

# helper for copy paste
tbl(con, 'meta') %>%
  filter(id == 'X') %>%
  show_query()
# <SQL>
#   SELECT *
#   FROM `meta`
# WHERE (`id` = 'X')
tbl(con, 'meta') %>%
  filter(meta %in% searchable) %>%
  filter(LIKE(id, 'X') | LIKE(info, 'X')) %>%
  show_query()
# SELECT *
# FROM (SELECT *
#         FROM `meta`
#       WHERE (`meta` IN ('Name', 'Alternative Name', 'Locus Tag', 'Alternative Locus Tag')))
# WHERE (LIKE(`id`, 'X') OR LIKE(`info`, 'X'))
tbl(con, 'meta') %>%
  filter(meta %in% searchable) %>%
  # filter(LIKE(id, '%sig%') | LIKE(info, '%sig%')) %>%
  filter(id == 'BSGatlas-gene-1234') %>%
  select(id) %>%
  left_join(
    tbl(con, 'meta') %>%
      filter(Resource == '1 BSGatlas'),
    'id'
  ) %>%
  show_query()
# <SQL>
#   SELECT `LHS`.`id` AS `id`, `RHS`.`meta` AS `meta`, `RHS`.`info` AS `info`, `RHS`.`Resource` AS `Resource`
# FROM (SELECT `id`
#       FROM (SELECT *
#               FROM `meta`
#             WHERE (`meta` IN ('Name', 'Alternative Name', 'Locus Tag', 'Alternative Locus Tag')))
#       WHERE (`id` = 'BSGatlas-gene-1234')) AS `LHS`
# LEFT JOIN (SELECT *
#              FROM `meta`
#            WHERE (`Resource` = '1 BSGatlas')) AS `RHS`
# ON (`LHS`.`id` = `RHS`.`id`)


# improved search concept
tbl(con, 'search') %>%
  filter((id == 'BSGatlas-gene-1234')            |
         LIKE(`Alternative Locus Tag`, '%BSU00010%') |
         LIKE(`Locus Tag`, '%BSU00010%')         |
         LIKE(`Name`, '%Sig%')                   |
         LIKE(`desc`, '%gyrase%')) %>%
  show_query()
# <SQL>
#   SELECT *
#   FROM `search`
# WHERE ((`id` = 'BSGatlas-gene-1234') OR
#        LIKE(`Alternative Locus Tag`, '%BSU00010%') OR
#        LIKE(`Locus Tag`, '%BSU00010%') OR
#        LIKE(`Name`, '%Sig%') OR
#        LIKE(`desc`, '%gyrase%'))

DBI::dbDisconnect(con)


