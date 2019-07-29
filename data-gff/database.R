# Convert annotation to a sqlite version for the front-end

source('analysis/00_load.R')

source('scripts/frame_helpers.R')


load('analysis/02_merging.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/07_isoforms.rda')
load('data-hub/03_meta.full.rda')

# Idea structure a file in the style of meta.full
# But with info for all elements

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
  rename(id = merged_id) -> genes


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

bsg.boundaries$TSS %>%
  transmute(id,
            Position = TSS,
            Strand = strand,
            `Sigma Factor` = sigma,
            `Resolution` = sprintf('&#177;%s bp', res.limit)) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(tss.src, tss.cite) %>%
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
  
bsg.boundaries$terminator %>%
  transmute(
    id, 
    Position = sprintf('%s..%s', start, end),
    Strand = strand,
    `Free Energy` = sprintf('%s [kcal/mol]', energy)
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(tts.src) %>%
  unique %>%
  arrange(id, desc(meta)) -> tts.meta


#########################################################################

UTRs %>%
  {.$internal_UTR %<>% select(- boundary) ; .} %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() %>%
  transmute(
    id,
    Position = sprintf('%s..%s', start, end),
    Type = type,
    Strand = strand,
    `Nearest Gene` = gene,
    `Associated TSS/TTS` = boundary
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  bind_rows(tts.src) %>%
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
  mutate(info = sprintf('%s (%s)', merged_name, features)) %>%
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
    Position = sprintf('%s..%s', utr.start, utr.end),
    Strand = strand,
    `List of Isoforms` = transcripts %>%
      str_replace_all(';', ', ')
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
    Position = sprintf('%s..%s', start, end),
    Strand = strand,
    `Transcribed elements` = features %>%
      str_replace_all(';', ', '),
    `Underlying TU` = TUs,
    `Promoter/TSS` = TSS,
    `Terminator/TTS` = Terminator
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  unique %>%
  mutate_at('info', replace_na, 'unknown') %>%
  arrange(id, meta) -> op.trans
#########################################################################
isoforms$tus %>%
  transmute(
    id = id,
    Position = sprintf('%s..%s', start, end),
    Strand = strand,
    `Based on` = src %>%
      str_replace_all(';', ', '),
    `Contained Genes` = genes %>%
      str_replace_all(';', ', ')
  ) %>%
  gather('meta', 'info', - id) %>%
  mutate(Resource = '1 BSGatlas') %>%
  unique %>%
  arrange(id, desc(meta)) -> op.tu
#########################################################################

all.meta <- bind_rows(genes, op.meta, op.trans, op.tu,
                      tss.meta, tts.meta, utr.meta)

# fix syntax on lon strings
all.meta %<>%
  mutate(info = ifelse(str_detect(meta, 'Expression '),
                        str_replace_all(info, ';', ', '),
                        info)) 

save(all.meta, file = 'data-gff/meta.rda')

#########################################################################

# Generated cleaned up maps for these gene set types
interest <- c('Category', 'Gene Ontology', 'Functions')
# ecs require some cleanup work
ec <- 'Enzyme Classifications'

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


# EC extra work
# get Ec nice names
ec.names <- read_tsv('data-gff/ec.tsv') %>%
  mutate_at('ID', str_remove, ' \\(transferred.*$')

all.meta %>%
  filter(meta == ec) %>%
  select(group = meta, gene = id, set = info) %>%
  separate_rows(set, sep = ';') %>%
  unique %>%
  mutate_at('set', str_remove, '^EC-') %>%
  left_join(ec.names, c('set' = 'ID')) %>%
  # filter(is.na(Name)) %>%
  # View
  drop_na %>%
  mutate(set = paste(set, Name)) %>%
  select(-Name) %>%
  arrange(group, set, gene) -> p2


bind_rows(p1, p2) %>%
  arrange(group, set, gene) -> genesets


#########################################################################
# Make an SQLite version

searchable <- c(
  "Name", "Alternative Name",
  "Locus Tag", "Alternative Locus Tag"
)

meta <- all.meta
con <- DBI::dbConnect(RSQLite::SQLite(), 'data-gff/meta.sqlite')
copy_to(con, meta, temporary = FALSE)
copy_to(con, genesets, temporary = FALSE)

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

DBI::dbDisconnect(con)


