# Convert annotation to a sqlite version for the front-end

load('analysis/00_load.R')

library(dplyr)


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
  rename(info = nice) -> tss.cite


#########################################################################


bsg.boundaries$terminator


#########################################################################

UTRs %>%
  {.$internal_UTR %<>% select(- boundary) ; .} %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() -> utr
