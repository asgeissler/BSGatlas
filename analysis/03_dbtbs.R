source('analysis/00_load.R')

source('scripts/frame_helpers.R')

################################################################################
# First: Parse the xml
library(xml2)

dbtbs.xml <- read_xml(
  file.path(rprojroot::find_rstudio_root_file(), 'data-raw', 'dbtbs.xml.gz')
)

# The top-most children types
top <- xml_children(dbtbs.xml) %>% xml_name() %>% unique
# the individual tags per top node
subtypes <- top %>%
  purrr::set_names(.) %>%
  map(function(i) {
    xml_find_all(dbtbs.xml, i) %>%
      xml_contents() %>%
      xml_name() %>%
      unique()
  })

# add references paths to promoter
subtypes$promoter %<>% append(
  c('reference/experiment', 'reference/pubmed')
) %>%
  # and remove the reference sub-node, the needed leaves were added above
  discard(. == 'reference')
# discard synonym as it is text only
subtypes %<>% discard(names(.) == 'synonym')

# similarly modify operon subtypes to reflect terminator and ref. sub-tree
subtypes$operon %<>%
  append(c('terminator/stemloop', 'terminator/energy')) %>%
  discard(. == 'terminator') %>%
  append(c('reference/pubmed')) %>%
  discard(. == 'reference')

subtypes$tfac %<>%
  append(c('reference/pubmed')) %>%
  discard(. == 'reference')

mod_paste <- function(i, sep = ';') {
  res <- replace_na(i, '') %>%
    discard(. == '') %>%
    unlist %>%
    vec_paste(sep = sep)
  if (identical(res, character(0))) {
    ''
  } else {
    res
  }
}

subtypes %>%
  map2(names(.), function(subs, name) {
    # name <- 'operon'; subs <- c('name', 'genes')
    # tibble(!! name := xml_find_all(dbtbs.xml, name) %>% head %>% as.list)
    nodes <- xml_find_all(dbtbs.xml, name)
    xs <- crossing(row = 1:length(nodes), label = subs)
    xs %>%
      rowwise %>%
      mutate(search = nodes[[row]] %>%
               xml_find_all(label) %>%
               map(xml_text) %>%
               mod_paste)
  }) -> tiblfied

tiblfied %<>% map(spread, label, search)

# bring synonyms in
tiblfied$synonym <- tibble(i = xml_find_all(dbtbs.xml, 'synonym') %>%
                             xml_text() %>%
                             as.list()
) %>%
  separate(i, ': ', into = c('a', 'b'))

# provide more reasonable names
DBTBS <- list(
  promoter =  set_names(
    tiblfied$promoter,
    c('row', 'comment', 'upstream_gene', 'sigma_binding_locations', 'promoter_name',
      'experimental_evidence', 'reference', 'mode', 'binding_sequence',
      'sigma_factor', 'transcription_factor')
  ),
  transcription_factor = set_names(
    tiblfied$tfac,
    c('row', 'comment', 'protein_domain', 'gene', 'reference', 'binding_motif',
      'factor_family')
  ),
  operon = set_names(
    tiblfied$operon,
    c('row', 'comment', 'experimental_evidence', 'genes', 'operon_name',
      'terminator_reference', 'terminator_energies', 'terminator_sequence')
  ),
  synonyms = set_names(tiblfied$synonym, c('A', 'B'))
)

# for improved readability spin out the terminator information
DBTBS$terminator <- DBTBS$operon %>%
  select(row, operon_name, starts_with('terminator_')) %>%
  set_names(c('row', 'operon_name', 'reference', 'energies', 'sequence')) %>%
  mutate(reference = ifelse(reference == '', NA_character_, reference))
  separate_rows(energies, sequence, sep = ';')

# and remove cols
DBTBS$operon %<>% select(-starts_with('terminator_'))

# store this as the 'unprocessed' version
dbtbs_xml <- DBTBS
save(dbtbs_xml, file = 'analysis/03_dbtbs_xml.rda')

################################################################################
# Second: Fill in the numbers for easier use