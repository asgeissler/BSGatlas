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
save(dbtbs_xml, file = 'data/03_dbtbs_xml.rda')

################################################################################
# Second: Fill in the numbers for easier use


### Find promoter binding sites

source('scripts/biostrings.R')
load('analysis/01_refseq.rda')
genome <- refseq$seq[[1]]

DBTBS$promoter %<>%
  mutate(clean.seq = str_remove_all(binding_sequence, '[^ACTG]'))

promoter.lookup <- DBTBS$promoter %>%
  select(clean.seq) %>%
  unique %>%
  filter(clean.seq != '') %>%
  pull(clean.seq) %>%
  find_pattern(., genome, 0)
  
# promoter.lookup$counts$sequence %>% unique %>% length
# 1262
# promoter.lookup$counts %>%
#   group_by(sequence) %>%
#   summarize(n = sum(count_matches)) %>%
#   mutate(mode = case_when(
#     n == 0 ~ 'none',
#     n == 1 ~ 'exact',
#     TRUE ~ 'ambiguous'
#   )) %>%
#   count(mode)
# mode      n
# <chr> <int>
#   1 exact  1236
#   2 none     26 (2.1%)


DBTBS$promoter %<>%
  left_join(promoter.lookup$positions, c('clean.seq' = 'sequence'))


### Find terminator positions
DBTBS$terminator %<>%
  mutate(clean.seq = str_remove_all(sequence, '[^ACTG]'))
term.lookup <- DBTBS$terminator %>%
  select(clean.seq) %>%
  unique %>%
  filter(clean.seq != '') %>%
  pull(clean.seq) %>%
  find_pattern(., genome, 0)

# term.lookup$counts$sequence %>% unique %>% length
# 1031
# term.lookup$counts %>%
#   group_by(sequence) %>%
#   summarize(n = sum(count_matches)) %>%
#   mutate(mode = case_when(
#     n == 0 ~ 'none',
#     n == 1 ~ 'exact',
#     TRUE ~ 'ambiguous'
#   )) %>%
#   count(mode)
# # mode          n
# # ambiguous     1
# # exact       926
# # none        104 (10.1%)

DBTBS$terminator %<>%
  left_join(term.lookup$positions, c('clean.seq' = 'sequence'))

### Separate transcription factors compute TSS
DBTBS$promoter %>%
  mutate(
    # binding posiitons relative to TSS
    left = str_split(sigma_binding_locations, '[; :]+') %>% 
      map(as.integer) %>%
      map(min) %>%
      unlist,
    right = str_split(sigma_binding_locations, '[; :]+') %>% 
      map(as.integer) %>%
      map(max) %>% 
      unlist,
    # cleanup mode, needed now
    mode = case_when(
      mode == 'Positive (SigA promoter); Negative (SigE promoter)' ~ 'Positive/Negative',
      mode == 'Pos/Neg' ~ 'Positive/Negative',
      startsWith(mode, 'Negative') ~ 'Negative',
      startsWith(mode, 'Positive') ~ 'Positive',
      startsWith(mode, 'Positive') ~ 'Positive',
      startsWith(mode, 'Promoter') ~ 'Promoter',
      TRUE ~ mode
    ),
    # compute TSS
    TSS = case_when(
      mode != 'Promoter' ~ NA_integer_,
      strand == '+'      ~ as.integer(end - right + 1),
      strand == '-'      ~ as.integer(start + right - 1),
      TRUE               ~ NA_integer_
    )
  ) -> intermediate


intermediate %>%
  select(TSS, strand, name = promoter_name, sigma = transcription_factor,
         reference) %>%
  mutate_if(is.character, function(i) ifelse(i == '', NA_character_, i)) %>%
  drop_na(TSS) %>%
  arrange(TSS)  -> TSSs

intermediate %>%
  filter(is.na(TSS), mode != 'Promoter') %>%
  arrange(start, desc(end)) %>%
  transmute(row = 1:n(),
            start, end, strand,
            factor = transcription_factor, mode, upstream_gene,
            reference) -> tfs

### Infer operons and determine IDs
load('analysis/02_merging.rda')
DBTBS$operon %>% 
  select(row, operon_name, gene = genes) %>%
  separate_rows(gene, sep = '[;,]') -> operon



operon %>%
  select(row, operon_name, name = gene) %>%
  left_join(
    merging$merged_src %>%
      select(merged_id, src, locus, name),
    'name'
  ) %>%
  # filter(str_detect(src, 'refseq')) %>%
  select(name, merged_id) %>%
  unique %>%
  drop_na -> direct.map

direct.map %>%
  count(name) %>%
  bind_rows(
    transmute(operon, name = gene, n = 0) %>%
      anti_join(direct.map, 'name') %>%
      unique
  ) %>%
  count(n) %>%
  rename('merged_ids matching DBTBS gene name' = n , 'count' = nn)
# `merged_ids matching DBTBS gene name` count
# 0             344
# 1            1848
# 2               9
# > 344 + 1848 + 9
# [1] 2201

load('analysis/03_subtiwiki.rda')
operon %>%
  select(name = gene) %>%
  unique %>%
  anti_join(direct.map, 'name') %>%
  left_join(
    # DBTBS$synonyms, c('name' = 'A')
    subtiwiki$genes %>%
      select(merged_id, name = names) %>%
      separate_rows(name, sep = ';'),
    'name'
  ) -> subti.synonym

subti.synonym %>%
  count(name) %>%
  bind_rows(
    transmute(operon, name = gene, n = 0) %>%
      anti_join(direct.map, 'name') %>%
      anti_join(subti.synonym, 'name') %>%
      unique
  ) %>%
  count(n) %>%
  rename('merged_ids matching via synonyms' = n , 'count' = nn)
# `merged_ids matching via synonyms` count
# 1   304
# 2    39
# 3     1


bind_rows(
  mutate(direct.map, mode = 'direct'),
  mutate(subti.synonym, mode = 'synonym')
) %>%
  unique -> gene.map

# gene.map %>%
#   left_join(count(gene.map, name)) %>%
#   left_join(merging$merged_src, 'merged_id') %>%
#   View
# not easy to resolve -> dump them away

gene.map %<>%
  select(name, merged_id) %>%
  drop_na %>%
  unique

gene.map %<>%
  anti_join(
    gene.map %>%
      count(name) %>%
      filter(n != 1),
    'name'
  )

# gene.map %>%
#   count(merged_id) %>%
#   count(n)
# n    nn
#  1  2165
#  2     6

operon %<>%
  left_join(gene.map, c('gene' = 'name'))
  # filter(!complete.cases(.)) %>%
  # View
# 24 troubled transcripts due to name look up

operon %<>%
  filter(is.na(merged_id)) %>%
  transmute(row, incomplete.flag = TRUE) %>%
  full_join(operon, 'row') %>%
  mutate_at('incomplete.flag', replace_na, FALSE) %>%
  drop_na %>%
  merge_variants('row') %>%
  select(operon = operon_name, merged_id, incomplete.flag) %>%
  mutate_at('incomplete.flag', as.logical)

# count(operon, incomplete.flag)
# # A tibble: 2 x 2
# incomplete.flag     n
# <lgl>           <int>
#   1 FALSE            1107
# 2 TRUE               16

dbtbs <- list(
  operon = operon,
  tss = TSSs,
  factor = tfs,
  term = DBTBS$terminator %>%
    select(operon_name, reference, energies, start, end, strand) %>%
    drop_na(strand) %>%
    arrange(start, desc(end))
)

save(dbtbs, file = 'data/03_dbtbs.rda')
