# generate meta page for merged genes
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


###########
# 1. collect all meta information in unified form

# --------
# RefSeq

refseq[c('coding', 'noncoding')] %>%
  bind_rows %>%
  transmute(
    'Resource' = 'RefSeq',
    'Locus Tag' = locus,
    'Name' = name,
    'Alternative Locus Tag' = old_locus,
    'Title' = title,
    'Functions' = fnc,
    'Description' = description,
    'Type' = type,
    'Enzyme Classifications' = ec
  ) %>%
  mutate_all(function(i) ifelse(i == '', NA_character_, i)) %>%
  left_join(
    merging$merged_src %>%
      filter(str_detect(src, 'refseq')) %>%
      select(`Locus Tag` = locus, merged_id),
    'Locus Tag'
  ) %>%
  gather('meta', 'info', - c(Resource, merged_id)) %>%
  drop_na %>%
  unique  -> refseq.meta

# --------
# BsubCyc

# prettify citation
bsubcyc[c('coding', 'noncoding')] %>%
  bind_rows %>%
  transmute(
    'Resource' = 'BsubCyc',
    'Locus Tag' = locus,
    'Type' = type,
    'Name' = name,
    'Alternative Name', synonyms,
    'Description' = description,
    'Molecular weight' = mw,
    'Comment' = comment,
    'Enzyme Classifications' = ec
  ) %>%
  mutate_all(function(i) ifelse(i == '', NA_character_, i)) %>%
  left_join(
    merging$merged_src %>%
      filter(str_detect(src, 'bsubcyc')) %>%
      select(`Locus Tag` = locus, merged_id),
    'Locus Tag'
  ) -> bsub.overview
    

# lookup GO terms
AnnotationDbi::select(
  GO.db::GO.db,
  AnnotationDbi::keys(GO.db::GO.db, "GOID"),
  c("TERM")) %>%
  as_tibble -> go.look

bsubcyc[c('coding', 'noncoding')] %>%
  bind_rows %>%
  transmute(
    'Resource' = 'BsubCyc',
    'Locus Tag' = locus,
    meta = 'Gene Ontology',
    go = go) %>%
  separate_rows(go, sep = ';') %>%
  left_join(go.look, c('go' = 'GOID')) %>%
  transmute(Resource, `Locus Tag`, meta, info = paste(go, TERM)) -> bsub.go
  

# --------
# Rfam/Dar/nicolas predictions

merging$merged_src %>%
  transmute(
    merged_id,
    Resource = case_when(
      startsWith(src, 'rfam') ~ 'Rfam',
      src == 'nicolas lower' ~ 'Nicolas et al. predictions',
      src == 'dar riboswitches' ~ 'Dar et al. riboswitches',
      TRUE ~ NA_character_
    ),
    'Locus Tag' = locus,
    'Name' = name,
    'Type' = type
  ) %>%
  drop_na(Resource) %>%
  gather('meta', 'info', - c(merged_id, Resource)) %>%
  unique %>%
  # remove selfmade loci
  filter(!(str_detect(info, 'row') & (meta == 'Locus Tag'))) %>%
  drop_na -> meta.other
  


#######################################################



library(ReportingTools)

mg <- 'BSGatlas-gene-5'




