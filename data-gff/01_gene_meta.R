# generate meta page for merged genes
# assumptions working dir is rproject dir
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

###########
# 1. collect all meta information in unified form

# --------
# RefSeq

refseq[c('coding', 'noncoding')] %>%
  bind_rows %>%
  transmute(
    'Resource' = '3 RefSeq',
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
      filter(str_detect(src, 'RefSeq')) %>%
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
    'Resource' = '4 BsubCyc',
    'Locus Tag' = locus,
    'Type' = type,
    'Name' = name,
    'Alternative Name' = synonyms,
    'Description' = description,
    'Molecular weight' = mw,
    'Comment' = comment,
    'Enzyme Classifications' = ec
  ) %>%
  mutate_all(function(i) ifelse(i == '', NA_character_, i)) %>%
  left_join(
    merging$merged_src %>%
      filter(str_detect(src, 'BsubCyc')) %>%
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
    'Locus Tag' = locus,
    go = go) %>%
  separate_rows(go, sep = ';') %>%
  drop_na %>%
  filter(go != '') %>%
  left_join(go.look, c('go' = 'GOID')) %>%
  transmute(`Locus Tag`, 'Gene Ontology' = paste(go, TERM)) -> bsub.go


bsubcyc[c('coding', 'noncoding')] %>%
  bind_rows %>%
  select(`Locus Tag` = locus, pid = cite) %>%
  separate_rows(pid, sep = ';') %>%
  left_join(bsubcyc$cite, 'pid') %>%
  filter(pid != '') %>%
  # filter(!complete.cases(.))
  transmute(
    `Locus Tag`,
    Citation = sprintf(
      
      '%s<br/>%s<br/>%s (%s)<br/>PUBMED: <a href="https://pubmed.ncbi.nlm.nih.gov/%s">%s</a>',
      authors, title, journal, year, pid, pid)) -> bsub.cite

bsub.overview %>%
  left_join(bsub.go, 'Locus Tag') %>%
  left_join(bsub.cite, 'Locus Tag') %>%
  mutate_all(function(i) ifelse(i == '', NA_character_, i)) %>%
  gather('meta', 'info', - c(Resource, merged_id)) %>%
  drop_na %>%
  unique  -> bsub.meta

# --------
# Rfam/Dar/nicolas predictions

merging$merged_src %>%
  filter(startsWith(src, 'Rfam')) %>%
  select(merged_id, id = src_locus) %>%
  left_join(
    bind_rows(
      rfam$conservative %>%
        mutate(id = paste0('Rfam conservative%', id)),
      rfam$medium %>%
        mutate(id = paste0('Rfam medium%', id) %>%
                 str_replace('-(\\d+)', ' \\1'))
    ),
    'id'
  ) %>%
  transmute(
    merged_id, Resource = '5 Rfam',
    Family = family,
    Name = name,
    Type = type
  ) %>%
  unique %>%
  gather('meta', 'info', - c(merged_id, Resource)) %>%
  unique %>%
  drop_na() -> meta.rfam


merging$merged_src %>%
  transmute(
    merged_id,
    Resource = case_when(
      src == 'Dar et al. riboswitches' ~ '6 Dar et al. riboswitches',
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
  filter(meta != 'Locus Tag') %>%
  drop_na -> meta.other
  

# --------
# SubtiWiki
subtiwiki$genes %>%
  transmute(
    merged_id, Resource = '2 SubtiWiki',
    'Locus Tag' = locus,
    'Name' = name,
    'Description' = description,
    'Molecular weight' = mw,
    'Isoelectric point' = pI,
    'Function' = `function`,
    'Product' = product,
    'Is essential?' = essential,
    'Enzyme Classifications' = ec,
    'Alternative Name' = names
  ) %>%
  mutate_all(function(i) ifelse(i == '', NA_character_, i)) %>%
  gather('meta', 'info', - c(merged_id, Resource)) %>%
  drop_na %>%
  unique -> subti.overview

subtiwiki$gene.categories %>%
  left_join(subtiwiki$categories, c('category' = 'id')) %>%
  transmute(Resource = '2 SubtiWiki', merged_id,
            meta = 'Category', info = paste(category, desc)) -> subti.path

# --------
# Correlations from Nicolas
merging$merged_src %>%
  inner_join(nicolas$all.features, 'locus') %>%
  transmute(merged_id,
            Resource = '7 Nicolas et al. predictions',
            Name = name.x,
            'Expression neg. correlated with' = negative_correlated_genes,
            'Expression pos. correlated with' = positive_correlated_genes
            ) %>%
  gather('meta', 'info', -c(merged_id, Resource)) -> meta.nic.cor

merging$merged_src %>%
  inner_join(nicolas$all.features, 'locus') %>%
  transmute(merged_id,
            Resource = '7 Nicolas et al. predictions',
            'Highly expressed condition' = highly_expressed_conditions,
            'Lowely expressed condition' = lowely_expressed_conditions
            ) %>%
  gather('meta', 'info', -c(merged_id, Resource)) %>%
  separate_rows(info, sep = '[ ,;]+') %>%
  mutate(info = case_when(
    info == 'M9-stat' ~ 'M9stat',
    info == 'M9-tran' ~ 'M9tran',
    info == 'M9-exp'  ~ 'M9exp',
    info == 'Lbstat' ~ 'LBstat',
    info == 'Lbtran' ~ 'LBtran',
    info == 'Lbexp'  ~ 'LBexp',
    info == 'M/G'    ~ 'M+G',
    TRUE ~ info
  )) %>%
  left_join(
    nicolas$conditions %>%
      fill(desc) %>%
      select(id, desc) %>%
      drop_na %>%
      mutate(id = id %>%
               str_split('_') %>%
               map(1) %>%
               unlist) %>%
      unique,
    c(info = 'id')
  ) %>%
  # no idea what 'G/S' or 'UNK1' are
  drop_na(desc) %>%
  mutate(
    info = sprintf('(%s) %s', info, desc)
  ) %>%
  select(-desc) -> meta.nic.cond

# --------
# And our conclusion
merging$merged_genes %>%
  transmute(merged_id, Resource = '1 BSGatlas',
            Type = type, Name = merged_name) %>%
  gather('meta', 'info', Type, Name) -> meta.bsg
#######################################################
# Load Kegg pathways

ko.names <- 'http://rest.kegg.jp/list/pathway/bsu' %>%
  read_tsv(col_names = FALSE) %>%
  set_names(c('ko', 'pathway_name'))

ko.genes <- 'http://rest.kegg.jp/link/bsu/pathway' %>%
  read_tsv(col_names = FALSE) %>%
  set_names(c('ko', 'gene'))

# match to BSGatlas
merging$merged_src %>%
  select(merged_id, locus) %>%
  left_join(
    ko.genes %>%
      mutate_at('gene', str_remove, '^bsu:'),
    c('locus' = 'gene')
  ) %>%
  drop_na() -> kegg.map

# Short investigation on mapping

# ko.genes %>% select(gene) %>% unique %>% nrow
# 1359
# kegg.map %>% select(locus) %>% unique %>% nrow
# 1359
# kegg.map %>% select(merged_id) %>% unique %>% nrow
# 1360
# -> all found with only one double enrty

# Prettify output and save

# Note: 'bsu02020' pathway is contained by the reference numbers
# 'map02020' and 'ko02020'. The latter seems what users are used to see.
# Thus we present them in that way.

kegg.map %>%
  left_join(ko.names, 'ko') %>%
  select(gene = merged_id, pathway = ko, pathway_name) %>%
  mutate_at('pathway', str_remove, '^path:bsu') %>%
  mutate_at('pathway', ~ paste0('ko', .x)) %>%
  mutate_at('pathway_name', str_remove,
            ' - Bacillus subtilis subsp. subtilis 168$') -> kegg.ids
write_tsv(kegg.ids, 'data-gff/kegg_mapping.tsv')

kegg.ids %>%
  transmute(
    merged_id = gene, Resource = '8 KEGG Pathways',
    meta = 'Pathway',
    info = sprintf(
      '<a href="https://www.genome.jp/kegg-bin/show_pathway?%s">%s (%s)</a>',
      pathway, pathway_name, pathway
    )
  ) -> kegg.meta

#######################################################

meta <- bind_rows(
  refseq.meta, 
  bsub.meta,
  meta.other,
  subti.overview,
  subti.path,
  kegg.meta,
  meta.bsg,
  meta.nic.cond,
  meta.nic.cor,
  meta.rfam
) %>%
  arrange(Resource, meta, info) %>%
  unique

# count(meta, meta)$meta
# [1] "Alternative Locus Tag"           "Alternative Name"                "Category"                       
# [4] "Citation"                        "Comment"                         "Description"                    
# [7] "Enzyme Classifications"          "Expression neg. correlated with" "Expression pos. correlated with"
# [10] "Family"                          "Function"                        "Functions"                      
# [13] "Gene Ontology"                   "Highly expressed condition"      "Is essential?"                  
# [16] "Isoelectric point"               "Locus Tag"                       "Lowely expressed condition"     
# [19] "Molecular weight"                "Name"                            "Pathway"                        
# [22] "Product"                         "Title"                           "Type"

to.split <- c('Alternative Name', 'Functions')

meta2 <- bind_rows(
  meta %>%
    filter(meta %in% to.split) %>%
    separate_rows(info, sep = ';'),
  meta %>%
    filter(! (meta %in% to.split))
)


#######################################################
# Add outlinks

  
rfam %>%
  map2(names(.), ~ mutate(.x, src_locus = sprintf(
    'rfam %s%%%s', 
    .y, id
  ))) %>%
  bind_rows() %>%
  select(family, src_locus) -> fams


tribble(
  ~db, ~db.name, ~url,
  'bsubcyc', 'BsubCyc',
  'https://bsubcyc.org/gene?orgid=BSUB&id=',
  'subti', 'SubtiWiki',
  'http://subtiwiki.uni-goettingen.de/v3/gene/search/exact/',
  'rfam', 'Rfam',
  'http://rfam.xfam.org/search?q='
) -> dbs

merging$merged_src %>%
  mutate(
    db = case_when(
      str_detect(src, 'BsubCyc') ~ 'bsubcyc',
      str_detect(src, 'Rfam') ~ 'rfam',
      str_detect(src, 'RefSeq') ~ 'subti',
      str_detect(src, 'Nicolas') ~ 'subti'
    )
  ) %>%
  select(merged_id, db, locus, src_locus) %>%
  drop_na %>%
  unique  %>%
  left_join(fams, 'src_locus') %>%
  mutate(key = ifelse(is.na(family), locus, family)) %>%
  mutate_at('key', str_remove, '_match_.*$') %>%
  select(merged_id, db, key) %>%
  mutate(
    key = ifelse(db == 'subti',
                 str_replace(key, 'BSU_', 'BSU'),
                 key)
  ) %>%
  unique %>%
  left_join(dbs, 'db') %>%
  mutate(
    txt = sprintf('<a href="%s%s" target="_blank">%s</a>', url, key, db.name),
    Resource = '1 BSGatlas',
    meta = 'Outside Links'
  ) %>%
  select(Resource, merged_id, meta, info = txt) -> out.links
  
bind_rows(meta2, out.links) %>%
  mutate_at('meta', fct_relevel,
            "Name", "Alternative Name", "Locus Tag", "Family",
            "Alternative Locus Tag", "Functions",
             "Title",  "Description",  "Type",
             "Product", "Category",  "Enzyme Classifications",
             "Outside Links", "Function", "Is essential?",
             "Isoelectric point",  "Molecular weight",
             "Citation", "Gene Ontology", "Comment",
             "Expression neg. correlated with",
            "Expression pos. correlated with", 
            "Highly expressed condition", "Lowely expressed condition"
            ) %>%
  arrange(merged_id, Resource, meta, info) %>%
  mutate_at('meta', as.character) -> meta.full

# meta.full %>%
#   filter(merged_id == 'BSGatlas-gene-3138') %>%
#   arrange(Resource, meta) %>%
#   View


searchable <- c(
  "Name", "Alternative Name",
  "Locus Tag", "Alternative Locus Tag"
)
meta.full %>%
  filter(meta %in% searchable) %>%
  select(merged_id, info) %>%
  mutate_at('info', stringi::stri_enc_toascii) %>%
  filter(!str_detect(info, '^[ ]+$')) %>%
  mutate_at('info', na_if, y = '') %>%
  drop_na %>%
  unique -> search

# stringi::stri_enc_isascii(search$info) %>% table

search %>%
  write_delim(path = 'data-hub/search.txt', delim = ' ', col_names = FALSE)

# and again manual execution
# ixIxx search.txt search.ix search.ixx

save(meta.full, file = 'data-gff/01_meta.full.rda')

