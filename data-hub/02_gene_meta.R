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
    'Locus Tag' = locus,
    go = go) %>%
  separate_rows(go, sep = ';') %>%
  left_join(go.look, c('go' = 'GOID')) %>%
  transmute(`Locus Tag`, 'Gene Ontology' = paste(go, TERM)) -> bsub.go


bsubcyc[c('coding', 'noncoding')] %>%
  bind_rows %>%
  select(`Locus Tag` = locus, pid = cite) %>%
  separate_rows(pid, sep = ';') %>%
  left_join(bsubcyc$cite, 'pid') %>%
  transmute(
    `Locus Tag`,
    Citation = sprintf(
      '%s<br/>%s<br/>%s (%s)<br/>PubMed: %s',
      authors, title, journal, year, pid)) -> bsub.cite

bsub.overview %>%
  left_join(bsub.go, 'Locus Tag') %>%
  left_join(bsub.cite, 'Locus Tag') %>%
  mutate_all(function(i) ifelse(i == '', NA_character_, i)) %>%
  gather('meta', 'info', - c(Resource, merged_id)) %>%
  separate_rows(info, sep = ';') %>%
  drop_na %>%
  unique  -> bsub.meta

# --------
# Rfam/Dar/nicolas predictions

merging$merged_src %>%
  transmute(
    merged_id,
    Resource = case_when(
      startsWith(src, 'rfam') ~ '5 Rfam',
      src == 'nicolas lower' ~ '7 Nicolas et al. predictions',
      src == 'dar riboswitches' ~ '6 Dar et al. riboswitches',
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
  separate_rows(info, sep = ';') %>%
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
            'Highly expressed condition' = highly_expressed_conditions,
            'Lowely expressed condition' = lowely_expressed_conditions,
            'Expression neg. correlated with' = negative_correlated_genes,
            'Expression pos. correlated with' = positive_correlated_genes
            ) %>%
  gather('meta', 'info', -c(merged_id, Resource)) %>%
  separate_rows(info, sep = '[ ,;]+') %>%
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
  mutate(
    info = ifelse(
      str_detect(meta, 'condition'),
      sprintf('(%s) %s', info, desc),
      info
    )
  ) -> meta.nic

meta.nic %<>% select(-desc)

# --------
# And our conclusion
merging$merged_genes %>%
  transmute(merged_id, Resource = '1 BSGatlas',
            Type = type, Name = merged_name) %>%
  gather('meta', 'info', Type, Name) -> meta.bsg

#######################################################

meta <- bind_rows(
  refseq.meta, 
  bsub.meta,
  meta.other,
  subti.overview,
  subti.path,
  meta.bsg,
  meta.nic
) %>%
  arrange(Resource, meta, info)



helper <- function(mg, out) {
  # mg <- 'BSGatlas-gene-5'
  
  meta %>%
    filter(merged_id == mg) %>%
    select(-merged_id) %>%
    mutate(
      # group same categories together
      meta = ifelse(meta == lag(meta, default = ''), '', meta)
    ) -> meta.tab
  
  
  meta.tab %>%
    mutate(row = 1:n()) %>%
    group_by(Resource) %>%
    summarise(from = min(row), to = max(row)) -> meta.groups
    
  
  meta.tab %>%
    select(meta, info) %>%
    kable('html', 
          caption = mg,
          escape = FALSE
          ) %>%
    kable_styling(bootstrap_options = c("striped", "hover"),
                  font_size = 14) -> tab
  meta.groups %>%
    rowwise %>%
    do(foo = {
      j <- .
      tab <<- pack_rows(tab, j$Resource, j$from, j$to,
                        label_row_css = "background-color: #666; color: #fff;") 
      j
    })
  
  to <- sprintf('%s/%s.html', out, mg)
  tab %>%
    save_kable(
      file = to,
      self_contained = FALSE,
      title = paste0('BSGatlas - ', mg)
    )
}
  
jobs <- merging$merged_genes$merged_id
worker <- safely(partial(helper, out = 'data-hub/meta/'))
cl <- makeForkCluster(detectCores())
jobs %>%
  set_names(., .) %>%
  parLapply(cl = cl, X = ., fun = worker) -> result

result %>%
  map('error') %>%
  map(negate(is.null)) %>%
  unlist %>%
  sum
jobs[result %>%
  map('error') %>%
  map(negate(is.null)) %>%
  unlist] %>%
  head
  
stopCluster(cl)


               