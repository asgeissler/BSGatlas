
source('analysis/00_load.R')

source('scripts/frame_helpers.R')


# 1. Loading raw data
path <- file.path(rprojroot::find_rstudio_root_file(),
          'data-raw', 'subtiwiki')
path %>%
  list.files() %>%
  # not the helper script
  discard(endsWith, 'sh') %>%
  # not the categories, they need special care
  discard(`==`, 'categories.csv.gz') %>%
  set_names(str_split(., '[.]') %>% map(1) %>% unlist) %>%
  map(function(i) file.path(path, i)) %>%
  map(read_csv) -> subti

# special treatment for categories
file.path(path, 'categories.csv.gz') %>%
  read_lines() %>%
  map(str_remove, '"$') %>%
  map(str_remove, '^[,]*\"') %>%
  unlist %>%
  tibble(raw = .) %>%
  separate(raw, c('id', 'desc'), sep = '(?<=SW [0-9.]{1,10})[.] ') -> categories
  
# unify separator
subti$export_fixed %<>%
  mutate_at('names', str_replace_all, ',[ ]*', ';') %>%
  mutate_at('names', str_replace_all, '[ ]*;[ ]*', ';') %>%
  # force primary name to be part of all synonyms
  rowwise %>%
  mutate(names = names %>%
           strsplit(';') %>%
           c(name) %>%
           clean_paste(sep = ';')) %>%
  ungroup %>%
  # and remove 11 entries that don't even have locus
  filter(!is.na(locus))

# 2. match their IDs to the merged gene set
load('analysis/01_refseq.rda')
load('analysis/02_merging.rda')

# loci frome merge_src + the alternatives of refseq
lookup.dat <- merging$merged_src %>%
  bind_rows(merging$merged_src %>%
              select(merged_id, src, alt = locus) %>%
              inner_join(
                refseq[c('noncoding', 'coding')] %>%
                  bind_rows %>%
                  select(alt = locus, locus = old_locus),
                'alt'
              ) %>%
              left_join(merging$merged_src,
                        c('merged_id', 'src', 'alt' = 'locus')) %>%
              select(-alt)
  ) %>%
  gather('q.key', 'q', locus, name) %>%
  separate_rows(q, sep = ';')
  
# make a putative map
putative.map <- subti$export_fixed %>%
  transmute(id = locus, locus, names, description) %>%
  gather('search.key', 'search', locus, names) %>%
  drop_na(search) %>%
  separate_rows(search, sep = ';') %>%
  left_join(lookup.dat, c('search' = 'q'))

# clear case, locus match
clear.map <- putative.map %>%
  filter(((q.key == 'locus') & (search.key == 'locus')) |
           # also accept matchings of nicolas specific indep transcripts
          (startsWith(search, 'S') & (search.key == 'names'))) %>%
  select(id = id, merged_id) %>%
  drop_na %>%
  unique

# rest <- putative.map %>%
#   anti_join(clear.map, 'id')
# 
# rest %>%
#   filter(!str_detect(description, 'UTR'), !str_detect(description, 'intergenic')) %>%
#   filter(!str_detect(description, '3\''), !str_detect(description, '5\'')) %>%
#   View
# # 16 cases that are either puative or we already have


# * how often is the matching ambiguous?
# count(clear.map, id) %>%
#   filter(n > 1) %>%
#   View
# only the 6 short peptides of unclear function and position

# count(clear.map, merged_id) %>%
#   filter(n > 1) %>%
#   View
# #none


# * how many genes are not described by them?
# merging$merged_genes %>%
#   anti_join(clear.map, 'merged_id') %>%
#   # count(type)
#   nrow
# #158
#   type                    n
#   <chr>               <int>
#   1 asRNA                   4
#   2 CDS                    70
#   3 intron                  3
#   4 putative-coding         9
#   5 putative-non-coding     2
#   6 riboswitch             52
#   7 sRNA                   18

# 3 prettify meta data, lookup proper ids

categories %<>% mutate(level = str_count(id, '[.]') + 1)
genes <- subti$export_fixed %>%
  left_join(clear.map, c('locus' = 'id'))
gene.categories <- subti$geneCategories %>%
  left_join(clear.map, c('gene' = 'id')) %>%
  # filter(is.na(merged_id))
  # only 19 weid cases
  drop_na(merged_id) %>%
  select(merged_id, category1, category2, category3, category4, category5) %>%
  gather(... = starts_with('category')) %>%
  left_join(categories, c('value' = 'desc')) %>%
  select(merged_id, category = id) %>%
  drop_na

interactions <- subti$interactions %>%
  select(partner1 = `prot1 locus tag`, partner2 = `prot2 locus tag`) %>%
  mutate(row = 1:n()) %>%
  gather('key', 'locus', partner1, partner2) %>%
  left_join(clear.map, c('locus' = 'id')) %>%
  select(row, key, merged_id) %>%
  filter(!(row == '144' & merged_id == 'BSGatlas-gene-2367')) %>%
  spread(key, merged_id) %>%
  select(- row)


regulations <- subti$regulations %>%
  select(regulon, old.regulator = `regulator locus tag`, mode, old.target = `locus tag`) %>%
  unique %>%
  left_join(clear.map, c('old.regulator' = 'id')) %>%
  rename(regulator = merged_id) %>%
  left_join(clear.map, c('old.target' = 'id')) %>%
  rename(target = merged_id) %>%
  select(-starts_with('old')) %>%
  drop_na %>%
  unique
  
  

# 4 lookup transcripts, complement
transcripts <- subti$operons

transcripts %<>%
  mutate(operon = ifelse(is.na(operon), genes, operon)) %>%
  drop_na %>%
  filter(operon != 'operon') %>%
  rename(id = operon) %>%
  unique

# count(transcripts, id) %>% filter(n > 1)
# okay, unique id

transcripts %<>%
  separate_rows(genes, sep = '(?<!rrn[ED])(->|([ ]*[-][ ]*)|→|–)') %>%
  mutate(
    id = ifelse(id == '<i>yitJ</i>', 'yitJ', id),
    genes = ifelse(genes == '<i>yitJ</i>', 'yitJ', genes),
    genes = ifelse(genes == 'yfkK]', 'yfkK', genes),
    genes = ifelse(genes == 'mneP', 'ydfM', genes)
  )

bad <- anti_join(transcripts, genes, by = c('genes' = 'name'))
# # now only 9 unusable cases
transcripts %<>% unique()

transcripts %<>%
  inner_join(genes, by = c('genes' = 'name')) %>%
  select(transcript = id, gene = merged_id)

transcripts %<>%
  # put a flag on troublesome cases, just 2
  left_join(
    bind_rows(
      bad %>%
        select(transcript = id) %>%
        unique,
      transcripts %>%
        filter(is.na(gene)) %>%
        select(transcript) %>%
        unique
    ) %>%
      unique %>%
      mutate(incomplete.flag = TRUE),
    'transcript'
  ) %>%
  mutate_at('incomplete.flag', replace_na, FALSE) %>%
  drop_na

# transcripts %>%
#   filter(incomplete.flag) %>%
#   select(transcript) %>%
#   unique %>%
#   nrow
# transcripts %>%
#   select(transcript) %>%
#   unique %>%
#   nrow
# 10 problamatic transcripts of 2259


subtiwiki <- list(
  genes = genes %>% drop_na(merged_id),
  categories = categories,
  gene.categories = gene.categories,
  interactions = interactions,
  regulations = regulations,
  transcripts = transcripts
)

save(subtiwiki, file = 'analysis/03_subtiwiki.rda')
