source('analysis/00_load.R')
source('scripts/overlap_matching.R')


list(
  'data-gff/BSGatlas_v1.1.gff' = 40,
  'data-gff/BSGatlas_v1.0.gff' = 39
) %>% map2(names(.), function(skip, path) {
  # path <- 'data-gff/BSGatlas_v1.1.gff' ; skip <- 40
  read_tsv(path, col_names = FALSE, skip = skip) %>%
    select(gff.type = X3,
           start = X4, end = X5, strand = X7,
           attrs = X9) %>%
    mutate(id = str_extract(attrs, '(?<=^ID=)BSGatlas-[^;]*-[^;]*(?=;)')) -> foo
  
  foo %>%
    filter(!str_detect(id, '_((exon)|(transcribed))$')) %>%
    drop_na(id) %>%
    mutate(type = str_extract(id, '(?<=-).*(?=-)')) %>%
    select(- gff.type, - attrs) -> res
  # query gene biotype
  foo %>%
    filter(str_detect(id, '_transcribed$')) %>%
    mutate_at('id', str_remove, '_transcribed$') %>%
    select(id, type = gff.type) %>%
    mutate(type = ifelse(type == 'CDS',
                         'gene - coding',
                         'gene - non-coding')) -> bio
  
  res %>%
    left_join(bio, 'id') %>%
    mutate(type = ifelse(is.na(type.y), type.x, type.y)) %>%
    select(id, type, start, end, strand)
}) %>%
  set_names(names(.) %>% str_extract('v1..')) -> annot


overlap_matching(annot$v1.1, annot$v1.0) %>%
  filter(!antisense) -> cmp

assertthat::are_equal(
  cmp %>%
    filter(x == y) %>%
    filter(x.type != y.type) %>%
    nrow,
  0
)
# Version did not change type

cmp %>%
  filter(x == y) %>%
  transmute(x, type = x.type,
            change = x5.dist + x3.dist,
            unchanged = mode == 'equal') -> kept

kept %>%
  count(type, unchanged) %>%
  mutate_at('unchanged', ~ ifelse(.x, 'same coordinates',
                                  'changed coordinates')) %>%
  spread(unchanged, n, fill = 0) %>%
  full_join(
    annot$v1.1 %>%
      anti_join(annot$v1.0, 'id') %>%
      count(type) %>%
      rename('newly annotated' = n),
    'type'
  ) %>%
  left_join(
    annot$v1.0 %>%
      anti_join(annot$v1.1, 'id') %>%
      count(type) %>%
      rename('removed annotation' = n),
    'type'
  ) %>%
  select(type, `same coordinates`, `changed coordinates`,
         `newly annotated`, `removed annotation`) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  mutate_if(is.numeric, prettyNum, big.mark = ',') %>%
  rename('annotation biotype' = type) %>%
  mutate(version = 'v1.1') -> stat

kept %>%
  filter(!unchanged)  %>%
  transmute(id = x, type, mode = 'changed coordinates') %>%
  bind_rows(
    anti_join(annot$v1.1, annot$v1.0, 'id') %>%
      transmute(id, type, mode = 'newly annotated'),
    anti_join(annot$v1.0, annot$v1.1, 'id') %>%
      transmute(id, type, mode = 'removed annotation')
  ) %>%
  mutate(version = 'v1.1') -> details
#################
con <- DBI::dbConnect(RSQLite::SQLite(), 'data-gff/04_changelog.sqlite')
copy_to(con, stat, temporary = FALSE)
copy_to(con, details, temporary = FALSE)

DBI::dbListTables(con)

# helper for copy paste
tbl(con, 'stat') %>%
  filter(version == 'v1.1') %>%
  select(- version) %>%
  show_query()
# SELECT `annotation biotype`, `same coordinates`, `changed coordinates`, `newly annotated`, `removed annotation`
# FROM `stat`
# WHERE (`version` = 'v1.1')

tbl(con, 'details') %>%
  filter(version == 'v1.1') %>%
  select(- version) %>%
  arrange(mode, type, id) %>%
  show_query()
# <SQL>
# SELECT `id`, `type`, `mode`
# FROM `details`
# WHERE (`version` = 'v1.1')
# ORDER BY `mode`, `type`, `id`

DBI::dbDisconnect(con)


