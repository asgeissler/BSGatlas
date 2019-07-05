# make summary pages for the gene types and categories

source('analysis/00_load.R')


load('analysis/02_merging.rda')
load('data-hub/03_meta.full.rda')

hubURL = 'https://127.0.0.1/browser/localHubs/BSGatlas_v1_ac/data-hub/hub.txt'
urlPattern = 'https://nextprod.dk/browser/cgi-bin/hgTracks?genome=BSGatlas&hubUrl=%s&position=basu168:%s-%s'

cats <- tribble(
  ~Resource, ~meta,
  '1 BSGatlas', 'Type',
  '5 Rfam', 'Family',
  
  '3 RefSeq', 'Enzyme Classifications',
  '4 BsubCyc', 'Enzyme Classifications',
  '2 SubtiWiki', 'Enzyme Classifications',
  
  '4 BsubCyc', 'Gene Ontology',
  '2 SubtiWiki', 'Category'
)

meta.full %>%
  semi_join(cats) %>%
  select(meta, info, merged_id) %>%
  arrange(meta, info, merged_id) %>%
  unique %>%
  left_join(merging$merged_genes, 'merged_id') %>%
  mutate(url = sprintf(urlPattern, hubURL, start, end),
         link = sprintf('<a href="%s" target="_blank">%s</a>', url, merged_name)
  ) %>%
  select(meta, info, link) -> dat


dat %>%
  group_by(meta, info) %>%
  mutate(
    i = 0:(n() - 1),
    row = i %/% 4,
    col = paste0('col', i %% 4)
  ) %>%
  ungroup %>%
  select(-i) %>%
  mutate(link = sprintf('<td>%s</td>', link)) %>%
  spread(col, link, fill = '') %>%
  mutate(row = paste('<tr>', col0, col1, col2, col3, '</tr>', sep = '\n')) %>%
  select(meta, info, row) %>%
  mutate(
    lagged = lag(info),
    info2 = ifelse(is.na(lagged) | (lagged != info),
                   info,
                   NA),
    info2.id = str_replace_all(info2, ' ', '_'),
    info2.a = ifelse(is.na(info2), 
                     NA,
                     sprintf('<a id="%s">%s</a>', info2.id, info2)),
    info2.toc = ifelse(is.na(info2), 
                     NA,
                     sprintf('<li><a href="#%s">%s</a></li>', info2.id, info2))
  ) -> dat2


dat2$meta %>%
  unique %>%
  map(function(i) {
    dat2 %>%
      filter(meta == i) -> tbl
    
    c('<h1>', i, '</h1>',
      '<ul>',
      tbl$info2.toc %>%
        discard(is.na),
      '</ul>',
      tbl %>%
        group_by(info) %>%
        summarize(
          html = c(first(info2.a), 
                   '<table>',
                   row,
                   '</table>') %>%
            invoke(.f = paste, sep = '\n')
        ) %>%
        pull(html)
    ) %>%
      invoke(.f = paste, sep = '\n') %>%
      strsplit('\n') %>%
      unlist %>%
      write_lines(sprintf('data-hub/summaries/%s.html', i))
  })
