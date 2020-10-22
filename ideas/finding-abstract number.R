isoforms$transcripts %>%
  select(trans = id, id = TUs) %>%
  left_join(isoforms$tus, 'id') %>%
  separate_rows(src, sep = ';') %>%
  group_by(src) %>%
  do(i = list(.$trans)) %>%
  ungroup %>%
  with(set_names(i, src)) %>%
  map(1) %>%
  venn::venn()


isoforms$tus %>%
  separate_rows(src, sep = ';') %>%
  group_by(src) %>%
  do(i = list(.$id)) %>%
  ungroup %>%
  with(set_names(i, src)) %>%
  map(1) %>%
  venn::venn()
