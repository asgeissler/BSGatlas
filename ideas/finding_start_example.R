# Find a ncRNA less well described in 
# SubtiWiki,
# how many are there, can we find an awesome bragging example?


source('analysis/00_load.R')

load('data-gff/03_meta.rda')

all.meta %>%
  filter(meta == 'Type') %>%
  select(Resource, id, type = info) %>%
  unique %>%
  group_by(id, Resource) %>%
  summarize(type = invoke(type, .f = str_c)) %>%
  spread(Resource, type) -> foo

foo %>%
  View
  filter()
