# Make sure liftover did not disrupt tiling-array resolution

library(tidyverse)
library(plyranges)

# the bed version of lited_genome/probes.bb
'data-hub/tiling-array/probes.bed' %>%
  rtracklayer::import() %>%
  mutate(rel = ifelse(strand == '+', start, end)) %>%
  as_tibble %>%
  filter(start < end) %>%
  select(name, strand, rel) -> foo

foo %>%
  arrange(strand, rel) %>%
  group_by(strand) %>%
  mutate(dist = lead(rel) - rel) %>%
  drop_na %>%
  count(strand, dist) -> bar

bar %>%
  filter(dist == 22) %>%
  pull(n) %>%
  sum() %>%
  `/`(sum(bar$n))

bar %>%
  filter(dplyr::between(dist, 20, 25)) %>%
  pull(n) %>%
  sum() %>%
  `/`(sum(bar$n))

bar  %>%
  ungroup %>%
  mutate(cl = cut(dist, c(0, 25, 50, 100, 200, Inf),
                  include.lowest = TRUE)) %>%
  group_by(cl) %>%
  summarize(nr = sum(n)) %>%
  mutate(rt = nr / sum(bar$n) * 100)

         
         