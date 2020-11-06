source('analysis/00_load.R')
source('scripts/distance_matching.R')

load('data/01_nicolas.rda')
load('analysis/02_merging.rda')


nicolas$upshifts %>%
  transmute(id, start = pos, end = pos, strand) %>%
  distance_matching(
    bind_rows(
      select(merging$merged_genes,
             id = merged_id,
             start, end, strand),
      nicolas$all.features %>%
        filter(is.na(type) | str_detect(type, "3'")) %>%
        mutate(five = ifelse(strand == '+', start, end)) %>%
        select(id = locus, start = five, end = five, strand)
    )
  ) %>%
  filter(!antisense) -> cmp


cmp %>%
  group_by(x) %>%
  summarize_at('distance', min) %>%
  filter(distance <= 500) %>%
  ggplot(aes(distance)) +
  geom_histogram()

cmp %>%
  filter(distance <= 50) %>%
  select(x) %>%
  unique -> foo

cmp %>%
  anti_join(foo) %>%
  select(x) %>%
  unique %>%
  nrow
