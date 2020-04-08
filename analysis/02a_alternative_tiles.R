# this script merged the gene annotation, as collected in `01_tomerge.rda`

source('analysis/00_load.R')

source('scripts/overlap_matching.R')

load('analysis/01_tomerge.rda')

data <- to.merge$todo

# Quarantee that the helper separator has no conflicts
assertthat::assert_that(!any(str_detect(data$id, '%')))
assertthat::assert_that(!any(str_detect(data$src, '%')))

# ids are unique
assertthat::are_equal(data %>% nrow,
                      data %>% unique %>% nrow)

# find overlapping pairs
search <- data %>%
  # ignore senseitve rfam
  filter(priority <= 4) %>%
  transmute(
    id = paste(src, id, sep = '%'),
    start, end, strand, type, name, priority
  )

over <- overlap_matching(search, search) %>%
  filter(
    # only same sense comparisons
    ! antisense,
    # don't compare identity
    x != y,
    # special case do not merge known CDS and riboswitch
    # ! (paste0(x.type, y.type) %in% c('CDSriboswitch', 'riboswitchCDS'))
  )

foo.help <- function(i) {
  fct_recode(
    i,
    'RefSeq\ncoding' = "refseq coding",
    'BsubCyc\ncoding' = "bsubcyc coding",
    'Nicolas et al.\nLit. review' = "nicolas lit review",
    'Nicolas et al.\nTrusted predictions' = "nicolas trusted", 
    'Rfam screen\nConservative' = "rfam conservative",
    'RefSeq\nnon-coding' = "refseq noncoding",
    'BsubCyc\nnon-coding' = "bsubcyc noncoding", 
    'Dar et al.\nriboswitches' = "dar riboswitches",
    'Rfam screen\nmedium' = "rfam medium",
    'Nicolas et al.\nRemaining predictions' = "nicolas lower"
  )
}

to.merge$todo %>%
  select(src, priority) %>%
  unique %>%
  arrange(priority) %>%
  pull(src) %>%
  foo.help %>%
  as.character() -> src.ord

over %>%
  mutate(
    `jaccard similarity` = cut(jaccard, seq(0, 1, 0.1), include.lowest = TRUE)
  ) %>%
  separate(x, c('src.x', 'id.x'), sep = '%') %>%
  separate(y, c('src.y', 'id.y'), sep = '%') %>%
  mutate_at('src.x', foo.help) %>%
  mutate_at('src.y', foo.help) %>%
  mutate_at('src.x', fct_relevel, src.ord) %>%
  mutate_at('src.y', fct_relevel, src.ord) %>%
  # ggplot(aes(x = `jaccard similarity`, fill = `jaccard similarity`)) +
  filter(as.integer(src.x) >= as.integer(src.y)) -> baz

foo <- c("RefSeq\ncoding", "BsubCyc\ncoding")

baz %>%
  filter(src.x %in% foo) %>%
  filter(src.y %in% foo) %>%
  ggplot(aes(x = `jaccard similarity`)) +
  geom_bar() +
  xlab('Jaccard Index') +
  theme_bw(18) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_grid(src.x ~ src.y,
             switch = 'both',
             scales = 'free_y', drop = TRUE) -> A

baz %>%
  filter(!(src.x %in% foo)) %>%
  filter(src.y %in% foo) %>%
  ggplot(aes(x = `jaccard similarity`)) +
  geom_bar() +
  xlab('Jaccard Index') +
  theme_bw(18) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_grid(src.y ~ src.x,
             switch = 'both',
             scales = 'free_y', drop = TRUE) -> B
baz %>%
  filter(!(src.x %in% foo)) %>%
  filter(!(src.y %in% foo)) %>%
  ggplot(aes(x = `jaccard similarity`)) +
  geom_bar() +
  xlab('Jaccard Index') +
  theme_bw(18) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_grid(src.x ~ src.y,
             switch = 'both',
             scales = 'free_y', drop = TRUE) -> C

cowplot::plot_grid(
  cowplot::plot_grid(A, B,
                     nrow = 1,
                     rel_widths = c(1, 2.5),
                     labels = c('(a)', '(b)'), label_size = 24),
  C,
  nrow = 2,
  rel_heights = c(1, 1.5),
  labels = c(NA, '(c)'), label_size = 24
)

ggsave(filename = 'analysis/02a_alternative_tiles.pdf',
       width = 4 * 4, height = 5 * 4,
       units = 'in')


