path <- '/home/projects/nextprod/results/CRISPR/gw_crispr_cas9_grnas/20190321_AEB284x_assemblies/runs/ASM904v1/'

library(tidyverse)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")

guides <- read_tsv('data-raw/CRISPRoff/00_guides.tsv.gz')
gids <- read_tsv('data-raw/CRISPRoff/01_gids.tsv.gz')
guides %>%
  count(guide) %>%
  mutate(multi.match = n > 1) %>%
  select(-n) %>%
  left_join(gids) %>%
  select(gid, multi.match) %>%
  with(set_names(multi.match, gid)) -> multi.flag


###########
# Estimate a cut-off for mismatches up-to which the full off-target list
# should be collected
'_rslurm_offtargets/results_*.RDS' %>%
  Sys.glob() %>%
  `[`(1) %>%
  readRDS() -> rds

rds %>%
  map('bindings') %>%
  map(count, mismatches)  %>%
  map(full_join, tibble(mismatches = 0:8), 'mismatches') %>%
  map(arrange, mismatches) -> foo

foo %>%
  map(pull, n) %>%
  pmap(function(...) {
    i <- c(...)
    i <- replace_na(i, 0L) 
    # list(
    #   m = mean(i),
    #   var = var(i)
    # )
    summary(i)
  }) -> bar
bar %>%
  map(as.list) %>%
  bind_rows %>%
  mutate(mis = 1:n() - 1) -> mis.stat
# > length(foo)
# [1] 7562
# length(foo) / nrow(gids) * 100
# ~ 2%

mis.stat %>%
  # ggplot(aes(as.character(mis), m,
  #            ymin = m - sqrt(var), ymax = m + sqrt(var))) +
  # geom_errorbar() +
  # geom_point(size = 5) +
  ggplot(aes(as.character(mis),
             ymin = `Min.`,
             lower = `1st Qu.`,
             middle = Median,
             upper= `3rd Qu.`,
             ymax = `Max.`)) +
  geom_boxplot(stat = 'identity') +
  xlab('gRNA mismatches') +
  ylab('Number of on-/off-target binding sites') +
  ggforce::facet_zoom(ylim = c(0, 50)) +
  theme_bw(18)

ggsave('data-raw/CRISPRoff/02_nr.pdf', width = 8, height = 6)

  

# up to 4 it is

####################
# helpers
mk.row <- function(xs) {
  paste0(
    '<tr><td>',
    str_c(xs, collapse = '</td><td>'),
    '</td></tr>'
  )
}
assertthat::are_equal(
  mk.row(c('a', 'b')),
  '<tr><td>a</td><td>b</td></tr>'
)

tbl.c <- function(...) {
  paste0(
    '<table>',
    str_c(..., collapse = ''),
    '</table>'
  )
  
}
assertthat::are_equal(
  tbl.c('a', 'b'),
  '<table>ab</table>'
)

tbl2tbl <- function(x) {
  c(
    mk.row(names(x)),
    x %>%
      rowwise() %>%
      do(i = mk.row(.)) %>%
      pull(i)
  ) %>%
    invoke(.f = tbl.c)
}
assertthat::are_equal(
  tbl2tbl(tribble(
    ~ A, ~B,
    1, 2
  )),
  '<table><tr><td>A</td><td>B</td></tr><tr><td>1</td><td>2</td></tr></table>'
)
####################

# Aggregate the various guide off-targets results into the following files

header <- paste(
  'Guide sequence',
  'Number of potential targets',
  'Potential on-targets',
  'Potential off-targets',
  sep = '\t'
) %>%
  paste0('\n')

con <- file('data-raw/CRISPRoff/02_targets.tab', 'w')
offset <- 0L
# keep track of byte offsets when writing
custom.writer <- function(x) {
  writeLines(x, con, sep = '')
  offset <<- offset + nchar(x)
  flush(con)
}

custom.writer(header)
# flush(con)
# close(con)


# The main worker to process all the info
guide.info <- function(x) {
  # x <- rds$AAAAAAAAAAAAAAAAACAAGGG
  i <- x$bindings$gid[1]
  g <- gids$guide[i]
  multi.flag[[i]]
  
  # Collect overall mismatch info
  mis.info <- tbl.c(
    mk.row(paste(0:8, 'mismatches')),
    x$bindings %>%
      count(mismatches) %>%
      right_join(tibble(mismatches = 0:8), 'mismatches') %>%
      arrange(mismatches) %>%
      mutate_at('n', replace_na, 0L) %>%
      pull(n) %>%
      mk.row
  )
  
  
  # short version of potential targets
  x$targets %>%
    group_by(cut.pos) %>%
    summarize(over = str_c(
      sprintf('%s: %s (%s)', type, name, ID),
      collapse = '<br/>'
    )) %>%
    ungroup -> target.over
  
  
  gids %>%
    filter(gid == i) %>%
    select(guide) %>%
    left_join(guides) -> on.guides
  on.guides %>%
    left_join(x$bindings, c('start', 'end', 'strand', 'cut.pos')) %>%
    left_join(target.over, 'cut.pos') %>%
    transmute(
      CRISPRoff = CRISPRspec.y %>% round(3),
      CRISPRspec = CRISPRspec.x %>% round(3),
      Azimuth = Azimuth %>% round(3),
      Position = sprintf('%s (%s)', cut.pos, strand),
      'Potential targets' = over
    ) %>%
    arrange(desc(CRISPRoff)) %>%
    tbl2tbl() -> on.target
  
  x$bindings %>%
    filter(mismatches <= 4) %>%
    anti_join(on.guides, 'cut.pos') %>%
    left_join(target.over, 'cut.pos') %>%
    transmute(
      CRISPRoff = CRISPRspec %>% round(3),
      mismatches,
      'Mismatched nucleotides' = cmp.pam,
      Position = sprintf('%s (%s)', cut.pos, strand),
      'Potential off-targets' = over
    ) %>%
    arrange(desc(CRISPRoff)) %>%
    tbl2tbl() -> off.target
  
  paste(
    paste('gRNA:',  str_replace(g, '(?=...$)', ' ')),
    mis.info,
    on.target,
    off.target,
    sep = '\t'
  ) %>%
    paste0('\n') %>%
    custom.writer()
  return(c(gid = i, offset = offset))
}

# do the work
rds %>%
  head %>%
  lapply(guide.info) %>%
  invoke(.f = bind_rows) -> offs

write_tsv(offs, 'data-raw/CRISPRoff/02_offsets.tsv.gz')

flush(con)
close(con)