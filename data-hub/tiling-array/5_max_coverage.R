# prepare the max coverage file, is used for annotation verification

library(tidyverse)
library(rtracklayer)
library(plyranges)

library(conflicted)

conflict_prefer('filter', 'dplyr')
conflict_prefer('first', 'dplyr')
conflict_prefer('lag', 'dplyr')
conflict_prefer('desc', 'dplyr')
conflict_prefer('rename', 'dplyr')
conflict_prefer("n", "dplyr")


# Sys.glob('lifted_genome/*_fwd.bw') %>%
#   head %>%
#   map(import) -> foo
# assertthat::are_equal(foo[[1]]$ranges, foo[[2]]$ranges)
# Gives: Ranges are indeed equal after lifting -> simplified aggregation

agg.helper <- function(suffix, s) {
  sprintf('data-hub/tiling-array/lifted_genome/*_%s.bw', suffix) %>%
    Sys.glob() %>%
    map(import) -> foo
  foo %>%
    map(score) %>%
    pmap(max) %>%
    unlist -> agg
  res <- foo[[1]]
  score(res) <- agg
  strand(res) <- s
  return(res)
}
max.c <- c(
  agg.helper('fwd', '+'),
  agg.helper('rev', '-')
)

export.bed(max.c, 'data-hub/tiling-array/max_coverage.bed')
