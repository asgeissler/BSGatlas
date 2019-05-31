# Small helper script to unify the library setup in each analysis
# step

# prevent ambiguous methods
library(conflicted)

# make R awesome
library(tidyverse)
# and some additional pipe symbols
library(magrittr)

# per default prettify plots
theme_set(theme_bw())

# Bioconductor packages
library(GenomicRanges)
library(Biostrings)

library(genbankr)


# helper for multicore processing
library(parallel)

# multi assignment
# c(a, b) %<-% f(x)
library(zeallot)

# default behaviour of knitr
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  cache = FALSE
)
# Prettify the display of numbers
knitr::knit_hooks$set(inline = function(x) {
  if(!is.numeric(x)){
    x
  } else {
    paste0('$', prettyNum(round(x, 2), big.mark=","), '$')
  }
})

# nice tables
library(kableExtra)
# don't add linesep automatically
kable <-partial(kable, linesep = "")

# Resolve conflicts

c(
  'clusterApply', 'clusterApplyLB', 'clusterCall', 'clusterEvalQ',
  'clusterExport', 'clusterMap', 'clusterSplit',
  'parApply', 'parCapply', 'parLapply', 'parLapplyLB', 'parRapply',
  'parSapply', 'parSapplyLB'
) %>%
  map(conflict_prefer, 'parallel')

c(
  'collapse', 'combine', 'desc', 'filter', 'first', 'intersect', 'lag',
  'rename', 'setdiff', 'setequal', 'slice', 'union', 'group_rows'
) %>%
  map(conflict_prefer, 'dplyr')


c(
  'expand', 'extract'
) %>%
  map(conflict_prefer, 'tidyr')

c(
  'compact', 'reduce'
) %>%
  map(conflict_prefer, 'purrr')

c(
  'Position'
) %>%
  map(conflict_prefer, 'ggplot2')

c(
  'strsplit'
) %>%
  map(conflict_prefer, 'base')


print("Remaining conflicts:")
x <- conflict_scout()
x[map(x, length) %>%
    keep(`>`, 1) %>%
    names()]
  

