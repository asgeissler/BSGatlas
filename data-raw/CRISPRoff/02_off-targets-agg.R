path <- '/home/projects/nextprod/results/CRISPR/gw_crispr_cas9_grnas/20190321_AEB284x_assemblies/runs/ASM904v1/'

library(tidyverse)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")

guides <- read_tsv('data-raw/CRISPRoff/00_guides.tsv.gz')
gids <- read_tsv('data-raw/CRISPRoff/01_gids.tsv.gz')

####
header <- paste(
  'Guide sequence',
  'Number of potential targets',
  'Note',
  'Potential on-targets',
  'Potential off-targets<br/>(With up to 4 mismatches)',
  sep = '\t'
)
con <- file('data-raw/CRISPRoff/02_targets.tab', 'w')
# keep track of byte offsets when writing
custom.writer <- function(x) {
  writeLines(x, con, sep = '')
  writeLines('\n', con, sep = '')
  flush(con)
  return(nchar(x) + 1)
}
header.bytes <- custom.writer(header)
####
offs <- tibble(gid = 0L, guide = "A", bytes = 0L, .rows = 0)
for (i in Sys.glob('_rslurm_offtargets/results_*.RDS')) {
  cat(i)
  rds <- readRDS(i)
  rds %>%
    map2(names(.), function(xs, g) {
      tibble(
        gid = first(xs$bindings$gid),
        guide = g,
        bytes = custom.writer(xs$meta)
      )
    }) %>%
    invoke(.f = bind_rows) -> sub
  offs <- bind_rows(offs, sub)
  # break
}
flush(con)
close(con)
####
offs %>%
  mutate(
    lagged = lag(bytes, 1, header.bytes),
    offset = cumsum(lagged)
  )  %>%
  select(gid, guide, offset) -> offsets
write_tsv(offsets, 'data-raw/CRISPRoff/02_offsets.tsv.gz')

