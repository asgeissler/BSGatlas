source('analysis/00_load.R')

source('scripts/biocyc_parser.R')

path <- file.path(rprojroot::find_rstudio_root_file(),
                  'data-raw', 'bsubcyc')

files <- list.files(path) %>%
  set_names(strsplit(., '[.]') %>% map(1) %>% unlist) %>%
  keep(endsWith, suffix = '.dat.gz')

dat <- file.path(path, files) %>%
  set_names(str_split(files, '[.]') %>% map(1) %>% unlist) %>%
  map(safely(parse_dat))


bsub_raw <- map(dat, 'result')

save(file = 'analysis/01_bsub_raw.rda', bsub_raw)
