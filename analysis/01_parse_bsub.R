source('analysis/00_load.R')

source('scripts/biocyc_parser.R')
source('scripts/frame_helpers.R')

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

# Some pre-processing to simplify annotation

pubs <- bsub_raw$pubs %>%
  transmute(
    pid = `PUBMED-ID`,
    journal = `SOURCE`,
    authors = `AUTHORS`,
    year = YEAR,
    title = TITLE
  ) %>%
  drop_na(pid)


clean.pid <- function(x) {
  x %>%
    strsplit(';') %>%
    map(function(i) {
      strsplit(i, ':') %>%
        map(1) %>%
        unlist
    }) %>%
    map(intersect, y = pubs$pid) %>%
    map(clean_paste) %>%
    unlist
}

bsub_raw$promoters %>%
  transmute(
    id = `UNIQUE-ID`,
    name = `COMMON-NAME`,
    tss = as.integer(`ABSOLUTE-PLUS-1-POS`),
    cite = clean.pid(CITATIONS),
    sigma.box = `PROMOTER-BOXES`
  ) %>%
  drop_na(sigma.box) %>%
  separate_rows(sigma.box, sep = ';') %>%
  left_join(
    bsub_raw$proteins %>%
      select(sigma = `COMMON-NAME`, sigma.box = `RECOGNIZED-PROMOTERS`) %>%
      drop_na(sigma.box) %>%
      separate_rows(sigma.box, sep = ';') %>%
      mutate(sigma = case_when(
sigma == 'RNA polymerase ECF (extracytoplasmic function)-type sigma factor (sigma-Y)' ~ 'Y',
sigma == 'RNA polymerase ECF (extracytoplasmic function)-type sigma factor (sigma(M))' ~ 'M',
sigma == 'RNA polymerase ECF-type sigma factor' ~ 'YlaC',
sigma == 'RNA polymerase ECF(extracytoplasmic function)-type sigma factor (sigma(V))' ~ 'V',
sigma == 'RNA polymerase ECF(extracytoplasmic function)-type sigma factor sigma(X)' ~ 'X',
sigma == 'RNA polymerase ECF(extracytoplasmic function)-type sigma factor W' ~ 'W',
sigma == 'RNA polymerase major sigma-43 factor (sigma-A)' ~ 'A',
sigma == 'RNA polymerase sigma factor' ~ 'I',
sigma == 'RNA polymerase sigma-28 factor (sigma-D)' ~ 'D',
sigma == 'RNA polymerase sigma-30 factor (sigma(H))' ~ 'H',
sigma == 'RNA polymerase sigma-37 factor (sigma(B))' ~ 'B',
sigma == 'RNA polymerase sigma-54 factor (sigma-L)' ~ 'L',
sigma == 'RNA polymerase sporulation-specific sigma factor (pro-&sigma;K)' ~ 'K',
sigma == 'RNA polymerase sporulation-specific sigma factor (sigma-F)' ~ 'G',
sigma == 'RNA polymerase sporulation-specific sigma factor (sigma-G)' ~ 'G',
sigma == 'RNA polymerase sporulation-specific sigma-29 factor (sigma-E)' ~ 'E',
sigma == 'two-subunit &sigma;-factor' ~ 'O'
      )),
    'sigma.box'
  ) %>%
  drop_na(sigma) %>%
  select(-sigma.box) -> TSSs

coords <- bsub_raw$genes %>%
  transmute(
    start = as.integer(`LEFT-END-POSITION`),
    end = as.integer(`RIGHT-END-POSITION`),
    strand = `TRANSCRIPTION-DIRECTION`,
    locus = `UNIQUE-ID`,
    name = `COMMON-NAME`,
    synonyms = SYNONYMS,
    # the id to look up
    fid = `PRODUCT`
  )


coding <- coords %>%
  inner_join(
    bsub_raw$proteins %>%
      transmute(
        fid = `UNIQUE-ID`,
        description = 'COMMON-NAME',
        mw = `MOLECULAR-WEIGHT-SEQ`,
        go = `GO-TERMS`,
        comment = COMMENT,
        cite = clean.pid(`CITATIONS`)
        ),
    'fid'
  )
