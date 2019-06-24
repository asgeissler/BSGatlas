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

dedup.go <- function(x) {
  x %>%
    strsplit(';') %>%
    map(clean_paste) %>%
    unlist
}


coding <- coords %>%
  separate_rows(fid, sep = ';') %>%
  inner_join(
    bsub_raw$proteins %>%
      transmute(
        fid = `UNIQUE-ID`,
        description = `COMMON-NAME`,
        mw = `MOLECULAR-WEIGHT-SEQ`,
        go = dedup.go(`GO-TERMS`),
        comment = COMMENT,
        cite = clean.pid(`CITATIONS`)
        ),
    'fid'
  )

# Add ec numbers
bsub_raw$enzrxns %>%
  select(fid = ENZYME, react = REACTION) %>%
  left_join(
    bsub_raw$reactions %>%
      select(ec = `EC-NUMBER`, react = `UNIQUE-ID`),
    'react'
  ) %>%
  select(-react) -> ec.lookup


coding %<>% left_join(ec.lookup, 'fid')

coding %<>% select(- fid)

coding %<>% merge_variants('locus')
  
coding %<>% mutate_at(c('start', 'end'), as.integer)

ncRNA <- coords %>%
  inner_join(
    bsub_raw$rnas %>%
      transmute(
        fid = `UNIQUE-ID`,
        description = `COMMON-NAME`,
        comment = COMMENT,
        type = TYPES,
        cite = clean.pid(`CITATIONS`)
        ),
    'fid'
  ) %>%
  mutate(
    type = case_when(
      str_detect(type, 'tRNA') ~ 'tRNA',
      str_detect(type, 'rRNA') ~ 'rRNA',
      str_detect(description, 'antisense') ~ 'asRNA',
      str_detect(description, 'scRNA') ~ 'SRP',
      str_detect(description, 'tmRNA') ~ 'tmRNA',
      str_detect(description, 'small') ~ 'small',
      description %in% c('T-box', 'SAM') ~ 'riboswitch',
      TRUE ~ 'unclear'
    )
  )

ncRNA %<>% select(- fid)

# coords %>%
#   anti_join(coding, 'locus') %>%
#   anti_join(ncRNA, 'locus')
# we got all!

bsub_raw$regulons %>%
  transmute(
    actor.id = `UNIQUE-ID`,
    name = `COMMON-NAME`,
    reg.ids = REGULATES,
    cite = clean.pid(CITATIONS)
  ) %>% 
  separate_rows(reg.ids, sep = ';') %>%
  left_join(
    bsub_raw$regulation %>%
      select(
        reg.ids = `UNIQUE-ID`,
        type = TYPES,
        affected.id = `REGULATED-ENTITY`
      ),
    'reg.ids'
  ) %>%
  # ids are mixed, look up trans id
  left_join(
    bsub_raw$transunits %>%
      select(trans.id = `UNIQUE-ID`, content = `COMPONENTS`) %>%
      separate_rows(content, sep = ';'),
    c('affected.id' = 'content')
  ) %>%
  mutate(trans.id = ifelse(is.na(trans.id), affected.id, trans.id)) %>%
  select(- reg.ids, - affected.id) %>%
  # now get content and keep only genes
  left_join(
    bsub_raw$transunits %>%
      select(trans.id = `UNIQUE-ID`, locus = `COMPONENTS`) %>%
      separate_rows(locus, sep = ';'),
    c('trans.id')
  ) %>%
  semi_join(bind_rows(
    select(coding, locus),
    select(ncRNA, locus)
  )) %>%
  select(- trans.id) %>%
  rename(affected = locus) %>%
  # look up actor name
  left_join(
    bsub_raw$proteins %>%
      select(actor.id = `UNIQUE-ID`, GENE, COMPONENTS) %>%
      gather('foo', 'actor', GENE, COMPONENTS),
    'actor.id'
  ) %>%
  # handle special cases for protein complexes
  mutate(
    actor = ifelse(foo == 'COMPONENTS', 
                   strsplit(actor, '-') %>%
                     map(1) %>%
                     unlist(),
                   actor)
  ) %>%
  select(actor, type, name, affected, cite) %>%
  drop_na %>%
  arrange(actor, type, name, affected) %>%
  unique -> regulations


terminator <- bsub_raw$terminators %>%
  transmute(
    id = `UNIQUE-ID`,
    cite = clean.pid(CITATIONS),
    start = `LEFT-END-POSITION` %>% as.numeric,
    end = `RIGHT-END-POSITION` %>% as.numeric,
    energy1 = str_replace(COMMENT,
                          '^.*Free energy \\(in kcal/mol\\): (-[:digit:]+\\.[:digit:]+).*$',
                          '\\1') %>% as.numeric,
    energy2 = str_replace(COMMENT,
                          '^.*Free energy (-[:digit:]+\\.[:digit:]+) kcal/mol.*$',
                          '\\1') %>% as.numeric,
    energy = case_when(
      !is.na(energy1) ~ energy1,
      !is.na(energy2) ~ energy2,
      T ~ NA_real_
    ),
    wtf = COMMENT,
    rho.independent = (TYPES == 'Rho-Independent-Terminators')
  ) %>%
  select(-wtf, -energy1, -energy2)


# Determine boundaries for transcripts, and get strand
bsub_raw$transunits %>%
  transmute(
    id = `UNIQUE-ID`,
    content = COMPONENTS,
    cite = clean.pid(CITATIONS)
  ) %>%
  separate_rows(content, sep = ';') %>%
  left_join(
    bind_rows(
      transmute(TSSs, content = id, start = tss, end = tss, strand = ''),
      transmute(coords, content = locus, start, end, strand),
      transmute(terminator, content = id, start, end, strand = '')
    ),
    'content'
  ) %>%
  drop_na(start, end) %>%
  group_by(id) %>%
  summarize(
    start = min(start),
    end = max(end),
    strand = clean_paste(strand)
  ) -> trans.bounds

# manual lookup
trans.bounds %<>%
  mutate(
    change = strand == '',
    start = ifelse(change, NA, start),
    end = ifelse(change, NA, end),
    strand = ifelse(change, '-', strand)
  ) %>%
  select(-change)

# attach strand to terminator and TSSs

TSSs %<>%
  drop_na(tss) %>%
  left_join(bsub_raw$promoters %>%
              select(id = `UNIQUE-ID`, trans = `COMPONENT-OF`),
            'id') %>%
  separate_rows(trans, sep = ';') %>%
  inner_join(select(trans.bounds, trans = id, strand), 'trans') %>%
  select(tss, strand, sigma, name, cite)

terminator %<>% 
  left_join(bsub_raw$terminator %>%
              select(id = `UNIQUE-ID`, trans = `COMPONENT-OF`),
            'id') %>%
  separate_rows(trans, sep = ';') %>%
  inner_join(select(trans.bounds, trans = id, strand), 'trans') %>%
  select(start, end, strand, energy, rho.independent, energy)

bsub_raw$transunits %>%
  transmute(
    id = `UNIQUE-ID`,
    genes = COMPONENTS,
    cite = clean.pid(CITATIONS)
  ) %>%
  separate_rows(genes, sep = ';') %>%
  semi_join(coords, c('genes' = 'locus')) %>%
  merge_variants('id') %>%
  left_join(trans.bounds, 'id') -> trans

trans %<>% select(genes, cite, start, end, strand)

bsubcyc <- list(
  coding = coding,
  noncoding = ncRNA,
  regulation = regulations,
  transunit = trans,
  terminator = terminator,
  TSS = TSSs,
  cite = pubs
)

save(file = 'analysis/01_bsubcyc.rda', bsubcyc)

# code on % of TUs without both TSS and terminator
# > bsub_raw$transunits %>%
#   +   transmute(row = 1:n(), part = COMPONENTS) %>%
#   +   separate_rows(part, sep = ';') %>%
#   +   left_join(
#     +     bsub_raw$promoters %>%
#       +       select(part = `UNIQUE-ID`, tss = `ABSOLUTE-PLUS-1-POS`),
#     +     'part') %>%
#   +   left_join(
#     +     bsub_raw$terminators %>%
#       +       select(part = `UNIQUE-ID`, left = `LEFT-END-POSITION`),
#     +     'part') %>%
#   +   gather('key', 'value', tss, left) %>%
#   +   drop_na %>%
#   +   select(row, key) %>%
#   +   unique %>%
#   +   mutate(foo = 'yes') %>%
#   +   unique %>%
#   +   spread(key, foo, fill = 'no') %>%
#   +   count(left, tss) %>%
#   +   arrange(n)
