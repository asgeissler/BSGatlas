source('analysis/00_load.R')

source('scripts/frame_helpers.R')

refseq.in <- readGenBank('data-raw/refseq.gb.gz')

refseq <- list(
  seq = getSeq(refseq.in),
  seqinfo = seqinfo(refseq.in),
  date =  genbankr::locus(refseq.in) %>%
    strsplit(' ') %>%
    unlist() %>%
    `[`(length(.)),
  coding = cds(refseq.in) %>%
    as.tibble %>%
    transmute(
      start, end, strand, name = gene, locus = locus_tag,
      old_locus = old_locus_tag,
      fnc = map(`function.`, clean_paste),
      title = product,
      description = note,
      ec = map(EC_number, clean_paste),
      type = ifelse(pseudo, 'putative', 'CDS') 
    ),
  noncoding = otherFeatures(refseq.in) %>%
    as.tibble %>%
    transmute(
      start, end, strand, name = gene, locus = locus_tag,
      old_locus = old_locus_tag,
      fnc = map(`function.`, clean_paste),
      title = product,
      type = case_when(
        str_detect(locus, 'rRNA') ~ 'rRNA',
        str_detect(locus, 'tRNA') ~ 'tRNA',
        str_detect(title, 'riboswitch') ~ 'riboswitch',
        str_detect(title, 'T-box') ~ 'riboswitch',
        str_detect(title, 'tmRNA') ~ 'tmRNA',
        name == 'scr' ~ 'SRP',
        str_detect(title, 'putative') ~ 'putative',
        str_detect(title, 'small') ~ 'small',
        str_detect(title, 'antisense') ~ 'asRNA',
        TRUE ~ 'unclear'
      )
    )
)
save(file = 'analysis/01_refseq.rda', refseq)
