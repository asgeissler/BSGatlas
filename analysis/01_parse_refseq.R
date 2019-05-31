source('analysis/00_load.R')

source('scripts/frame_helpers.R')

refseq <- readGenBank('data-raw/refseq.gb.gz')

refseq <- list(
  seq = getSeq(refseq),
  seqinfo = seqinfo(refseq),
  date =  genbankr::locus(refseq) %>%
    strsplit(' ') %>%
    unlist() %>%
    `[`(length(.)),
  coding = cds(refseq) %>%
    as.tibble %>%
    transmute(
      start, end, strand, name = gene, locus = locus_tag,
      fnc = map(`function.`, clean_paste),
      title = product,
      description = note,
      ec = map(EC_number, clean_paste),
      type = ifelse(pseudo, 'putative', 'CDS') 
    ),
  noncoding = otherFeatures(refseq) %>%
    as.tibble %>%
    transmute(
      start, end, strand, name = gene, locus = locus_tag,
      fnc = map(`function.`, clean_paste),
      title = product,
      type = case_when(
        str_detect(locus, 'rRNA') ~ 'rRNA',
        str_detect(locus, 'tRNA') ~ 'tRNA',
        str_detect(title, 'riboswitch') ~ 'riboswitch',
        str_detect(title, 'T-box') ~ 'riboswitch',
        str_detect(title, 'tmRNA') ~ 'tmRNA',
        name == 'scr' ~ 'SRP',
        # parsing description not clear
        TRUE ~ 'ncRNA'
      )
    )
)
save(file = '01_refseq.rda', refseq)
