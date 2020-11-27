# helpers
mk.row <- function(xs) {
  paste0(
    '<tr><td style="font-weight: normal; width: auto;">',
    str_c(xs, collapse = '</td><td style="font-weight: normal; width: auto;">'),
    '</td></tr>'
  )
}
tbl.c <- function(...) {
  paste0(
    '<table class="bedExtraTbl">',
    str_c(..., collapse = ''),
    '</table>'
  )
  
}
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
####################
guide.meta <- function(bindings, targets, guide, gid) {
  
  # Collect overall mismatch info
  mis.info <- tbl.c(
    mk.row(paste(0:8, 'mismatches')),
    bindings %>%
      count(mismatches) %>%
      right_join(tibble(mismatches = 0:8), 'mismatches') %>%
      arrange(mismatches) %>%
      mutate_at('n', replace_na, 0L) %>%
      pull(n) %>%
      mk.row
  )
  
  
  # short version of potential targets
  targets %>%
    group_by(cut.pos) %>%
    summarize(over = str_c(
      sprintf('%s: %s (<a href="https://rth.dk/resources/bsgatlas/details.php?id=%s">%s</a>)', type, name, ID, ID),
      collapse = '<br/>'
    )) %>%
    ungroup -> target.over
  
  
  tibble(guide = guide) %>%
    left_join(guides, 'guide') -> on.guides
  on.guides %>%
    left_join(bindings, c('start', 'end', 'strand', 'cut.pos')) %>%
    left_join(target.over, 'cut.pos') %>%
    arrange(desc(CRISPRoff)) %>%
    transmute(
      'CRISPR off-target score<br/>(CRISPRoff)' = CRISPRoff %>% round(3),
      CRISPRspec = CRISPRspec %>% round(3),
      Azimuth = Azimuth %>% round(3),
      'Cut-Position' = sprintf(
        '<a href=hgTracks?position=basu168:%s-%s>%s (%s)',
        start, end, cut.pos, strand),
      'Potential gene on-targets' = over
    ) %>%
    tbl2tbl() -> on.target
  
  bindings %>%
    filter(mismatches <= 4) %>%
    anti_join(on.guides, 'cut.pos') %>%
    left_join(target.over, 'cut.pos') %>%
    arrange(desc(CRISPRoff)) %>%
    transmute(
      'CRISPR off-target score<br/>(CRISPRoff)' = CRISPRoff %>% round(3),
      mismatches,
      # Brutal way of fixing width
      # replace 'auto;">' with '15em;">
      # Backspace must be manually evaluated below
      '\b\b\b\b\b\b\b15em;">Mismatched nucleotides' = cmp.pam %>%
        str_replace(' ', '&nbsp;') %>%
        paste0('<tt>', ., '</tt>'),
      'Cut-Position' = sprintf(
        '<a href=hgTracks?position=basu168:%s-%s>%s (%s)',
        start, end, cut.pos, strand),
      'Potential gene off-targets' = over
    ) %>%
    mutate_all(replace_na, '') %>%
    tbl2tbl() %>%
    # 7 times for 7 backspaces....
    str_remove('[^\b]\b') %>%
    str_remove('[^\b]\b') %>%
    str_remove('[^\b]\b') %>%
    str_remove('[^\b]\b') %>%
    str_remove('[^\b]\b') %>%
    str_remove('[^\b]\b') %>%
    str_remove('[^\b]\b') -> off.target
  
  paste(
    paste('gRNA:',  str_replace(guide, '(?=...$)', ' ')),
    mis.info,
    ifelse(multi.flag[[gid]],
           'This guide had potentially multiple on-target sites (exact sequence matches)',
           'This guide has only one potential on-target site'),
    on.target,
    off.target,
    sep = '\t'
  )
}
