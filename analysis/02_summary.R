source('analysis/00_load.R')

# as used during merging, needed to compute nr of changed types
type.helper <- function(i) {
  if (identical(sort(i), c('asRNA', 'sRNA'))) {
    'asRNA'
  } else if (identical(sort(i), c('riboswitch', 'sRNA'))) {
    'sRNA'
  } else {
    stopifnot(
      assertthat::are_equal(1, length(i)) |
        {print(sprintf('unresolved type: %s', i)) ; FALSE}
    )
    i
  }
}

load('analysis/02_merging.rda')

prots <- c('CDS', 'putative-coding')

# A quick helper to show where the ncRNAs are from
# merging$merged_genes %>%
#   filter(type != 'putative-coding', type != 'CDS') %>%
#   select(merged_id) %>%
#   left_join(merging$merged_src) %>%
#   select(merged_id, src) %>%
#   mutate(src = case_when(
#     startsWith(src, 'rfam') ~ 'Rfam',
#     src == 'refseq noncoding' ~ 'RefSeq',
#     src == 'bsubcyc noncoding' ~ 'BsubCyc',
#     src %in% c('dar riboswitches', 'nicolas lit review',
#                'nicolas trusted') ~ 'Literature',
#     src == 'nicolas lower' ~ 'Nicolas et al'
#   )) %>%
#   unique %>%
#   # pull(merged_id) %>% unique %>% length
#   # 441 (check)
#   # remove those that refseq also desribes
#   anti_join(merging$merged_src %>%
#               filter(startsWith(src, 'refseq')),
#             'merged_id')  %>%
#   # select(merged_id) %>% unique %>% nrow
#   # [1] 229 > check
#   # group_by(src) %>%
#   # do(i = list(.$merged_id)) %>%
#   # with(set_names(i, src)) %>%
#   # map(1) %>%
#   # map(unique) %>%
#   # venn::venn(cexil = 1.3, cexsn = 1.3)
#   left_join(merging$merged_genes) %>%
#   filter(type != 'putative-non-coding') %>%
#   # select(merged_id) %>% unique %>% nrow
#   # 61 ?
#   group_by(src) %>%
#   do(i = list(.$merged_id)) %>%
#   with(set_names(i, src)) %>%
#   map(1) %>%
#   map(unique) %>%
#   venn::venn(cexil = 1.3, cexsn = 1.3)


#  Make statistics

# fairly long helper functions to easily compute merging stat
# for both the individual resources, and category wise version
summarize.helper <- function(input, prefix) {
  # input <- merging
  # prefix <- 'detailed'
  
  # for which gene has type/coord changed?
  stat <- input$merged_src %>%
    left_join(input$merged_genes, 'merged_id', suffix = c('', '.merged')) %>%
    mutate(type.equal  = type == type.merged,
           # indicate of change from putative-coding to RNA classification
           type.reclassRNA = (type %in% prots) & !(type.merged %in% prots),
           # and the other way
           type.reclassProt = !(type %in% prots) & (type.merged %in% prots),
           coord.equal = (start == start.merged) & (end == end.merged)) %>%
    transmute(merged_id, src, priority, type.equal, coord.equal,
              type.reclassRNA, type.reclassProt,
              # is.coding = type.merged %in% prots,
              # this should not be according to the merged type!
              # otherwise the stat for re-class cases is messed up on
              # categorization
              is.coding = type %in% prots,
             `gene type` = ifelse(is.coding, 'coding', 'non-coding')) %>%
    unique
  
  # how often is a gene specific in a source
  stat.src_uniq <- stat %>%
    count(merged_id, `gene type`) %>%
    filter(n == 1) %>%
    left_join(stat) %>%
    count(`gene type`, src) %>%
    rename(uniq = n)
  
  # how often has source type the same as merged
  stat.src_type <- stat %>%
    count(src, `gene type`, type.equal) %>%
    mutate(type.equal = ifelse(type.equal, 'same type', 'diff type')) %>%
    spread(type.equal, n, fill = 0)
  
  # total number of genes per resource
  stat.src_total <- stat %>%
    count(src, `gene type`) %>%
    rename(Total = n)
    
  
  # assert: both reclass can't happen at the same time
  assertthat::are_equal(
    0, stat %>% filter(type.reclassProt, type.reclassRNA) %>% nrow
  )
  stat.src_reclass <- stat %>%
    count(src, type.reclassRNA, type.reclassProt) %>%
    filter(type.reclassProt | type.reclassRNA) %>%
    gather('key', 'value', c(type.reclassProt, type.reclassRNA)) %>%
    filter(value) %>%
    select(- value) %>%
    mutate(key = ifelse(key == 'type.reclassProt',
                        'non-coding%putative ncRNA reclassified as coding',
                        'coding%putative coding reclassified as non-coding')) %>%
    separate(key, c('gene type', 'desc'), sep = '%') %>%
    spread(desc, n)
  
  
  # how often has source coord the same as merged
  stat.src_coord <- stat %>%
    count(src, `gene type`, coord.equal) %>%
    mutate(coord.equal = ifelse(coord.equal, 'same coord', 'diff coord')) %>%
    spread(coord.equal, n, fill = 0)
  
  # general providence stat
  stat.src_genes <- input$merged_src %>%
    select(merged_id, src, type) %>%
    unique
    
  # as here is a multi source aggregation, I need to work on the type
  foo <- stat.src_genes %>%
    select(merged_id, src, type) %>%
    filter(!str_detect(type, 'putative')) %>%
    unique %>%
    group_by(merged_id, src) %>%
    summarize_at('type', type.helper) %>%
    ungroup
    
  # now the actual stat
  stat.src_genes %<>%
    left_join(foo, c('merged_id', 'src')) %>%
    mutate(type = ifelse(is.na(type.y), type.x, type.y)) %>%
    count(type, src) %>%
    rename(count = n, `gene type` = type)
  
  
  # summary <-
  tmp <- stat.src_coord %>%
    left_join(stat.src_type,  c('src', 'gene type')) %>%
    left_join(stat.src_total,  c('src', 'gene type')) %>%
    left_join(stat.src_uniq,  c('src', 'gene type')) %>%
    left_join(stat.src_reclass,  c('src', 'gene type'))
  
  
  # make sure numbers add up before continuing
  assertthat::assert_that(with(tmp, `diff coord` + `same coord` == Total) %>% all)
  assertthat::assert_that(with(tmp, `diff type` + `same type` == Total) %>% all)
  
  # The coding part
  coding.part <- tmp %>%
    filter(`gene type` == 'coding') %>%
    left_join(
      # stat on who is putative
      stat.src_genes %>%
        transmute(src, count, specific = `gene type`,
                  `gene type` = ifelse(`gene type` %in% prots, 'coding', 'non-coding')) %>%
        filter(`gene type` == 'coding') %>%
        spread(specific, count)
    )
  # if(prefix == 'detailed') {
  #   coding.part %<>% filter(src != 'nicolas lower')
  # }
  # assertthat::assert_that(with(coding.part, `CDS` + `putative-coding` == Total) %>%
  #                           all(na.rm = TRUE))
  # Because of that weird reclass not given
    
  pct.helper <-function(i, total) {
    sprintf('%s (%s %%)',
            i,
            round(i / total * 100) %>% as.character
    )
  }
  coding.part2 <- coding.part %>%
    replace_na(list(
      'putative-coding' = 0
    )) %>%
    # special case in which a putative ncRNA has been re.classified as coding
    filter(src != 'Nicolas et al predictions') %>%
    transmute(
      src,
      'Protein Coding Genes' = Total,
      'of these hypothetical' = pct.helper(`putative-coding`, Total),
      'Merging refined coordinates of' = pct.helper(`diff coord`, Total),
      'Merging removed hypothetical status for' = pct.helper(`diff type`, Total),
      'Genes uniq to this resource' = ifelse(is.na(uniq),
                                             '-',
                                             pct.helper(`uniq`, Total))
    ) %>%
    # transpose
    gather('key', 'value', - src) %>%
    spread(src, value) %>%
    arrange(desc(1:n()))
    
    
    
    
    
  
  # The non-coding part
  nc.part <- tmp %>%
    filter(`gene type` == 'non-coding') %>%
    left_join(
      # stat on who is putative
      stat.src_genes %>%
        transmute(src, count, specific = `gene type`,
                  `gene type` = ifelse(`gene type` %in% prots, 'coding', 'non-coding')) %>%
        filter(`gene type` == 'non-coding') %>%
        spread(specific, count)
    )
  assertthat::assert_that(with(
    nc.part, pmap(list(
      `putative-non-coding`, `putative ncRNA reclassified as coding`,
      `asRNA`, `intron`, `putative-non-coding`, `riboswitch`, `ribozyme`,
      `rRNA`, `sRNA`, `SRP`, `tmRNA`, `tRNA` 
    ),  sum
    ) == Total) %>%
      all(na.rm = TRUE)
  )
    
  nc.part2 <- nc.part %>%
    mutate_at(c(
      'putative-non-coding', 'asRNA', 'intron', 'riboswitch', 'ribozyme', 'rRNA',
      'sRNA', 'SRP', 'tmRNA', 'tRNA'
      ),
      replace_na, 0
    ) %>%
    transmute(
      src,
      'Non-Coding Structures and Genes' = Total,
      'of these predecited/hypothetical' = pct.helper(`putative-non-coding`, Total),
      'specific ncRNA type known' = pmap(list(
          `asRNA`, `intron`, `riboswitch`, `ribozyme`,
          `rRNA`, `sRNA`, `SRP`, `tmRNA`, `tRNA` 
        ),  sum) %>%
        unlist %>%
        pct.helper(Total),
      'ribosomal RNA (rRNA)' = pct.helper(`rRNA`, Total),
      'transfer RNA (tRNA)' = pct.helper(`tRNA`, Total),
      'small regulatory RNA (sRNA)' = pct.helper(`sRNA`, Total),
      'regulatory antisense RNA (asRNA)' = pct.helper(`asRNA`, Total),
      'riboswitch' = pct.helper(`riboswitch`, Total),
      'self-splicing intron' = pct.helper(intron, Total),
      'other (ribozyme, SRP, tmRNA)' = pct.helper(ribozyme + SRP + tmRNA,
                                                  Total),
      
      
      'Merging refined coordinates of' = pct.helper(`diff coord`, Total),
      'Merging removed hypothetical status for' = pct.helper(`diff type`, Total),
      'Merging reclassified putative ncRNA as coding' = ifelse(
       is.na(`putative ncRNA reclassified as coding`),
       '-',
       `putative ncRNA reclassified as coding`
      ),
      'Genes uniq to this resource' = ifelse(is.na(uniq),
                                             '-',
                                             pct.helper(`uniq`, Total))
    )
  
  # transpose and sort columns correctly
  nc.part3 <- nc.part2 %>%
    gather('key', 'value', - src) %>%
    spread(src, value) %>%
    mutate_at('key', fct_relevel,
              names(nc.part2) %>%
                discard(`==`, 'src') %>%
                unlist) %>%
    arrange(key) %>%
    mutate_at('key', as.character)
    
  
  # stat for the final product
  foo <- input$merged_genes %>%
    select(type) %>%
    mutate(gene.type = ifelse(type %in% prots, 'total.coding', 'total.noncoding'))
  
  stat.bsgatlas <- foo %>%
    count(gene.type) %>%
    set_names(c('type', 'n')) %>%
    bind_rows(count(foo, type)) %>%
    spread(type, n) %>%
    transmute(
      src = 'BSGatlas',
      'Protein Coding Genes' = total.coding,
      'of these hypothetical' = pct.helper(`putative-non-coding`, total.coding),
      
      'Non-Coding Structures and Genes' = total.noncoding,
      'of these predecited/hypothetical' = pct.helper(`putative-non-coding`, total.noncoding),
      'ribosomal RNA (rRNA)' = pct.helper(`rRNA`, total.noncoding),
      'transfer RNA (tRNA)' = pct.helper(`tRNA`, total.noncoding),
      'small regulatory RNA (sRNA)' = pct.helper(`sRNA`, total.noncoding),
      'regulatory antisense RNA (asRNA)' = pct.helper(`asRNA`, total.noncoding),
      'riboswitch' = pct.helper(`riboswitch`, total.noncoding),
      'self-splicing intron' = pct.helper(intron, total.noncoding),
      'other (ribozyme, SRP, tmRNA)' = pct.helper(ribozyme + SRP + tmRNA,
                                                  total.noncoding)
    ) %>%
    gather('key', 'value', -src) %>%
    spread(src, value)
  
  
  res <- bind_rows(coding.part2, nc.part3) %>%
    left_join(stat.bsgatlas, 'key') %>%
    mutate_all(replace_na, '-')
  
  return(res)
}

############################################################################


stat.all <- summarize.helper(merging, 'detailed')

stat.all <- bind_rows(
  merging$merged_src %>%
    select(src, priority) %>%
    unique %>%
    spread(src, priority) %>%
    mutate_all(as.character) %>%
    mutate(key = 'Resource Priority'),
  stat.all
)

# Nice formatting full table
stat.all %<>%
  # names
  select(
    description = key,
    BSGatlas,
    `RefSeq Coding` = `refseq coding`,
    `BsubCyc Coding` = `bsubcyc coding`,
    `RefSeq Non-Coding` = `refseq noncoding`,
    `Nicolas et al's literature review` = `nicolas lit review`,
    `Nicolas et al trusted predictions` = `nicolas trusted`,
    `rfam (conservative)` = `rfam conservative`,
    `BsubCyc Non-Coding` = `bsubcyc noncoding`,
    `Dar et al term-seq` = `dar riboswitches`,
    `rfam (medium)` = `rfam medium`,
    `Nicolas et al predictions` = `nicolas lower`
  )


############################################################################

categorized <- merging
# manually rename sources
categorized$merged_src %<>%
  mutate(src = case_when(
    startsWith(src, 'bsubcyc') ~ 'BsubCyc',
    startsWith(src, 'refseq') ~ 'RefSeq',
    startsWith(src, 'rfam') ~ 'rfam',
    src %in% c("nicolas lit review", 'dar riboswitches') ~ 'literature review',
    src %in% c('nicolas trusted', 'nicolas lower') ~ 'Nicolas et al predictions'
  ))
# reduce entries but only keep positions from highest priority and unify
# if there are more then one, type most specific type
categorized$merged_src %<>%
  group_by(merged_id, src) %>%
  filter(priority == min(priority)) %>%
  summarize(
    priority = first(priority),
    # most cmoprising coordinates
    start = min(start), end = max(end), strand = first(strand),
    # most specific type, or the putative one
    type = {
      specific <- type %>%
        unique %>%
        discard(str_detect, 'putative') %>%
        unlist
      other <- type %>%
        unique
      if (length(specific) == 0) {
        other
      } else {
        type.helper(specific)
      }
    }
    ) %>%
  ungroup

categorized.stat <- summarize.helper(categorized, 'categorized') %>%
  select(
    description = key,
    BSGatlas,
    `RefSeq`,
    `BsubCyc`,
    `Various Literature resources` = `literature review`,
    `rfam (various sensetivity levels)` = `rfam`,
    `Nicolas et al predictions`
  )

merging_stat <- list(
  stat.all = stat.all,
  categorized.stat = categorized.stat
)

save(merging_stat, file = 'analysis/02_merging_stat.rda')

############################################################################
# Make nice tex version of the jsut generated tables


# individual
merging_stat$stat.all %>%
  # provide shorter descriptions with linebreaks where needed
  mutate(description = c(
    "Resource Priority",
    "Protein Coding Genes",
    "\\hspace{1em}putative/predictions",
    "Hypothetical status removed",
    "Coordinates refined",
    "Resource specific genes",
    "Non-Coding RNAs",
    "\\hspace{1em}putative/predictions",
    "known specific types",
    "\\hspace{1em}ribosomal RNA",
    "\\hspace{1em}transfer RNA",
    "\\hspace{1em}small regulatory RNA",
    "\\hspace{1em}regulatory antisense RNA",
    "\\hspace{1em}riboswitch",
    "\\hspace{1em}self-splicing intron",
    "\\hspace{1em}ribozyme, SRP, tmRNA",
    "Coordinate refined",
    "Hypothetical status removed",
    "Reclassifed as coding",
    "Resource specific genes"
)) %>% 
  # mutate(row = 1:n()) %>%
  # select(description, row) %>% View
  slice(1, 2, 3, 4, 6,
        7, 10:16, 8, 18, 19, 20 ) %>%
  mutate_all(replace_na, '--') %>%
  set_names(c(
    " ",
    "Merged",
    "RefSeq Coding",
    "BsubCyc Coding",
    "RefSeq Non-Coding",
    "\\specialcell{Nicolas \\emph{et al.}'s\\\\ literature review}",
    "\\specialcell{Nicolas \\emph{et al.}\\\\ trusted predictions}",
    "\\specialcell{\\emph{Rfam}\\\\ (conservative)}",
    "BsubCyc Non-Coding",
    "\\specialcell{Dar \\emph{et al.}\\\\ term-seq}",
    "\\specialcell{\\emph{Rfam}\\\\ (medium)}",
    "\\specialcell{Nicolas \\emph{et al.}\\\\ predictions}"
  )) %>%
  kableExtra::kable('latex', booktabs = TRUE,
    caption = 'Comparison of the individual resources with the merged result',
    escape = FALSE,
    linesep = '') %>%
  kable_styling(latex_options = 'scale_down') %>%
  row_spec(5, hline_after = TRUE) %>%
  str_split('\\n') %>%
  unlist %>%
  # escape percent
  str_replace_all('%\\)', '\\\\%)') %>%
  # thousand digit mark
  str_replace_all('(\\d)(\\d{3})', '\\1,\\2') %>%
  # without environment
  `[`(4:(length(.) - 1)) %>%
  write_lines(path = 'analysis/02_stat_all.tex')



# categorized
merging_stat$categorized.stat %>%
  # provide shorter descriptions with linebreaks where needed
  mutate(description = c(
    "Protein Coding Genes",
    "\\hspace{1em}putative/predictions",
    "Hypothetical status removed",
    "Coordinates refined",
    "Resource specific genes",
    "Non-Coding RNAs",
    "\\hspace{1em}putative/predictions",
    "known specific types",
    "\\hspace{1em}ribosomal RNA",
    "\\hspace{1em}transfer RNA",
    "\\hspace{1em}small regulatory RNA",
    "\\hspace{1em}regulatory antisense RNA",
    "\\hspace{1em}riboswitch",
    "\\hspace{1em}self-splicing intron",
    "\\hspace{1em}ribozyme, SRP, tmRNA",
    "Coordinates refined",
    "Hypothetical status removed",
    "Reclassifed as coding",
    "Resource specific genes"
)) %>% 
  # mutate(row = 1:n()) %>%
  # select(description, row) %>% View
  slice(1, 2, 3, 5,
        6, 9:15, 7, 17:19) %>%
  mutate_all(replace_na, '--') %>%
  set_names(c(
    " ",
    "Merged",
    "RefSeq",
    "BsubCyc",
    "\\specialcell{Literature\\\\resources}",
    "\\emph{Rfam}",
    "\\specialcell{Nicolas \\emph{et al.}\\\\ predictions}"
  )) %>%
  kableExtra::kable('latex', booktabs = TRUE,
caption = 'Comparison of the resources with the merged results, aggregated by category. The shown comparison statistics are relative to the gene annotation with the highest priority within each category',
    escape = FALSE, linesep = '') %>%
  kable_styling(latex_options = 'scale_down') %>%
  row_spec(4, hline_after = TRUE) %>%
  str_split('\\n') %>%
  unlist %>%
  # escape percent
  str_replace_all('%\\)', '\\\\%)') %>%
  # thousand digit mark
  str_replace_all('(\\d)(\\d{3})', '\\1,\\2') %>%
  # without environment
  `[`(4:(length(.) - 1)) %>%
  write_lines(path = 'analysis/02_stat_categorized.tex')


##############################################################################

# overlaps ofter mreging
x <- merging$merged_genes %>%
  rename(id = merged_id)

source('scripts/overlap_matching.R')
final.over <- overlap_matching(x, x) %>%
  filter(!antisense) %>%
  filter(x < y) %>%
  drop_na(jaccard)

final.over %>%
  # pull(jaccard) %>% summary
  mutate(
    class = case_when(
       (x.type %in% prots) &  (y.type %in% prots) ~ 'protein-protein overlap',
      !(x.type %in% prots) &  (y.type %in% prots) ~ 'RNA-protein overlap',
       (x.type %in% prots) & !(y.type %in% prots) ~ 'RNA-protein overlap',
      !(x.type %in% prots) & !(y.type %in% prots) ~ 'RNA-RNA overlap'
    ),
    `jaccard similarity` = cut(jaccard, seq(0, 1, 0.1), include.lowest = TRUE)
  ) %>%
  # select(x, y, jaccard, x.type, y.type, class, overlap, mode, x.length, y.length) %>%
  # filter(jaccard > 0.1) %>%
  # View
  count(`jaccard similarity`, class) -> foo
expand(foo, `jaccard similarity`, class) %>%
  mutate(n = 0) -> bar
foo %>%
  bind_rows(bar %>% anti_join(foo, c('jaccard similarity', 'class'))) %>%
  ggplot(aes(x = `jaccard similarity`, y = n, fill = `jaccard similarity`)) +
  scale_fill_brewer(palette = 'RdYlBu', direction = -1) +
  geom_bar(stat='identity') +
  xlab(NULL) +
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  facet_wrap( ~ class, scales = 'free') +
  ggtitle('Comparison overlaps within the merged gene set')

ggsave(filename = 'analysis/02_overlaps_merged.pdf',
       width = 20, height = 10, units = 'cm')


# investigation of loci
merging$merged_src %>%
  transmute(merged_id, coding = type %in% prots,
            locus = str_remove_all(locus, '_')) %>%
  unique %>%
  count(merged_id, coding) %>%
  rename(`count.loci` = n) %>%
  # count(coding, count.loci)
  # head
  filter(count.loci > 1, coding) %>%
  select(merged_id) %>%
  left_join(merging$merged_src) %>%
  View
# nothing: that should not have been merged
  
