source('analysis/00_load.R')


# 7. Make statistics
prots <- c('CDS', 'putative-coding')
stat <- merged_src %>%
  left_join(merged_genes, 'merged_id', suffix = c('', '.merged')) %>%
  mutate(type.equal  = type == type.merged,
         # indicate of change from putative-coding to RNA classification
         type.reclassRNA = (type %in% prots) & !(type.merged %in% prots),
         # and the other way
         type.reclassProt = !(type %in% prots) & (type.merged %in% prots),
         coord.equal = (start == start.merged) & (end == end.merged)) %>%
  transmute(merged_id, src, priority, type.equal, coord.equal,
            type.reclassRNA, type.reclassProt,
            is.coding = type.merged %in% prots,
           `gene type` = ifelse(is.coding, 'coding', 'non-coding')) %>%
  mutate(src = case_when(
    startsWith(src, 'bsubcyc') ~ 'BsubCyc',
    startsWith(src, 'refseq') ~ 'RefSeq',
    startsWith(src, 'rfam') ~ 'rfam',
    src %in% c("nicolas trusted", "nicolas lit review") ~ 'literature review',
    src == 'nicolas lower' ~ 'Nicolas et al predictions'
  )) %>%
  unique

p.src_count <- stat %>%
  count(merged_id, `gene type`) %>%
  count(n, `gene type`) %>%
  mutate_at(c('gene type', 'n'), as.factor) %>%
  right_join(expand(., `gene type`, `n`)) %>% 
  replace_na(list(nn = 0)) %>%
  ggplot(aes(x = n, y = nn,
             fill =  `gene type`,
             group = `gene type`)) +
  ggsci::scale_fill_jama(name = 'Gene Type') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  xlab('Number of sources contributing to a merged gene') +
  ylab('Count of occurences')

p.type_count <- merged_genes %>%
  mutate(type = case_when(
    type %in% c('riboswitch', 'CDS', 'putative-coding',
              'putative-non-coding') ~ type,
    TRUE ~ 'non-coding with known type (tRNA, sRNA, ...)'
  ) %>%
    fct_infreq) %>%
  pull(type) %>%
  fct_count(prop = TRUE) %>% mutate(cs = rev(cumsum(rev(p)))) %>%
  rename(type = f) %>%
  mutate(
    lab = round(p * 100),
    lab = ifelse(lab >= 0,
                 sprintf('%s %%', as.character(lab)),
                 ''),
    at = cs * sum(n) - n / 2
  ) %>%
  ggplot(aes(x = '', y = n, fill = type)) +
  geom_bar(stat = 'identity') +
  ggsci::scale_fill_jama(name = NULL) +
  coord_polar(theta = 'y') +
  geom_label(aes(y = at,# + c(0, cumsum(n)[-length(n)]), 
                x = c(1, 1, 1.3, 1.05, 0.8),
                label = lab), size=3,
            color = 'white') +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text.x=element_blank()
  ) +
  guides(fill = guide_legend(ncol = 2))
p.type_count


cowplot::plot_grid(
  p.type_count, p.src_count,
  labels = 'auto'
) %>% ggsave(filename = 'analysis/02_merge_grid.pdf',
             width = 30, height = 10, units = 'cm')

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
  gather(... = c(type.reclassProt, type.reclassRNA)) %>%
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
stat.src_genes <- merged_src %>%
  mutate(src = case_when(
    startsWith(src, 'bsubcyc') ~ 'BsubCyc',
    startsWith(src, 'refseq') ~ 'RefSeq',
    startsWith(src, 'rfam') ~ 'rfam',
    src %in% c("nicolas trusted", "nicolas lit review") ~ 'literature review',
    src == 'nicolas lower' ~ 'Nicolas et al predictions'
  )) %>%
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
assertthat::assert_that(with(coding.part, `CDS` + `putative-coding` == Total) %>%
                          all(na.rm = TRUE))
  
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
  gather(... = - src) %>%
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
  gather(... = - src) %>%
  spread(src, value) %>%
  mutate_at('key', fct_relevel,
            names(nc.part2) %>%
              discard(`==`, 'src') %>%
              unlist) %>%
  arrange(key) %>%
  mutate_at('key', as.character)
  

# stat for the final product
foo <- merged_genes %>%
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
  gather(... = -src) %>%
  spread(src, value) %>%
  rename(description = key)

# Show full table
bind_rows(coding.part2, nc.part3) %>%
  mutate_all(replace_na, '-') %>%
  # names
  select(
    description = key,
    RefSeq, BsubCyc, 
    `Literature Review` = `literature review`,
    `rfam (various sensetivity levels)` = `rfam`,
    `Nicolas et al predictions`
  ) %>%
  left_join(stat.bsgatlas) %>%
  replace_na(list(BSGatlas = '-')) %>%
  kable('latex', booktabs = TRUE) %>%
  add_indent(c(2, 7, 9:15)) %>%
  row_spec(5, hline_after = TRUE) %>%
  kable_styling() -> tbl_tex

as_image(tbl_tex)

tbl_tex %>%
  str_split('\\n') %>%
  unlist %>%
  discard(str_detect, pattern = '\\{table\\}') %>%
  write_lines(path = 'analysis/02_stat.tex')


# overlaps ofter mreging
x <- merging$merged_genes %>%
  rename(id = merged_id)
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
    )
  ) %>%
  ggplot(aes(x = jaccard)) +
  geom_histogram(bins = 10) +
  facet_wrap(~class, scales = 'free')


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
  left_join(merging$merged_src)
# nothing: that should not have been merged
  
merging$merged_src %>%
  transmute(merged_id, coding = type %in% prots,
            locus = str_remove_all(locus, '_')) %>%
  unique %>%
  count(locus, coding) %>%
  rename(`mergedgenes.per.loci` = n) %>%
  count(coding, mergedgenes.per.loci)
# 11 coding that might/should have been merged but weren't
  # filter(mergedgenes.per.loci > 1, coding) %>%
  # left_join(merging$merged_src %>%
  #             mutate_at('locus', str_remove_all, '_'),
  #           'locus') %>%
  # select(locus, merged_id, src_locus) %>%
  # arrange(locus) %>%
  # left_join(over, c('src_locus' = 'x')) %>%
  # select(locus, merged_id)
  # View

