# Cook together a stupid excel file
library(tidyverse)

load('analysis/02_merging.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/07_isoforms.rda')

library(xlsx)
library(conflicted)
conflict_scout()
conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")

##############################################################################
out <- createWorkbook()

# bold font: 44546A
# header row: 335BB7
# row highlight: B4C6E7
createFonts <- function(wb) {
  list(
    data = Font(wb, heightInPoints = 16, name='Calibri'),
    head = Font(wb, heightInPoints = 18, name='Calibri',
                color = '#EEEEEE', isBold = TRUE),
    bold = Font(wb, heightInPoints = 18, name='Calibri',
                 isBold = TRUE, color = '#44546A'),
    title = Font(wb, heightInPoints = 22, name='Calibri',
                 isBold = TRUE, color = '#44546A'),
    subtitle = Font(wb, heightInPoints = 18, name='Calibri',
                    isBold = FALSE, isItalic = TRUE)
  )
}
fonts <- createFonts(out)

highlight <- Fill('#B4C6E7', '#B4C6E7', 'SOLID_FOREGROUND')
header <- Fill('#335BB7', '#335BB7', 'SOLID_FOREGROUND')
##############################################################################
# I am a legend

title <- 'BSGatlas annotation Tables (version 1.1)'
sub <- "This file contains multiple worksheets with the annotations provided in the BSGatlas. 
All cooridnates are relative to the RefSeq reference assembly genome ASM904v1.
The contained columns provide the following information:
(Most of the BSGatlas identifiers are clickable hyperlinks directly to the BSGatlas webpage)"									

more <- tribble(
  ~Sheet, ~Content,
  'Genes',
  'Table of the coding/non-coding genes after the gene mergin procedure and indication from which resources the annotaiton originated.',
  'Operons',
  'Table of the operons and their contained genes and the genome coordinates of the genes. The rows are colored alternatingly between operons.',
  'TU',
  'Transcriptional Units (TU). The set of co-transcribed genes without indication of TSS/Terminator.',
  'Transcripts',
  "List of transcripts. The contained genes are indicated with the list of names separated by an underscore '_'. If available, the associated TSS and Terminator are shown. The length of the 5' and 3' UTR and the number internal UTRs are indicates, as well as which annotaiton resources provided evidence for co-transcription. The BSGatlas does not provide UTR lengths shorter than 15 bp, such short UTR lengths are generally stated with '<15'.",
  "5'/3' UTRs",
  "Table of 5'/3' UTRs with the genome coordinates and the nearest up-/down-stream associated gene and TSS/Terminator.",
  'Internal UTRs',
  'Table of internal UTRs and list of which resources gave rise to the annotaiton by indiating co-transcriptional evidence for the neighboring genes.',
  'TSSs',
  'List of TSS positions with indication from which resource they originated. All TSSs have a resource dependent esolution limit. If available, the binding sigma factor and related PubMed Id citations are shown.',
  'Terminators',
  'Table of Terminator annotations which resource annotated them and what free binding energy the terminating RNA structure has.'
)

###
sheet <- createSheet(out, 'Legend')

headerrow <- createRow(sheet, c(1, 3))
headercell <- createCell(headerrow, 1:2)

addMergedRegion(sheet, 1, 1, 1, 2)
lapply(headercell[1,],function(cell) {
  setCellValue(cell, title)
  setCellStyle(cell, CellStyle(out) + fonts$title)
})

addMergedRegion(sheet, 3, 3, 1, 2)
lapply(headercell[2,],function(cell) {
  setCellValue(cell, sub)
  setCellStyle(cell, CellStyle(out) + fonts$subtitle +
                 Alignment(wrapText = TRUE))
})
setRowHeight(getRows(sheet, 3), multiplier = 10)
# first entry bold, other normal
cslist <- list(
  `1` = CellStyle(out) + fonts$data + fonts$bold,
  `2` = CellStyle(out) + fonts$data + Alignment(wrapText = TRUE)
)
addDataFrame(
  as.data.frame(more),
  sheet,
  col.names = TRUE, row.names = FALSE,
  startRow = 5,
  colStyle = cslist,
  colnamesStyle = CellStyle(out) + header + fonts$head
)

setRowHeight(getRows(sheet, 5 + (1:nrow(more))), multiplier = 4.5)

setRowHeight(getRows(sheet, 9), multiplier = 7)

setColumnWidth(sheet, 1, 50)
setColumnWidth(sheet, 2, 100)

##############################################################################
# Genes
merging$merged_src %>%
  select(merged_id, src) %>%
  group_by(merged_id) %>%
  summarize(src = str_c(src, collapse = ', ')) -> foo

merging$merged_genes %>%
  left_join(foo, 'merged_id') %>%
  select(gene = merged_id, name = merged_name, type,
         start, end, strand,
         'annotation origin' = src) -> foo

load('data/03_subtiwiki.rda')


foo %>%
  arrange(start) %>%
  left_join(
    select(subtiwiki$genes, merged_id, locus),
    c('gene' = 'merged_id')
  ) %>%
  select(gene, locus, name, everything()) %>%
  mutate(highlight = 1:n() %% 2 == 0) -> dat
##################################
# some amazing helper
myAdd <- function(xs, lab) {
  sheet <- createSheet(out, lab)
  
  xs %>%
    # also takes care of chr conversion
    mutate_all(replace_na, '') %>%
    select(-highlight) %>%
    as.matrix -> xs2
  
  cb <- CellBlock(sheet, 1, 1, nrow(xs2) + 1, ncol(xs2))
  
  # header
  CB.setRowData(
    cb,
    colnames(xs2),
    1,
    rowStyle = CellStyle(out) + fonts$head + header
  )
  
  # The data
  CB.setMatrixData(
    cb,
    xs2,
    2, 1,
    cellStyle = CellStyle(out) + fonts$data
  )
  #The row-wise highlight
  xs3 <- matrix(FALSE, nrow(xs2), ncol(xs2))
  xs3[xs$highlight, ] <- TRUE
  xs3 <- which(xs3, arr.ind = TRUE)
  CB.setFill(cb, highlight, xs3[, 1] + 1, xs3[, 2])
  
  autoSizeColumn(sheet, 1:ncol(xs2))
  
  return(sheet)
}

sheet <- myAdd(dat, 'Genes')
##############################################################################
# Operons
isoforms$operons %>%
  select(id, transcripts) %>%
  separate_rows(transcripts, sep = ';') %>%
  left_join(isoforms$transcripts, c('transcripts' = 'id')) %>%
  select(operon = id, id = features) %>%
  separate_rows(id, sep = ';') %>%
  inner_join(merging$merged_genes, c('id' = 'merged_id')) %>%
  select(operon, gene = id, name = merged_name, type, start, end, strand) %>%
  mutate(tmp = operon %>%
           strsplit('-operon-') %>%
           map(2) %>%
           unlist %>%
           as.integer) %>%
  arrange(tmp, start, end) %>%
  select(-tmp) %>%
  unique  %>%
  mutate(
    toggle = replace_na(operon != lag(operon), FALSE),
    highlight = accumulate(toggle, xor)
  ) %>%
  select(-toggle) %>%
  myAdd("Operons")
##############################################################################
# TU 
isoforms$tus %>%
  select(id, src, start, genes) %>%
  separate_rows(genes, sep = ';')  %>%
  inner_join(merging$merged_genes, c('genes' = 'merged_id')) %>%
  mutate_at('merged_name', replace_na, 'unnamed_gene') %>%
  select(id, src, start.x, genes, name = merged_name, start.y) %>%
  unique %>%
  group_by(id, src, start.x) %>%
  arrange(id, start.y) %>%
  summarize(cotranscribed_genes = name %>%
              invoke(.f = paste, sep = ', '),
            cotranscribed_genes_id = str_c(genes, collapse = ', ')) %>%
  ungroup %>%
  left_join(
    isoforms$operons %>%
      select(operon = id, id = TUs, o.start = utr.start) %>%
      separate_rows(id, sep = ';'),
    'id'
  ) %>%
  arrange(o.start, start.x) %>%
  mutate_at('src', str_replace_all, ';', ', ') %>%
  select(operon, TU = id, 'annotation origin' = src,
         'co-transcribed genes' = cotranscribed_genes,
         'ids' = cotranscribed_genes_id) %>%
  mutate(
    toggle = replace_na(operon != lag(operon), FALSE),
    highlight = accumulate(toggle, xor)
  ) %>%
  select(-toggle) -> tu.dat

myAdd(tu.dat, 'TU')
##############################################################################
# Trnascripts incl. UTR stat

isoforms$transcripts %>%
  select(id, features) %>%
  separate_rows(features, sep = ';')  %>%
  filter(str_detect(features, 'UTR')) %>%
  left_join(bind_rows(UTRs), c('features' = 'id')) %>%
  mutate(len = end - start + 1) -> foo

foo %>%
  filter(type == 'internal_UTR') %>%
  count(id) %>%
  rename('nr. of internal UTRs' = n) %>%
  left_join(
    foo %>%
      filter(type == '5\'UTR') %>%
      select(id, "5' UTR length" = len),
    'id'
  ) %>%
  left_join(
    foo %>%
      filter(type == '3\'UTR') %>%
      select(id, "3' UTR length" = len),
    'id'
  ) -> utr.stat



# combine together with the TSS/terminators info
isoforms$transcripts %>%
  rename(TU = TUs) %>%
  left_join(tu.dat, 'TU') %>%
  left_join(utr.stat, 'id') %>%
  mutate_at('nr. of internal UTRs', replace_na, 0) %>%
  mutate_at('5\' UTR length', replace_na, '<15') %>%
  mutate_at('3\' UTR length', replace_na, '<15') %>%
  select(
    transcript = id,
    start, end, strand, 
    `nr. of internal UTRs`, `5\' UTR length`, `3\' UTR length`,
    'transcribed genes' = `co-transcribed genes`, 
    'underlying TU' = TU,
    TSS, Terminator,
    operon
  ) %>%
  left_join(
    isoforms$operon %>%
      select(operon = id, o.start = utr.start)
  ) %>%
  arrange(o.start, start) %>%
  select(- o.start) %>%
  mutate(
    toggle = replace_na(operon != lag(operon), FALSE),
    highlight = accumulate(toggle, xor)
  ) %>%
  select(-toggle) %>%
  myAdd('Transcripts')
##############################################################################
UTRs %>%
  bind_rows %>%
  select(-boundary)  %>%
  left_join(
    isoforms$transcripts %>%
      select(transcript = id, id = features) %>%
      separate_rows(id, sep = ';'),
    'id'
  ) %>%
  group_by(type, id, start, end, strand) %>%
  summarize_at('transcript', str_c, collapse = ', ') %>%
  arrange(type, start, desc(end))  %>%
  group_by(type) %>%
  do(i = list(.)) %>%
  with(set_names(i, type)) %>%
  map(1) %>%
  map(select, - type) %>%
  map(mutate, highlight = 1:n() %% 2 == 0) %>%
  map2(names(.), ~ myAdd(.x, .y))
    
##############################################################################
# TSS
bsg.boundaries$TSS %>%
  mutate_at(c('pubmed', 'src', 'sigma'), str_replace_all, ';', ', ') %>%
  rename(TSS = start) %>%
  arrange(TSS) %>%
  select(TSS = id, position = TSS, strand,
         'resolution' = res.limit,
         'binding sigma factors' = sigma,
         'related PubMed entries' = pubmed,
         'annotation origin' = src) %>%
  mutate(highlight = 1:n() %% 2 == 0) %>%
  myAdd('TSS')
##############################################################################
# Terminator

bsg.boundaries$terminator %>%
  arrange(start, desc(end)) %>%
  mutate(energy = sprintf('%.2f', energy)) %>%
  select(Terminator = id, start, end, strand,
         'Free energy [kcal/mol]' = energy,
         'annotation origin' = src) %>%
  mutate(highlight = 1:n() %% 2 == 0) %>%
  myAdd('Terminators')
##############################################################################
##############################################################################
saveWorkbook(out, 'xlsx-table/BSGatlas-v1.1.xlsx')
##############################################################################





