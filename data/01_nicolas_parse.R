library(readxl)

source('analysis/00_load.R')

source('scripts/frame_helpers.R')


# load Nicolas data
path <- file.path(rprojroot::find_rstudio_root_file(),
                  'data-raw', 'nicolas')
files <- list.files(path)
paths <- file.path(path, files) %>%
  as.list %>%
  set_names(str_replace(files, '^Table(S[0-9]+)\\.xlsx', '\\1'))

paths %>%
  map(excel_sheets)

# from http://genome.jouy.inra.fr/basysbio/bsubtranscriptome/

# Table	Location	Description
# Table S1		List of the growth conditions used for the 269 mRNA profiling experiments
conditions <- read_excel(paths$S1, sheet = 'Data') %>%
  select(id = `Condition key *`, desc = `Condition description`,
         collected = `Cells collected`)
# Table S2		Expression levels for all features in the 269 mRNA profiling experiments
features <- read_excel(paths$S2, sheet = 'Data') %>%
  select(name = Name, locus = Locus_tag, start = StartV3, end = EndV3,
         strand = Strand, type = classif) %>%
  mutate(
    strand = ifelse(strand > 0, '+', '-'),
    type = ifelse(type == '-', NA, type)
  )
# Table S4		List of promoter up-shifts with cluster information and TU definition
upshifts <- read_excel(paths$S4, sheet = 'Data') %>%
  transmute(id, strand, pos=posV3, sigma = sig,
            without.tu = (beginTU == -1),
            strand = ifelse(strand > 0, '+', '-'))
# Table S5		Information summary for each feature
feature_names_correlation <- read_excel(paths$S5, sheet = 'Data') %>%
  select(locus = Locus_tag, name = Name,
         strand = Strand, start = StartV3, end = EndV3,
         highly_expressed_conditions = HighExp,
         lowely_expressed_conditions = LowExp,
         positive_correlated_genes = CorName,
         negative_correlated_genes = CorNegName) %>%
  mutate(strand = ifelse(strand > 0, '+', '-'))


# Table S7		List of high confidence down-shifts and classification
downshifts <- read_excel(paths$S7,
                         sheet = 'Data') %>%
  transmute(id, strand, pos = posV3, energy,
            without.tu = (TU == "NA"),
            strand = ifelse(strand > 0, '+', '-'))
# Table S11 	List of all sense-antisense RNA pairs predicted
asRNA.predicted <- read_excel(paths$S11, sheet = 'Data') %>%
  select(locus = Locus_tag...2, expression_correlation = corrcoef,
         pvalue = `pvalue (<0,05)`, target = Locus_tag...28)

# reformat features to contain correlated genes with loci
correlated_target <- features %>% left_join(
  feature_names_correlation,
  by = c('locus', 'name', 'strand', 'start', 'end')
) %>%
  gather(positive_correlated_genes, negative_correlated_genes,
         key = 'correlation', value = 'correlated_name') %>%
  separate_rows(correlated_name, sep = ', ') %>%
  left_join(
    select(features, correlated_locus = locus, name),
    by = c('correlated_name' = 'name')
  ) %>%
  select(- correlated_name) %>%
  group_by(
    # everything except for the correlated_locus column
    name, locus, start, end, strand, type,
    highly_expressed_conditions, lowely_expressed_conditions,
    correlation
  ) %>%
  summarize(
    correlated_locus = clean_paste(correlated_locus, sep = ';')
  ) %>%
  spread(correlation, correlated_locus) %>%
  ungroup


## Divide into condifence levels

# high
list(
  all.features = correlated_target,
  conditions = conditions,
  upshifts = upshifts,
  downshifts = downshifts,
  predicted.asRNA.targets = asRNA.predicted
) -> nicolas

save(file = 'data/01_nicolas.rda', nicolas)
