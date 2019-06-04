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
  select(id, strand, pos=posV3, sigma = sig) %>%
  mutate(strand = ifelse(strand > 0, '+', '-'))
# Table S5		Information summary for each feature
feature_names_correlation <- read_excel(paths$S5, sheet = 'Data') %>%
  select(locus = Locus_tag, name = Name,
         strand = Strand, start = StartV3, end = EndV3,
         highly_expressed_conditions = HighExp,
         lowely_expressed_conditions = LowExp,
         positive_correlated_genes = CorName,
         negative_correlated_genes = CorNegName) %>%
  mutate(strand = ifelse(strand > 0, '+', '-'))
# Table S6		Analysis of newly annotated RNA features

# they found these ncRNA in their literature search
ncRNA.lit <- read_excel(paths$S6,
                        sheet = 'Table S6.5',
                        skip = 2,
                        n_max = 25) %>%
  drop_na %>%
  transmute(name = Name, start = StartV3, end = EndV3,
            promoter = str_remove(Promoter, '\\..*'),
            cite = PubMed) %>%
  left_join(select(upshifts, promoter = id, strand)) %>%
  select(- promoter)

# their RNA features validated in at least two other studies
ncRNA.trusted <- read_excel(paths$S6,
                            sheet = 'Table S6.3',
                            skip = 2) %>%
  transmute(name = Name,
            start = StartV3,
            end = EndV3,
            strand = ifelse(Strand > 0, '+', '-'),
            Irnov, 
            Rasm = `Rasm.`) %>%
  gather('cite', 'fid', Irnov, Rasm) %>%
  drop_na() %>%
  # substitute by pubmed ID
  mutate(cite = ifelse(cite == 'Irnov', '20525796', '19682248')) %>%
  select(-fid) %>%
  group_by(name) %>%
  summarize(
    start = min(start), 
    end = max(start),
    strand = first(strand),
    cite = clean_paste(cite)
  ) %>%
  mutate_at(c('start', 'end'), as.integer)

# Table S7		List of high confidence down-shifts and classification
downshifts <- read_excel(paths$S7,
                         sheet = 'Data') %>%
  select(id, strand, pos = posV3, energy) %>%
  mutate(strand = ifelse(strand > 0, '+', '-'))
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
  high = list(
    their.lit_review = ncRNA.lit,
    their.trusted = ncRNA.trusted
  ),
  lower =  features %>%
    drop_na() %>%
    anti_join(
      ncRNA.trusted %>%
        transmute(name = str_replace(name, ';', ';S')) %>%
        separate_rows(name, sep = ';'),
      'name'
    ),
  all.features = correlated_target,
  conditions = conditions,
  upshifts = upshifts,
  downshifts = downshifts,
  predicted.asRNA.targets = asRNA.predicted
) -> nicolas

save(file = 'analysis/01_nicolas.rda', nicolas)
