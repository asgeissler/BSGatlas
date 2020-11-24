# The full CRISPRspec pipeline out put path (too big for repo)
path <- '/home/projects/nextprod/results/CRISPR/gw_crispr_cas9_grnas/20190321_AEB284x_assemblies/runs/ASM904v1/'

# Folder structure
# > tree $path
# ├── AAAA_energy_parameters.tsv
# ├── AAAA_GCF_000009045.1_ASM904v1_genomic.suff
# |   ├── AAAAAAAAAAAAAAAAACAAGGG.CRISPRoff.tsv  (contains off-targets per guide)
# |   ├── AAAAAAAAACAAGGGCTCTTAGG.CRISPRoff.tsv
# |   ├── AAAAAAAAATGAAACCATGGCGG.CRISPRoff.tsv
# |   ├── AAAAAAAACAAGGGCTCTTAGGG.CRISPRoff.tsv
# |   ├── AAAAAAAACCGATTTCCGGAAGG.CRISPRoff.tsv
# |   ├── AAAAAAAACGAATTAGGGAAAGG.CRISPRoff.tsv
# |   ├── AAAAAAAACGCTATGCCGAGTGG.CRISPRoff.tsv
# |   ├──  ....
# ├── AAAA_specificity_efficiency_report.tsv     (contains the guides)
# ├── AAAC_energy_parameters.tsv
# ...


# Step 01: load the guides (*_report.tsv)
library(tidyverse)

parse.guides <- function(i) {
  read_tsv(i) %>%
    select(guide = Guide_ID,
           pos = Genomic_position,
           CRISPRspec =  CRISPRspec_specificity_score, 
           Azimuth = Azimuth_ontarget_score) %>%
    drop_na %>%
    separate(pos, c('chr', 'start', 'end', 'strand'), sep = '\\|')
}

file.path(path, '*specificity_efficiency_report.tsv') %>%
  Sys.glob() %>%
  map(parse.guides) %>%
  bind_rows %>%
  select(-chr) %>%
  mutate_at(c('start', 'end'), as.integer) %>%
  mutate(cut.pos = ifelse(strand == '+', end - 6 - 1, start + 6)) -> guides

write_tsv(guides, 'data-raw/CRISPRoff/00_guides.tsv.gz')
