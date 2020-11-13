# create BSGatlas annotation GFF file

# please note, that the GFF standard distingushies genes
# in to separate part, the biological region on the DNA, and the 
# biological region on the transcript

# UTRs and riboswitches transcript regions

source('analysis/00_load.R')

source('scripts/frame_helpers.R')

load('analysis/02_merging.rda')
load('analysis/05_bsg_boundaries.rda')
load('analysis/06_utrs.rda')
load('analysis/07_isoforms.rda')

# use same color scheme as in UCSC hub
color.scheme <- read_tsv('data-hub/color_scheme.tsv')

# The meta information used in the gene details view
load('data-gff/01_meta.full.rda')
# > meta.full %>% pull(meta) %>% unique
# [1] "Name"                            "Type"                            "Outside Links"                  
# [4] "Alternative Name"                "Locus Tag"                       "Description"                    
# [7] "Product"                         "Expression neg. correlated with" "Expression pos. correlated with"
# [10] "Highly expressed condition"      "Lowely expressed condition"      "Category"                       
# [13] "Enzyme Classifications"          "Function"                        "Is essential?"                  
# [16] "Isoelectric point"               "Molecular weight"                "Alternative Locus Tag"          
# [19] "Functions"                       "Title"                           "Citation"                       
# [22] "Gene Ontology"                   "Comment"                         "Pathway"                        
# [25] "Family"

####################
# 1. the operons

isoforms$operons %>%
  select(ID = id, start = utr.start, end = utr.end, strand) %>%
  mutate(type = 'operon') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  select(-rgb) %>%
  rename(color = html) -> oper

####################
# 2. the DNA regions


meta.full %>%
  filter(meta == 'Locus Tag') %>%
  select(ID = merged_id, locus_tag = info)  %>%
  unique %>%
  merge_variants('ID', sep=',') -> loci

isoforms$operons %>%
  select(Parent = id, TUs) %>%
  separate_rows(TUs, sep = ';') %>%
  left_join(isoforms$tus, c('TUs' = 'id')) %>%
  select(ID = genes, Parent) %>%
  separate_rows(ID, sep = ';') %>%
  unique -> gene.parent
  

merging$merged_genes %>%
  rename(ID = merged_id, Name = merged_name) %>%
  mutate(rgb.type = case_when(
    type %in% c('putative-coding', 'CDS') ~ 'protein',
    type %in% c('tRNA', 'rRNA', 'riboswitch') ~ type,
    type %in% c('sRNA', 'asRNA') ~ 'shortRNA',
    TRUE ~ 'other'
  )) %>%
  left_join(color.scheme, c('rgb.type' = 'type', 'strand')) %>%
  select(-rgb.type, - rgb) %>%
  rename(color = html) %>%
  left_join(loci, 'ID') %>%
  left_join(gene.parent, 'ID') %>%
  mutate(type = 'gene') -> dna_region
  

####################
# 3. transcripts

bind_rows(
  isoforms$operons %>%
    select(Parent = id, ID = transcripts) %>%
    separate_rows(ID, sep = ';'),
  isoforms$transcripts %>%
    select(ID = id, Parent = features) %>%
    separate_rows(Parent, sep = ';') %>%
    filter(str_detect(Parent, 'gene'))
) -> trans.parent

isoforms$transcripts %>%
  select(ID = id, TUs) %>%
  left_join(isoforms$tus, c('TUs' = 'id')) %>%
  transmute(
    ID,
    comment = sprintf('Based on: %s',
                      str_replace_all(src, ';', ','))
  ) -> trans.origin

isoforms$transcripts %>%
  select(transcribed_from = TSS, ID = id) %>%
  drop_na -> trans.tss

isoforms$transcripts %>%
  transmute(ID = id, start, end, strand,
         type = 'transcript') %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  select(- rgb) %>%
  rename(color = html) %>%
  left_join(trans.parent, 'ID') %>%
  merge_variants('ID', sep=',') %>%
  left_join(trans.origin, 'ID') %>%
  left_join(trans.tss, 'ID') %>%
  mutate_at(c('start', 'end'), as.integer) -> trans


####################
# 4. TSS


bsg.boundaries$TSS %>%
  transmute(ID = id, start, end, strand,
            type = 'TSS',
            comment = sprintf(
              'Sigma: %s, Resolution: %s, Based on: %s, PubMed: %s',
              sigma, res.limit,
              str_replace_all(src, ';', ','),
              str_replace_all(pubmed, ';', ','))) %>%
  left_join(color.scheme, c('type', 'strand')) %>%
  select(- rgb) %>%
  rename(color = html) -> tss


####################
# 5. RNA regions incl UTRs and Terminators

isoforms$transcripts %>%
  transmute(Parent = id, Terminator, features) %>%
  gather('key', 'ID', Terminator, features) %>%
  select(-key) %>%
  drop_na %>%
  separate_rows(ID, sep = ';') %>%
  unique %>%
  merge_variants('ID', sep = ',') -> region.parent

meta.full %>%
  filter(meta %in% c('Alternative Name', 'Category', 'Gene Ontology',
                     'Enzyme Classifications', 'Title')) %>%
  mutate_at('meta', fct_recode,
            'synonyms' = 'Alternative Name',
            'subtiwiki_category' = 'Category',
            'go' = 'Gene Ontology',
            'ec' = 'Enzyme Classifications',
            'comment' = 'Title') %>%
  select(-Resource) %>%
  rename(ID = merged_id) %>%
  bind_rows(
    read_tsv('data-gff/kegg_mapping.tsv') %>%
      transmute(ID = gene, meta = 'kegg_pathways',
                info = sprintf('%s (%s)', pathway_name, pathway)) %>%
      mutate_at('info', str_remove_all, ',')
  ) %>%
  unique %>%
  arrange(ID, meta, info) %>%
  group_by(ID, meta) %>%
  summarize_at('info', clean_paste, sep = ',') -> rna.meta


rna.meta %<>% spread(meta, info)

merging$merged_genes %>%
  rename(ID = merged_id, Name = merged_name) %>%
  left_join(
    dna_region %>%
      select(ID, color),
    'ID'
  ) %>% 
  mutate(
    type = case_when(
      type == 'putative-coding' ~ 'CDS',
      type == 'putative-non-coding' ~ 'ncRNA',
      TRUE ~ type
    )
  ) %>%
  left_join(rna.meta, 'ID') %>%
  mutate(next.ID = paste0(ID, '_transcribed')) -> rna.genes

UTRs %>%
  bind_rows %>%
  # miniscule adjustment on boundary such that there is no single nt 
  # overlap on 5'/3' UTRs
  mutate(
    start = case_when(
      (type == "5'UTR") & (strand == '+') ~ start,
      (type == "5'UTR") & (strand == '-') ~ start + 1L,
      (type == "3'UTR") & (strand == '+') ~ start + 1L,
      (type == "3'UTR") & (strand == '-') ~ start,
      TRUE ~ start
    ),
    end = case_when(
      (type == "5'UTR") & (strand == '+') ~ end - 1L,
      (type == "5'UTR") & (strand == '-') ~ end,
      (type == "3'UTR") & (strand == '+') ~ end,
      (type == "3'UTR") & (strand == '-') ~ end - 1L,
      TRUE ~ end
    )
  ) %>%
  transmute(ID = id, type, start, end, strand,
            comment = type,
            type = 'UTR') -> rna.utrs

bsg.boundaries$terminator %>%
  transmute(
    ID = id, start, end, strand,
    type = 'terminator',
    comment = sprintf(
      'Based on: %s, Energy: %.2f [kcal/mol]',
      str_replace_all(src, ';', ','),
      energy
    )
  ) -> rna.term

bind_rows(
  rna.genes,
  rna.term,
  rna.utrs
) %>%
  left_join(region.parent, 'ID') %>%
  mutate(
    derives_from = ifelse(!is.na(next.ID), ID, NA),
    ID = ifelse(is.na(next.ID), ID, next.ID)
  ) %>%
  select(- next.ID) -> rna_region

####################
# finally. collect all

bind_rows(
  oper,
  dna_region,
  trans,
  tss,
  rna_region,
  # Helper to get the transcript vizualization done
  trans %>%
    transmute(
      Parent = ID,
      ID = paste0(ID, '_exon'),
      type = 'exon',
      start, end, strand
    )
) %>%
  arrange(start, desc(end)) %>%
  mutate_at(c('comment', 'synonyms', 'ec', 'go', 'subtiwiki_category',
              'kegg_pathways'),
            ~ ifelse(is.na(.x), NA, sprintf('"%s"', .x))) -> dat.full

# (optional)
dat.full %>%
  filter(!(type %in% c('exon', 'TSS', 'terminator', 'operon'))) %>%
  select(- transcribed_from) %>%
  mutate_at('Parent', str_remove_all, 'BSGatlas-operon-[0-9]*') %>%
  mutate_at('Parent', str_replace_all, ',,*', ',') %>%
  mutate_at('Parent', str_remove_all, ',$') %>%
  mutate_at('Parent', str_remove_all, '^,') %>%
  mutate_at('Parent', ~ ifelse(.x == '', NA_character_, .x)) -> dat.simple

list(
  'data-gff/BSGatlas_v1.1.gff' = dat.full,
  'data-gff/BSGatlas_v1.1-simplified.gff' = dat.simple
) %>%
  map2(names(.), function(dat, path){
    lines <- bind_rows(
      tibble(line = c(
        '##gff-version 3',
        '#!genome-build ASM904v1',
        '#!genome-build-accession NCBI_Assembly:GCF_000009045.1',
        '##sequence-region NC_000964.3 1 4215606',
        '##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=224308',
        '#',
        '# BSGatlas version 1.1',
        '# A fine grained atlas of the Bacillus Subtillis genome',
        '#',
        '# This is a gff file representation of the BSGatlas for Bacillus subtilis',
        '# according to the NCBI assembly ASM904v1. The columns are:',
        '#',
        '# 1. Chromosome: Here `basu168` for coherence with the browser hub',
        '# 2. Source: Is BSGatlas for this annotation.',
        '# 3. Type: The type for an annotation feature',
        '# 4. Start: Start and ...',
        '# 5. End: ... ending position of a feature.',
        '# 6. Score: Here not used, constant dot as value',
        '# 7. Strand: Strand of a feature. Possible values are {"+", "-", "."}',
        '# 8. Phase: Translation offset phase for coding sequences. Is either "0" or ".".',
        '# 9. List of attributes',
        '# ',
        '# Attibutes are a semicolon separated list of <key>=<value> pairs.',
        '# Alternative values are a comma separated list.',
        '# Example: Name=GeneA;Synonyms=GeneX,GeneY',
        '# ',
        '# The following keys are used:',
        '# ',
        '# ID:                  Unique identifier for each row',
        '# Parent:              List of hierarchical parent relationships',
        '# Derives_from:        Link to gene from which a transcript feature originates',
        '# Name:                The name of a feature',
        '# locus_tag:           A locus identifier for genes; this often is the used identifier in external databases.',
        '# Description:         Free descriptive text of a feature',
        '# comment:             more ellaborate descriptive text of a feature',
        '# synonyms:            List of alternative names',
        '# go:                  List of gene function terms accodting to the Gene Ontology',
        '# ec:                  List of Enzyme Classifications',
        '# subtiwiki_category:  A gene classification system from SubtiWiki',
        '# kegg_pathways:       Information about potentially associated KEGG pathways'
      )),
      dat %>% rowwise() %>% do(line = {
        row <- .
        row <- row %>%
          as.list() %>%
          discard(is.na) %>%
          unlist()
        
        # entries without an explicit column
        rest <- row %>%
          discard(names(.) %in% c('type', 'start', 'end', 'strand'))
        attributes <- paste(
          names(rest), rest, sep = '='
        ) %>%
          clean_paste(sep = ';')
        
        row <- as.list(row)
        paste(
          'basu168',
          'BSGatlas',
          row$type,
          row$start,
          row$end,
          '.',
          row$strand,
          ifelse(row$type == 'CDS', '0', '.'),
          attributes,
          sep = '\t'
        )
      }) %>%
        unnest
    )
    
    pull(lines, line) %>% write_lines(path)
})
 
