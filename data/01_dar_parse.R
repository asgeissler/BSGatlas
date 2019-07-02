library(readxl)

source('analysis/00_load.R')


path <- 'data-raw/Dar_etal_S2-termseq-switches.xlsx'

# excel_sheets(path)
# [1] "term-seq regulators"      "meta-term-seq regulators"

dat <- read_excel(path, 'term-seq regulators')

# > names(dat)
# [1] "organism"                                                
# [2] "Regulatory element"                                      
# [3] "Mechanism of regulation"                                 
# [4] "Gene"                                                    
# [5] "Gene from"                                               
# [6] "Gene to"                                                 
# [7] "Strand"                                                  
# [8] "Gene annotation"                                         
# [9] "TSS"                                                     
# [10] "Term-seq position"                                       
# [11] "Regulator length"                                        
# [12] "Putative Intrinsic terminator upstream of Term-seq site?"
# [13] "Manually corrected? (original site)" 

dar_riboswitches <- dat %>%
  filter(organism == 'B. subtilis') %>%
  transmute(id = paste0('row-', 1:n()),
            name = `Regulatory element`,
            strand = Strand,
            start = ifelse(strand == '+', TSS, `Term-seq position`),
            end = ifelse(strand == '+', `Term-seq position`, TSS),
            upstream.gene = Gene,
            type = 'riboswitch')

save(dar_riboswitches, file = 'data/01_dar-riboswitches.rda')
