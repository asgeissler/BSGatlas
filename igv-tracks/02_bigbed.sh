# assumptions:
#  - conda enviornment already loaded
#  - working dir is data-hub
# 
# convert the various bed files to bigbed


bedToBigBed -type=bed9+1 -tab -as=genes.as \
  genes.bed                       \
  genome.info                              \
  -extraIndex=name,ID                      \
  genes.bb


  
bedToBigBed -type=bed9   -tab    \
  operons.bed           \
  genome.info                    \
  -extraIndex=name               \
  operons.bb
  
  
for i in {tss,terminator} ; do
bedToBigBed -type=bed9+2 -tab -as=src_extra.as \
  $i.bed                                       \
  genome.info                                  \
  -extraIndex=name                             \
  $i.bb
done

bedToBigBed -type=bed9 -tab \
  utrs.bed         \
  genome.info               \
  -extraIndex=name          \
  utrs.bb
  
bedToBigBed -type=bed9+2 -tab -as=src_extra.as \
  transcripts.bed                     \
  genome.info                                  \
  -extraIndex=name                             \
  transcripts.bb