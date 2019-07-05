# assumptions:
#  - conda enviornment already loaded
#  - working dir is data-hub
# 
# convert the various bed files to bigbed


bedToBigBed -type=bed9+1 -tab -as=genes.as \
  BSGatlas_genes.bed                       \
  genome.info                              \
  -extraIndex=name,ID                      \
  BSGatlas_genes.bb


bedToBigBed -type=bed9+1 -tab -as=genes.as \
  rfam-extra.bed                           \
  genome.info                              \
  -extraIndex=name,ID                      \
  rfam-extra.bb


for i in individual/*.bed ; do
suff=${i%%[.]bed*}
echo $suff
bedToBigBed -type=bed9+1 -tab -as=genes.as \
  $i                                       \
  genome.info                              \
  -extraIndex=name,ID                      \
  $suff.bb
done


bedToBigBed -type=bed9+1 -tab -as=transunit.as \
  BSGatlas_tus.bed                             \
  genome.info                                  \
  -extraIndex=name                             \
  BSGatlas_tus.bb
  
bedToBigBed -type=bed9   -tab    \
  BSGatlas_operons.bed           \
  genome.info                    \
  -extraIndex=name               \
  BSGatlas_operons.bb
  
  
for i in transbounds/*.bed ; do
suff=${i%%[.]bed*}
echo $suff
bedToBigBed -type=bed9+1 -tab -as=extra.as \
  $i                                       \
  genome.info                              \
  -extraIndex=name                         \
  $suff.bb
done


bedToBigBed -type=bed9         \
  transbounds/nicolas-utrs.bed \
  genome.info                  \
  transbounds/nicolas-utrs.bb
  
for i in BSGatlas_{tss,terminator} ; do
bedToBigBed -type=bed9+2 -tab -as=src_extra.as \
  $i.bed                                       \
  genome.info                                  \
  -extraIndex=name                             \
  $i.bb
done

bedToBigBed -type=bed9 -tab \
  BSGatlas_UTRs.bed         \
  genome.info               \
  -extraIndex=name          \
  BSGatlas_UTRs.bb
  
bedToBigBed -type=bed9+2 -tab -as=src_extra.as \
  BSGatlas_transcripts.bed                     \
  genome.info                                  \
  -extraIndex=name                             \
  BSGatlas_transcripts.bb