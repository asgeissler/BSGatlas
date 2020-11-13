P=$(pwd -L)
to="genome.fna.gz"
src="old_genome.fasta.gz"

tmp=$(mktemp -d)
cp $P/$src $tmp/from.fna.gz
cp $P/$to $tmp/to.fna.gz

cd $tmp

for i in from to ; do
    gunzip -c $i.fna.gz > $i.fna
    faToTwoBit $i.fna $i.2bit
    twoBitInfo $i.2bit $i.sizes
done


# align
lastz                             \
    `# output as lav`             \
    --format=lav                  \
    `# alignment options`         \
    --chain --gapped --gfextend   \
    --allocate:traceback=1G       \
    from.fna to.fna               \
    > alignment.lav
# finding the chains
lavToPsl alignment.lav alignment.psl
axtChain -linearGap=medium -psl alignment.psl    \
    from.2bit to.2bit chain.chain
chainAntiRepeat from.2bit to.2bit    \
    chain.chain chain.as
chainMergeSort chain.as > chain.as.sorted
# Netting
chainFilter -noHap chain.as.sorted > chain.as.filtered
chainPreNet chain.as.filtered                \
    from.sizes to.sizes chain.as.prenet
chainNet chain.as.prenet                     \
    from.sizes to.sizes                          \
    stdout /dev/null | netSyntenic stdin net.net
# the actual liftOver file
netChainSubset net.net chain.as.prenet net.subnet
chainStitchId net.subnet result.lift

gzip result.lift

# store output and cleanup
cp $tmp/result.lift.gz $P/liftover.gz
cp $tmp/from.2bit $P/old_genome.2bit
cp $tmp/from.sizes $P/old_genome.sizes
cp $tmp/to.2bit $P/genome.2bit
cp $tmp/to.sizes $P/genome.sizes

rm -rf $tmp
cd $P
