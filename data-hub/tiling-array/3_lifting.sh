# Do liftover of everything from old_genome to the one used in the BSGatlas
# Assumption: working dir is data-hub/tiling-array

for i in old_genome/GSM*.bed ; do
    x=$(basename $i)
    echo $x
    liftOver        \
        -minMatch=1 \
        $i liftover.gz lifted_genome/$x.pre /dev/null

    # filter out probes that changed strand
    os="\\$(head -n 1 $i | cut -f 6)"
    grep $os lifted_genome/$x.pre > lifted_genome/$x
    # to bigwig
    cut -f 1,2,3,7 lifted_genome/$x | \
        wigToBigWig stdin genome.sizes lifted_genome/${x%%.bed*}.bw
done


# similar for probes but no filtering and bigbed
echo Probe list
i=old_genome/probes.bed
x=$(basename $i)
liftOver        \
    -minMatch=1 \
    $i liftover.gz lifted_genome/$x /dev/null
bedSort lifted_genome/$x foo
bedToBigBed -type=bed6+ foo genome.sizes \
    lifted_genome/${x%%.bed*}.bb
rm foo

echo "DONE."
