table transunits
"my trans units"
    (
    string chrom;          "Reference sequence genome"
    uint   chromStart;     "Start position in genome"
    uint   chromEnd;       "End position in genome"
    string name;           "BSGatlas ID"
    uint score;            "Score, no used"
    char[1] strand;        "+ or - for strand"
    uint thickStart;       "Start of where display should be thick"
    uint thickEnd;         "End of where display should be thick"
    uint reserved;       "User friendly color"
    string src;             "Origin of this transcript annotation"
    string incomplete;             "Transcript might be incomplete?"
    )
