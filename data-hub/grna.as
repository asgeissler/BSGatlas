table crispr
"CRISPR guides incl. off-targets"
(
  string chrom;       "Reference sequence chromosome or scaffold"
  uint   chromStart;  "Start position in chromosome"
  uint   chromEnd;    "End position in chromosome"
  string name;        "Name or ID of item, ideally both human readable and unique"
  uint score;         "Score (CRISPRspec rounded)"
  char[1] strand;     "+ or - for strand"
  uint thickStart;    "Start of where display should be thick (start codon)"
  uint thickEnd;      "End of where display should be thick (stop codon)"
  uint reserved;       "color"
  string CRISPRspec;  "CRISPR specificity score (CRISPRspec)"
  string Azimuth;      "CRISPR efficiency score (Azimuth)"
  bigint _offset;        "Offset into tab-sep file for details page"
)
