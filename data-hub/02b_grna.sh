#!/bin/sh
module load hg
bedSort grna.bed grna.bed
bedToBigBed -tab \
  -type=bed9+3 -as=grna.as \
  grna.bed genome.info grna.bb