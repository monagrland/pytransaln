#!/bin/bash

set -e

case $1 in
  "build" )
    pip install .
    python -m unittest -v pytransaln.testmodule
    ;;
  "each" )
    pytransaln \
      --input testdata/ASVs_nonchimeras.fasta --code 5 \
      align --how each --threads 24
      # --input transAlign/benchMark_data/RBP3_unaligned.fas --code 1 \
          ;;
  "each_hmm" )
    pytransaln \
      --input testdata/ASVs_nonchimeras.fasta --code 5 --hmm data/coi_arthropod.hmm \
      align --how cons --threads 24
          ;;
  "user" )
    pytransaln \
      --input testdata/ASVs_nonchimeras.fasta --code 5 \
      align --how user --frame 1 --threads 24
          ;;
  "user_hmm" )
    pytransaln \
      --input testdata/ASVs_nonchimeras.fasta --code 5 --hmm data/coi_arthropod.hmm \
      align --how user --frame 1 --threads 24
          ;;
  "stats" )
    pytransaln \
      --input testdata/ASVs_nonchimeras.fasta --code 5 \
      stats 
          ;;
  "stats_hmm" )
    pytransaln \
      --input testdata/ASVs_nonchimeras.fasta --code 5 --hmm data/coi_arthropod.hmm \
      stats 
          ;;
  "clean" )
    rm test.*
    ;;
  *)
    echo build each cons stats clean
    ;;
esac
