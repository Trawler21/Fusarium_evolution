#!/bin/bash
cd /home/nn/papaya/wzy/eVO/20251028FIG1/BUSCO20251203/proteins/OrthoFinder/Results_Dec16/SinglecopyOG_fasta

mkdir -p OG_aln

echo "Running MAFFT on *.faa ..."
for f in *.faa; do
    echo "Processing $f ..."
    mafft --auto "$f" > OG_aln/$(basename "$f" .faa).aln
done

echo "All alignments done! Results in OG_aln/"

