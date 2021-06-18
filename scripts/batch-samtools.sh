#!/bin/bash

# loop over all .sam files with samtools to produce .bam files
# Author = graham thomas

# Run from data/seqs/trimmed/

for i in *.sam;
do 
samtools view -bS --threads 8 \
-T /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
$i > ${i%.sam}.bam;
done
