#!/bin/bash

# loop over all .sam files with samtools to produce .bam files
# Author = graham thomas

# Run from data/seqs/trimmed/

for i in *.bam;
do
samtools sort --threads 8 $i -o ${i%.bam}.sorted.bam;
done
