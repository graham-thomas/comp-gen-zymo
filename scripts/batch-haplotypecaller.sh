#!/bin/bash

# Shell script to loop through HaplotypeCaller.
# run from /data/seqs/trimmed/marked-dup-bam/
# Needed to employ script as nohup wont work on bash for loops
# Author = graham thomas

# Always crosscheck each step with echo statement.
# do echo $i; done

for i in $(ls *.bam);
do \
    gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -I $i \
    -O ${i/_mark_dup.RG.bam/.g.vcf.gz} \
    --emit-ref-confidence GVCF;
done
