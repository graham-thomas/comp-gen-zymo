#!/bin/bash

# Shell script to loop through AddOrReplaceReadGroups.
# run from /data/seqs/trimmed/marked-dup-bam/
# Needed to employ script as nohup wont work on bash for loops
# Author = graham thomas

# Always crosscheck each step with echo statement.
# do echo $i; done

for i in $(ls *.bam); 
do \
    gatk AddOrReplaceReadGroups \
    -I $i \
    -O ${i/.bam/.RG.bam} \
    -RGID ${i%_mark_dup.bam} \
    -RGLB wgs \
    -RGPL illumina \
    -RGPU ${i%_mark_dup.bam} \
    -RGSM ${i%_mark_dup.bam} \
    -CREATE_INDEX True;
done