#!/bin/bash

# Shell script to batch identify PCR duplicates in .sam files with picard
# Author = graham thomas
# Source = adapted from http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

# run from ...

# Always crosscheck each step with echo statement.

for i in $(ls data/seqs/trimmed/*.sorted.bam);
#do echo $i
do
java -jar /home/graham/Bin/picard.jar MarkDuplicates \
    -I $i \
    -O ${i%.sorted.bam}_marked_duplicates.bam \
    -M ${i%.sorted.bam}_marked_dup_metrics.txt;
done