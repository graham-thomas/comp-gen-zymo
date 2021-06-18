#!/bin/bash

# Shell script to batch align trimmed .fastq files with bowtie
# Author = cpad0112
# Source = https://www.biostars.org/p/467896/

# Always crosscheck each step with echo statement.

#This is what I actually used for individual step processing

for i in $(ls ../seqs/trimmed/new/*_R1.fastq); 
do 
    bowtie2 -p 4 --very-sensitive-local -x Ztritici \
    -1 $i \
    -2 ${i/_R1.fastq/_R2.fastq} \
    -S ${i%_R1.fastq}.sam; 
done

# This is what I should have used i.e. a pipeline to go straight to .sorted.bam output.
#for i in $(ls ../seqs/trimmed/new/*_R1.fastq); 
#do 
#    bowtie2 -p 4 --very-sensitive-local -x Ztritici \
#    -1 $i \
#    -2 ${i/_R1.fastq/_R2.fastq} \
#    -S ${i%_R1.fastq}.sam \
#    | samtools sort -O bam -o ${i%_R1.fastq}.sorted.bam; 
#done

# The above does not generate correct bam files, but the first script does generate
# correct .sam files.