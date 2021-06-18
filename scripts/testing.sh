#!/bin/bash

# Shell script demonstrating code testing with echo
# Author = cpad0112
# Source = https://www.biostars.org/p/467896/

# Always crosscheck each step with echo statement.

for i in $(ls /mnt/path/DEC2017/*1.fastq.gz); \
do echo $i ${i/1.fastq/2.fastq} ${i%1.fastq.gz}.bam; \
done

# This should print R1, R2, destination directory for bam and name of the bam file. 
# Then dry-run with appropriate commands as below:

for i in $(ls /mnt/path/DEC2017/*_1.fastq.gz); \
do echo "bowtie2 -p 4 -x testbuild -1 $i -2 ${i/_1.fastq/_2.fastq} \
| samtools view -b -o ${i%_1.fastq.gz}.bam"; \
done

# if the output looks correct, delete echo " " and run for real.

# both of these options work
for i in $(ls analysis/*.stderr); do echo "mv $i ${i/.stderr/.out}"; done
for i in $(ls analysis/*.stderr); do echo "mv $i ${i%.stderr}.out"; done

# real version without echo test
for i in $(ls analysis/*.stderr); do mv $i ${i/.stderr/.out}; done

# From; https://www.biostars.org/p/91736/

for i in *.bam;
 do y=$(echo $i|sed 's/bam/list/g');
 java -Xmx32g -jar $gatk -T IndelRealigner -I $i -R $reference -targetIntervals $y -o ${i%.list}.realignedBam.bam; 
done;

# adapted for my combineGVCFs.sh

for i in $(ls *.g.vcf.gz); 
    do y=$(echo --variant $i); 
    gatk CombineGVCFs \
    -R ../../../../zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    $y
    -O cohort.g.vcf.gz
done

#this doesn't work, throws error "A USER ERROR has occurred: Argument output was missing: Argument 'output' is required"
# just do it the long way and put every sample id in the script