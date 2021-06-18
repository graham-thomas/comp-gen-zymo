#!/bin/bash

# Shell script for gatk CombineGVCFs, looping through all .g.vcf.gz files.
# Author = graham thomas
# Date = 2021-06-15

# From; https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs
# example usage;

#gatk CombineGVCFs \
#   -R reference.fasta \
#   --variant sample1.g.vcf.gz \
#   --variant sample2.g.vcf.gz \
#   -O cohort.g.vcf.gz

# I have over 50 samples so want to use a loop;

gatk CombineGVCFs \
-R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
--variant SRR3740249.g.vcf.gz \
--variant SRR3740250.g.vcf.gz \
--variant SRR3740251.g.vcf.gz \
--variant SRR3740252.g.vcf.gz \
--variant SRR3740253.g.vcf.gz \
--variant SRR3740254.g.vcf.gz \
--variant SRR3740255.g.vcf.gz \
--variant SRR3740256.g.vcf.gz \
--variant SRR3740257.g.vcf.gz \
--variant SRR3740258.g.vcf.gz \
--variant SRR3740259.g.vcf.gz \
--variant SRR3740260.g.vcf.gz \
--variant SRR3740261.g.vcf.gz \
--variant SRR3740262.g.vcf.gz \
--variant SRR3740263.g.vcf.gz \
--variant SRR3740264.g.vcf.gz \
--variant SRR3740265.g.vcf.gz \
--variant SRR3740266.g.vcf.gz \
--variant SRR3740267.g.vcf.gz \
--variant SRR3740268.g.vcf.gz \
--variant SRR3740269.g.vcf.gz \
--variant SRR3740270.g.vcf.gz \
--variant SRR3740271.g.vcf.gz \
--variant SRR3740272.g.vcf.gz \
--variant SRR3740273.g.vcf.gz \
--variant SRR3740274.g.vcf.gz \
--variant SRR3740275.g.vcf.gz \
--variant SRR3740276.g.vcf.gz \
--variant SRR3740277.g.vcf.gz \
--variant SRR3740278.g.vcf.gz \
--variant SRR3740279.g.vcf.gz \
--variant SRR3740280.g.vcf.gz \
--variant SRR3740281.g.vcf.gz \
--variant SRR3740282.g.vcf.gz \
--variant SRR3740283.g.vcf.gz \
--variant SRR3740284.g.vcf.gz \
--variant SRR3740285.g.vcf.gz \
--variant SRR3740286.g.vcf.gz \
--variant SRR3740287.g.vcf.gz \
--variant SRR3740288.g.vcf.gz \
--variant SRR3740299.g.vcf.gz \
--variant SRR3740300.g.vcf.gz \
--variant SRR3740301.g.vcf.gz \
--variant SRR3740303.g.vcf.gz \
--variant SRR3740305.g.vcf.gz \
--variant SRR3740335.g.vcf.gz \
--variant SRR3740336.g.vcf.gz \
--variant SRR3740337.g.vcf.gz \
--variant SRR3740338.g.vcf.gz \
--variant SRR3740342.g.vcf.gz \
--variant SRR3740345.g.vcf.gz \
--variant SRR3740346.g.vcf.gz \
--variant SRR3740347.g.vcf.gz \
--variant SRR3740349.g.vcf.gz \
-O cohort.g.vcf.gz
