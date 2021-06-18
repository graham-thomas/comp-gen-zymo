#!/bin/bash

# Shell script to run variantFiltration with nohup
# Author = graham thomas
# date = 2021-06-15

gatk VariantFiltration \
-R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
-V genotype.vcf.gz \
-O genotype-varFilt-one.vcf.gz \
--filter-expression "BaseQRankSum < -2.0" \
--filter-name "BaseQRankSumNeg" \
--filter-expression "BaseQRankSum > 2.0" \
--filter-name "BaseQRankSumPos" \
--filter-expression "FS > 0.1" \
--filter-name "FSFilter" \
--filter-expression "Low quality" \
--filter-name "LowQual" \
--filter-expression "MQ < 30.0" \
--filter-name "MQ" \
--filter-expression "MQRankSum < -2.0" \
--filter-name "MQRankSumNeg" \
--filter-expression "MQRankSum > 2.0" \
--filter-name "MQRankSumPos" \
--filter-expression "QD < 20.0" \
--filter-name "QDFilter" \
--filter-expression "QUAL < 250" \
--filter-name "QualFilter" \
--filter-expression "ReadPosRankSum < -2.0" \
--filter-name "ReadPosRankSumNeg" \
--filter-expression "ReadPosRankSum > 2.0" \
--filter-name "ReadPosRankSumPos"

