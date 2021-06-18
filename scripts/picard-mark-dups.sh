#!/bin/bash

java -jar /home/graham/Bin/picard.jar MarkDuplicates \
-I SRR3740249.sam \
-O SRR3740249_marked_duplicates.sam \
-M SRR3740249_marked_dup_metrics.txt \

