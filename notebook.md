# comp-gen-zymo

## Comparative genomic analysis between strains of *Zymoseptoria tritici*

Notebook on reproducing the first steps of analysis from Hartmann et al., 2017.

### 2021-05-14 

Set up env for project with packages required for analysis;

    conda create --name comp-gen-zymo trimmomatic bowtie2 gatk picard
    conda activate comp-gen-zymo

Generate yaml file for project config;

    conda env export --name comp-gen-zymo > comp-gen-zymo.yml

----------------------------------------------------------------------------------

Starting from the raw sequencing reads (fastq files) from the paper which have been deposited in the NCBI Short Read Archive; BioProject PRJNA327615.

Downloaded SraRunInfo.csv

Am following https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/

Need to create sra list from .csv. 

Try 

    awk -F"," '{print $1}' SraRunInfo.csv > sraAccList.txt

kept header, but otherwise success. After removing header file has 128 entries.

Now run

    prefetch --option-file sraAccList.txt

had to configure sratools (followed quick configure guide from website)

Move all .sra files to comp-gen-zymo directory

    mv **/*.sra ./

Delete redundant SRR directories and try fasterq-dump command again.

    for f in *.sra; do fasterq-dump $f; done

That worked. It created dir 'fasterq.tmp.deathstar.5186', and is producing all the fastq files. Once finished confirmed completion without errors.

    echo $?
    0

Check disk usage.

    du sra/
    107209852
    du fastq/
    544424340

Delete sra files.

    rm -r sra/

### 2021-05-17

Next step in the methods; trim reads using trimmomatic. From trimmomatic mannual (online), under paired-end reads;

    java -jar /home/graham/Bin/trimmomatic-master/bin/trimmomatic.jar 
    PE -threads 8 -basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> 
    | <unpaired output 1> <paired output 2> <unpaired output 2> <step 1> ...

From paper;

    illuminaclip=TruSeq3-PE.fa:2:30:10 leading=10 trailing=10 slidingwindow= 5:10 minlen=50

Downloaded batch trimmomatic processing script to 'scripts/auto_trim' from https://github.com/lakhujanivijay/I-Shell-Do-This/blob/master/automation_of_trimmomatic/auto_trim.sh

Add scripts dir to path;

    export PATH="$PATH:/media/Data/gt293/comp-gen-zymo/scripts"

Need to edit file names to read R1 and R2;

    for file in SRR*; do mv -i "${file}" "${file/_1/_R1}"; done
    for file in SRR*; do mv -i "${file}" "${file/_2/_R2}"; done

Edit auto-trim script to point to my trimmomatic.jar and use options from paper, then execute. For some reason cant just type the command - says command not found. Have to source.

    source ../scripts/auto_trim *.fastq

### 2021-05-19

Having googled TruSeq3-PE.fa the file is found on github, contains four lines. Have copy pasted into new file of same name in 'fastq/'.

Try again;

    cd ../seqs/
    nohup ../scripts/auto_trim.sh *.fastq 2> auto_trim.stderr # 2> auto_trim.stderr redirects stderr to file rather than stdout.

This is working, with every run giving 'TrimmomaticPE: Completed successfully'.

Have tidied project, deleting empty directories. 

Next step is aligning with bowtie2.

### 2021-05-20 

- Consult bowtie2 documentation http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner

### 2021-06-01

Usage e.g.

    bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> \
    | --sra-acc <acc> | b <bam>} -S [<sam>]

Script to run bowtie on multiple samples from https://www.biostars.org/p/467896/ modified for my use and stored in /scripts/batch-bowtie.sh

change permission & run script from /data/zymo-reference;

    chmod +x /scripts/batch-bowtie.sh
    nohup ../../scripts/batch-bowtie.sh &

Success. Using nohup stores sterr in file 'nohup'. Renamed 'bowtie2.out' and moved to /analysis.

### 2021-06-02

Next step of analysis from paper;
- PCR duplicates were identified using Picard tools version 1.118 (http://broadinstitute.github.io/ picard).

Cloned git repo;

    cd /home/graham/Bin/
    git clone https://github.com/broadinstitute/picard.git

instruction's on git page;

#### Building Picard

* First, clone the repo:
```
    git clone https://github.com/broadinstitute/picard.git
    cd picard/
```

* Picard is now built using [gradle](http://gradle.org/). A wrapper script (`gradlew`) is included which will download the appropriate version of gradle on the first invocation.
    
* To build a fully-packaged, runnable Picard jar with all dependencies included, run:
```
    ./gradlew shadowJar
```

* The resulting jar will be in `build/libs`. To run it, the command is:
```
    java -jar build/libs/picard.jar
    
    or
    
    java -jar build/libs/picard-<VERSION>-all.jar 
```    
Encountered error when running
    ./gradlew shadowJar
Tried;
    ./gradlew --debug shadowJar
which threw a more detailed error message. cp paste to google and found; 
https://stackoverflow.com/questions/61328936/gradle-error-could-not-initialize-class-org-codehaus-groovy-runtime-invokerhelp

gives recommendation to install gradle version 6.3
    ./gradlew -version

    ------------------------------------------------------------
    Gradle 5.6
    ------------------------------------------------------------

    Build time:   2019-08-14 21:05:25 UTC
    Revision:     f0b9d60906c7b8c42cd6c61a39ae7b74767bb012

    Kotlin:       1.3.41
    Groovy:       2.5.4
    Ant:          Apache Ant(TM) version 1.9.14 compiled on March 12 2019
    JVM:          14 (Oracle Corporation 14+36-1461)
    OS:           Linux 4.15.0-142-generic amd64

    sudo apt install gradle
    gradle -v

    WARNING: An illegal reflective access operation has occurred
    WARNING: Illegal reflective access by org.codehaus.groovy.reflection.CachedClass (file:/usr/share/java/groovy-all.jar) to method java.lang.Object.finalize()
    WARNING: Please consider reporting this to the maintainers of org.codehaus.groovy.reflection.CachedClass
    WARNING: Use --illegal-access=warn to enable warnings of further illegal reflective access operations
    WARNING: All illegal access operations will be denied in a future release

    ------------------------------------------------------------
    Gradle 4.4.1
    ------------------------------------------------------------

    Build time:   2012-12-21 00:00:00 UTC
    Revision:     none

    Groovy:       2.4.16
    Ant:          Apache Ant(TM) version 1.10.5 compiled on March 28 2019
    JVM:          14 (Oracle Corporation 14+36-1461)
    OS:           Linux 4.15.0-142-generic amd64

    sudo apt update
    gradle -v

    WARNING: An illegal reflective access operation has occurred
    WARNING: Illegal reflective access by org.codehaus.groovy.reflection.CachedClass (file:/usr/share/java/groovy-all.jar) to method java.lang.Object.finalize()
    WARNING: Please consider reporting this to the maintainers of org.codehaus.groovy.reflection.CachedClass
    WARNING: Use --illegal-access=warn to enable warnings of further illegal reflective access operations
    WARNING: All illegal access operations will be denied in a future release

    ------------------------------------------------------------
    Gradle 4.4.1
    ------------------------------------------------------------

    Build time:   2012-12-21 00:00:00 UTC
    Revision:     none

    Groovy:       2.4.16
    Ant:          Apache Ant(TM) version 1.10.5 compiled on March 28 2019
    JVM:          14 (Oracle Corporation 14+36-1461)
    OS:           Linux 4.15.0-142-generic amd64

So it seems I haven't upgraded/updated gradle?? And I'm getting a warning for the `org.codehaus.groovy.reflection.CachedClass` which is what threw the error in the first place.

Going back to the picard install page,
http://broadinstitute.github.io/picard/
it says download the picard.jar file, so uninstall the git clone and try 
    wget https://github.com/broadinstitute/picard/releases/download/2.25.5/picard.jar

#### Test Installation
To test that you can run Picard tools, run the following command in your terminal application, providing either the full path to the picard.jar file:

    java -jar /path/to/picard.jar -h 
Success.

Read through documentation http://broadinstitute.github.io/picard/
identified tool 'MarkDuplicates', go to specific documentation at 
http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

States - 'The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.'

Gives example usage;

    java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt

Next step of analysis from paper;

- Single nucleotide polymorphism (SNP) calling was performed using the Genome Analysis Toolkit (GATK) version 3.3-0.
Go to;

https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4.

    
    cd /home/graham/Bin/
    wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
    unzip *.zip
    rm gatk-4.2.0.0.zip
    export PATH="/home/graham/Bin/gatk-4.2.0.0/:$PATH"
    alias gatk='/home/graham/Bin/gatk-4.2.0.0/gatk'

### Test that it works

    gatk --help

Success.

Installing anaconda3 to /home/graham/anaconda3/.
Installation completed successfully but gives notification message;

If you'd prefer that conda's base environment not be activated on startup, 
set the auto_activate_base parameter to false: 

conda config --set auto_activate_base false

### Run GATK and Picard commands

From https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4

Available tools are listed and described in some detail in the Tool Index section, along with available options. The basic syntax for invoking any GATK or Picard tool is the following:

    gatk [--java-options "jvm args like -Xmx4G go here"] ToolName [GATK args go here]

So for example, a simple GATK command would look like:

    gatk --java-options "-Xmx8G" \
    HaplotypeCaller -R reference.fasta -I input.bam -O output.vcf

You can find more information about GATK command-line syntax here.

### Syntax for Picard tools

When used from within GATK, all Picard tools use the same syntax as GATK. The conversion relative to the "Picard-style" syntax is very straightforward; wherever you used to do e.g. I=input.bam, you now do -I input.bam. So for example, a simple Picard command would look like:

    gatk ValidateSamFile  -I input.bam -MODE SUMMARY

So lets try this with 'MarkDuplicates'.

### 2021-06-03

MarkDuplicates not working with sam files. Try using bam instead of sam.
Implement samtools;

    samtools view -bS --threads 8 -T /media/Data/gt293/zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa SRR3740249.sam > SRR3740249.bam

Success.

- run same command through loop over all .sam files. `batch-samtools.sh`.
- run `sort-bam.sh` to generate *.sorted.bam in data/seqs/trimmed/.
- edit `batch-picard.sh` and execute.

Success.

Back to next stage of analysis. From above;
So for example, a simple GATK command would look like:

    gatk --java-options "-Xmx8G" HaplotypeCaller -R reference.fasta 
    -I input.bam -O output.vcf

And from paper;

    gatk --java-options "-Xmx8G" HaplotypeCaller -emitRefConfidence GVCF
    –variant_index_type LINEAR –variant_index_parameter 128000 –sample_ploidy 1

Combine for my use;

    nohup gatk --java-options "-Xmx8G" HaplotypeCaller -emitRefConfidence GVCF \
    –variant_index_type LINEAR –variant_index_parameter 128000 –sample_ploidy 1 \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.\
    toplevel.fa -I SRR3740249_marked_duplicates.bam -O SRR3740249.vcf &

Failed.

    A USER ERROR has occurred: -emitRefConfidence is not a recognized option

Output saved `/analysis/HaplotypCaller.out`

### 2021-06-05

Check the bam file;

    samtools view SRR3740249_marked_duplicates.bam | head
    
### 2021-06-07

Tried searching GATK, turns out the syntax for that option is --emit-ref-confidence. Check all options correct syntax. and re-run.
Just found this; https://gatk.broadinstitute.org/hc/en-us/community/posts/360076574871-HaplotypeCaller-ERC-GVCF-variant-index-type-LINEAR-

    gatk HaplotypeCaller \
    --reference ~/hs37d5.fa \
    --input ~/GATK/test_joint_genotype/$id”Al_sort_md_sN_RG_bqsr.bam” \
    --output ~/GATK/test_joint_genotype/$id”Al_sort_md_sN_RG_bqsr_HCT.g.vcf” \
    -ERC GVCF \
    --variant_index_type LINEAR \
    --variant_index_parameter 128000

You are using a very old version of GATK which we do not support anymore. Please upgrade to the latest version of GATK v 4.1.9.0. The latest version of GATK does not require "-variant_index_type" argument anymore. That was a temporary workaround for an issue with index compression in a very old version of GATK.

Adjust accordingly;

    nohup gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -I SRR3740249_marked_duplicates.bam \
    -O SRR3740249..g.vcf.gz \
    --emit-ref-confidence GVCF \
    –-variant_index_parameter 128000 \
    –-sample_ploidy 1 &

Error msg;

    A USER ERROR has occurred: Illegal argument value: \
    Positional arguments were provided \
    ',–-variant_index_parameter{128000{–-sample_ploidy{1}' \
    but no positional argument is defined for this tool.

Try excluding `–-variant_index_parameter`

    nohup gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -I SRR3740249_marked_duplicates.bam \
    -O SRR3740249..g.vcf.gz \
    --emit-ref-confidence GVCF \
    –-sample_ploidy 1

Error msg;

    A USER ERROR has occurred: Illegal argument value: \
    Positional arguments were provided ',–-sample_ploidy{1}' \
    but no positional argument is defined for this tool.

Try excluding `–-sample_ploidy`

    nohup gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -I SRR3740249_marked_duplicates.bam \
    -O SRR3740249.g.vcf.gz \
    --emit-ref-confidence GVCF \

Error msg;

    A USER ERROR has occurred: Fasta dict file file:///media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.dict for reference file:///media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa does not exist. Please see http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference for help creating it.

This is progress. Follow guidelines and create dict file.

From; https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format

*Creating the FASTA sequence dictionary file*

    gatk-launch CreateSequenceDictionary -R ref.fasta

My usage;

    gatk CreateSequenceDictionary -R /media/Data/gt293/zymo-reference\
    /zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa

Success.

Retry;

    nohup gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -I SRR3740249_marked_duplicates.bam \
    -O SRR3740249.g.vcf.gz \
    --emit-ref-confidence GVCF

Error msg;

    A USER ERROR has occurred: Argument emit-ref-confidence has a bad value: Can only be used in single sample mode currently. Use the --sample-name argument to run on a single sample out of a multi-sample BAM file.

From; https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups

    samtools view -H sample.bam | grep '^@RG'

My usage;

    samtools view -H SRR3740249_marked_duplicates.bam | grep '^@RG'

Returns nothing i.e. my bam files do not have read groups (@RG) assigned.

From; https://gatk.broadinstitute.org/hc/en-us/articles/360035532352

If your sample collection's BAM files lack required fields, or do not differentiate pertinent factors within the fields, use Picard's AddOrReplaceReadGroups to add or appropriately rename the read group fields in your input file.

Here's an example:

    # throws an error
    gatk HaplotypeCaller \
        -R reference.fasta \
        -I reads_without_RG.bam \
        -o output.vcf

    # fix the read groups
    java -jar picard.jar AddOrReplaceReadGroups \
        I=reads_without_RG.bam \
        O=reads_with_RG.bam \
        SORT_ORDER=coordinate \
        RGID=foo \
        RGLB=bar \
        RGPL=illumina \
        RGSM=Sample1 \
        CREATE_INDEX=True

    # runs without error
    gatk HaplotypeCaller \
        -R reference.fasta \
        -I reads_with_RG.bam \
        -o output.vcf

Note - read group `ID`s are composed using the flowcell name and lane number, making them a globally unique identifier across all sequencing data in the world.

I have emailed Jahcub Trew and Paul O'Neall for comment. Jahcub came back with;
The two main things before haplotype caller are marking duplicates (which you've done) and adding read groups, so that you know which samples are which if you pool them in later analyses for joint SNP calls.

Command is something like this:

    java -jar picard.jar AddOrReplaceReadGroups \
    I=(name of sample) \
    O=(name of sample) \
    RGID=(name of sample) \
    RGLB=(Sequence library used) \
    RGPL=(sequencing platform used) \
    RGPU=(sequencing run) \
    RGSM=(name of sample) \
    CREATE_INDEX=True

They're not all "name of sample" officially, but they amount to the same thing. It doesn't make a  difference to the analysis as long as you use IDs that mean something to you.

My usage;

    gatk AddOrReplaceReadGroups \
    -I SRR3740249_mark_dup.bam \
    -O SRR3740249_mark_dup.RG.bam \
    -RGID SRR3740249 \
    -RGLB wgs \
    -RGPL illumina \
    -RGPU SRR3740249 \
    -RGSM SRR3740249 \
    -CREATE_INDEX True

Success.

Loop through all files.

    nohup \
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

Failed, nohup doesn't work with bash for-loops. Write it to a script and run;

    nohup add-read-groups-to-bam.sh

Success. Saved stderr to `/analysis/addReadGroup.out`

Retry;

    samtools view -H SRR3740249_mark_dup.RG.bam | grep '^@RG'

Returns;

    @RG     ID:SRR3740249   LB:wgs  PL:illumina     SM:SRR3740249   PU:SRR3740249

Success.

Now re-run HaplotypeCaller;

    nohup gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -I SRR3740249_mark_dup.RG.bam \
    -O SRR3740249.g.vcf.gz \
    --emit-ref-confidence GVCF

Success. stderr saved to `/analysis/haplotypeCaller-error3.out` but may need to rename as, even though there are many warnings in the output, it did generate the expected files and when inspecting them with `zcat SRR3740249.g.vcf.gz | less` there is some sensible info listing e.g.

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR3740249
    1       1       .       C       <NON_REF>       .       .       END=87  GT:DP:GQ:MIN_DP:PL      0/0:1:0:0:0,0,0

Runtime = 92.26 minutes.

Temporarily move processed files to done dir;

    mkdir done
    mv SRR3740249* done/

Write shell script to loop over remaining files;

    nohup ../../../../scripts/batch-haplotypecaller.sh 

### 2021-06-08

That is running, but the stderr for each run is ~10,000 lines. I should have included in the loop to 2>stderr for each run. Currently all stderr is appended to nohup.out.

Create seperate files for haplotypecaller stderr in `/analysis/haplotypeCaller/`

    touch haplotypeCaller_SRR37402{49..88}.out

This creates all the files in one step.

### 2021-06-10

Realised I downloaded the first 40 samples without checking their geographic source. Having now checked this I've learnt the samples I am working with are primarily from Israel and Australia (see file srr-strain-geo-old.txt). Ideally I would move forward with 10 strains from each of the 4 geographic locations. Therefore acquire extra samples in file srr-strain-geo-new.txt.

    cd /data/seqs/raw-seqs/
    prefetch --option-file ../../sraAccList2.txt
    mv **/*.sra ./
    mkdir new
    mv *.sra new/
    cd new/
    for f in *.sra; do fasterq-dump $f; done
    for file in SRR*; do mv -i "${file}" "${file/_1/_R1}"; done
    for file in SRR*; do mv -i "${file}" "${file/_2/_R2}"; done
    for i in *_pair.fastq ; do mv $i ${i/_pair.fastq/.trim.fastq}; done
    
Repeat previous steps to get these new files to same place as originals with trimmed, mark_dups, add @RG and process with `gatk haplotypecaller`. Once all the new files are in the /data/seqs/trimmed/marked-dup-bam/ dir then proceed to next step.

## 2021-06-15

Next step of analysis;

Joint variant calls were performed using GenotypeGVCFs on a merged gvcf variant file using the option -maxAltAlleles 2.

So create the merged gvcf file. From; https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs

    cd /data/seqs/trimmed/marked-dup-bam/

    gatk CombineGVCFs \
    -R reference.fasta \
    --variant sample1.g.vcf.gz \
    --variant sample2.g.vcf.gz \
    -O cohort.g.vcf.gz

Can I loop through my sample files e.g.

    gatk CombineGVCFs \
    -R ../../../../zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    for i in *.g.vcf.gz; do \
        --variant $i; done \
    -O cohort.g.vcf.gz

Dry run;

    echo "gatk CombineGVCFs \
    -R ../../../../zymo-reference/zymo-dna-fasta\
    /Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    for i in $(ls *.g.vcf.gz); do \
        --variant $i; done \
    -O cohort.g.vcf.gz"

Seems to have worked, try for real with nohup

    nohup gatk CombineGVCFs \
    -R ../../../../zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
        for i in $(ls *.g.vcf.gz); do \
            --variant $i; done \
    -O cohort.g.vcf.gz

Not working. Try putting it into a script /scripts/combineGVCFs.sh

Still not working. Current work around, create sample file list;

    for i in *.g.vcf.gz; do echo --variant $i '\'; done; 

Cut and paste the output to /scripts/combineGVCFs.sh

    nohup ../../../../scripts/combineGVCFs.sh

Success. Saved stderr to /analysis/CombineGVCFs.out


    nohup gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -V cohort.g.vcf.gz \
    -O genotype.vcf.gz \
    --max-alternate-alleles 2

Success. Save stderr to /analysis/genotypeGVCFs.out.

Move these 4 files (cohort.g.vcf.gz, cohort.g.vcf.gz.tbi, genotype_vcf.gz and genotype_vcf.gz.tbi to new dir `combined-vcf`.

## 2021-06-16

(I missed this step, completed the next step and now realising the mistake have come back to complete this step).
    
Next step of analysis from supplementary methods;

"The joint variant calls were subset to include only SNPs."

From; https://vcftools.github.io/man_latest.html

Output a new vcf file from the input vcf file that removes any indel sites.

    vcftools --gzvcf input_file.vcf.gz --remove-indels --recode --recode-INFO-all --out SNPs_only.vcf.gz

My usage;

    cd combined-vcf/
    vcftools --gzvcf genotype_vcf.gz --remove-indels --recode --recode-INFO-all --out genotype_SNPs_only.vcf.gz

vcftools not installed.

    conda install vcftools

Failed to install, try;

    conda install -c bioconda vcftools

Threw up lots of package conflicts and appears not to have installed vcftools.
Clone github repo into /home/graham/Bin/vcftools/

Follow install guide from github page; https://github.com/vcftools/vcftools

Now retry;

    cd combined-vcf/
    vcftools --gzvcf genotype_vcf.gz --remove-indels --recode --recode-INFO-all --out genotype_SNPs_only.vcf.gz

Throws error;

```
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
    --gzvcf genotype_vcf.gz
    --recode-INFO-all
    --out genotype_SNPs_only.vcf.gz
    --recode
    --remove-indels

stat error: No such file or directory
Error: Can't determine file type of genotype_vcf.gz
```

Mis-typed filename as genotype_vcf.gz instead of genotype.vcf.gz

Retry;

    vcftools --gzvcf genotype.vcf.gz --remove-indels --recode --recode-INFO-all --out genotype_SNPs_only.vcf.gz

Success. Stderr saved to genotype_SNPs_only.vcf.gz.log and moved to /analysis/

I now tried to re-run the variantFiltration.sh on the genotype_SNPs_only.vcf file but it seems to have not worked properly. It did generate the output files but the genotype_SNPs_only-filt.vcf file only has 1 line of data, after listing all the filtering parameters?

This script did work previously when I used a compressed file as input.

For some reason the previous vcftools cmd did not generate the output file as a compressed file.
Do;

    bgzip -c genotype_SNPs_only.vcf > genotype_SNPs_only.vcf.gz
    tabix -p vcf genotype_SNPs_only.vcf.gz

tabix not installed;

    conda install tabix

Failed.

## 2021-06-17

Alternative work around omitting the vcftools step; from https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants

Can select SNPs only with this tool. E.g. usage;

    gatk SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V input.vcf \
    --select-type-to-include SNP \
    -O output.vcf

My usage;

    gatk SelectVariants \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -V genotype.vcf.gz \
    --select-type-to-include SNP \
    -O genotype-selectedSNPs.vcf.gz

Success. 

    wc -l genotype.vcf.gz

2663555 genotype.vcf.gz

    wc -l genotype-selectedSNPs.vcf.gz

2280927 genotype-selectedSNPs.vcf.gz

Now edit and run `variantFiltration.sh`.
Success, but the resulting file is very small.

    wc -l genotype-selectedSNPs-filt.vcf.gz

16 genotype-selectedSNPs-filt.vcf.gz

From the papers vcf file, I am aiming for 

    wc -l GWAS_SNP_dataset_106_Zymoseptoria_tritici_isolates.vcf

**779249** GWAS_SNP_dataset_106_Zymoseptoria_tritici_isolates.vcf

Edited variantFiltration.sh again to reflect the paramaters used in /ena/GWAS_SNP_dataset_106_Zymoseptoria_tritici_isolates.vcf. Then run;

    gatk SelectVariants \
    -R /media/Data/gt293/zymo-reference/zymo-dna-fasta/Zymoseptoria_tritici.MG2.dna.toplevel.fa \
    -V genotype-varFilt-one.vcf.gz \
    --select-type-to-include SNP \
    -O genotype-varFilt-one-SNPs.vcf.gz

    wc -l genotype-varFilt-one-SNPs.vcf.gz
    2268775 genotype-varFilt-one-SNPs.vcf.gz

---------------------------------------------------------------------------------

Next step of analysis from supplementary methods; "SNPs were hard filtered using the GATK VariantFiltration and SelectVariants tools. We chose quality cutoffs following the GATK Best Practices recommendations (DePristo et al., 2011; Van der Auwera et al., 2013) and using empirical evaluations of filter value distributions in our own joint variant file. The final filtering was performed using the filtering cutoffs: QUAL < 250; QD < 20.0; MQ < 30.0; -2 > BaseQRankSum > 2; -2 > MQRankSum > 2; -2 > ReadPosRankSum > 2; FS > 0.1."

    gatk VariantFiltration

example usage from; https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration

    gatk VariantFiltration \
        -R reference.fasta \
        -V input.vcf.gz \
        -O output.vcf.gz \
        --filterExpression "AB < 0.2 || MQ0 > 50" \
        --filterName "my_filters"

Also see; https://gatk.broadinstitute.org/hc/en-us/articles/360035891011?id=1255

My usage;

    nohup ../../../../scripts/variantFiltration.sh 

Failed. Did not produce output file genotype_filt.vcf.gz. Stderr saved to /analysis/variantFiltration-error.out

I hadn't included the `--filter-name`. Have now added 1 name per filter listed as per filter names used in original analysis.
Also had to reformat entries like `-2 > BaseQRankSum > 2` to `BaseQRankSum > 2` & `BaseQRankSum < -2`.

Success. Moved stderr to /analysis/variantFiltration.out

--------------------------------------------------------------------------

Next step from supplementary methods;

"After hard-filtering, we excluded SNPs with a genotyping rate <90% and a minor allele frequency (MAF) <5%. We also excluded SNPs on accessory chromosomes (i.e. chromosomes not found in all isolates of the species). These filtration steps were performed with vcftools version 0.1.12b (Danecek et al., 2011) and plink 1.9 version (https://www.cog- genomics.org/plink2; (Chang et al., 2015). Finally, we converted tri-allelic SNPs to bi-allelic SNPs by recoding the least frequent allele as a missing genotype. Conversion of tri-allelic SNPs was performed to satisfy requirements of association mapping software."


-------------------------------------------------------------------------------------------------------

Generate .pdf from this notebook;

    pandoc --from markdown --to latex notebook.md --output output.pdf

-------------------------------------------------------------------------------------------
Notes;

## Tajima's D statistic

The purpose of Tajima's test is to identify sequences which do not fit the neutral theory model at equilibrium between mutation and genetic drift.

S = number of variable sites.
Tajima's D = normalised version of pairwise difference.
$\pi$ = average pairwise difference (observed).
$\theta$ = expectation of $\pi$ (expected).

If Tajima's D is negative, the inference is the population experience selection reducing the variation, or experienced a recent population expansion.

If Tajima's D is positive, the inference is the population is under selection maintaining the variation, or experienced a recent contraction.