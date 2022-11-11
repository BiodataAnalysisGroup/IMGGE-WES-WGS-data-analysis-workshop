### 0. Linux basic commands ###


### 1. Introduction on Variant calling ###

# Different types of sequencing for different purposes. For example: 
# - Targeted Sequencing (e.g., COVID-19 and lineagespot?)
# - Whole-Exome Sequencing (WES)
# - Whole-Genome Sequencing (WGS)

# Goal of variant calling, i.e. identification of variants that are associated with the phenotype under study (e.g., genetic disorders, cancer)

### 2. FastQC - Quality control ###

# Overview of FastQ format

# The need for quality control of raw data, following sequencing

# Running FastQC
fastqc --version
fastqc --help
fastqc *.fq.gz

# 3. Trim galore / Trimmomatic - Trimming for low-quality, adapter and duplicated (e.g., optical or PCR duplicates, PCR bias) sequences
# Trim_galore
trim_galore --version
trim_galore --help
mkdir trim_galore_data
# set the option --path_to_cutadapt to the Cutadapt executable if Cutadapt is not in the PATH (default)
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired father_R1.fq.gz father_R2.fq.gz
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired mother_R1.fq.gz mother_R1.fq.gz
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired proband_R1.fq.gz proband_R2.fq.gz
# OR in a single command
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired \
father_R1.fq.gz father_R2.fq.gz \
mother_R1.fq.gz mother_R1.fq.gz \
proband_R1.fq.gz proband_R2.fq.gz

# Trimmomatic
mkdir trimmomatic_data
# paired-end mode
#
java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 father_R1.fq.gz father_R2.fq.gz \
trimmomatic_data/father_R1_paired.fq.gz trimmomatic_data/father_R1_unpaired.fq.gz \
trimmomatic_data/father_R2_paired.fq.gz trimmomatic_data/father_R2_unpaired.fq.gz \
ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
SLIDINGWINDOW:4:15 LEADING:5 TRAILING:3 MINLEN:36
#
java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 mother_R1.fq.gz mother_R2.fq.gz \
trimmomatic_data/mother_R1_paired.fq.gz trimmomatic_data/mother_R1_unpaired.fq.gz \
trimmomatic_data/mother_R2_paired.fq.gz trimmomatic_data/mother_R2_unpaired.fq.gz \
ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
SLIDINGWINDOW:4:15 LEADING:5 TRAILING:3 MINLEN:36
#
java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 proband_R1.fq.gz proband_R2.fq.gz \
trimmomatic_data/proband_R1_paired.fq.gz trimmomatic_data/proband_R1_unpaired.fq.gz \
trimmomatic_data/proband_R2_paired.fq.gz trimmomatic_data/proband_R2_unpaired.fq.gz \
ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
SLIDINGWINDOW:4:15 LEADING:5 TRAILING:3 MINLEN:36

### Mapping reads to reference genome - BWA-MEM ###

# Create index for reference genome
bwa index hg19_chr8.fa

# Map pre-processed reads with BWA-MEM
# Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments. (-a)
# Mark shorter split hits as secondary (for Picard compatibility). (-M)
bwa mem -t 16 -a -M hg19_chr8.fa trimmomatic_data/father_R1_paired.fq.gz trimmomatic_data/father_R2_paired.fq.gz > father.sam
bwa mem -t 16 -a -M hg19_chr8.fa trimmomatic_data/mother_R1_paired.fq.gz trimmomatic_data/mother_R2_paired.fq.gz > mother.sam
bwa mem -t 16 -a -M hg19_chr8.fa trimmomatic_data/proband_R1_paired.fq.gz trimmomatic_data/proband_R2_paired.fq.gz > proband.sam

### Pre-processing of mapped reads before variant calling ###

# Explain SAM & BAM format

# Convert SAM-to-BAM file
samtools view --threads 16 -b father.sam > father.bam
samtools view --threads 16 -b mother.sam > mother.bam
samtools view --threads 16 -b proband.sam > proband.bam

# Sort reads based on genomic coordinates
samtools sort --threads 16 -o father_sorted.bam father.bam
samtools sort --threads 16 -o mother_sorted.bam mother.bam
samtools sort --threads 16 -o proband_sorted.bam proband.bam

# Generate mapping statistics on BAM files
# 1. FASTQC on sorted BAM files
fastqc *_sorted.bam

# 2. samtools flagstat
samtools flagstat father_sorted.bam > father_sorted.flagstat
samtools flagstat mother_sorted.bam > mother_sorted.flagstat
samtools flagstat proband_sorted.bam > proband_sorted.flagstat
# OR (in case of large number of samples, use loops)
for i in $(ls *_sorted.bam)
do
    # extract sample name (remove suffix)
    sample=${i%_sorted*}
    # run command
    samtools flagstat $i > ${sample}.flagstat
done

# 3. samtools stats & plot-bamstat
for i in $(ls *_sorted.bam)
do
    # extract sample name (remove suffix)
    sample=${i%_sorted*} # keep the *_sorted* part of the filename
    # run samtools stat
    samtools stats $i > ${sample}.bc
    # extract summary statistics
    cat ${sample}.bc | grep ^SN | cut -f 2- > ${sample}.bc.summary
    # produce plot-bamstat plot
    plot-bamstats -p ./bamstat_plots/${sample} ${sample}.bc
done

# 4. custom script (e.g., in R) OR ready-to-use software solution (e.g., BAMStats (https://bamstats.sourceforge.net/))
Rscript summary_plot_generation.R

# Filter the paired-end reads of all samples to retain only those read pairs, 
# for which both the forward and the reverse read have been mapped to the reference successfully
# https://broadinstitute.github.io/picard/explain-flags.html

# exclude flags with -F
# require flags with -f
# read mapped in proper pair (0x2) - 2
# not primary alignment (0x100) - 256
samtools view --threads 16 -b -f 2 -o father_filtered.bam father_sorted.bam 
samtools view --threads 16 -b -f 2 -o mother_filtered.bam mother_sorted.bam
samtools view --threads 16 -b -f 2 -o proband_filtered.bam proband_sorted.bam

# Set read groups information (Picard)
java -jar /usr/local/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups \
    I=father_filtered.bam \
    O=father_aor.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=father
#
java -jar /usr/local/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups \
    I=mother_filtered.bam \
    O=mother_aor.bam \
    RGID=2 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=mother
#
java -jar /usr/local/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups \
    I=proband_filtered.bam \
    O=proband_aor.bam \
    RGID=3 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=proband

# MarkDuplicates (Picard)
java -jar /usr/local/picard-tools-2.5.0/picard.jar MarkDuplicates \
      I=father_aor.bam \
      O=father_md.bam \
      M=father_md_metrics.txt
#
java -jar /usr/local/picard-tools-2.5.0/picard.jar MarkDuplicates \
      I=mother_aor.bam \
      O=mother_md.bam \
      M=mother_md_metrics.txt
#
java -jar /usr/local/picard-tools-2.5.0/picard.jar MarkDuplicates \
      I=proband_aor.bam \
      O=proband_md.bam \
      M=proband_md_metrics.txt

# Generate index (.fai) for the reference fasta
samtools faidx hg19_chr8.fa

### Variant calling and filtering - FreeBayes/GATK HaplotypeCaller/bcftools ###
# default settings (for low to high depth sequencing in haploid and diploid samples)
freebayes -f hg19_chr8.fa -C 2 -F 0.2 father_md.bam > father.vcf
freebayes -f hg19_chr8.fa -C 2 -F 0.2 mother_md.bam > mother.vcf
freebayes -f hg19_chr8.fa -C 2 -F 0.2 proband_md.bam > proband.vcf

# Variant filtering (bcftools) and index generation (tabix)
bcftools filter -s QUAL100 -e '%QUAL<100' -Oz father.vcf | bcftools filter -s DP5 -e 'DP<5' -m+ -Oz -o father_filtered.vcf.gz
tabix -p vcf father_filtered.vcf.gz
bcftools filter -s QUAL100 -e '%QUAL<100' -Oz mother.vcf | bcftools filter -s DP5 -e 'DP<5' -m+ -Oz -o mother_filtered.vcf.gz
tabix -p vcf mother_filtered.vcf.gz
bcftools filter -s QUAL100 -e '%QUAL<100' -Oz proband.vcf | bcftools filter -s DP5 -e 'DP<5' -m+ -Oz -o proband_filtered.vcf.gz
tabix -p vcf proband_filtered.vcf.gz

# Create list of VCF files
ls *_filtered.vcf.gz > tmp_vcf_file_list.txt

# Merge VCF 
bcftools merge --threads 16 --filter-logic + --missing-to-ref -Oz --file-list tmp_vcf_file_list.txt | bcftools view --threads 16 -f PASS -Ov -o merged.vcf

# Post-processing
bcftools norm --threads 16 -f hg19_chr8.fa --multiallelics -both -Ov -o merged_norm.vcf merged.vcf

# Variant annotation - SnpEff & SnpSift
# check available databases for Homo sapiens
java -jar snpEff/snpEff.jar databases | grep -i "Homo_sapiens"
# Build UCSC hg19 database
java -jar snpEff/snpEff.jar download -v hg19
# Annotate with snpEff hg19 genome and dbSNP138 - hg19 - chr8
java -jar snpEff/snpEff.jar ann -v -c snpEff/snpEff.config -noStats hg19 merged_norm.vcf | \
java -jar snpEff/SnpSift.jar Annotate -v -id dbsnp_138.hg19.chr8.vcf > merged_norm_anno.vcf

# GEMINI
# Run Galaxy Europe tools: GEMINI load and GEMINI inheritance pattern with the tutorial-suggested parameters
