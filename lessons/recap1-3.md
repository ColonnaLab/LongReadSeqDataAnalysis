[back to course home page ](../README.md)

# Day 4 - Recap and Consolidation Practical

This page contains a guided recap of Days 1, 2, and 3. We will use it to consolidate the practical skills needed to run a short-read bioinformatics workflow from remote login to filtered variant calls.

The lesson is planned as a **4-hour practical session**.

```diff
+ Main tools: ssh, bash, FastQC, Trimmomatic, bwa/bowtie2, samtools, bcftools
+ Main formats: FASTA, FASTQ, SAM/BAM, VCF
+ Main workflow: QC -> trimming -> alignment -> variant calling -> filtering
```

## Learning objectives

By the end of this recap, you should be able to:

- connect to a remote machine with `ssh`
- navigate the Unix filesystem and organize an analysis directory
- recognize FASTA, FASTQ, SAM/BAM, and VCF files
- run and interpret basic quality control with `FastQC`
- trim paired-end reads with `Trimmomatic`
- align reads to a reference genome with `bwa` or `bowtie2`
- convert, sort, index, and inspect alignments with `samtools`
- produce and filter variant calls with `bcftools`
- explain the input and output of each step in the short-read workflow

---

## 1. Recap: connection to the remote machine

We start from the local computer and connect to the remote server.

```bash
ssh user1@212.189.205.193
```

After login, commands are executed on the remote machine.

```diff
! EXERCISE:
! 1. Connect to the remote machine with ssh
! 2. Check where you are with pwd
! 3. List the files in your home directory
! 4. Move to your analysis directory
```

Useful commands:

```bash
pwd
ls
ls -lh
cd
cd ~
mkdir
mkdir -p
```

```diff
+ Remember:
+ / is the root of the filesystem
+ ~ is your home directory
+ . means current directory
+ .. means parent directory
```

## 2. Recap: shell syntax and good working habits

A shell command usually has this structure:

```bash
command -option argument
```

Examples:

```bash
ls -lh
head -n 4 file.fastq
samtools view -h alignment.bam
```

Good habits for bioinformatics analyses:

- keep raw data unchanged
- make a dedicated analysis directory
- use clear names for directories and files
- save new output files in separate folders
- write down commands and parameters

```diff
! EXERCISE:
! Create or check the following directory structure inside ~/seq-analysis
```

```bash
cd ~
mkdir -p seq-analysis
cd seq-analysis
mkdir -p fastqc-res trimmed refgenome recap-results/sam recap-results/bam recap-results/bcf recap-results/vcf
ls -R
```

---

## 3. Recap: main sequence file formats

### Where is the data 
We have already downloaded the data, and it lives here (make sure you type the right user number ...  ): 

```bash 
/data/user1/genomic-lesson-data
```

To make the commands shorter, create a symbolic link to the data inside the directory where you are working:

```bash
cd ~/seq-analysis
ln -s /data/user1/genomic-lesson-data genomic-lesson-data
ls -lh
```

```diff
+ COMMAND ln -s creates a symbolic link
+ The data stay in /data/user1/genomic-lesson-data
+ The link genomic-lesson-data lets us access them from ~/seq-analysis
```

```diff
+ genomic-lesson-data/results contains example results that were already prepared
+ recap-results will contain the files that you produce during this recap
```

### FASTA

FASTA stores one or more sequences. It is often used for reference genomes.

```text
>CP000819.1 Escherichia coli B str. REL606
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGT
```

```diff
+ Line beginning with >: sequence identifier
+ Following lines: nucleotide or protein sequence
```

Inspect a FASTA file:

```bash
head genomic-lesson-data/refgenome/ecoli_rel606.fasta
grep ">" genomic-lesson-data/refgenome/ecoli_rel606.fasta
```

### FASTQ

FASTQ stores sequencing reads and their base quality scores. Each read has four lines.

```text
@SRR2584863.1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTT
+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIF
```

```diff
+ Line 1: read identifier
+ Line 2: sequence
+ Line 3: separator
+ Line 4: quality scores
```

Inspect a compressed FASTQ file:

```bash
zcat genomic-lesson-data/untrimmed_fastq/SRR2584863_1.fastq.gz | head -n 4
```

### SAM/BAM

SAM and BAM store alignments of reads against a reference genome.

- SAM is text
- BAM is the compressed binary version
- sorted and indexed BAM files are used by many downstream tools

Inspect a BAM file:

```bash
samtools view genomic-lesson-data/results/bam/SRR2584863.aligned.sorted.bam | head
samtools view genomic-lesson-data/results/bam/SRR2584863.aligned.sorted.bam | less -S
samtools flagstat genomic-lesson-data/results/bam/SRR2584863.aligned.sorted.bam
```

### VCF

VCF stores genetic variants.

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
CP000819.1  9972  .  T  G  38.4  PASS  DP=24
```

Important fields:

- `CHROM`: reference sequence
- `POS`: position
- `REF`: reference allele
- `ALT`: alternative allele
- `QUAL`: variant quality
- `FILTER`: filtering status

Inspect a VCF file:

```bash
grep -v "^##" genomic-lesson-data/results/vcf/myfirst.filtered.vcf | head
```

```diff
! EXERCISE:
! For each file type, answer:
! 1. Is it raw data, reference data, alignment data, or variant data?
! 2. Is it human-readable text or binary/compressed?
! 3. Which tool produced or used it?
```

---

## 4. Recap: short-read workflow overview

```diff
+ Raw FASTQ
+   |
+   | FastQC
+   v
+ Quality report
+
+ Raw FASTQ
+   |
+   | Trimmomatic
+   v
+ Trimmed FASTQ
+   |
+   | bwa or bowtie2
+   v
+ SAM
+   |
+   | samtools view/sort/index
+   v
+ sorted BAM + BAI
+   |
+   | bcftools mpileup/call
+   v
+ VCF
+   |
+   | bcftools filter
+   v
+ filtered VCF
```

Questions to keep asking during the workflow:

- What is the input?
- What is the output?
- Which file format is used?
- Which software is responsible for the transformation?
- What quality check tells us if the step worked?

---

## 5. Guided consolidation exercise

In this exercise, we repeat the complete short-read workflow using the same organization introduced in Days 1-3.

```diff
+ In Day 3 we used freebayes for variant calling and bcftools for filtering.
+ In this recap we use bcftools for both variant calling and filtering.
```

Start in your analysis directory:

```bash
cd /data/user1/seq-analysis
pwd
ls
```

### Step 1 - Check the raw reads

```bash
ls genomic-lesson-data/untrimmed_fastq/
zcat genomic-lesson-data/untrimmed_fastq/SRR2584863_1.fastq.gz | head -n 4
zcat genomic-lesson-data/untrimmed_fastq/SRR2584863_2.fastq.gz | head -n 4
```

```diff
! EXERCISE:
! What tells you these files are paired-end FASTQ reads?
! What happens to the quality scores near the end of the read?
```

### Step 2 - Run FastQC on raw reads

```bash
mkdir -p fastqc-res
fastqc genomic-lesson-data/untrimmed_fastq/SRR2584863_1.fastq.gz -o fastqc-res
fastqc genomic-lesson-data/untrimmed_fastq/SRR2584863_2.fastq.gz -o fastqc-res
ls fastqc-res
```

Transfer one HTML report to your laptop with `scp` and open it in a browser:

```bash
scp user1@212.189.205.193:/home/user1/seq-analysis/fastqc-res/SRR2584863_1_fastqc.html .
```

```diff
! EXERCISE:
! Identify at least one FastQC module that passes, warns, or fails.
! Why do we inspect raw reads before trimming?
```

### Step 3 - Trim paired-end reads with Trimmomatic

```bash
mkdir -p trimmed
cd trimmed
```

```bash
java -jar /data/Trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -threads 4 \
  ../genomic-lesson-data/untrimmed_fastq/SRR2584863_1.fastq.gz \
  ../genomic-lesson-data/untrimmed_fastq/SRR2584863_2.fastq.gz \
  trimmed/SRR2584863_1.trimmed.fastq trimmed/SRR2584863_1.un.trimmed.fastq \
  trimmed/SRR2584863_2.trimmed.fastq trimmed/SRR2584863_2.un.trimmed.fastq \
  ILLUMINACLIP:/data/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:20
```

Check the output:

```bash
ls -lh
wc -l SRR2584863_1.trimmed.fastq
wc -l SRR2584863_2.trimmed.fastq
cd ..
```

```diff
! EXERCISE:
! Which files contain paired reads?
! Which files contain reads where the mate was discarded?
! Why do we keep the paired files for paired-end alignment?
```

Run FastQC again on the trimmed reads:

```bash
fastqc trimmed/SRR2584863_1.trimmed.fastq -o fastqc-res
fastqc trimmed/SRR2584863_2.trimmed.fastq -o fastqc-res
```

```diff
! EXERCISE:
! Compare raw and trimmed FastQC reports.
! Did trimming improve the end-of-read quality?
```

### Step 4 - Prepare the reference genome

```bash
mkdir -p refgenome
cd refgenome
```

If the reference genome is not already present, download and uncompress it:

```bash
curl -L -o ecoli_rel606.fasta.gz \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

gunzip ecoli_rel606.fasta.gz
```

Inspect the FASTA file:

```bash
head ecoli_rel606.fasta
grep ">" ecoli_rel606.fasta
```

Index the reference for `bwa` and create the FASTA index used by `samtools` and `bcftools`:

```bash
bwa index ecoli_rel606.fasta
samtools faidx ecoli_rel606.fasta
ls
cd ..
```

```diff
! EXERCISE:
! Why does the aligner need an index of the reference genome?
! Why do samtools and bcftools need the .fai index?
```

### Step 5 - Align reads to the reference genome

Create the result directories for the files produced during this recap:

```bash
mkdir -p recap-results/sam recap-results/bam recap-results/bcf recap-results/vcf
```

Align the trimmed reads with `bwa mem`:

```bash
bwa mem \
  refgenome/ecoli_rel606.fasta \
  trimmed/SRR2584863_1.trimmed.fastq \
  trimmed/SRR2584863_2.trimmed.fastq \
  > recap-results/sam/SRR2584863.aligned.sam
```

Inspect the SAM file:

```bash
head -n 3 recap-results/sam/SRR2584863.aligned.sam
less -S recap-results/sam/SRR2584863.aligned.sam
```

```diff
! EXERCISE:
! In the SAM file, identify:
! QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQ, QUAL
```

Alternative aligner reminder:

```diff
+ bowtie2 can also be used for short-read alignment.
+ Like bwa, it requires an indexed reference genome before alignment.
```

### Step 6 - Convert SAM to sorted BAM

```bash
samtools view -S -b recap-results/sam/SRR2584863.aligned.sam \
  > recap-results/bam/SRR2584863.aligned.bam

samtools sort \
  -o recap-results/bam/SRR2584863.aligned.sorted.bam \
  recap-results/bam/SRR2584863.aligned.bam

samtools index recap-results/bam/SRR2584863.aligned.sorted.bam
```

Check alignment summary:

```bash
ls -lh recap-results/bam
samtools flagstat recap-results/bam/SRR2584863.aligned.sorted.bam
samtools idxstats recap-results/bam/SRR2584863.aligned.sorted.bam
```

Visualize the alignment:

```bash
samtools tview recap-results/bam/SRR2584863.aligned.sorted.bam refgenome/ecoli_rel606.fasta
```

```diff
! EXERCISE:
! Why do we sort and index BAM files before variant calling?
```

### Step 7 - Variant calling with bcftools

Create genotype likelihoods in BCF format:

```bash
bcftools mpileup \
  -Ou \
  -f refgenome/ecoli_rel606.fasta \
  recap-results/bam/SRR2584863.aligned.sorted.bam \
  > recap-results/bcf/SRR2584863.bcf
```

Call variants:

```bash
bcftools call \
  -mv \
  -Ov \
  recap-results/bcf/SRR2584863.bcf \
  > recap-results/vcf/SRR2584863.vcf
```

Inspect the VCF:

```bash
grep -v "^##" recap-results/vcf/SRR2584863.vcf | head
grep -v "^#" recap-results/vcf/SRR2584863.vcf | wc -l
```

Filter low-quality variants:

```bash
bcftools filter \
  -i 'QUAL>20' \
  recap-results/vcf/SRR2584863.vcf \
  > recap-results/vcf/SRR2584863.filtered.vcf
```

Compare before and after filtering:

```bash
grep -v "^#" recap-results/vcf/SRR2584863.vcf | wc -l
grep -v "^#" recap-results/vcf/SRR2584863.filtered.vcf | wc -l
grep -v "^##" recap-results/vcf/SRR2584863.filtered.vcf | head
```

```diff
! EXERCISE:
! How many variants were present before filtering?
! How many variants remain after filtering?
! Which VCF column did we use for filtering?
```

---

## 6. Final consolidation checklist

At the end of the practical, your `seq-analysis` directory should contain:

```text
seq-analysis/
|-- fastqc-res/
|-- trimmed/
|-- refgenome/
`-- recap-results/
    |-- sam/
    |-- bam/
    |-- bcf/
    `-- vcf/
```

Main expected files:

- raw and trimmed FastQC reports
- paired trimmed FASTQ files
- indexed reference FASTA
- SAM alignment file
- sorted and indexed BAM file
- unfiltered and filtered VCF files

```diff
! FINAL EXERCISE:
! Explain the complete workflow using file formats:
! FASTQ -> FASTQ -> SAM -> BAM -> VCF
!
! Then explain the same workflow using tools:
! FastQC/Trimmomatic -> bwa or bowtie2 -> samtools -> bcftools
```

## 7. Short glossary

- `ssh`: connect to a remote machine
- `bash`: Unix shell used to run commands
- `FastQC`: quality control of sequencing reads
- `Trimmomatic`: adapter and quality trimming for Illumina reads
- `bwa`: short-read aligner
- `bowtie2`: short-read aligner
- `samtools`: tools for SAM/BAM files
- `bcftools`: tools for BCF/VCF files and variant calling/filtering
- `FASTA`: reference or assembled sequence format
- `FASTQ`: reads plus base quality scores
- `SAM`: text alignment format
- `BAM`: compressed binary alignment format
- `VCF`: variant call format
