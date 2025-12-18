[back to course home page ](../README.md)

## Variant Calling Workflow
This page contains a short and adapted version of the Software Carpentry lesson [Variant Calling Workflow](https://datacarpentry.github.io/wrangling-genomics/04-variant_calling.html). We will use it as notes to key concepts we will discuss during our lesson

This workflow demonstrates how to identify genetic variants from aligned sequencing data. Each step transforms the data from raw alignments to filtered variant calls ready for downstream analysis. The pipeline uses standard tools like BWA for alignment, SAMtools for processing, and bcftools for variant calling.

![variant-calling](../img/variant-calling-workflow.png)

```diff
+ We will continue learning bash commands while performing variant calling: samtools, bcftools, grep, awk...

```

### **1. Understanding SAM/BAM Format**

[SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map) format stores read alignments to a reference genome:

```
@HD     VN:1.6  SO:coordinate
@SQ     SN:CP000819.1   LN:4558953
SRR2584863.1    83      CP000819.1      9972    60      151M    =       9972    0       GAACGCGAGCAGCGCAGCAGCAGGAGGAACAGGAGGATGGGAATTCCGCGGTGCGGGATGGGA
```
```diff
+ Column 1: Read name
+ Column 2: FLAG (alignment information encoded as number)
+ Column 3: Reference sequence name  
+ Column 4: Position (1-based)
+ Column 5: Mapping quality
+ Column 6: CIGAR string (alignment details)

+ BAM is the binary (compressed) version of SAM
```

```diff
! EXERCISE: Navigate to your seq-analysis directory and create a new folder called 'variant-calling'
! Move into this new directory - this is where we'll perform our analysis
```

```diff
+ COMMAND: samtools - tools for alignments in SAM format
+ COMMAND: bcftools - tools for variant calling and VCF/BCF files
```

### **2. Alignment to Reference Genome**

Before calling variants, we need to align our reads to a reference genome. 

##### **Indexing the Reference**
The reference genome must be indexed before alignment:

```bash
user1@vm-corso-colonna:~/seq-analysis/variant-calling$ bwa index ../data/ref_genome/ecoli_rel606.fasta
```

```diff
! Why do we index? Indexing creates a searchable data structure that speeds up alignment
```

##### **Performing Alignment**
Align paired-end reads using BWA-MEM:

```bash
bwa mem ../data/ref_genome/ecoli_rel606.fasta \
        ../trimmed_fastq/SRR2584863_1.trim.sub.fastq \
        ../trimmed_fastq/SRR2584863_2.trim.sub.fastq > results/sam/SRR2584863.aligned.sam
```

```diff
+ BWA MEM is optimized for reads 70bp-1Mbp
+ Outputs SAM format to stdout (we redirect to file with >)
+ Uses paired-end mode when given two input files
```

### **3. SAM to BAM Conversion**

Convert SAM to BAM format for efficiency:

```bash
samtools view -S -b results/sam/SRR2584863.aligned.sam > results/bam/SRR2584863.aligned.bam
```

```diff
+ -S: input is SAM format
+ -b: output as BAM format
! BAM files are ~4x smaller than SAM files
```

### **4. Sorting and Indexing BAM Files**

##### **Why Sort BAM Files?**
- Many tools require coordinate-sorted BAM files
- Improves performance for variant calling
- Required for BAM indexing

```bash
samtools sort -o results/bam/SRR2584863.aligned.sorted.bam \
              results/bam/SRR2584863.aligned.bam
```

##### **Creating BAM Index**
```bash
samtools index results/bam/SRR2584863.aligned.sorted.bam
```

```diff
! This creates a .bai file that allows random access to the BAM file
! Always required before variant calling
```

### **5. Variant Calling with bcftools**

##### **Calculate Genotype Likelihoods**
```bash
bcftools mpileup -O b -o results/bcf/SRR2584863_raw.bcf \
                 -f ../data/ref_genome/ecoli_rel606.fasta \
                 results/bam/SRR2584863.aligned.sorted.bam
```

```diff
+ -O b: output BCF format (binary VCF)
+ -f: reference genome FASTA file
+ mpileup: calculates genotype likelihoods at each position
```

##### **Call Variants**
```bash
bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584863_variants.vcf \
              results/bcf/SRR2584863_raw.bcf
```

```diff
+ --ploidy 1: haploid genome (bacteria)
+ -m: multiallelic caller
+ -v: output variant sites only
! For diploid organisms (like humans), use --ploidy 2
```

### **6. Understanding VCF Format**

[VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (Variant Call Format) stores genetic variants:

```
##fileformat=VCFv4.2
##reference=ecoli_rel606.fasta
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR2584863
CP000819.1      9972    .       T       G       225     .       DP=10   GT:PL   1:255,0
```

```diff
+ CHROM: Chromosome/contig name
+ POS: Position (1-based)
+ REF: Reference allele
+ ALT: Alternative allele(s)
+ QUAL: Phred-scaled quality score
+ INFO: Additional information (DP=depth)
+ FORMAT/Sample: Genotype information
```

### **7. Variant Filtering**

##### **Filter by Quality**
```bash
bcftools filter -s LowQual -e 'QUAL<20' results/vcf/SRR2584863_variants.vcf > results/vcf/SRR2584863_filtered.vcf
```

```diff
+ -s: add FILTER column annotation
+ -e: exclude sites matching expression
! Q20 = 99% confidence in variant call
```

##### **Filter by Read Depth**
```bash
bcftools filter -s LowDepth -e 'DP<10' -o results/vcf/SRR2584863_final.vcf \
                results/vcf/SRR2584863_filtered.vcf
```

```diff
! Low coverage sites are more prone to false calls
! Typical minimum depth: 10x-30x depending on application
```

### **8. Exploring Variant Statistics**

##### **Count Variants**
```bash
grep -v "^#" results/vcf/SRR2584863_final.vcf | wc -l
```

##### **Extract Variant Types**
```bash
bcftools stats results/vcf/SRR2584863_final.vcf | grep "TSTV"
```

```diff
+ Transition/Transversion ratio (Ts/Tv) is a quality metric
+ Expected Ts/Tv: ~2.0-2.1 for whole genome, ~3.0 for exomes
```

### **9. Automation with Shell Scripts**

Create a script to process multiple samples:

```bash
#!/bin/bash
set -e

# Variables
genome="../data/ref_genome/ecoli_rel606.fasta"
trimmed_dir="../trimmed_fastq"
results_dir="results"

# Create output directories
mkdir -p ${results_dir}/{sam,bam,bcf,vcf}

# Process each sample
for sample in ${trimmed_dir}/*_1.trim.sub.fastq
do
    base=$(basename ${sample} _1.trim.sub.fastq)
    
    echo "Processing sample: ${base}"
    
    # Alignment
    bwa mem ${genome} ${trimmed_dir}/${base}_1.trim.sub.fastq \
                      ${trimmed_dir}/${base}_2.trim.sub.fastq \
                      > ${results_dir}/sam/${base}.aligned.sam
    
    # Convert to BAM and sort
    samtools view -S -b ${results_dir}/sam/${base}.aligned.sam | \
    samtools sort -o ${results_dir}/bam/${base}.aligned.sorted.bam
    
    # Index BAM
    samtools index ${results_dir}/bam/${base}.aligned.sorted.bam
    
    # Variant calling
    bcftools mpileup -O b -o ${results_dir}/bcf/${base}_raw.bcf \
                     -f ${genome} ${results_dir}/bam/${base}.aligned.sorted.bam
    
    bcftools call --ploidy 1 -m -v -o ${results_dir}/vcf/${base}_variants.vcf \
                  ${results_dir}/bcf/${base}_raw.bcf
    
    # Filter variants
    bcftools filter -s LowQual -e 'QUAL<20 || DP<10' \
                    -o ${results_dir}/vcf/${base}_final.vcf \
                    ${results_dir}/vcf/${base}_variants.vcf
done
```

```diff
! EXERCISE: Save this script as 'variant_calling.sh' and make it executable with chmod +x
! Run it with: ./variant_calling.sh
```

### **Key Quality Considerations**

```diff
+ Mapping Quality (MAPQ): Confidence in read placement (0-60 scale)
+ Base Quality: Confidence in individual base calls
+ Variant Quality (QUAL): Confidence in variant call
+ Read Depth (DP): Number of reads supporting position
! Always check multiple metrics before trusting variants
```

### **Visualizing Results**

Transfer VCF files to your local computer for viewing:

```diff
! On your local computer:
scp user1@212.189.205.193:/home/user1/seq-analysis/variant-calling/results/vcf/*_final.vcf .

+ VCF files can be viewed in:
+ - Text editors (small files)
+ - IGV (Integrative Genomics Viewer)
+ - Excel/spreadsheet programs (after conversion)
```

### **Summary Commands Learned**

```diff
+ bwa index - Index reference genome
+ bwa mem - Align reads to reference
+ samtools view - Convert between SAM/BAM
+ samtools sort - Sort alignments
+ samtools index - Index BAM files
+ bcftools mpileup - Calculate genotype likelihoods
+ bcftools call - Call variants
+ bcftools filter - Filter variants
+ bcftools stats - Generate statistics
```