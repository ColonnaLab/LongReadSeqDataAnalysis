[back to course home page ](../README.md)

## Long reads variant calling exercise

In this practical we will use the aligned Oxford Nanopore reads from MRSA strain KUN1163 to call and inspect variants relative to the reference chromosome and plasmid.

The input files are:

- reference genome: `../data-longreads/reference/KUN1163_reference.fasta`
- aligned reads: `bam/DRR187567.KUN1163.sorted.bam`
- BAM index: `bam/DRR187567.KUN1163.sorted.bam.bai`

```diff
+ We will work inside ~/lr-working
+ We will use the BAM file produced in lesson 5
+ We will create VCF files for variant calls
+ We will inspect and filter the variant calls
```

### **1. Prepare the Working Directory**

Move to the long reads working folder and check that the alignment files exist:

```bash
user1@vm-corso-colonna:~$ cd ~/lr-working
user1@vm-corso-colonna:~/lr-working$ ls -lh bam/
user1@vm-corso-colonna:~/lr-working$ ls -lh ../data-longreads/reference/
```

Create folders for variant calling results:

```bash
user1@vm-corso-colonna:~/lr-working$ mkdir -p variants variants/logs
```

```diff
! EXERCISE: Confirm that the sorted BAM and BAI files are present
! EXERCISE: Confirm that the reference FASTA and FAI files are present
```

### **2. Check the BAM and Reference Indexes**

Variant calling tools need indexed input files.

If the BAM index is missing, create it:

```bash
user1@vm-corso-colonna:~/lr-working$ samtools index bam/DRR187567.KUN1163.sorted.bam
```

If the reference FASTA index is missing, create it:

```bash
user1@vm-corso-colonna:~/lr-working$ samtools faidx ../data-longreads/reference/KUN1163_reference.fasta
```

Check the indexed files:

```bash
user1@vm-corso-colonna:~/lr-working$ ls -lh bam/DRR187567.KUN1163.sorted.bam*
user1@vm-corso-colonna:~/lr-working$ ls -lh ../data-longreads/reference/KUN1163_reference.fasta*
```

### **3. Inspect Coverage Before Variant Calling**

Coverage affects variant calling reliability. Positions with very low coverage are more difficult to call confidently.

```bash
user1@vm-corso-colonna:~/lr-working$ samtools coverage \
  bam/DRR187567.KUN1163.sorted.bam \
  > variants/DRR187567.coverage.before_variant_calling.txt
```

```bash
user1@vm-corso-colonna:~/lr-working$ cat variants/DRR187567.coverage.before_variant_calling.txt
```

```diff
! EXERCISE: Is the chromosome coverage high enough for variant calling?
! EXERCISE: Is the plasmid coverage similar to the chromosome coverage?
```

### **4. Small Variant Calling with FreeBayes**

`freebayes` can call small variants such as SNVs and small indels from aligned reads.

Check the help page:

```bash
user1@vm-corso-colonna:~/lr-working$ freebayes -h | less
```

Run variant calling:

```bash
user1@vm-corso-colonna:~/lr-working$ freebayes \
  -f ../data-longreads/reference/KUN1163_reference.fasta \
  bam/DRR187567.KUN1163.sorted.bam \
  > variants/DRR187567.freebayes.vcf
```

Check the output:

```bash
user1@vm-corso-colonna:~/lr-working$ ls -lh variants/DRR187567.freebayes.vcf
```

```diff
+ The VCF file contains variant calls relative to the KUN1163 reference
+ Header lines start with ##
+ The column header line starts with #CHROM
+ Variant rows contain chromosome, position, REF, ALT, QUAL, FILTER, and INFO
```

### **5. Inspect the VCF File**

View the header:

```bash
user1@vm-corso-colonna:~/lr-working$ grep "^##" variants/DRR187567.freebayes.vcf | less
```

View the VCF column names:

```bash
user1@vm-corso-colonna:~/lr-working$ grep "^#CHROM" variants/DRR187567.freebayes.vcf
```

View the first variant records:

```bash
user1@vm-corso-colonna:~/lr-working$ grep -v "^#" variants/DRR187567.freebayes.vcf | head
```

Count the number of variant records:

```bash
user1@vm-corso-colonna:~/lr-working$ grep -vc "^#" variants/DRR187567.freebayes.vcf
```

```diff
! EXERCISE: How many variants were called?
! EXERCISE: Which reference sequence contains most variants?
! EXERCISE: Look at the REF and ALT columns. Do you see SNVs, indels, or both?
```

### **6. Compress and Index the VCF**

Many tools work better with compressed and indexed VCF files.

Compress the VCF:

```bash
user1@vm-corso-colonna:~/lr-working$ bgzip \
  -c variants/DRR187567.freebayes.vcf \
  > variants/DRR187567.freebayes.vcf.gz
```

Index the compressed VCF:

```bash
user1@vm-corso-colonna:~/lr-working$ tabix \
  -p vcf \
  variants/DRR187567.freebayes.vcf.gz
```

Check the output:

```bash
user1@vm-corso-colonna:~/lr-working$ ls -lh variants/DRR187567.freebayes.vcf.gz*
```

### **7. Summarize Variants with bcftools**

Use `bcftools stats` to summarize the VCF:

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools stats \
  variants/DRR187567.freebayes.vcf.gz \
  > variants/DRR187567.freebayes.stats.txt
```

Inspect the summary:

```bash
user1@vm-corso-colonna:~/lr-working$ less variants/DRR187567.freebayes.stats.txt
```

Count variants by reference sequence:

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools query \
  -f '%CHROM\n' \
  variants/DRR187567.freebayes.vcf.gz \
  | sort \
  | uniq -c
```

```diff
! EXERCISE: How many calls are on the chromosome?
! EXERCISE: How many calls are on the plasmid?
```

### **8. Filter Variants**

Filtering removes low-confidence calls or creates a smaller set of variants for inspection.

Filter by quality:

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools filter \
  -i 'QUAL>20' \
  variants/DRR187567.freebayes.vcf.gz \
  -Oz \
  -o variants/DRR187567.freebayes.QUAL20.vcf.gz
```

Index the filtered VCF:

```bash
user1@vm-corso-colonna:~/lr-working$ tabix \
  -p vcf \
  variants/DRR187567.freebayes.QUAL20.vcf.gz
```

Count variants before and after filtering:

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools view \
  -H \
  variants/DRR187567.freebayes.vcf.gz \
  | wc -l
```

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools view \
  -H \
  variants/DRR187567.freebayes.QUAL20.vcf.gz \
  | wc -l
```

```diff
! EXERCISE: How many variants remain after QUAL>20 filtering?
! EXERCISE: Did filtering remove many calls or only a few?
! EXERCISE: Why might low-quality variant calls be unreliable?
```

### **9. Inspect Variants in One Region**

Use `bcftools view` to inspect variants in a specific region.

Example for the chromosome:

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools view \
  -r AP020324.1:1-100000 \
  variants/DRR187567.freebayes.QUAL20.vcf.gz \
  | less
```

Example for the plasmid:

```bash
user1@vm-corso-colonna:~/lr-working$ bcftools view \
  -r AP020325.1:1-30220 \
  variants/DRR187567.freebayes.QUAL20.vcf.gz \
  | less
```

```diff
! EXERCISE: Are variants present on the plasmid?
! EXERCISE: Compare the QUAL values of variants in different regions
```

### **10. Structural Variant Calling Concepts**

Small variant calling is not enough to describe all long-read variation. Long reads are especially useful for structural variants.

Signals that can support structural variants include:

- split alignments
- supplementary alignments
- long insertions in CIGAR strings
- long deletions in CIGAR strings
- changes in coverage
- reads that span both sides of a rearrangement

```diff
! EXERCISE: Go back to lesson 5 section 6b
! How many supplementary alignments did we observe?
! Why are supplementary alignments useful for structural variant discovery?
```

### **11. Save Results for Later**

Files that should now be present:

- `variants/DRR187567.freebayes.vcf`
- `variants/DRR187567.freebayes.vcf.gz`
- `variants/DRR187567.freebayes.vcf.gz.tbi`
- `variants/DRR187567.freebayes.stats.txt`
- `variants/DRR187567.freebayes.QUAL20.vcf.gz`
- `variants/DRR187567.freebayes.QUAL20.vcf.gz.tbi`

Copy selected results to your local machine from the `lr-local-work` folder:

```bash
lr-local-work$ scp user1@212.189.205.193:/home/user1/lr-working/variants/*.txt .
lr-local-work$ scp user1@212.189.205.193:/home/user1/lr-working/variants/*.vcf .
```

```diff
! EXERCISE: Keep the filtered VCF and the variant statistics file
! These are the main outputs of this practical
```

