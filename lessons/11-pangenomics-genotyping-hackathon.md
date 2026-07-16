[back to course home page ](../README.md)

## Advanced 5: Pangenomics genotyping hackathon

This practical is adapted to the course machine and workspace from the [PanGenie workshop](https://pangenie-workshop.readthedocs.io/en/latest/index.html) by Jana Ebler.

Follow the original tutorial for the biological background and interpretation. This page records the course-specific paths and commands.

Dataset and outputs assumed on the server:

- downloaded dataset: `/scratch/user1/data-pangenie`
- working directory: `/scratch/user1/pangenie-working`
- PanGenie index output: `indexes/genotypes`
- genotype output folder: `genotypes`

```diff
+ You should follow the PanGenie workshop tutorial
+ You should keep input data separate from working outputs
+ You should build the PanGenie index before genotyping samples
+ You should inspect both PanGenie bubble genotypes and converted biallelic genotypes
```

### **1. Prepare the Working Directory**

Create the working directory:

```bash
user1@vm-corso-colonna:~$ mkdir -p /scratch/user1/pangenie-working
user1@vm-corso-colonna:~$ cd /scratch/user1/pangenie-working
```

Create a link to the PanGenie data:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ln -s /scratch/user1/data-pangenie data-pangenie
```

Check that the expected files are present:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh data-pangenie
```

Create output folders:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ mkdir -p indexes genotypes
```

```diff
! TASK: Open the PanGenie workshop page
! TASK: Confirm that the course dataset is in data-pangenie
! TASK: Identify the panel VCF, reference FASTA, and sample FASTA files
```

### **2. Inspect the Input Data**

Before running PanGenie, inspect the prepared dataset and connect each file to the workflow described in the tutorial.

List files with sizes:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ find data-pangenie -maxdepth 2 -type f | sort
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh data-pangenie
```

Inspect the reference and sample FASTA files:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ seqkit stats data-pangenie/*.fasta
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ seqkit seq -n data-pangenie/reference_chr5:50200000-50400000.fasta
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ seqkit seq -n data-pangenie/NA19189_chr5:50200000-50400000.fasta
```

If `seqkit` is not available, use basic shell commands:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ grep '^>' data-pangenie/*.fasta
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ head -n 4 data-pangenie/reference_chr5:50200000-50400000.fasta
```

Inspect the multi-allelic panel VCF used for indexing:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ less -S data-pangenie/panel_multi_chr5:50200000-50400000.vcf
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ grep -v '^##' data-pangenie/panel_multi_chr5:50200000-50400000.vcf | head
```

Inspect the biallelic panel VCF used during conversion:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools view -h data-pangenie/panel_bi_chr5:50200000-50400000.vcf.gz | less
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools view data-pangenie/panel_bi_chr5:50200000-50400000.vcf.gz | less -S
```

Inspect the helper script used later for converting bubble genotypes:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh data-pangenie/software
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ python3 data-pangenie/software/convert-to-biallelic.py --help
```

```diff
! TASK: How many FASTA records are in the reference and sample files?
! TASK: Which genomic region is represented by the file names?
! TASK: Which file is the multi-allelic panel for PanGenie-index?
! TASK: Which file is the biallelic panel for conversion?
! TASK: What samples are represented in the panel VCF header?
```

### **3. Build the PanGenie Index**

Build a PanGenie index from the multi-allelic panel VCF and reference FASTA:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ PanGenie-index \
  -v data-pangenie/panel_multi_chr5:50200000-50400000.vcf \
  -r data-pangenie/reference_chr5:50200000-50400000.fasta \
  -o indexes/genotypes \
  -t 1 \
  -e 100000 \
  > indexing.log 2>&1
```

Command parameters:

- `-v`: input pangenome panel VCF
- `-r`: reference FASTA used with the panel
- `-o`: output prefix for PanGenie index files
- `-t`: number of threads
- `-e`: maximum nesting level parameter used by the workshop example

Inspect the index output:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh indexes
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ less indexing.log
```

```diff
! TASK: Did PanGenie-index finish successfully?
! TASK: Which files were created under indexes?
! TASK: What input VCF and reference FASTA were used?
```

### **4. Genotype the Samples**

Run PanGenie for sample `NA19189`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ PanGenie \
  -f indexes/genotypes \
  -i data-pangenie/NA19189_chr5:50200000-50400000.fasta \
  -o genotypes/NA19189_pangenie_multi \
  -t 1 \
  -j 1 \
  -s NA19189 \
  -e 10
```

Run PanGenie for sample `NA19190`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ PanGenie \
  -f indexes/genotypes \
  -i data-pangenie/NA19190_chr5:50200000-50400000.fasta \
  -o genotypes/NA19190_pangenie_multi \
  -t 1 \
  -j 1 \
  -s NA19190 \
  -e 10
```

Run PanGenie for sample `NA19191`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ PanGenie \
  -f indexes/genotypes \
  -i data-pangenie/NA19191_chr5:50200000-50400000.fasta \
  -o genotypes/NA19191_pangenie_multi \
  -t 1 \
  -j 1 \
  -s NA19191 \
  -e 10
```

Command parameters:

- `-f indexes/genotypes`: prefix of the index files produced earlier by `PanGenie-index -o indexes/genotypes`
- `-i data-pangenie/NA19189_chr5:50200000-50400000.fasta`: input sequence file for the sample being genotyped; replace the sample name for `NA19190` and `NA19191`
- `-o genotypes/NA19189_pangenie_multi`: output prefix; PanGenie will write files using this prefix, including `NA19189_pangenie_multi_genotyping.vcf`
- `-t 1`: number of threads for the main genotyping algorithm
- `-j 1`: number of threads for k-mer counting
- `-s NA19189`: sample name written into the output VCF; change this to match each sample
- `-e 10`: Jellyfish hash size setting used for this small workshop example

In this section, the command is run once per sample. The index from section 3 is reused, but the `-i`, `-o`, and `-s` values change for each sample.

Inspect the genotype outputs:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh genotypes
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ less -S genotypes/NA19189_pangenie_multi_genotyping.vcf
```

The `*_histogram.histo` file is a k-mer abundance histogram for the sample sequence data. It is not another genotype or variant file. PanGenie creates it while counting sample k-mers and uses it to estimate k-mer coverage before calculating genotype probabilities.

```diff
! TASK: Which output files were produced for each sample?
! TASK: Which sample name appears in each VCF header?
! TASK: What genotype fields are reported by PanGenie?
```

### **5. Convert Bubble Genotypes to Biallelic Genotypes**

Convert the PanGenie bubble genotypes for `NA19189`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ cat genotypes/NA19189_pangenie_multi_genotyping.vcf \
  | python3 data-pangenie/software/convert-to-biallelic.py \
      data-pangenie/panel_bi_chr5:50200000-50400000.vcf.gz \
  | bgzip \
  > genotypes/NA19189_pangenie_bi_genotyping.vcf.gz

user1@vm-corso-colonna:/scratch/user1/pangenie-working$ tabix -p vcf genotypes/NA19189_pangenie_bi_genotyping.vcf.gz
```

Convert the PanGenie bubble genotypes for `NA19190`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ cat genotypes/NA19190_pangenie_multi_genotyping.vcf \
  | python3 data-pangenie/software/convert-to-biallelic.py \
      data-pangenie/panel_bi_chr5:50200000-50400000.vcf.gz \
  | bgzip \
  > genotypes/NA19190_pangenie_bi_genotyping.vcf.gz

user1@vm-corso-colonna:/scratch/user1/pangenie-working$ tabix -p vcf genotypes/NA19190_pangenie_bi_genotyping.vcf.gz
```

Convert the PanGenie bubble genotypes for `NA19191`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ cat genotypes/NA19191_pangenie_multi_genotyping.vcf \
  | python3 data-pangenie/software/convert-to-biallelic.py \
      data-pangenie/panel_bi_chr5:50200000-50400000.vcf.gz \
  | bgzip \
  > genotypes/NA19191_pangenie_bi_genotyping.vcf.gz

user1@vm-corso-colonna:/scratch/user1/pangenie-working$ tabix -p vcf genotypes/NA19191_pangenie_bi_genotyping.vcf.gz
```

Inspect the converted VCF files:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh genotypes/*_pangenie_bi_genotyping.vcf.gz*
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools view -h genotypes/NA19189_pangenie_bi_genotyping.vcf.gz | less
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools view genotypes/NA19189_pangenie_bi_genotyping.vcf.gz | less -S
```

```diff
! TASK: Confirm that each converted VCF has a tabix index
! TASK: Compare the multi-allelic PanGenie output with the converted biallelic output
! TASK: Record one genotype example from each sample
```

### **6. Merge the Per-Sample VCF Files**

After converting each sample to a compressed and indexed biallelic VCF, merge the three samples into one multi-sample VCF:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools merge \
  genotypes/NA19189_pangenie_bi_genotyping.vcf.gz \
  genotypes/NA19190_pangenie_bi_genotyping.vcf.gz \
  genotypes/NA19191_pangenie_bi_genotyping.vcf.gz \
  -Oz \
  -o genotypes/pangenie_bi_merged.vcf.gz

user1@vm-corso-colonna:/scratch/user1/pangenie-working$ tabix -p vcf genotypes/pangenie_bi_merged.vcf.gz
```

Command parameters:

- `bcftools merge`: combine VCF files that contain the same variant records but different samples
- `genotypes/NA19189_pangenie_bi_genotyping.vcf.gz`: one compressed and indexed single-sample VCF
- `-Oz`: write compressed VCF output
- `-o genotypes/pangenie_bi_merged.vcf.gz`: output merged multi-sample VCF
- `tabix -p vcf`: create an index for the merged compressed VCF

Inspect the merged VCF:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ ls -lh genotypes/pangenie_bi_merged.vcf.gz*
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools query -l genotypes/pangenie_bi_merged.vcf.gz
user1@vm-corso-colonna:/scratch/user1/pangenie-working$ bcftools view genotypes/pangenie_bi_merged.vcf.gz | less -S
```

```diff
! TASK: Which samples are present in the merged VCF?
! TASK: How many variant records are in the merged VCF?
! TASK: Choose one variant and compare the genotypes across NA19189, NA19190, and NA19191
```

### **7. Report**

Prepare a short report:

```markdown
# PanGenie genotyping report

## Source tutorial
- tutorial:
- author:

## Input data
- panel VCF:
- reference FASTA:
- samples:

## Indexing
- command:
- output files:

## Genotyping
- sample:
- output VCF:
- number of records:

## Biallelic conversion
- converted VCF:
- index file:
- example genotype:

## Merged VCF
- merged VCF:
- samples:
- example variant across samples:
```

### **8. Credit**

This practical is adapted from the [PanGenie workshop](https://pangenie-workshop.readthedocs.io/en/latest/index.html) by Jana Ebler.
