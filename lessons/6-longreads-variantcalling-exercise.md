[back to course home page ](../README.md)

## Long reads variant calling exercise

In this practical we will use the aligned Oxford Nanopore reads from MRSA strain KUN1163 to call and inspect small variants relative to the reference chromosome and plasmid.

The input files are:

- reference genome: `data-longreads/reference/KUN1163_reference.fasta`
- aligned reads: `bam/DRR187567.KUN1163.sorted.bam`
- BAM index: `bam/DRR187567.KUN1163.sorted.bam.bai`

```diff
+ We will work inside /data/user1/lr-working
+ We will use the BAM file produced in lesson 5
+ We will call variants with Clair3 and DeepVariant
+ We will inspect, filter, normalize, and compare VCF files
```

### **1. Prepare the Working Directory**

Move to the long reads working folder and check that the alignment files exist:

```bash
user1@vm-corso-colonna:~$ cd /data/user1/lr-working
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh bam/
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh data-longreads/reference/
```

Create folders for variant calling results:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ mkdir -p variants variants/clair3 variants/deepvariant variants/compare variants/logs
```

```diff
! EXERCISE: Confirm that the sorted BAM and BAI files are present
! EXERCISE: Confirm that the reference FASTA and FAI files are present
```

### **2. Check the BAM and Reference Indexes**

Variant calling tools need indexed input files.

If the BAM index is missing, create it:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ samtools index bam/DRR187567.KUN1163.sorted.bam
```

If the reference FASTA index is missing, create it:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ samtools faidx data-longreads/reference/KUN1163_reference.fasta
```

Check the indexed files:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh bam/DRR187567.KUN1163.sorted.bam*
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh data-longreads/reference/KUN1163_reference.fasta*
```

### **3. Inspect Coverage Before Variant Calling**

Coverage affects variant calling reliability. Positions with very low coverage are more difficult to call confidently.

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ samtools coverage \
  bam/DRR187567.KUN1163.sorted.bam \
  > variants/DRR187567.coverage.before_variant_calling.txt
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ cat variants/DRR187567.coverage.before_variant_calling.txt
```

```diff
! EXERCISE: Is the chromosome coverage high enough for variant calling?
! EXERCISE: Is the plasmid coverage similar to the chromosome coverage?
! EXERCISE: Which low-coverage regions might be difficult to call confidently?
```

### **4. Small Variant Calling with Clair3**

[Clair3](https://github.com/HKU-BAL/Clair3) is a long-read small-variant caller. It combines a fast pileup model with a more detailed full-alignment model.

For Oxford Nanopore data, Clair3 needs:

- the BAM file
- the reference FASTA
- `--platform=ont`
- a pre-trained model path
- an output directory

#### Activate the Clair3 conda environment

`conda` is a tool that manages software environments. An environment is a prepared set of programs and libraries. `run_clair3.sh` must be run from the Clair3 conda environment, otherwise Clair3 cannot find all the files it needs.

Activate the environment and check that the command is available:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ source $(conda info --base)/etc/profile.d/conda.sh
user1@vm-corso-colonna:/data/user1/lr-working$ conda activate clair3
user1@vm-corso-colonna:/data/user1/lr-working$ which run_clair3.sh
user1@vm-corso-colonna:/data/user1/lr-working$ echo $CONDA_PREFIX
```

The `source` command loads conda activation support for the current terminal. This is needed if `conda activate clair3` gives an error asking you to run `conda init`.

`which run_clair3.sh` shows where the Clair3 command is installed. `echo $CONDA_PREFIX` shows which conda environment is currently active.

Check the available Clair3 models:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls $CONDA_PREFIX/bin/models
```

`$CONDA_PREFIX` is the path to the active conda environment. After `conda activate clair3`, `$CONDA_PREFIX/bin/models` is the folder where Clair3 models are installed in this environment.

#### Choose the Clair3 model

Choose the model that best matches how the reads were generated. For this lesson the reads are Oxford Nanopore reads, so the model must be an ONT model.

Useful clues in model names:

- `r941`: Oxford Nanopore R9.4.1 chemistry
- `r1041`: Oxford Nanopore R10.4.1 chemistry
- `sup`: super-accurate basecalling
- `hac`: high-accuracy basecalling
- `hifi`: PacBio HiFi, not for ONT reads
- `ilmn`: Illumina, not for ONT reads

If you do not know the exact chemistry or basecaller, use the ONT model recommended by the instructor or the one that best matches the dataset documentation.

Choose the ONT model name from the list printed by `ls`. The model name in `--model_path` must match one folder name exactly.

#### Run Clair3

Before running Clair3, check that the model folder exists. Replace `EXACT_MODEL_NAME_FROM_LS` with one of the model names printed by the previous `ls` command:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls -ld $CONDA_PREFIX/bin/models/EXACT_MODEL_NAME_FROM_LS
```

Run Clair3 using that exact model name:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ run_clair3.sh \
  --bam_fn=bam/DRR187567.KUN1163.sorted.bam \
  --ref_fn=data-longreads/reference/KUN1163_reference.fasta \
  --threads=4 \
  --platform=ont \
  --model_path=$CONDA_PREFIX/bin/models/r941_prom_sup_g5014 \
  --output=/data/user1/lr-working/variants/clair3 \
  --include_all_ctgs \
  --no_phasing_for_fa \
  --haploid_precise
```

```diff
+ --platform=ont tells Clair3 that the reads are Oxford Nanopore reads
+ --haploid_precise is appropriate for a haploid bacterial genome
+ --include_all_ctgs asks Clair3 to call all reference sequences, including the plasmid
+ --output uses an absolute path to avoid output path warnings
```

Check the Clair3 output:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh variants/clair3
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh variants/clair3/merge_output.vcf.gz*
```

The main small-variant output is:

```diff
+ variants/clair3/merge_output.vcf.gz
```

```diff
! EXERCISE: Did Clair3 create a compressed VCF file?
! EXERCISE: Did Clair3 create a VCF index file?
```

### **5. Inspect the Clair3 VCF**

View the header:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -h \
  variants/clair3/merge_output.vcf.gz \
  | less
```

View the first variant records:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -H \
  variants/clair3/merge_output.vcf.gz \
  | head
```

Count the number of Clair3 variant records:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -H \
  variants/clair3/merge_output.vcf.gz \
  | wc -l
```

Count variants by reference sequence:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools query \
  -f '%CHROM\n' \
  variants/clair3/merge_output.vcf.gz \
  | sort \
  | uniq -c
```

```diff
! EXERCISE: How many variants did Clair3 call?
! EXERCISE: Which reference sequence contains most variants?
! EXERCISE: Look at the REF and ALT columns. Do you see SNVs, indels, or both?
```

### **6. Small Variant Calling with DeepVariant**

[DeepVariant](https://github.com/google/deepvariant) is a deep-learning variant calling pipeline. It uses `--model_type` to choose the model appropriate for the sequencing data.

For ONT R10.4 data, DeepVariant uses:

```diff
+ --model_type=ONT_R104
```

DeepVariant's released germline models are mainly trained for human, diploid variant calling. The MRSA reads in this lesson are bacterial and may also have been produced with an older ONT chemistry. Treat the DeepVariant result as a teaching comparison rather than a validated production call set.

#### Activate the DeepVariant conda environment

Like Clair3, DeepVariant must be run from the conda environment where it is installed.

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ source $(conda info --base)/etc/profile.d/conda.sh
user1@vm-corso-colonna:/data/user1/lr-working$ conda activate deepvariant
user1@vm-corso-colonna:/data/user1/lr-working$ echo $CONDA_PREFIX
user1@vm-corso-colonna:/data/user1/lr-working$ ls $CONDA_PREFIX/bin/run_deepvariant*
```

The last command checks which DeepVariant run scripts are installed in the active environment.

#### Run DeepVariant

Run DeepVariant with the Python wrapper installed in the conda environment:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ python $CONDA_PREFIX/bin/run_deepvariant.py \
  --model_type=ONT_R104 \
  --ref=data-longreads/reference/KUN1163_reference.fasta \
  --reads=bam/DRR187567.KUN1163.sorted.bam \
  --output_vcf=variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  --output_gvcf=variants/deepvariant/DRR187567.deepvariant.g.vcf.gz \
  --num_shards=4 \
  --haploid_contigs=AP020324.1,AP020325.1
```

Check the DeepVariant output:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh variants/deepvariant
```

If the VCF index is missing, create it:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ tabix \
  -p vcf \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz
```

```diff
! EXERCISE: Which DeepVariant model type did we use?
! EXERCISE: Why should model choice be treated carefully?
! EXERCISE: Why did we mark both MRSA reference sequences as haploid contigs?
```

### **7. Inspect the DeepVariant VCF**

View the header:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -h \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  | less
```

View the first variant records:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -H \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  | head
```

Count the number of DeepVariant variant records:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -H \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  | wc -l
```

Count variants by reference sequence:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools query \
  -f '%CHROM\n' \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  | sort \
  | uniq -c
```

```diff
! EXERCISE: How many variants did DeepVariant call?
! EXERCISE: Are calls present on both the chromosome and plasmid?
! EXERCISE: Does the DeepVariant VCF contain SNVs, indels, or both?
```

### **8. Summarize the VCF Files with bcftools**

Create summary files for both callers:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools stats \
  variants/clair3/merge_output.vcf.gz \
  > variants/clair3/DRR187567.clair3.stats.txt
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools stats \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  > variants/deepvariant/DRR187567.deepvariant.stats.txt
```

Inspect the summaries:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ less variants/clair3/DRR187567.clair3.stats.txt
user1@vm-corso-colonna:/data/user1/lr-working$ less variants/deepvariant/DRR187567.deepvariant.stats.txt
```

Extract the main summary lines:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ grep "^SN" variants/clair3/DRR187567.clair3.stats.txt
user1@vm-corso-colonna:/data/user1/lr-working$ grep "^SN" variants/deepvariant/DRR187567.deepvariant.stats.txt
```

```diff
! EXERCISE: Which caller produced more total variant records?
! EXERCISE: Which caller produced more SNPs?
! EXERCISE: Which caller produced more indels?
```

### **9. Filter Variants**

Filtering removes low-confidence calls or creates a smaller set of variants for inspection.

Filter Clair3 calls by quality:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools filter \
  -i 'QUAL>20' \
  variants/clair3/merge_output.vcf.gz \
  -Oz \
  -o variants/clair3/DRR187567.clair3.QUAL20.vcf.gz
```

Filter DeepVariant calls by quality:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools filter \
  -i 'QUAL>20' \
  variants/deepvariant/DRR187567.deepvariant.vcf.gz \
  -Oz \
  -o variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz
```

Index the filtered VCF files:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ tabix -p vcf variants/clair3/DRR187567.clair3.QUAL20.vcf.gz
user1@vm-corso-colonna:/data/user1/lr-working$ tabix -p vcf variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz
```

Count variants before and after filtering:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view -H variants/clair3/merge_output.vcf.gz | wc -l
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view -H variants/clair3/DRR187567.clair3.QUAL20.vcf.gz | wc -l
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view -H variants/deepvariant/DRR187567.deepvariant.vcf.gz | wc -l
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view -H variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz | wc -l
```

```diff
! EXERCISE: How many Clair3 variants remain after QUAL>20 filtering?
! EXERCISE: How many DeepVariant variants remain after QUAL>20 filtering?
! EXERCISE: Did filtering affect the two callers in the same way?
```

### **10. Inspect Variants in One Region**

Use `bcftools view` to inspect variants in a specific region.

Example for the chromosome:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -r AP020324.1:1-100000 \
  variants/clair3/DRR187567.clair3.QUAL20.vcf.gz \
  | less
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -r AP020324.1:1-100000 \
  variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz \
  | less
```

Example for the plasmid:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -r AP020325.1:1-30220 \
  variants/clair3/DRR187567.clair3.QUAL20.vcf.gz \
  | less
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools view \
  -r AP020325.1:1-30220 \
  variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz \
  | less
```

```diff
! EXERCISE: Are variants present on the plasmid?
! EXERCISE: Do Clair3 and DeepVariant report the same calls in these regions?
! EXERCISE: Compare the QUAL values of variants in different regions
```

### **11. Normalize VCF Files Before Comparison**

Indels can be represented in more than one equivalent way. Normalization makes VCF comparison more reliable.

Normalize the filtered Clair3 VCF:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools norm \
  -f data-longreads/reference/KUN1163_reference.fasta \
  -Oz \
  -o variants/compare/DRR187567.clair3.QUAL20.norm.vcf.gz \
  variants/clair3/DRR187567.clair3.QUAL20.vcf.gz
```

Normalize the filtered DeepVariant VCF:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools norm \
  -f data-longreads/reference/KUN1163_reference.fasta \
  -Oz \
  -o variants/compare/DRR187567.deepvariant.QUAL20.norm.vcf.gz \
  variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz
```

Index the normalized files:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ tabix -p vcf variants/compare/DRR187567.clair3.QUAL20.norm.vcf.gz
user1@vm-corso-colonna:/data/user1/lr-working$ tabix -p vcf variants/compare/DRR187567.deepvariant.QUAL20.norm.vcf.gz
```

### **12. Compare Clair3 and DeepVariant Calls**

Use `bcftools isec` to compare the two normalized VCF files.

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ mkdir -p variants/compare/isec
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ bcftools isec \
  -p variants/compare/isec \
  variants/compare/DRR187567.clair3.QUAL20.norm.vcf.gz \
  variants/compare/DRR187567.deepvariant.QUAL20.norm.vcf.gz
```

The output files include:

- `0000.vcf`: records unique to the first file, Clair3
- `0001.vcf`: records unique to the second file, DeepVariant
- `0002.vcf`: records from the first file that are shared
- `0003.vcf`: records from the second file that are shared

Count each group:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ grep -vc "^#" variants/compare/isec/0000.vcf
user1@vm-corso-colonna:/data/user1/lr-working$ grep -vc "^#" variants/compare/isec/0001.vcf
user1@vm-corso-colonna:/data/user1/lr-working$ grep -vc "^#" variants/compare/isec/0002.vcf
user1@vm-corso-colonna:/data/user1/lr-working$ grep -vc "^#" variants/compare/isec/0003.vcf
```

```diff
! EXERCISE: How many filtered variants are shared by both callers?
! EXERCISE: How many filtered variants are unique to Clair3?
! EXERCISE: How many filtered variants are unique to DeepVariant?
! EXERCISE: Why might two high-quality callers disagree?
```

### **13. Structural Variant Calling Concepts**

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

### **14. Save Results for Later**

Files that should now be present:

- `variants/DRR187567.coverage.before_variant_calling.txt`
- `variants/clair3/merge_output.vcf.gz`
- `variants/clair3/merge_output.vcf.gz.tbi`
- `variants/clair3/DRR187567.clair3.stats.txt`
- `variants/clair3/DRR187567.clair3.QUAL20.vcf.gz`
- `variants/deepvariant/DRR187567.deepvariant.vcf.gz`
- `variants/deepvariant/DRR187567.deepvariant.vcf.gz.tbi`
- `variants/deepvariant/DRR187567.deepvariant.stats.txt`
- `variants/deepvariant/DRR187567.deepvariant.QUAL20.vcf.gz`
- `variants/compare/DRR187567.clair3.QUAL20.norm.vcf.gz`
- `variants/compare/DRR187567.deepvariant.QUAL20.norm.vcf.gz`
- `variants/compare/isec/`

Copy selected results to your local machine from the `lr-local-work` folder:

```bash
lr-local-work$ scp user1@212.189.205.193:/data/user1/lr-working/variants/*/*.txt .
lr-local-work$ scp user1@212.189.205.193:/data/user1/lr-working/variants/clair3/*.vcf.gz .
lr-local-work$ scp user1@212.189.205.193:/data/user1/lr-working/variants/deepvariant/*.vcf.gz .
```

```diff
! EXERCISE: Keep the filtered VCF files and the variant statistics files
! These are the main outputs of this practical
```
