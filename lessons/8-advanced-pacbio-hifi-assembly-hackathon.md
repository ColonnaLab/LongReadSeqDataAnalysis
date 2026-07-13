[back to course home page ](../README.md)

## Advanced 2: PacBio HiFi assembly hackathon

This is a small-group challenge using PacBio HiFi reads from *Escherichia coli* K-12.

Dataset:

- SRA accession: `SRR10971019`
- data type: PacBio HiFi / CCS reads
- organism: *Escherichia coli* K-12
- expected genome size: approximately 4.64 Mb
- expected structure: one circular chromosome
- reference: `NC_000913.3`
- reference assembly: `GCF_000005845.2`
- full download size: approximately 3 GB
- full coverage: approximately 290x
- server data folder: `/data/user1/data-pacbio`

```diff
+ This is an advanced hackathon
+ You should decide whether to use the full dataset or a subsample
+ You should justify assembler choice and parameters
+ You should compare assembly statistics with reference-based evaluation
```

### **1. Prepare the Working Directory**

Keep downloaded data in `/data/user1/data-pacbio` and write analysis results in a separate working folder.

```bash
user1@vm-corso-colonna:~$ mkdir -p /scratch/user1/pacbio-working
user1@vm-corso-colonna:~$ cd /scratch/user1/pacbio-working
```

Create output folders:

```bash
user1@vm-corso-colonna:/data/user1/pacbio-working$ mkdir -p assembly assembly/hifiasm assembly/flye assembly/qc assembly/reference_comparison assembly/report
```


```diff
! TASK: Identify the HiFi reads file
! TASK: Identify the reference FASTA file
! TASK: Check file sizes
! TASK: Check free disk space before running assembly
```

Expected file names used below:

- reads: `/scratch/user1/pacbio-working/data-pacbio/ecoli_hifi.fastq
- reference: `/data/user1/data-pacbio/NC_000913.3.fasta`

If your downloaded files have different names or are inside subfolders, use the correct paths from the `find` command.

### **2. Inspect the HiFi Reads**

Summarize the reads:

```bash
user1@vm-corso-colonna:/data/user1/pacbio-working$ seqkit stats  -a \
 data-pacbio/ecoli_hifi.fastq \
 > assembly/qc/ecoli.seqkit.stats.txt

```

```bash
user1@vm-corso-colonna:/data/user1/pacbio-working$ cat assembly/qc/SRR10971019.seqkit.stats.txt
```

```diff
! TASK: Record the number of reads
! TASK: Record total bases
! TASK: Record mean read length
! TASK: Record N50 read length
! TASK: Estimate coverage using a 4.64 Mb genome size
```

### **3. Optional: Subsample the Reads**

The full dataset is approximately 290x coverage. A smaller subset may assemble faster during class.

Example: keep about 25% of the reads.

```bash
user1@vm-corso-colonna:/data/user1/pacbio-working$ seqkit sample \
  -p 0.25 \
  -s 11 \
   data-pacbio/ecoli_hifi.fastq \
  -o assembly/qc/ecoli_hifi.subsample25.fastq.gz
```

```diff
! OPTIONAL TASK: Estimate the coverage of the subsampled reads
! OPTIONAL TASK: Decide whether to assemble full data or subsampled data
```

The commands below use the full dataset. If you use the subsample, replace the reads path with:

```diff
+ assembly/qc/ecoli_hifi.subsample25.fastq.gz
```

### **4. Run a First HiFi Assembly with hifiasm**

Run `hifiasm`:

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ hifiasm \
  -o assembly/hifiasm/ecoli_k12 \
  -t 4 \
  data-pacbio/ecoli_hifi.fastq \
  2> assembly/hifiasm/hifiasm.log
```

`hifiasm` writes graph files in GFA format. Convert the primary contig graph to FASTA:

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ awk \
  '/^S/{print ">"$2"\n"$3}' \
  assembly/hifiasm/ecoli_k12.bp.p_ctg.gfa \
  > assembly/hifiasm/ecoli_k12.hifiasm.p_ctg.fasta
```

```diff
! TASK: Did hifiasm finish successfully?
! TASK: Which GFA files were produced?
! TASK: How many contigs are in the converted FASTA?
```

### **5. Optional: Run Flye for Comparison**

Flye can also assemble PacBio HiFi reads.

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ flye \
  --pacbio-hifi data-pacbio/ecoli_hifi.fastq \
  --genome-size 4.64m \
  --threads 4 \
  --out-dir assembly/flye
```

```diff
! OPTIONAL TASK: Compare hifiasm and Flye results
! OPTIONAL TASK: Which assembly is more contiguous?
! OPTIONAL TASK: Which assembly is closer to the expected genome size?
```

### **6. Summarize the Assembly**

Summarize the hifiasm FASTA:

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ seqkit stats \
  assembly/hifiasm/ecoli_k12.hifiasm.p_ctg.fasta \
  > assembly/qc/ecoli_k12.hifiasm.seqkit.stats.txt
```

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ cat assembly/qc/ecoli_k12.hifiasm.seqkit.stats.txt
```

```diff
! TASK: How many contigs were produced?
! TASK: What is the total assembly length?
! TASK: What is the longest contig?
! TASK: What is the N50?
! TASK: Is the assembly close to 4.64 Mb?
```

### **7. Compare the Assembly to the Reference**

Use QUAST if available:

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ quast.py \
  assembly/hifiasm/ecoli_k12.hifiasm.p_ctg.fasta \
  -r /data/user1/data-pacbio/NC_000913.3.fasta \
  -o assembly/reference_comparison/quast_hifiasm
```

```diff
! TASK: Record genome fraction
! TASK: Record number of contigs
! TASK: Record major warnings or misassemblies, if reported
```

If QUAST is not available, align the assembly to the reference:

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ minimap2 \
  -ax asm5 \
  /data/user1/data-pacbio/NC_000913.3.fasta \
  assembly/hifiasm/ecoli_k12.hifiasm.p_ctg.fasta \
  > assembly/reference_comparison/hifiasm_vs_reference.sam
```

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ samtools sort \
  -o assembly/reference_comparison/hifiasm_vs_reference.sorted.bam \
  assembly/reference_comparison/hifiasm_vs_reference.sam
```

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ samtools index \
  assembly/reference_comparison/hifiasm_vs_reference.sorted.bam
```

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ samtools coverage \
  assembly/reference_comparison/hifiasm_vs_reference.sorted.bam \
  > assembly/reference_comparison/hifiasm_vs_reference.coverage.txt
```

```bash
user1@vm-corso-colonna:/scratch/user1/pacbio-working$ cat assembly/reference_comparison/hifiasm_vs_reference.coverage.txt
```

```diff
! TASK: Does one contig cover most or all of NC_000913.3?
! TASK: Are there uncovered reference regions?
! TASK: Are there signs of fragmentation or duplicated sequence?
```

### **8. Visual Inspection**

Copy these files to your local machine for IGV inspection:

- `/data/user1/data-pacbio/NC_000913.3.fasta`
- `/data/user1/data-pacbio/NC_000913.3.fasta.fai`
- `assembly/reference_comparison/hifiasm_vs_reference.sorted.bam`
- `assembly/reference_comparison/hifiasm_vs_reference.sorted.bam.bai`

```diff
! TASK: Load the E. coli reference in IGV
! TASK: Load the assembly-reference BAM
! TASK: Inspect whether the assembly covers the chromosome continuously
```

### **9. Final Report**

Prepare a group report in:

```diff
+ assembly/report/pacbio_hifi_assembly_report.md
```

Use this template:

```markdown
# PacBio HiFi assembly hackathon report

## Team
- names:

## Input data
- accession:
- reads file:
- total bases:
- estimated coverage:
- full or subsampled data:

## Assembly command
- assembler:
- command:
- parameters:

## Assembly summary
- number of contigs:
- total assembly length:
- longest contig:
- N50:

## Reference comparison
- reference:
- chromosome recovered?
- major differences or warnings:

## Interpretation
- what worked?
- what was uncertain?
- what would we try next?
```

```diff
! FINAL TASK: Compare the HiFi assembly experience with the ONT assembly experience
```
