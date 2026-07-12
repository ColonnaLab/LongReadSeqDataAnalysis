[back to course home page ](../README.md)

## Advanced 1: ONT long-read assembly hackathon

This is a small-group challenge. Work in teams, make decisions, keep notes, and prepare a short report.

The goal is to assemble the MRSA KUN1163 Oxford Nanopore Technologies (ONT) reads and evaluate whether the assembly recovered the chromosome and plasmid.

Input files:

- reads: `data-longreads/fastq/DRR187567.fastq.gz`
- reference genome for comparison: `data-longreads/reference/KUN1163_reference.fasta`
- data type: Oxford Nanopore reads

```diff
+ This is not a copy-paste beginner exercise
+ You should decide parameters, inspect outputs, and justify choices
+ Keep track of every command that produces an important result
```

### **1. Prepare the Working Directory**

Start from the long-read working directory:

```bash
user1@vm-corso-colonna:~$ cd /data/user1/lr-working
```

Create output folders:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ mkdir -p assembly assembly/flye assembly/qc assembly/reference_comparison assembly/report
```

```diff
! TASK: Check that the reads and reference FASTA are present
! TASK: Check how much free disk space is available
! TASK: Decide where large temporary and output files should be written
```

Useful commands:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh data-longreads/fastq/
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh data-longreads/reference/
user1@vm-corso-colonna:/data/user1/lr-working$ df -h .
```

### **2. Inspect the Input Reads**

Use the QC results from lesson 4 or regenerate a compact summary.

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ seqkit stats \
  data-longreads/fastq/DRR187567.fastq.gz \
  > assembly/qc/DRR187567.seqkit.stats.txt
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ cat assembly/qc/DRR187567.seqkit.stats.txt
```

```diff
! TASK: Record the number of reads
! TASK: Record total bases
! TASK: Record mean read length
! TASK: Record N50 read length
! TASK: Estimate whether the dataset has enough coverage for bacterial assembly
```

### **3. Choose Assembly Parameters**

Use the published KUN1163 reference as a guide for expected genome size. The reference contains a chromosome and plasmid, so the expected assembly size should be close to the combined reference length.

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ cat data-longreads/reference/KUN1163_reference.fasta.fai
```

```diff
! TASK: Estimate the total reference length
! TASK: Choose an approximate genome size for the assembler
! TASK: Decide how many threads to use
```

For Flye, useful options include:

```diff
+ --nano-raw: Oxford Nanopore reads
+ --genome-size: approximate genome size
+ --threads: number of CPU threads
+ --out-dir: output directory
```

### **4. Run a First Assembly**

Recommended first attempt with Flye:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ flye \
  --nano-raw data-longreads/fastq/DRR187567.fastq.gz \
  --genome-size 3m \
  --threads 4 \
  --out-dir assembly/flye
```

```diff
! TASK: Did Flye finish successfully?
! TASK: Which output file contains the assembled contigs?
! TASK: Did Flye report any circular contigs?
```

Check the output:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ ls -lh assembly/flye
```

### **5. Summarize the Assembly**

Use `seqkit stats` on the assembly FASTA:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ seqkit stats \
  assembly/flye/assembly.fasta \
  > assembly/qc/flye.assembly.seqkit.stats.txt
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ cat assembly/qc/flye.assembly.seqkit.stats.txt
```

```diff
! TASK: How many contigs were produced?
! TASK: What is the total assembly length?
! TASK: What is the longest contig?
! TASK: What is the assembly N50?
! TASK: Is the total length close to the expected genome size?
```

### **6. Compare the Assembly to the Reference**

Use QUAST if available:

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ quast.py \
  assembly/flye/assembly.fasta \
  -r data-longreads/reference/KUN1163_reference.fasta \
  -o assembly/reference_comparison/quast_flye
```

```diff
! TASK: Open the QUAST report
! TASK: Record genome fraction
! TASK: Record number of contigs
! TASK: Record major warnings or misassemblies, if reported
```

If QUAST is not available, use assembly-to-reference alignment.

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ minimap2 \
  -ax asm5 \
  data-longreads/reference/KUN1163_reference.fasta \
  assembly/flye/assembly.fasta \
  > assembly/reference_comparison/flye_vs_reference.sam
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ samtools sort \
  -o assembly/reference_comparison/flye_vs_reference.sorted.bam \
  assembly/reference_comparison/flye_vs_reference.sam
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ samtools index \
  assembly/reference_comparison/flye_vs_reference.sorted.bam
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ samtools coverage \
  assembly/reference_comparison/flye_vs_reference.sorted.bam \
  > assembly/reference_comparison/flye_vs_reference.coverage.txt
```

```bash
user1@vm-corso-colonna:/data/user1/lr-working$ cat assembly/reference_comparison/flye_vs_reference.coverage.txt
```

```diff
! TASK: Does the assembly align to the chromosome?
! TASK: Does the assembly align to the plasmid?
! TASK: Are any reference sequences poorly covered?
```

### **7. Visual Inspection**

Copy the assembly-reference BAM and reference files to your local machine for IGV inspection.

Files to copy:

- `data-longreads/reference/KUN1163_reference.fasta`
- `data-longreads/reference/KUN1163_reference.fasta.fai`
- `assembly/reference_comparison/flye_vs_reference.sorted.bam`
- `assembly/reference_comparison/flye_vs_reference.sorted.bam.bai`

```diff
! TASK: Inspect the chromosome in IGV
! TASK: Inspect the plasmid in IGV
! TASK: Look for breaks, missing regions, or suspicious alignments
```

### **8. Optional Extensions**

Choose one extension if time allows:

- Run a second assembler and compare results
- Change Flye parameters and evaluate whether the result improves
- Polish the assembly
- Inspect the assembly graph, if available
- Investigate unmapped reads from lesson 5
- Compare assembly-based differences with variant calls from lesson 6

```diff
! OPTIONAL TASK: Explain what changed and whether it improved the assembly
```

### **9. Final Report**

Prepare a short group report in `assembly/report/assembly_report.md`.

Use this template:

```markdown
# Long-read assembly hackathon report

## Team
- names:

## Input data
- reads:
- number of reads:
- total bases:
- estimated genome size:
- estimated coverage:

## Assembly command
- assembler:
- command:
- parameters we chose:

## Assembly summary
- number of contigs:
- total assembly length:
- longest contig:
- N50:
- circular contigs, if reported:

## Reference comparison
- chromosome recovered?
- plasmid recovered?
- major differences or warnings:

## Interpretation
- what worked?
- what was uncertain?
- what would we try next?
```

```diff
! FINAL TASK: Be ready to explain your assembly decisions to the group
```
