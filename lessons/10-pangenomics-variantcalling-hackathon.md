[back to course home page ](../README.md)

## Advanced 4: Pangenomics variant-calling hackathon

This practical uses a small graph-first toy dataset to study graph-aware read mapping and graph-based variant interpretation with VG Giraffe. The dataset starts from an existing GFA graph, not from a VCF.

Dataset and outputs assumed on the server:

- downloaded dataset: `/scratch/user1/data-pangen-vc`
- input graph: `/scratch/user1/data-pangen-vc/toy_graph.gfa`
- paired reads: `/scratch/user1/data-pangen-vc/toy_reads_R1.fastq.gz`
- paired reads: `/scratch/user1/data-pangen-vc/toy_reads_R2.fastq.gz`
- working directory: `/scratch/user1/pangenome-vc-working`

The source toy dataset used to prepare this practical is:

```text
/home/enza/iusatoProtocol/teaching/testlongdir/pangen/giraffe/data-pangen-vc
```

```diff
+ You should inspect a supplied pangenome graph before mapping reads
+ You should build Giraffe indexes directly from a GFA graph
+ You should call read-supported graph variants after mapping reads
+ You should compare graph paths, graph-derived read alignments, and expected truth tables
```

### **1. Prepare the Working Directory**

Create a new working directory for the variant-calling practical:

```bash
user1@vm-corso-colonna:~$ mkdir -p /scratch/user1/pangenome-vc-working
user1@vm-corso-colonna:~$ cd /scratch/user1/pangenome-vc-working
```

Create a link to the variant-calling dataset:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ ln -s /scratch/user1/data-pangen-vc data-pangen-vc
```

Create output folders:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ mkdir -p graph/vg graph/qc graph/report
```

Confirm that the expected toy files are present:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ ls -lh data-pangen-vc
```

Important files:

- `toy_graph.gfa`: supplied GFA 1.1 graph
- `toy_reads_R1.fastq.gz` and `toy_reads_R2.fastq.gz`: paired-end reads
- `graph_bubbles.tsv`: expected variable sites
- `path_truth.tsv`: expected path alleles and path lengths
- `read_truth.tsv`: expected read-pair origin and bubble overlap
- `run_practical.sh`: helper script for indexing, mapping, and text outputs

```diff
! TASK: Confirm that this practical is using data-pangen-vc, not the lesson 9 graph
! TASK: Confirm that toy_graph.gfa and the paired read files exist
! TASK: Inspect metadata.json before running VG
```

### **2. Inspect the Supplied Graph**

The input graph is already built:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ cat  data-pangen-vc/toy_graph.gfa | less -S 
```

Count graph segments, links, and walks:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep '^S' data-pangen-vc/toy_graph.gfa | wc -l
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep '^L' data-pangen-vc/toy_graph.gfa | wc -l
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep '^W' data-pangen-vc/toy_graph.gfa
```

The graph contains:

- 11 sequence segments
- 14 links
- three W-line walks
- two SNPs, one 3 bp insertion, and one 3 bp deletion
- reference sample `REF`, marked in the GFA header as `RS:Z:REF`

Expected graph paths:

- `REF#0#toy_chr1`, 1,200 bp
- `TOY_DIPLOID#1#toy_chr1`, 1,203 bp
- `TOY_DIPLOID#2#toy_chr1`, 1,197 bp

Inspect the truth tables:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ column -t data-pangen-vc/graph_bubbles.tsv | less -S
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ column -t data-pangen-vc/path_truth.tsv | less -S
```

```diff
! TASK: Which sample is marked as the reference?
! TASK: Which graph nodes encode alternate alleles?
! TASK: Which haplotype contains the 3 bp insertion?
! TASK: Which haplotype skips the deletion-reference node?
```

### **3. Build VG Giraffe Indexes**

Build Giraffe indexes directly from the supplied GFA graph:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg autoindex \
  --workflow giraffe \
  --prefix graph/vg/toy_graph_index \
  --gfa data-pangen-vc/toy_graph.gfa \
  --threads 2
```

Command parameters:

- `--workflow giraffe`: build indexes for short-read mapping with VG Giraffe
- `--prefix graph/vg/toy_graph_index`: write index files under `graph/vg` with this filename prefix
- `--gfa data-pangen-vc/toy_graph.gfa`: use the supplied graph-first GFA as the input graph
- `--threads 2`: use two CPU threads

If your VG version does not recognize `giraffe` as the workflow name, try:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg autoindex \
  --workflow sr-giraffe \
  --prefix graph/vg/toy_graph_index \
  --gfa data-pangen-vc/toy_graph.gfa \
  --threads 2
```

The fallback uses the same inputs and outputs. Only the workflow name changes from `giraffe` to `sr-giraffe`.

Check generated index files:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ ls -lh graph/vg/toy_graph_index*
```

Expected index types include:

- GBZ graph index
- distance index
- minimizer index
- zipcodes file, depending on VG version

```diff
! TASK: Which exact index filenames did your VG version produce?
! TASK: Why is the GFA graph enough for this autoindex run?
! TASK: Why is no VCF needed for this dataset?
```

### **4. Map the Paired Reads**

The data folder is read-only input. Do not run `run_practical.sh` inside `data-pangen-vc`, because it writes output files in its current directory. Instead, run `vg giraffe` from the working directory and write all mapping outputs to `graph/vg` and `graph/qc`.

Map the paired reads to GAM:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg giraffe \
  -t 2 \
  -Z graph/vg/toy_graph_index.giraffe.gbz \
  -d graph/vg/toy_graph_index.dist \
  -m graph/vg/toy_graph_index.shortread.withzip.min \
  -z graph/vg/toy_graph_index.shortread.zipcodes \
  --fragment-mean 260 \
  --fragment-stdev 10 \
  -N toy_reads \
  -f data-pangen-vc/toy_reads_R1.fastq.gz \
  -f data-pangen-vc/toy_reads_R2.fastq.gz \
  > graph/vg/toy_mapped.gam
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg stats -a graph/vg/toy_mapped.gam > graph/qc/toy_mapping_stats.txt
```

Command parameters:

- `-t 2`: use two CPU threads
- `-Z graph/vg/toy_graph_index.giraffe.gbz`: indexed graph used by Giraffe
- `-d graph/vg/toy_graph_index.dist`: distance index used to evaluate graph distances between seeds
- `-m graph/vg/toy_graph_index.shortread.withzip.min`: minimizer index used to seed read mapping
- `-z graph/vg/toy_graph_index.shortread.zipcodes`: zipcode index paired with the minimizer index
- `--fragment-mean 260`: expected paired-end fragment length for this toy dataset
- `--fragment-stdev 10`: expected fragment length standard deviation
- `-N toy_reads`: sample/read-group name prefix for the output alignments
- `-f data-pangen-vc/toy_reads_R1.fastq.gz`: read 1 FASTQ input
- `-f data-pangen-vc/toy_reads_R2.fastq.gz`: read 2 FASTQ input
- `>`: write GAM alignments to `graph/vg/toy_mapped.gam`
- `vg stats -a`: summarize graph alignments from the GAM file

Create a text GAF alignment file:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg giraffe \
  -t 2 \
  -Z graph/vg/toy_graph_index.giraffe.gbz \
  -d graph/vg/toy_graph_index.dist \
  -m graph/vg/toy_graph_index.shortread.withzip.min \
  -z graph/vg/toy_graph_index.shortread.zipcodes \
  --fragment-mean 260 \
  --fragment-stdev 10 \
  -N toy_reads \
  -f data-pangen-vc/toy_reads_R1.fastq.gz \
  -f data-pangen-vc/toy_reads_R2.fastq.gz \
  -o GAF \
  > graph/vg/toy_mapped.gaf
```

Command parameters:

- `-o GAF`: write alignments in GAF text format
- `>`: write the GAF output to `graph/vg/toy_mapped.gaf`
- the graph index, fragment, sample, and read parameters are the same as in the GAM command above

Create a reference-projected BAM alignment file:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg giraffe \
  -t 2 \
  -Z graph/vg/toy_graph_index.giraffe.gbz \
  -d graph/vg/toy_graph_index.dist \
  -m graph/vg/toy_graph_index.shortread.withzip.min \
  -z graph/vg/toy_graph_index.shortread.zipcodes \
  --fragment-mean 260 \
  --fragment-stdev 10 \
  -N toy_reads \
  -f data-pangen-vc/toy_reads_R1.fastq.gz \
  -f data-pangen-vc/toy_reads_R2.fastq.gz \
  -o BAM \
  > graph/vg/toy_mapped.bam
```

Command parameters:

- `-o BAM`: write alignments projected to the graph reference path as BAM
- `>`: write the BAM output to `graph/vg/toy_mapped.bam`
- the graph index, fragment, sample, and read parameters are the same as in the GAM command above

Main outputs:

- `graph/vg/toy_mapped.gam`: graph alignment file
- `graph/vg/toy_mapped.gaf`: text alignment file useful for inspection
- `graph/vg/toy_mapped.bam`: alignment projected to the reference path
- `graph/qc/toy_mapping_stats.txt`: mapping summary

Inspect mapping statistics and GAF alignments:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ cat graph/qc/toy_mapping_stats.txt
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ head graph/vg/toy_mapped.gaf
```

```diff
! TASK: How many read pairs were in the dataset?
! TASK: Do the contaminant reads appear unmapped or low confidence?
! TASK: What graph paths appear in the GAF alignments?
```

### **5. Compute Read Support and Call Variants**

`vg giraffe` maps reads to the graph. To call read-supported variants, first convert the GAM alignments into graph support with `vg pack`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg pack \
  -x graph/vg/toy_graph_index.giraffe.gbz \
  -g graph/vg/toy_mapped.gam \
  -Q 5 \
  -o graph/vg/toy_mapped.pack
```

Command parameters:

- `-x graph/vg/toy_graph_index.giraffe.gbz`: graph on which read support is accumulated
- `-g graph/vg/toy_mapped.gam`: graph alignments produced by `vg giraffe`
- `-Q 5`: minimum mapping quality for reads to contribute support
- `-o graph/vg/toy_mapped.pack`: compressed read-support file for `vg call`

Call variants using the packed read support:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg call \
  graph/vg/toy_graph_index.giraffe.gbz \
  -k graph/vg/toy_mapped.pack \
  -a \
  -p 'REF#0#toy_chr1' \
  -s toy_reads \
  > graph/vg/toy_read_supported_calls.vcf
```

Command parameters:

- `graph/vg/toy_graph_index.giraffe.gbz`: graph to call variants from
- `-k graph/vg/toy_mapped.pack`: read support generated by `vg pack`
- `-a`: call all sites, including reference calls
- `-p 'REF#0#toy_chr1'`: graph path used as the VCF coordinate system
- `-s toy_reads`: sample name written in the VCF
- `>`: write read-supported variant calls to `graph/vg/toy_read_supported_calls.vcf`

Inspect the read-supported VCF:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep -v '^##' graph/vg/toy_read_supported_calls.vcf
```

```diff
! TASK: How many read-supported variant records were produced?
! TASK: Which reference path is used for VCF coordinates?
! TASK: Which calls are supported by reads from the toy diploid sample?
```

### **6. Compare Alignments with Truth Tables**

Use the read truth table to find reads that overlap each graph bubble:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ column -t data-pangen-vc/read_truth.tsv | less -S
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ column -t data-pangen-vc/graph_bubbles.tsv | less -S
```

Search for reads from one bubble class:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep 'bubble_ins1' data-pangen-vc/read_truth.tsv
```

Then look for those read names in the GAF output:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep 'READ_NAME' graph/vg/toy_mapped.gaf
```

Replace `READ_NAME` with an actual read name from `read_truth.tsv`.

```diff
! TASK: Choose one SNP, insertion, or deletion bubble
! TASK: Identify read pairs expected to cross that bubble
! TASK: Compare the expected truth table with the mapped GAF path
! TASK: Record the mapping quality for those reads
```

### **7. Deconstruct Graph Variation**

Deconstruct graph variation relative to the reference path:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg paths \
  -x graph/vg/toy_graph_index.giraffe.gbz \
  -L
```

Use the reference path `REF#0#toy_chr1`:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ vg deconstruct \
  graph/vg/toy_graph_index.giraffe.gbz \
  -p 'REF#0#toy_chr1' \
  > graph/vg/toy_path_variants.vcf
```

Inspect the VCF:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ grep -v '^##' graph/vg/toy_path_variants.vcf
```

```diff
! TASK: How many path-derived variant records were produced?
! TASK: Do they correspond to the four expected variable sites?
! TASK: How do path-derived calls differ from read-supported evidence?
```

### **8. Optional Graph Experiment**

Make a working copy of the GFA and remove the insertion branch:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-vc-working$ cp data-pangen-vc/toy_graph.gfa graph/vg/toy_graph.noins.gfa
```

Edit only `graph/vg/toy_graph.noins.gfa`:

- delete segment line `S 5 ...`
- delete links `4 -> 5` and `5 -> 6`
- remove `>5` from the affected W-line

Rebuild indexes from the edited working-copy graph and remap reads. Focus on read pairs labelled `bubble_ins1` in `data-pangen-vc/read_truth.tsv`.

```diff
! OPTIONAL TASK: Do insertion-spanning reads become clipped, mismatched, lower MAPQ, or differently placed?
! OPTIONAL TASK: What does this show about graph completeness and read mapping?
```

### **9. Final Report**

Prepare a short group report:

```diff
+ graph/report/pangenomics_variantcalling_report.md
```

Use this template:

```markdown
# Pangenomics variant-calling hackathon report

## Team
- names:

## Input graph
- GFA:
- reference sample:
- reference path:
- number of segments:
- number of links:
- number of W-line paths:

## Giraffe indexing
- command:
- GBZ:
- distance index:
- minimizer index:
- zipcodes file, if present:

## Read mapping
- read files:
- number of read pairs:
- mapping statistics:
- contaminant-read behavior:

## Read-supported variant calling
- pack file:
- read-supported VCF:
- number of records:
- coordinate/reference path:

## Graph variant interpretation
- expected variable sites:
- selected bubble:
- reads expected to overlap the bubble:
- observed GAF path and MAPQ:

## Graph deconstruction
- VCF file:
- number of records:
- coordinate/reference path:
- comparison with graph_bubbles.tsv:

## Interpretation
- what worked?
- what was uncertain?
- how did graph structure affect read placement?
- how is graph-based interpretation different from linear-reference variant calling?
```

```diff
! FINAL TASK: Explain how the graph paths, read alignments, and truth tables together support variant interpretation
```
