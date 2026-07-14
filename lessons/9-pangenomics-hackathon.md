[back to course home page ](../README.md)

## Advanced 3: Pangenomics hackathon

This practical is a three-hour small-group session on pangenome graph construction, graph inspection, and graph-based read mapping for variant calling.

The data used for this tutorial come from the [PGGB GitHub repository](https://github.com/pangenome/pggb/), specifically the `data/LPA` example dataset.

Dataset assumed on the server:

- source repository: `https://github.com/pangenome/pggb.git`
- source dataset inside the repository: `data/LPA`
- downloaded dataset: `/scratch/user1/data-pangenome/LPA`
- optional paired reads: `/scratch/user1/data-pangenome/reads/sample_R1.fastq.gz`
- optional paired reads: `/scratch/user1/data-pangenome/reads/sample_R2.fastq.gz`
- working directory: `/scratch/user1/pangenome-working`

```diff
+ You should document graph-building parameters
+ You should inspect graph structure before variant interpretation
+ You should distinguish graph construction, read mapping, and variant calling
```

### **1. Prepare the Working Directory**

Create the working directory for analysis outputs:

```bash
user1@vm-corso-colonna:~$ mkdir -p /scratch/user1/pangenome-working
user1@vm-corso-colonna:~$ cd /scratch/user1/pangenome-working
```

Create link to data for pangenomic lesson 

```bash 

user1@vm-corso-colonna-11:/scratch/user1/pangenome-working$ ln -s /scratch/user1/data-pangenome/ data-pangenome 

```

Check that the dataset is present:

```bash
user1@vm-corso-colonna:~$ find data-pangenome/LPA -maxdepth 2 -type f | sort
```



Create output folders:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ mkdir -p graph/pggb graph/odgi graph/vg graph/qc graph/report
```

Check the input data:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ ls -lh /scratch/user1/data-pangenome/LPA
```

```diff
! TASK: Confirm that the assembly FASTA exists
! TASK: Confirm whether optional read files are available
! TASK: Check free disk space before graph construction
```

### **2. Inspect the Input Assemblies**

List the LPA files:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ ls -lh data-pangenome/LPA
```

The pangenome input FASTA is:

```diff
+ data-pangenome/LPA/LPA.fa.gz
```

Summarize the assembly FASTA:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ seqkit stats \
  data-pangenome/LPA/LPA.fa.gz \
  > graph/qc/assemblies.seqkit.stats.txt
```

Inspect sequence names. PGGB-compatible pangenome FASTA headers commonly follow the [PanSN naming specification](https://github.com/pangenome/PanSN-spec):

```text
sample#haplotype#contig
```

This naming helps graph tools distinguish samples, haplotypes, and contigs when they embed paths in the graph.

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ seqkit seq \
  -n \
  data-pangenome/LPA/LPA.fa.gz \
  > graph/qc/assembly_paths.txt

user1@vm-corso-colonna:/scratch/user1/pangenome-working$ head graph/qc/assembly_paths.txt
```

Count input sequences:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ wc -l graph/qc/assembly_paths.txt
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ wc -l data-pangenome/LPA/LPA.sample.list
```

```diff
! TASK: How many sequences are in the input FASTA?
! TASK: Do the names identify sample, haplotype, and contig?
! TASK: Which path should be treated as the reference for VCF interpretation?
```

### **3. Build a Pangenome Graph with PGGB**

Reference paper: [PGGB: the PanGenome Graph Builder](https://www.nature.com/articles/s41592-024-02430-3).

For a first classroom run, use conservative settings and a modest thread count.

Set `-n` to the number of input genomes or haplotypes. For the LPA classroom dataset, use the number of records in `data-pangenome/LPA/LPA.sample.list`.

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ pggb \
  -i data-pangenome/LPA/LPA.fa.gz \
  -o graph/pggb \
  -n NUMBER_OF_SAMPLES \
  -t 4 \
  -p 95 \
  -s 5k
```

Command parameters:

- `-i data-pangenome/LPA/LPA.fa.gz`: input FASTA containing the assembled sequences to combine into a pangenome graph
- `-o graph/pggb`: output directory for PGGB results
- `-n NUMBER_OF_SAMPLES`: number of input genomes or haplotypes; for this dataset, use the count from `LPA.sample.list`
- `-t 4`: number of CPU threads
- `-p 95`: minimum pairwise identity percentage used during all-vs-all alignment
- `-s 5k`: segment length used for graph construction; smaller values can represent finer variation but may increase graph complexity

Inspect the output folder 
```bash
user1@vm-corso-colonna-11:/scratch/user1/pangenome-working$ ls graph/pggb/
```

The PGGB output includes `.paf` files from `wfmash`. PAF is the [Pairwise mApping Format](https://github.com/lh3/miniasm/blob/master/PAF.md), a tab-delimited format for sequence alignments.


```diff
! TASK: Did PGGB finish successfully?
! TASK: What GFA file was produced?
! TASK: Which parameters would you change if the graph looked too fragmented or too collapsed?
```

### **4. Build and Inspect an ODGI Graph**

Documentation: [ODGI documentation](https://odgi.readthedocs.io/en/latest/index.html).

Convert the GFA graph into ODGI format:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ odgi build \
  -g graph/pggb/*.seqwish.gfa \
  -o graph/odgi/pangenome.og \
  -t 4 \
  -P
```

Command parameters:

- `-g graph/pggb/*.seqwish.gfa`: input GFA graph produced by PGGB
- `-o graph/odgi/pangenome.og`: output ODGI-format graph
- `-t 4`: number of CPU threads
- `-P`: show progress while building the ODGI graph

Summarize graph dimensions:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ odgi stats \
  -i graph/odgi/pangenome.og \
  -S \
  > graph/qc/pangenome.odgi.stats.txt

user1@vm-corso-colonna:/scratch/user1/pangenome-working$ cat graph/qc/pangenome.odgi.stats.txt
```

Command parameters:

- `-i graph/odgi/pangenome.og`: input ODGI graph
- `-S`: print a summary of graph size, including nucleotides, nodes, edges, paths, and steps
- `>`: redirect the summary into a text file

List graph paths:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ odgi paths \
  -i graph/odgi/pangenome.og \
  -L \
  > graph/qc/pangenome.paths.txt

user1@vm-corso-colonna:/scratch/user1/pangenome-working$ head graph/qc/pangenome.paths.txt
```

Command parameters:

- `-i graph/odgi/pangenome.og`: input ODGI graph
- `-L`: list all embedded paths in the graph

```diff
! TASK: How many graph paths are present?
! TASK: Do graph paths match the input assembly names?
! TASK: Which path will you use as a reference coordinate system?
```

In the variant-calling commands below, replace `ACTUAL_PATH_NAME` with one path from the graph path list.

### **5. Sort and Visualize the Graph**

Sort the graph to make visualization easier:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ odgi sort \
  -i graph/odgi/pangenome.og \
  -o graph/odgi/pangenome.sorted.og \
  -P
```

Command parameters:

- `-i graph/odgi/pangenome.og`: input ODGI graph
- `-o graph/odgi/pangenome.sorted.og`: output sorted ODGI graph
- `-P`: show progress while sorting

Make a 1D visualization:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ odgi viz \
  -i graph/odgi/pangenome.sorted.og \
  -o graph/odgi/pangenome.sorted.png \
  -x 1600 \
  -y 600
```

Command parameters:

- `-i graph/odgi/pangenome.sorted.og`: input sorted ODGI graph
- `-o graph/odgi/pangenome.sorted.png`: output PNG image
- `-x 1600`: image width in pixels
- `-y 600`: image height in pixels

```diff
! TASK: Copy or open the PNG visualization
! TASK: Identify regions where paths diverge
! TASK: Decide whether the graph appears mostly simple or highly complex
```

### **6. Build VG Giraffe Indexes**

Use `vg autoindex` to prepare indexes for graph read mapping:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg autoindex \
  --workflow giraffe \
  -g graph/pggb/*.seqwish.gfa \
  -p graph/vg/pangenome \
  -t 4
```

Command parameters:

- `--workflow giraffe`: build the graph indexes needed by `vg giraffe`
- `-g graph/pggb/*.seqwish.gfa`: input GFA graph
- `-p graph/vg/pangenome`: prefix for output index files
- `-t 4`: number of CPU threads

Check generated files:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ ls -lh graph/vg
```

```diff
! TASK: Was a GBZ graph produced?
! TASK: Were minimizer and distance indexes produced?
! TASK: Why does vg giraffe need graph indexes before mapping?
```

### **7. Deconstruct Graph Variation from Paths**

Infer variants directly from graph paths using the selected reference path.

First list the path names available in the VG graph:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg paths \
  -x graph/vg/pangenome.giraffe.gbz \
  -L \
  > graph/qc/vg.paths.txt

user1@vm-corso-colonna:/scratch/user1/pangenome-working$ head graph/qc/vg.paths.txt
```

Choose one path from `graph/qc/vg.paths.txt` and use it after `-p`. Do not copy `ACTUAL_PATH_NAME` literally.

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg deconstruct \
  graph/vg/pangenome.giraffe.gbz \
  -p 'ACTUAL_PATH_NAME' \
  > graph/vg/pangenome.path_variants.vcf
```

Command parameters:

- `graph/vg/pangenome.giraffe.gbz`: graph containing embedded assembly paths
- `-p 'ACTUAL_PATH_NAME'`: graph path to use as the reference coordinate system
- `>`: write path-derived variants to a VCF file

```diff
! TASK: How many path-derived variant records were produced?
! TASK: Which reference path or coordinate system appears in the VCF?
! TASK: These calls come from assembled genomes. How is that different from read-supported calls?
```

### **8. Final Report**

Prepare a short group report:

```diff
+ graph/report/pangenomics_report.md
```

Use this template:

```markdown
# Pangenomics hackathon report

## Team
- names:

## Input data
- assembly FASTA:
- number of input sequences:
- chosen reference path:

## PGGB graph build
- command:
- number of genomes or haplotypes:
- identity threshold:
- segment length:
- GFA:

## ODGI graph inspection
- graph statistics:
- number of paths:
- main observations from visualization:

## Graph/path-derived variant calling
- VCF file:
- number of records:
- coordinate/reference path:
- interpretation limits:

## Optional read-supported analysis
- read files:
- VG Giraffe mapping summary:
- read-supported VCF:
- read-supported variant count:

## Interpretation
- what worked?
- what was uncertain?
- what would we change in the graph build?
```

```diff
! FINAL TASK: Explain how pangenome graph variant calling differs from linear-reference variant calling
```

### **9. Optional: Map Reads to the Pangenome Graph**

The LPA dataset provides assemblies for graph construction. If matching read files are available for the course, use their paths in `vg giraffe`.

Single-end read example:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg giraffe \
  -Z graph/vg/pangenome.giraffe.gbz \
  -m graph/vg/pangenome.min \
  -d graph/vg/pangenome.dist \
  -f data-pangenome/reads/sample.fastq.gz \
  -t 4 \
  > graph/vg/sample.gam
```

Command parameters:

- `-Z graph/vg/pangenome.giraffe.gbz`: indexed graph used for mapping
- `-m graph/vg/pangenome.min`: minimizer index used to seed read mapping
- `-d graph/vg/pangenome.dist`: distance index used to evaluate graph distances between seeds
- `-f data-pangenome/reads/sample.fastq.gz`: single-end input FASTQ
- `-t 4`: number of CPU threads
- `>`: write alignments to a GAM file

Paired-end read example:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg giraffe \
  -Z graph/vg/pangenome.giraffe.gbz \
  -m graph/vg/pangenome.min \
  -d graph/vg/pangenome.dist \
  -f data-pangenome/reads/sample_R1.fastq.gz \
  -f data-pangenome/reads/sample_R2.fastq.gz \
  -t 4 \
  > graph/vg/sample.gam
```

For paired-end data, use two `-f` options: the first for read 1 and the second for read 2.

Basic alignment QC:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg stats \
  -a \
  graph/vg/sample.gam \
  > graph/qc/sample.gam.stats.txt

user1@vm-corso-colonna:/scratch/user1/pangenome-working$ cat graph/qc/sample.gam.stats.txt
```

Command parameters:

- `-a`: report alignment statistics for the GAM file
- `graph/vg/sample.gam`: graph alignment file produced by `vg giraffe`

```diff
! OPTIONAL TASK: Did reads map to the graph?
! OPTIONAL TASK: Are there unmapped reads?
! OPTIONAL TASK: What could cause poor graph mapping?
```

### **10. Optional: Compute Read Support and Call Variants**

Compute graph read support:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg pack \
  -x graph/vg/pangenome.giraffe.gbz \
  -g graph/vg/sample.gam \
  -Q 5 \
  -o graph/vg/sample.pack
```

Command parameters:

- `-x graph/vg/pangenome.giraffe.gbz`: graph to accumulate read support on
- `-g graph/vg/sample.gam`: read alignments from `vg giraffe`
- `-Q 5`: minimum mapping quality for reads to contribute support
- `-o graph/vg/sample.pack`: output support file used by `vg call`

Call variants from read support:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ vg call \
  graph/vg/pangenome.giraffe.gbz \
  -k graph/vg/sample.pack \
  -a \
  -p 'ACTUAL_PATH_NAME' \
  -s sample \
  > graph/vg/sample.graph_calls.vcf
```

Command parameters:

- `graph/vg/pangenome.giraffe.gbz`: graph to call variants from
- `-k graph/vg/sample.pack`: read support generated by `vg pack`
- `-a`: call all snarls, including reference calls
- `-p 'ACTUAL_PATH_NAME'`: graph path to use as the VCF reference coordinate system
- `-s sample`: sample name to write in the VCF
- `>`: write variant calls to a VCF file

Inspect the VCF:

```bash
user1@vm-corso-colonna:/scratch/user1/pangenome-working$ grep -v '^##' graph/vg/sample.graph_calls.vcf | head
```

```diff
! OPTIONAL TASK: How many read-supported variant records were produced?
! OPTIONAL TASK: Which reference path or coordinate system appears in the VCF?
! OPTIONAL TASK: Compare read-supported variants with path-derived variants
```
