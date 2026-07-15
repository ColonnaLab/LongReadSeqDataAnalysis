[back to course home page ](../README.md)

## Advanced 3: Pangenomics hackathon

This practical is a three-hour small-group session on pangenome graph construction and graph inspection.

Graph-based variant calling continues in [the pangenomics variant-calling hackathon brief](10-pangenomics-variantcalling-hackathon.md).

The data used for this tutorial come from the [PGGB GitHub repository](https://github.com/pangenome/pggb/), specifically the `data/LPA` example dataset.

Dataset assumed on the server:

- source repository: `https://github.com/pangenome/pggb.git`
- source dataset inside the repository: `data/LPA`
- downloaded dataset: `/scratch/user1/data-pangenome/LPA`
- working directory: `/scratch/user1/pangenome-working`

```diff
+ You should document graph-building parameters
+ You should inspect graph structure before downstream interpretation
+ You should identify graph paths that could be used as coordinate systems later
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
! TASK: Which path could be treated as a reference coordinate system later?
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

Record one candidate reference path for the variant-calling lesson.

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

### **6. Final Report**

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

## Interpretation
- what worked?
- what was uncertain?
- what would we change in the graph build?
- which graph path could be used as a reference coordinate system later?
```

```diff
! FINAL TASK: Explain how graph construction choices could affect downstream variant interpretation
```
