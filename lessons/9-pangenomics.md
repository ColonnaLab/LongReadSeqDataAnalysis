[back to course home page ](../README.md)

## Advanced 3: Pangenomics

This three-hour session introduces pangenome graph analysis. The goal is to build a small pangenome graph from genome assemblies, then inspect and manipulate it with graph-aware tools.

The practical challenge is in [the pangenomics hackathon brief](9-pangenomics-hackathon.md).
Graph-based variant calling is covered in [the pangenomics variant-calling lesson](10-pangenomics-variantcalling.md).

Dataset assumed for the practical:

- source repository: `https://github.com/pangenome/pggb.git`
- source dataset inside the repository: `data/LPA`
- server data folder after download: `/scratch/user1/data-pangenome/LPA`
- input assembly FASTA: selected from the LPA folder during the practical
- reference path name for downstream variant interpretation: decided from the graph path list

```diff
+ We will build a pangenome graph from multiple assemblies
+ We will inspect graph paths, size, and topology with odgi
+ We will decide which graph path could be used as a reference coordinate system later
```

### **1. What Is a Pangenome?**

A single reference genome represents one linear sequence. A pangenome represents the sequence diversity of multiple genomes from the same species or population.

A pangenome can describe:

- core sequence shared by all genomes
- accessory sequence present in only some genomes
- SNPs and small indels
- structural variants
- inversions, duplications, insertions, and deletions
- strain-specific or haplotype-specific regions

For bacterial isolates, a pangenome can help represent genes or regions that are absent from the reference. For diploid organisms, a pangenome can help represent haplotype diversity and reduce reference bias.

```diff
! DISCUSSION: What information is lost when all samples are compared only to one linear reference?
! DISCUSSION: Why might pangenome graphs be useful after we have assembled multiple genomes?
```

### **2. What Is a Pangenome Graph?**

A pangenome graph stores sequence as nodes and relationships between sequences as edges. Genome assemblies are embedded as paths through the graph.

In simple terms:

- nodes contain DNA sequence
- edges connect nodes
- paths represent genomes, contigs, haplotypes, or references
- bubbles often represent alternative alleles or structural differences

```diff
+ shared sequence:  ACGT
+ alternative 1:    ACGT-GGA-TT
+ alternative 2:    ACGT-CTA-TT
+ graph idea:       shared nodes with a bubble for GGA versus CTA
```

Graph coordinates are more complex than linear reference coordinates. For variant reports, we usually choose one path as the reference path and project graph variation relative to that path.

### **3. Building Pangenome Graphs with PGGB**

[`pggb`](https://github.com/pangenome/pggb) builds pangenome graphs from a FASTA file containing multiple assemblies or haplotypes.

The input FASTA should contain related sequences. Sequence names should be informative. A common convention is PanSN-style naming:

```text
sample#haplotype#contig
```

Examples:

```text
K12#1#chromosome
isolateA#1#chromosome
isolateB#1#chromosome
```

Important PGGB parameters:

- `-i`: input FASTA
- `-o`: output directory
- `-n`: number of input genomes or haplotypes, optional when names allow automatic detection
- `-t`: number of threads
- `-p`: minimum pairwise identity used during graph construction
- `-s`: segment length used during graph construction
- `-V`: optional VCF generation relative to a named reference path

Example command:

```bash
pggb \
  -i data-pangenome/LPA/LPA.fa.gz \
  -o graph/pggb \
  -n NUMBER_OF_SAMPLES \
  -t 8 \
  -p 95 \
  -s 5k
```

For the LPA run, the graph output used below is the `seqwish` GFA file:

```text
graph/pggb/*.seqwish.gfa
```

### **4. Inspecting and Manipulating Graphs with ODGI**

[`odgi`](https://github.com/pangenome/odgi) is used to build, inspect, sort, subset, and visualize pangenome graphs.

Common ODGI commands:

- `odgi build`: convert GFA into ODGI format
- `odgi stats`: summarize graph size and structure
- `odgi paths`: list or extract embedded graph paths
- `odgi sort`: reorder graph nodes for easier visualization and analysis
- `odgi viz`: make a 1D PNG visualization
- `odgi view`: convert ODGI back to GFA

Example:

```bash
odgi build \
  -g graph/pggb/*.seqwish.gfa \
  -o graph/odgi/pangenome.og \
  -t 4 \
  -P
```

Summarize the graph:

```bash
odgi stats \
  -i graph/odgi/pangenome.og \
  -S
```

List graph paths:

```bash
odgi paths \
  -i graph/odgi/pangenome.og \
  -L \
  > graph/odgi/paths.txt
```

Sort and visualize:

```bash
odgi sort \
  -i graph/odgi/pangenome.og \
  -o graph/odgi/pangenome.sorted.og \
  -P

odgi viz \
  -i graph/odgi/pangenome.sorted.og \
  -o graph/odgi/pangenome.sorted.png \
  -x 1600 \
  -y 600
```

### **5. Key Points**

- A pangenome represents diversity across multiple genomes, not only one reference sequence
- A graph stores shared and variable sequence in one structure
- PGGB builds pangenome graphs from assembled sequences
- ODGI is useful for graph inspection, manipulation, and visualization
- Downstream variant interpretation depends on the chosen reference path and graph construction choices
