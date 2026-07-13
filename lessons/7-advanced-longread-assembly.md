[back to course home page ](../README.md)

## Advanced 1: ONT long-read genome assembly

This advanced session is organized as a small-group hackathon. The goal is not to follow a fixed recipe, but to make analysis decisions, document them, and evaluate whether the results make biological and technical sense.

The practical challenge is in [the ONT assembly hackathon brief](7-advanced-longread-assembly-hackathon.md).

This session focuses on **Oxford Nanopore Technologies (ONT)** assembly. The dataset is the same MRSA KUN1163 ONT dataset used in lessons 4-6:

- reads: `data-longreads/fastq/DRR187567.fastq.gz`
- reference genome for comparison: `data-longreads/reference/KUN1163_reference.fasta`
- reference chromosome: `AP020324.1`
- reference plasmid: `AP020325.1`

```diff
+ We will assemble long reads into contigs
+ We will evaluate assembly size, contiguity, and completeness
+ We will compare the assembly to the published KUN1163 reference
+ We will decide whether the chromosome and plasmid were recovered
```

### **1. What Is Genome Assembly?**

Genome assembly reconstructs longer genomic sequences from sequencing reads.

- Reads are short or long fragments sampled from the genome
- An assembler looks for overlaps or graph relationships between reads
- The output is a set of contigs
- A complete bacterial assembly may contain:
  - one chromosome contig
  - one or more plasmid contigs
  - sometimes extra fragmented or contaminant contigs

```diff
+ Reads -> assembler -> contigs -> assembly evaluation
```

Unlike reference alignment, assembly does not start by placing every read on a reference genome. This makes assembly useful when the sample differs from the reference, contains plasmids, or has structural differences.

### **2. Why Long Reads Help Assembly**

Long reads are useful for assembly because they can span:

- repetitive regions
- insertion sequences
- duplicated genes
- structural variants
- plasmid repeats
- low-complexity regions

Short reads often break assemblies at repeats because the reads are too short to connect unique sequence on both sides. Long reads can sometimes bridge the repeat and produce a more complete contig.

```diff
! DISCUSSION: Why can a long read spanning a repeat help connect two unique regions?
! DISCUSSION: Why can low-quality long reads still create assembly errors?
```

### **3. Haploid, Diploid, Collapsed, and Phased Assemblies**

The KUN1163 dataset in this session is bacterial, so it is treated as a **haploid** genome: each genomic region is expected to be represented once, except for repeats, plasmids, duplicated genes, or contamination.

Many animal and plant genomes are **diploid**. A diploid individual carries two copies of each chromosome, one inherited from each parent. These two chromosome copies are called **haplotypes**. The haplotypes are usually similar, but they can differ by:

- single-nucleotide variants
- small insertions and deletions
- structural variants
- copy-number differences
- larger rearrangements

These differences matter during assembly because the assembler must decide whether two similar sequences are:

- allelic versions of the same locus, one on each haplotype
- duplicated regions that should both be present in the final assembly
- sequencing errors that should be corrected into one consensus sequence

#### Collapsed assemblies

A **collapsed assembly** represents the two haplotypes as one consensus sequence for most regions.

```diff
+ haplotype A:  ATGCTACGTTAG
+ haplotype B:  ATGCTATGTTAG
+ collapsed:    ATGCTA?GTTAG
```

In practice, the consensus base is chosen from the read evidence, so the output does not usually contain `?`. The important idea is that the assembly has one contig for a region where the organism actually has two haplotypic copies.

Collapsed assemblies are common when:

- the two haplotypes are very similar
- read length or coverage is not sufficient to separate them
- the assembler is configured to produce a primary consensus assembly

Collapsed assemblies can be useful because they are simpler to analyze, but they can hide heterozygous variation and make some regions look artificially homozygous.

#### Phased haplotype assemblies

A **phased assembly** separates the two haplotypes and attempts to reconstruct each chromosome copy independently.

```diff
+ haplotype A assembly:  ATGCTACGTTAG
+ haplotype B assembly:  ATGCTATGTTAG
```

This is also called a **haplotype-resolved assembly**. It is especially useful when the goal is to study heterozygosity, allele-specific genes, structural variation, or inherited chromosome blocks.

Phased assemblies usually need stronger evidence than collapsed assemblies, such as:

- long reads that span multiple heterozygous variants
- high coverage
- low sequencing error
- parental data, Hi-C, linked reads, or other long-range information
- assemblers designed for haplotype resolution, such as `hifiasm` for PacBio HiFi data

#### Why this matters for evaluation

For a haploid bacterial genome, a good assembly size should usually be close to the expected genome size. For a diploid genome, interpretation depends on the assembly type:

- a collapsed primary assembly is expected to be close to one haploid genome size
- a fully duplicated haplotype-resolved assembly can approach two haploid genome copies
- unresolved repeats or duplicated haplotypes can inflate assembly size
- over-collapsing can reduce assembly size and hide real variation

```diff
! DISCUSSION: If a diploid genome has an expected haploid size of 3 Gb, should a 6 Gb assembly always be considered wrong?
! DISCUSSION: Why might a collapsed assembly miss biologically important variants?
```

### **4. Assembly Tools**

Common long-read assembly tools include:

- [`Flye`](https://github.com/mikolmogorov/Flye)
  - widely used for noisy long reads
  - supports Oxford Nanopore and PacBio reads
  - reports assembly information and possible circular contigs
- [`Raven`](https://github.com/lbcb-sci/raven)
  - fast long-read assembler
  - useful for comparison
- [`miniasm`](https://github.com/lh3/miniasm)
  - very fast overlap-based assembler
  - usually requires polishing because output is not consensus-corrected

In this hackathon, `Flye` is the recommended first assembler.

```diff
+ For ONT reads, a typical Flye mode is --nano-raw
+ For a bacterial genome, the expected genome size is approximately a few Mb
```

### **5. Assembly Evaluation**

Assembly is not finished when the assembler produces a FASTA file. We need to evaluate it.

Useful questions:

- How many contigs were produced?
- What is the total assembly length?
- What is the longest contig?
- What is the N50?
- Is the assembly close to the expected genome size?
- Was the plasmid recovered?
- Are contigs circular or linear?
- Does the assembly align cleanly to the reference?
- Are there duplicated, missing, or rearranged regions?

Useful tools:

- `seqkit stats`: basic FASTA statistics
- `QUAST`: assembly quality and reference comparison
- `minimap2`: assembly-to-reference alignment
- `samtools`: alignment sorting and indexing
- `IGV`: visual inspection of assembly-reference alignments
- `Bandage`: graph visualization, if an assembly graph is available

### **6. Reference Comparison**

Because this dataset has a published KUN1163 reference, we can compare the assembly to:

- chromosome `AP020324.1`
- plasmid `AP020325.1`

This comparison helps answer:

- Did the assembly recover both reference sequences?
- Is the chromosome assembled as one contig or multiple contigs?
- Is the plasmid present?
- Are there large disagreements between assembly and reference?

Reference comparison is useful, but it should not be treated as the only truth. If the sample and reference differ, a real biological difference can look like a disagreement.

### **7. Polishing**

Polishing attempts to correct base-level errors in an assembly.

- Long-read-only polishing uses the same long reads
- Hybrid polishing can use short reads if available
- Polishing can improve small-variant accuracy
- Polishing can also introduce problems if reads are misaligned or contaminated

For this session, polishing is optional. The core challenge is to produce and evaluate a first long-read assembly.

### **8. Expected Deliverable**

Each group should produce a short report:

```markdown
# Long-read assembly hackathon report

## Input data
- reads:
- estimated genome size:
- estimated coverage:

## Assembly command
- assembler:
- parameters:

## Assembly summary
- number of contigs:
- total assembly length:
- longest contig:
- N50:

## Reference comparison
- chromosome recovered?
- plasmid recovered?
- major differences:

## Interpretation
- what worked?
- what was uncertain?
- what would you try next?
```

### **9. Key Points**

- Genome assembly reconstructs contigs from reads
- Long reads can resolve repeats and structural regions better than short reads
- Haploid, collapsed diploid, and phased haplotype assemblies represent genome copies differently
- Assembly quality must be evaluated after assembly
- A good bacterial assembly should be close to the expected genome size
- Reference comparison helps interpret completeness and structural agreement
- A hackathon workflow should document decisions, not only commands
