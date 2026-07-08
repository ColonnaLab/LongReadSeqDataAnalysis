[back to course home page ](../README.md)

## Long reads alignment

This page contains the conceptual notes for the long reads alignment lesson. The practical commands using the MRSA KUN1163 dataset are in [the alignment exercise](5-longreads-alignment-exercise.md).

The dataset is the same used in lesson 4:

- Oxford Nanopore reads from MRSA strain KUN1163
- run accession: `DRR187567`
- reference chromosome: `AP020324.1`
- reference plasmid: `AP020325.1`

```diff
+ We will align long reads to a reference genome
+ We will use minimap2 with an Oxford Nanopore preset
+ We will produce sorted and indexed BAM files
+ We will evaluate alignment and coverage statistics
```

### **1. What Does Alignment Mean?**

- Alignment means finding where each sequencing read fits best on a reference genome
- The read is the query sequence
- The genome is the reference sequence
- The aligner compares the read to the reference and reports:
  - reference sequence name
  - alignment position
  - mapping quality
  - matches and mismatches
  - insertions and deletions
  - soft-clipped or split parts of the read



### **2. Why Long Reads Need Specific Aligners**

- Short-read aligners such as `bwa` are optimized for short, accurate Illumina reads
- Long reads are much longer and can contain different error profiles
- Oxford Nanopore reads may include:
  - more variable quality
  - insertions and deletions
  - very different read lengths
  - long reads that span repeats or structural variants
- Long-read aligners must handle:
  - long sequences
  - higher error rates
  - large gaps
  - split alignments
  - structural rearrangements

```diff
+ Short reads: many small pieces
+ Long reads: fewer but much longer pieces
+ The alignment strategy must change
```

### **3. minimap2**

- `minimap2` is a widely used aligner for long reads
- It supports different presets for different sequencing technologies
- Common presets:
  - `map-ont`: Oxford Nanopore genomic reads
  - `map-pb`: PacBio continuous long reads
  - `map-hifi`: PacBio HiFi reads
- In this lesson we use Oxford Nanopore reads, so the practical command uses:

```diff
+ minimap2 -ax map-ont
```

- Meaning of the main options:
  - `-a`: output SAM format
  - `-x map-ont`: use settings appropriate for Oxford Nanopore genomic reads

### **4. SAM and BAM Files**

- Alignment output is commonly stored as SAM or BAM
- SAM:
  - text format
  - human-readable
  - large file size
  - useful for learning and inspection
- BAM:
  - binary compressed version of SAM
  - smaller and faster to process
  - used by most downstream tools
- A sorted BAM file is ordered by reference genome coordinates
- A BAM index file (`.bai`) allows tools to access specific genome regions quickly

```diff
+ Long reads FASTQ -> minimap2 -> SAM/BAM -> sorted BAM -> indexed BAM
```

### **5. Important SAM/BAM Concepts**

- Reference name:
  - tells us whether a read aligned to the chromosome or plasmid
- Position:
  - tells us where the read starts on the reference
- MAPQ:
  - mapping quality
  - estimates how confidently the read was placed
- CIGAR string:
  - compact description of matches, insertions, deletions, clipping, and skipped regions
- FLAG:
  - encoded information about the read and its alignment status
- Primary and secondary alignments:
  - long reads may align well to more than one location
- Supplementary alignments:
  - one read may align in pieces, which is important for structural variants

```diff
! DISCUSSION: Why might one long read align in two different pieces?
! DISCUSSION: Why are supplementary alignments useful for structural variant detection?
```

### **6. Alignment Summary Metrics**

- Number of total reads
- Number of mapped reads
- Percentage of mapped reads
- Number of unmapped reads
- Reads aligned to each reference sequence
- Mean coverage of each reference sequence
- Fraction of the reference covered
- Coverage differences between chromosome and plasmid

Useful `samtools` commands:

```diff
+ samtools flagstat: summarizes mapped and unmapped reads
+ samtools idxstats: reports read counts per reference sequence
+ samtools coverage: reports coverage statistics
+ samtools tview: visualizes alignments in the terminal
```

### **7. Long Reads Alignment Compared with Short Reads**

- Long reads often produce fewer total reads than short-read datasets
- Each long read covers a larger region of the genome
- Long reads can span repetitive regions that short reads cannot resolve
- Long reads can reveal large insertions, deletions, inversions, and rearrangements
- Long-read alignments may contain more indels
- The CIGAR strings can be longer and more complex
- Coverage can be uneven depending on DNA quality, library preparation, and sequencing run

```diff
! EXERCISE: Compare short-read and long-read alignment outputs
! Which dataset has more reads?
! Which dataset has longer alignments?
! Which one is more informative for large structural variants?
```

### **8. Practical Workflow Overview**

- Prepare a working folder
- Check that the reads and reference are available
- Inspect the reference genome
- Run `minimap2` with the `map-ont` preset
- Convert SAM to BAM
- Sort and index BAM
- Generate alignment summaries with `samtools`
- Estimate coverage
- Inspect alignments in the terminal
- Save the summary files for the variant calling lesson

```diff
+ The practical commands are in: lessons/5-longreads-alignment-exercise.md
```

### **9. Key Points**

- Long reads need aligners designed for long, error-prone sequences
- `minimap2` is commonly used for Oxford Nanopore and PacBio data
- The correct preset matters
- SAM is readable but large
- BAM is compressed and suitable for downstream analyses
- Sorted and indexed BAM files are required by many tools
- Alignment and coverage summaries help decide whether the data are suitable for variant calling

