[back to course home page ](../README.md)
## 4. Long-read sequencing: molecule context

Long-read technologies help recover information that short reads fragment away: structural variation, repeats, haplotypes, isoforms, and epigenetic marks.

### PacBio HiFi sequencing

[PacBio HiFi sequencing](https://www.pacb.com/technology/hifi-sequencing/) combines long-read sequencing with high per-read accuracy. During library preparation, adapters are added to both ends of a double-stranded DNA fragment, creating a circular SMRTbell template. A DNA polymerase can then move around the same template multiple times during single-molecule real-time sequencing.

Each pass produces a subread of the insert. Because the passes observe the same original DNA molecule repeatedly, their base calls can be compared and combined into a circular consensus sequence (CCS). Random errors are less likely to occur at the same position in every pass, so consensus generation produces one long, highly accurate HiFi read.

The number of complete passes depends partly on insert length and polymerase read length. Shorter inserts may be read more times, while longer inserts may receive fewer passes. HiFi library design therefore balances read length, consensus accuracy, DNA quality, and the biological structures that need to be spanned.

```text
1. Linear double-stranded DNA fragment

   5' ---------------- DNA insert ---------------- 3'
   3' ---------------- DNA insert ---------------- 5'

2. Hairpin adapters create a circular SMRTbell template

                     DNA insert
                /------------------\
         adapter                    adapter
                \------------------/
                     DNA insert

3. Polymerase reads around the same molecule repeatedly

                pass 1  --->
             /------------------\
            |                    |
             \------------------/
                       <---  pass 2
                pass 3  --->   ...

4. Multiple subreads are combined

   pass 1:  A C G T T G C A C T G A
   pass 2:  A C G T T G C A C T G A
   pass 3:  A C G T C G C A C T G A
   pass 4:  A C G T T G C A C T G A
            | | | | | | | | | | | |
   HiFi:    A C G T T G C A C T G A

            repeated observations -> accurate consensus read
```

*A HiFi read is a consensus generated from repeated observations of the same circularized DNA molecule, not a consensus across different molecules or individuals.*

| | |
|---|---|
| **Strengths** | &bull; Long reads with high per-read accuracy<br>&bull; Strong for de novo assembly and structural variant discovery<br>&bull; Useful for kilobase-scale phasing<br>&bull; Can support DNA methylation detection<br>&bull; Well suited to haplotype-resolved assembly |
| **Limitations for pangenomics** | &bull; Requires high-molecular-weight DNA<br>&bull; Throughput and cost depend on instrument, chemistry, genome size, and coverage<br>&bull; Reads are usually shorter than ultra-long nanopore reads<br>&bull; Library preparation and DNA quality can be limiting |
| **Pangenomics role** | &bull; Build high-quality assemblies for graph construction<br>&bull; Discover structural variants precisely<br>&bull; Resolve paralogs, duplications, and complex loci<br>&bull; Produce haplotype-resolved graph paths |

: {.platform-summary tbl-colwidths="[24,76]"}

### Oxford Nanopore sequencing

[Oxford Nanopore sequencing](https://nanoporetech.com/platform/technology) measures changes in electrical current as a native DNA or RNA molecule passes through a protein nanopore embedded in a membrane. A motor protein controls movement through the pore so that successive groups of bases influence the current signal over time.

The instrument records this changing current as raw signal rather than nucleotide letters. A basecalling model interprets patterns in the signal and converts them into a sequence. Because neighboring bases affect the signal together, basecalling is a context-dependent inference rather than a direct observation of one isolated base at a time.

The molecule can be read continuously until it exits the pore or sequencing is interrupted. Read length is therefore strongly influenced by the length and integrity of the input molecule. Native-molecule sequencing can also retain signal associated with base modifications, allowing sequence and modification information to be inferred from the same experiment.

```text
1. A motor protein guides a native molecule through a nanopore

   DNA or RNA -----> [ motor ] -----> | nanopore | ----->
                                         membrane

2. Bases in and near the pore alter the electrical current

   molecule:   A C G T T G C A C T G A
               |---|---|---|---|---|---|
   current:    _/\___/\_/\____/\_/\____
                 signal changes over time

3. A basecaller converts signal patterns into sequence

   raw signal  ->  basecalling model  ->  A C G T T G C A C T G A

4. One molecule produces one continuous long read

   molecule:   |----------------------------------------------|
   read:       ACGTTGCA...----------------------------------->
```

*A nanopore read is inferred from the electrical signal generated by one molecule passing through a pore. Unlike HiFi sequencing, high accuracy is not produced by repeated passes around a circular template.*

| | |
|---|---|
| **Strengths** | &bull; Very long and ultra-long reads<br>&bull; Real-time sequencing and analysis<br>&bull; Portable and flexible instruments<br>&bull; Direct DNA or RNA sequencing with modification information<br>&bull; Spans large repeats and complex variants |
| **Limitations for pangenomics** | &bull; Accuracy depends on chemistry, basecaller, model, and sample quality<br>&bull; Ultra-high-molecular-weight DNA extraction is technically demanding<br>&bull; Raw signal data can be large<br>&bull; Calling and polishing require suitable tools and QC |
| **Pangenomics role** | &bull; Provide ultra-long-range continuity<br>&bull; Resolve large repeats and variants<br>&bull; Support rapid sequencing near the sample<br>&bull; Detect DNA or RNA modifications directly<br>&bull; Complement HiFi when length is the priority |

: {.platform-summary tbl-colwidths="[24,76]"}

### Approximate span of one sequence read

```text
                 100 bp       1 kb        10 kb       100 kb       1 Mb
                   |-----------|-----------|-----------|-----------|

Short read         |--|
                   tens to hundreds of bases

PacBio HiFi                    |-----------|
                               roughly 10-25 kb is common

Nanopore                       |-----------------------|
                               often tens of kb; highly variable

Ultra-long nanopore                            |-------------------|
                                               reads can exceed 100 kb

Longer reads preserve more context across:

   variants ---- repeats ---- genes ---- structural variants ---- haplotypes
```

*Sequence technologies differ not only in accuracy and throughput but also in the amount of continuous molecule context retained in each read. Short reads provide deep local evidence, HiFi reads combine long-range context with high accuracy, and nanopore sequencing can extend into ultra-long molecule ranges. The horizontal scale is approximate and logarithmic; actual read-length distributions depend on platform, library preparation, DNA quality, size selection, chemistry, and run settings.*

---

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

```diff
+ Reference:   ...ATCGATCGATCGATCG...
+ Read:           ATCGATCGTTTGATCG
+ Alignment:      matches + mismatches + possible insertions/deletions
```

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

