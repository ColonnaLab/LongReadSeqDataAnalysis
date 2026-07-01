[back to course home page ](../README.md)

## Long reads technologies and quality control

This page contains a bullet-point outline for the first long reads lesson. We will use it as notes to discuss the principles of PacBio HiFi and Oxford Nanopore sequencing, the main differences with short reads, the most common long reads data formats, and basic quality control.

```diff
+ We will move from short reads workflows to long reads workflows
+ Focus of this lesson: technologies, file formats, and quality control
+ Full practical commands will be added after revising the outline
```

### **1. Why Long Reads?**

- Long reads are sequencing reads that are much longer than standard Illumina short reads
- Typical read lengths:
  - Short reads: usually 100-300 bp
  - PacBio HiFi reads: often 10-25 kbp, with high per-base accuracy
  - Oxford Nanopore reads: often 10-100 kbp, with possible ultra-long reads above 100 kbp
- Long reads help resolve genomic regions that are difficult for short reads:
  - Repetitive regions
  - Structural variants
  - Large insertions and deletions
  - Complex rearrangements
  - Full-length transcripts
  - Haplotype phasing
- Main tradeoff:
  - Short reads: very accurate, inexpensive, but short
  - Long reads: much longer, more informative for structure, but technology-specific error profiles and data handling

```diff
! DISCUSSION: Why can a 150 bp read be hard to place uniquely in a repetitive genome?
! DISCUSSION: Which variant types are difficult to detect using only short reads?
```

### **2. PacBio HiFi Sequencing**

- PacBio sequencing reads DNA molecules using Single Molecule Real-Time sequencing
- DNA polymerase incorporates bases while the instrument records signal in real time
- A circular DNA molecule is read multiple times to produce a consensus sequence
- HiFi reads are also called CCS reads: Circular Consensus Sequencing reads
- Key properties of PacBio HiFi reads:
  - Long reads with high accuracy
  - Relatively random error profile
  - Useful for small variants, structural variants, assemblies, and phasing
  - Usually delivered as FASTQ or BAM files
- Important concepts:
  - Subreads: raw passes from the sequencing molecule
  - CCS/HiFi reads: consensus reads generated from multiple passes
  - Read length distribution: not all reads have the same length
  - Per-read quality: each read has an estimated quality value

```diff
+ PacBio HiFi = long reads + high consensus accuracy
+ We usually work with processed HiFi reads, not raw instrument signals
```

### **3. Oxford Nanopore Sequencing**

- Oxford Nanopore sequencing reads DNA or RNA molecules as they pass through a nanopore
- Changes in electrical current are converted into bases by basecalling software
- Sequencing can produce very long reads because native DNA molecules are read directly
- Key properties of Oxford Nanopore reads:
  - Very flexible read length, including ultra-long reads
  - Portable instruments are available
  - Signal-level data can be retained for re-basecalling or modified base analysis
  - Accuracy depends on chemistry, flowcell, basecaller, model, and filtering
- Important concepts:
  - Raw signal: original electrical signal from the pore
  - Basecalling: conversion from signal to nucleotide sequence
  - Duplex reads: reads where both strands are used to improve accuracy
  - Modified bases: methylation and other DNA modifications may be detected from signal

```diff
+ Oxford Nanopore = direct molecule sequencing + signal data
+ Read length depends strongly on DNA quality and library preparation
```

### **4. Long Reads vs Short Reads**

- Read length:
  - Short reads produce many small fragments
  - Long reads can span repeats and large variants
- Error profile:
  - Short reads usually have low per-base error rates
  - Long reads have technology-specific errors and quality distributions
- Alignment:
  - Short reads are usually aligned with tools such as BWA
  - Long reads are usually aligned with tools such as minimap2
- Variant calling:
  - Short reads are strong for SNVs and small indels in accessible regions
  - Long reads are strong for structural variants and phasing
  - HiFi reads can also perform very well for SNVs and indels
- Assembly:
  - Short reads often produce fragmented assemblies
  - Long reads make genome assembly more contiguous
- File size and compute:
  - Long reads workflows may require large files and substantial compute resources

```diff
! EXERCISE: Make a table comparing short reads, PacBio HiFi, and Oxford Nanopore
! Include: read length, accuracy, file formats, strengths, limitations
```

### **5. Long Reads File Formats**

- FASTQ:
  - Contains sequence and quality scores
  - Same basic four-line structure as short reads FASTQ
  - Long reads FASTQ entries can be thousands of bases long
- BAM:
  - Binary alignment format
  - Can store aligned reads
  - Can also store unaligned reads, especially in PacBio workflows
- uBAM:
  - Unaligned BAM
  - Contains reads and metadata before alignment
  - Useful for preserving sequencing and sample information
- POD5:
  - Oxford Nanopore signal-level file format
  - Stores raw signal data
  - Used for basecalling, re-basecalling, and signal-level analyses
- FAST5:
  - Older Oxford Nanopore signal-level format
  - Still encountered in older datasets and workflows

```diff
+ FASTQ: sequence + quality
+ BAM/uBAM: reads plus structured metadata
+ POD5/FAST5: Oxford Nanopore signal data
```

### **6. Quality Control Questions for Long Reads**

- How many reads do we have?
- What is the total yield in bases?
- What is the read length distribution?
- What is the N50 read length?
- What is the quality score distribution?
- Are there very short reads that should be filtered?
- Are there very low-quality reads that should be filtered?
- Is coverage sufficient for the planned analysis?
- Are there differences between samples or sequencing runs?
- For Oxford Nanopore:
  - Which basecaller and model were used?
  - Are reads simplex or duplex?
  - Are POD5/FAST5 files available if signal-level analysis is needed?
- For PacBio HiFi:
  - Are the reads already CCS/HiFi reads?
  - What minimum predicted accuracy was used?
  - Are reads provided as FASTQ or BAM?

```diff
! DISCUSSION: Why is read length distribution more important for long reads than for short reads?
! DISCUSSION: Why can filtering too aggressively remove useful long reads?
```

### **7. Common QC Metrics**

- Number of reads:
  - Total number of sequencing reads in the dataset
- Total bases:
  - Sum of all bases across all reads
- Mean read length:
  - Average read length
  - Can be affected by a small number of very long reads
- Median read length:
  - The middle read length after sorting reads by size
  - Often more robust than the mean
- N50:
  - Length such that 50% of all sequenced bases are in reads of that length or longer
  - Common metric for long reads datasets and assemblies
- Quality score:
  - Estimate of base-level or read-level accuracy
  - Interpretation depends on the sequencing technology and processing software
- Coverage:
  - Approximate number of times the genome is represented in the reads
  - Coverage = total bases / genome size

```diff
+ Example: 500 Mb of reads for a 5 Mb bacterial genome = about 100x coverage
```

### **8. Tools for Long Reads QC**

- `NanoPlot`
  - Commonly used for Oxford Nanopore read QC
  - Summarizes read length, quality, and yield
- `pycoQC`
  - Used with Oxford Nanopore sequencing summary files
  - Useful for run-level QC
- `seqkit`
  - General-purpose FASTA/FASTQ toolkit
  - Can summarize read counts, lengths, and basic statistics
- `samtools`
  - Useful for BAM files
  - Can inspect headers, count reads, and summarize alignments
- PacBio-specific tools:
  - Used depending on whether data are subreads, CCS reads, HiFi reads, FASTQ, or BAM

```diff
+ COMMAND: seqkit stats - summarize FASTA/FASTQ files
+ COMMAND: NanoPlot - generate long reads QC plots
+ COMMAND: samtools view - inspect BAM files
```

### **9. Practical QC Workflow Outline**

- Start from the provided long reads dataset
- Identify the file format:
  - FASTQ
  - BAM/uBAM
  - POD5/FAST5
- Create a dedicated analysis directory
- Inspect file size and file names
- Count reads and summarize sequence lengths
- Generate QC plots
- Interpret read length and quality distributions
- Decide whether filtering is needed
- Save QC reports in a dedicated results folder
- Keep raw data unchanged

```diff
! EXERCISE: Locate the long reads input files
! EXERCISE: Identify the file format from the file extension
! EXERCISE: Predict which QC tool is appropriate for each file type
```





### **10. Key Points**

- Long reads are useful because they can span genomic regions that short reads cannot resolve
- PacBio HiFi and Oxford Nanopore are both long reads technologies, but they generate different data types
- FASTQ still matters, but long reads workflows may also use BAM/uBAM and signal-level files
- Quality control for long reads must consider both quality and read length distribution
- The best QC decisions depend on the downstream task:
  - Alignment
  - Assembly
  - SNV/indel calling
  - Structural variant calling
  - Methylation or signal-level analysis

```diff
+ Poor QC decisions at this stage can affect alignment, assembly, and variant calling
+ We will use this QC information in the next lesson on long reads alignment
```
