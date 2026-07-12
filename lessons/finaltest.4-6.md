[back to course home page ](../README.md)

### Basic Bioinformatics and Long-Read Sequencing Data Analysis
##### Final Test - Lessons 4-6

Name: ________________Date: ________________

*Instructions: Answer all questions. Write commands exactly as you would type them.*
---

#### Part 1: Long-Read Sequencing and Quality Control

**1.** Name one important difference between Oxford Nanopore sequencing and PacBio HiFi sequencing.
```diff
+Answer: ____________________________________________________________________
```

**2.** In Oxford Nanopore sequencing, what is basecalling?
```diff
+Answer: ____________________________________________________________________
```

**3.** What does a Phred quality score estimate?
```diff
+Answer: ____________________________________________________________________
```

**4.** Which command did we use to summarize basic statistics from a FASTQ file?
```diff
+Answer: ____________________________________________________________________
```

**5.** Which tool did we use to generate long-read QC plots such as read length and quality distributions?
```diff
+Answer: ____________________________________________________________________
```

**6.** What does N50 read length mean?
```diff
+Answer: ____________________________________________________________________
```

---

#### Part 2: Long-Read Alignment

**7.** Which aligner did we use for Oxford Nanopore reads?
```diff
+Answer: ____________________________________________________________________
```

**8.** In the command below, what does `-x map-ont` mean?
```bash
minimap2 -ax map-ont reference.fasta reads.fastq.gz > alignment.sam
```
```diff
+Answer: ____________________________________________________________________
```

**9.** Put these files in the correct order for the long-read alignment workflow:
___ sorted BAM    ___ FASTQ    ___ SAM    ___ BAM

**10.** Why do we sort and index a BAM file?
```diff
+Answer: ____________________________________________________________________
```

**11.** Match each `samtools` command with its function:
- `samtools flagstat` ___    a) reports coverage by reference sequence
- `samtools idxstats` ___    b) summarizes mapped and unmapped reads
- `samtools coverage` ___    c) counts reads assigned to each reference sequence
- `samtools tview` ___       d) visualizes alignments in the terminal

**12.** What is a supplementary alignment, and why can it be useful for long-read analysis?
```diff
+Answer: ____________________________________________________________________
```

**13.** To inspect the BAM file in IGV, which files must be copied to the local machine?
```diff
+Answer: ____________________________________________________________________
```

---

#### Part 3: Long-Read Variant Calling

**14.** What is variant calling?
```diff
+Answer: ____________________________________________________________________
```

**15.** In a VCF file, what do the `REF` and `ALT` columns represent?
```diff
+Answer: ____________________________________________________________________
```
