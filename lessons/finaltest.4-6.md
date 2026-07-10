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

**4.** Write the formula for expected average read depth using genome size, number of reads, and average read length.
```diff
+Answer: ____________________________________________________________________
```

**5.** Which command did we use to summarize basic statistics from a FASTQ file?
```diff
+Answer: ____________________________________________________________________
```

**6.** Which tool did we use to generate long-read QC plots such as read length and quality distributions?
```diff
+Answer: ____________________________________________________________________
```

**7.** What does N50 read length mean?
```diff
+Answer: ____________________________________________________________________
```

---

#### Part 2: Long-Read Alignment

**8.** Which aligner did we use for Oxford Nanopore reads?
```diff
+Answer: ____________________________________________________________________
```

**9.** In the command below, what does `-x map-ont` mean?
```bash
minimap2 -ax map-ont reference.fasta reads.fastq.gz > alignment.sam
```
```diff
+Answer: ____________________________________________________________________
```

**10.** Put these files in the correct order for the long-read alignment workflow:
___ sorted BAM    ___ FASTQ    ___ SAM    ___ BAM

**11.** Why do we sort and index a BAM file?
```diff
+Answer: ____________________________________________________________________
```

**12.** Match each `samtools` command with its function:
- `samtools flagstat` ___    a) reports coverage by reference sequence
- `samtools idxstats` ___    b) summarizes mapped and unmapped reads
- `samtools coverage` ___    c) counts reads assigned to each reference sequence
- `samtools tview` ___       d) visualizes alignments in the terminal

**13.** What is a supplementary alignment, and why can it be useful for long-read analysis?
```diff
+Answer: ____________________________________________________________________
```

**14.** To inspect the BAM file in IGV, which files must be copied to the local machine?
```diff
+Answer: ____________________________________________________________________
```

---

#### Part 3: Long-Read Variant Calling

**15.** What is variant calling?
```diff
+Answer: ____________________________________________________________________
```

**16.** Give one example of a small variant and one example of a structural variant.
```diff
+Answer: ____________________________________________________________________
```

**17.** Why is it important to choose a Clair3 model that matches the sequencing technology and basecalling model?
```diff
+Answer: ____________________________________________________________________
```

**18.** What does this command do?
```bash
ls $CONDA_PREFIX/bin/models
```
```diff
+Answer: ____________________________________________________________________
```

**19.** Why did we use haploid options for the MRSA dataset?
```diff
+Answer: ____________________________________________________________________
```

**20.** In a VCF file, what do the `REF` and `ALT` columns represent?
```diff
+Answer: ____________________________________________________________________
```

**21.** What does this command count?
```bash
bcftools view -H variants/clair3/merge_output.vcf.gz | wc -l
```
```diff
+Answer: ____________________________________________________________________
```

**22.** Why do we normalize VCF files before comparing calls from Clair3 and DeepVariant?
```diff
+Answer: ____________________________________________________________________
```

**23.** What does `bcftools isec` help us compare?
```diff
+Answer: ____________________________________________________________________
```

---

**Bonus Question**

Long reads can help detect structural variants. Explain why one long read spanning a repetitive or rearranged region can be more informative than many short reads.
```diff
+Answer: ____________________________________________________________________
```

---
