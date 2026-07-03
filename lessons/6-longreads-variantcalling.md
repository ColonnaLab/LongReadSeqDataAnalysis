[back to course home page ](../README.md)

## Long reads variant calling

This page contains the conceptual notes for the long reads variant calling lesson. The practical commands using the MRSA KUN1163 dataset are in [the variant calling exercise](6-longreads-variantcalling-exercise.md).

The input for this lesson is the sorted and indexed BAM file produced in lesson 5:

- reference genome: `../data-longreads/reference/KUN1163_reference.fasta`
- aligned reads: `bam/DRR187567.KUN1163.sorted.bam`
- BAM index: `bam/DRR187567.KUN1163.sorted.bam.bai`

```diff
+ We will use aligned long reads to identify sequence differences from the reference
+ We will distinguish small variants from structural variants
+ We will produce VCF files
+ We will filter and inspect variant calls
```

### **1. What Is Variant Calling?**

- Variant calling means identifying differences between the sample reads and the reference genome
- These differences can be small or large
- Small variants include:
  - SNVs: single nucleotide variants
  - small insertions
  - small deletions
  - small complex substitutions
- Structural variants include:
  - large deletions (`DEL`)
  - insertions (`INS`)
  - duplications (`DUP`)
  - inversions (`INV`)
  - translocations (`TRA`)

```diff
+ Reference:  ATCGATCGATCG
+ Sample:     ATCGATTGATCG
+ Variant:         C -> T
```

### **2. Why Long Reads Are Useful for Variant Calling**

- Long reads can span repetitive regions
- Long reads can connect both sides of a large insertion or deletion
- Long reads can reveal split alignments caused by structural variants
- Long reads can support haplotype phasing
- Long reads can detect variants that are difficult to resolve with short reads

```diff
! DISCUSSION: Why is it useful for one read to span an entire deletion?
! DISCUSSION: Why are structural variants difficult for short reads?
```

### **3. Small Variants vs Structural Variants**

- Small variant calling:
  - focuses on SNVs and small indels
  - usually uses aligned reads in BAM format
  - output is usually VCF
- Structural variant calling:
  - focuses on larger genome changes
  - uses long alignments, split reads, and abnormal alignment patterns
  - output is also usually VCF, but with `SVTYPE` annotations

```diff
+ SNV: one base changes
+ Small indel: a few bases inserted or deleted
+ Structural variant: a larger genome segment changes
```

### **4. VCF Format**

- VCF stands for Variant Call Format
- VCF files contain:
  - header lines beginning with `##`
  - one column header line beginning with `#CHROM`
  - one row per variant
- Important VCF columns:
  - `CHROM`: reference sequence
  - `POS`: variant position
  - `REF`: reference allele
  - `ALT`: alternative allele
  - `QUAL`: variant quality
  - `FILTER`: whether the variant passed filters
  - `INFO`: extra annotations
  - `FORMAT`: sample-level fields

```diff
+ VCF is the standard text format for variant calls
+ It can store both small variants and structural variants
```

### **5. Variant Quality and Filtering**

- Variant callers assign quality values and annotations
- Not every called variant is reliable
- Low-confidence calls can be caused by:
  - low coverage
  - low base quality
  - mapping errors
  - repetitive sequence
  - strand bias
  - systematic sequencing errors
- Filtering removes or labels low-confidence calls

```diff
! DISCUSSION: Why can low coverage create false negative variants?
! DISCUSSION: Why can mapping errors create false positive variants?
```

### **6. Tools**

- `bcftools`
  - manipulates VCF/BCF files
  - can index, filter, count, view, and summarize variants
- `freebayes`
  - small variant caller
  - already used in the short-read lesson
- Long-read small variant callers:
  - tools may differ depending on ONT, PacBio HiFi, and analysis goal
- Structural variant callers:
  - use long-read alignment patterns to identify larger variants

```diff
+ In this draft practical, we focus on VCF generation, inspection, and filtering
+ Structural variant calling concepts will be introduced, then expanded as needed
```

### **7. Read-Based vs Assembly-Based Structural Variant Calling**

- Read-based SV calling:
  - uses alignments of sequencing reads to the reference
  - detects split reads, long gaps, and abnormal alignment patterns
  - works directly from BAM files
- Assembly-vs-reference SV calling:
  - first assembles the sample genome
  - aligns the assembly to the reference
  - compares assembled genome structure with the reference
  - can be useful for larger and more complex events

```diff
+ Read-based: reads -> reference -> SV calls
+ Assembly-based: reads -> assembly -> reference comparison -> SV calls
```

### **8. Practical Workflow Overview**

- Check that the sorted BAM and BAM index from lesson 5 are present
- Check that the reference genome is indexed
- Create output folders for VCF files
- Run a small variant caller
- Compress and index VCF files
- Inspect VCF header and variant rows
- Count variants
- Filter variants by quality or depth
- Discuss structural variant calling and how long reads make it possible

```diff
+ The practical commands are in: lessons/6-longreads-variantcalling-exercise.md
```

### **9. Key Points**

- Variant calling identifies differences between a sample and a reference genome
- Long reads are useful for both small variants and structural variants
- Structural variants include deletions, insertions, duplications, inversions, and translocations
- VCF is the common output format for variant calls
- Filtering is necessary because not every call is reliable
- The quality of variant calls depends on read quality, alignment quality, and coverage

