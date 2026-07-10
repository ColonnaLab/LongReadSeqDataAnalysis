[back to course home page ](../README.md)

## Long reads variant calling

This page contains the conceptual notes for the long reads variant calling lesson. The practical commands using the MRSA KUN1163 dataset are in [the variant calling exercise](6-longreads-variantcalling-exercise.md).

The input for this lesson is the sorted and indexed BAM file produced in lesson 5:

- reference genome: `data-longreads/reference/KUN1163_reference.fasta`
- aligned reads: `bam/DRR187567.KUN1163.sorted.bam`
- BAM index: `bam/DRR187567.KUN1163.sorted.bam.bai`

```diff
+ We will use aligned Oxford Nanopore reads to identify sequence differences from the reference
+ We will focus on modern long-read small-variant callers
+ We will compare Clair3 and DeepVariant outputs
+ We will inspect, count, filter, and compare VCF files
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
  - output is usually VCF or compressed VCF
- Structural variant calling:
  - focuses on larger genome changes
  - uses long alignments, split reads, and abnormal alignment patterns
  - output is also usually VCF, but with `SVTYPE` annotations

```diff
+ SNV: one base changes
+ Small indel: a few bases inserted or deleted
+ Structural variant: a larger genome segment changes
```

In this lesson we focus on small variant calling with two long-read-aware tools:

- [Clair3](https://github.com/HKU-BAL/Clair3)
- [DeepVariant](https://github.com/google/deepvariant)

Structural variant calling is introduced conceptually, but it is not the main practical workflow here.

### **4. Why We Use Specialized Long-Read Callers**

Long-read variant callers are designed for the error profiles of long-read sequencing platforms.

- Oxford Nanopore reads can contain more indel errors than short Illumina reads
- Raw read accuracy depends on pore chemistry, basecaller, and basecalling model
- Long reads can be phased across longer distances
- A caller trained for one platform or chemistry may not be appropriate for another

For this reason, long-read variant calling tools usually require:

- the sequencing platform, such as ONT or PacBio HiFi
- a trained model
- a sorted and indexed BAM file
- an indexed reference FASTA
- enough read coverage for reliable calls

```diff
! DISCUSSION: Why should an ONT dataset not automatically be analyzed with a PacBio HiFi model?
! DISCUSSION: Why might a deep-learning caller need a trained model?
```

### **5. Clair3**

[Clair3](https://github.com/HKU-BAL/Clair3) is a germline small-variant caller for long-read sequencing.

The Clair3 method is described in Zheng et al., [*Symphonizing pileup and full-alignment for deep learning-based long-read variant calling*](https://www.nature.com/articles/s43588-022-00387-x), published in *Nature Computational Science* in 2022.

Important ideas:

- It combines a fast pileup model with a more detailed full-alignment model
- It supports Oxford Nanopore, PacBio HiFi, and Illumina modes
- It produces VCF output for small variants
- It can use phasing information
- It includes options useful for haploid or non-diploid organisms

In human germline analysis, variant callers often assume diploidy. Bacteria are usually haploid. For this MRSA exercise, we therefore pay attention to haploid calling options such as:

```diff
+ --haploid_precise
+ --include_all_ctgs
+ --no_phasing_for_fa
```

The exact model matters. For ONT data, the model should match the pore chemistry and basecaller as closely as possible. In a training course, the model may already be installed on the server by the instructor.

### **6. DeepVariant**

[DeepVariant](https://github.com/google/deepvariant) is a deep-learning variant calling pipeline.

Important ideas:

- It converts read evidence around candidate variants into examples for a neural network
- It supports multiple sequencing data types through `--model_type`
- For long-read data, model choice is essential
- It produces VCF and optionally gVCF output
- The released DeepVariant germline models are mainly designed for diploid human variant calling

Common model types include:

```diff
+ WGS: Illumina whole-genome sequencing
+ PACBIO: PacBio data
+ ONT_R104: Oxford Nanopore R10.4 data
+ HYBRID_PACBIO_ILLUMINA: combined PacBio and Illumina data
```

For this MRSA lesson, DeepVariant is useful as a teaching comparison: it shows that different callers can use the same BAM and reference but may produce different call sets. Because the sample is bacterial and haploid, DeepVariant output should not be treated as a validated production call set without additional benchmarking.

```diff
! DISCUSSION: What could cause Clair3 and DeepVariant to disagree at a position?
! DISCUSSION: Why is a human-trained diploid caller not automatically ideal for bacterial data?
```

### **7. VCF Format**

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
+ It can store small variants, reference confidence records, and structural variants
```

### **8. Variant Quality and Filtering**

- Variant callers assign quality values and annotations
- Not every called variant is reliable
- Low-confidence calls can be caused by:
  - low coverage
  - low base quality
  - mapping errors
  - repetitive sequence
  - systematic sequencing errors
  - model mismatch
- Filtering removes or labels low-confidence calls

Filtering should be interpreted carefully because different callers use different annotations. A simple `QUAL` threshold is useful for teaching, but production workflows usually use caller-specific recommendations and validation data.

```diff
! DISCUSSION: Why can low coverage create false negative variants?
! DISCUSSION: Why can mapping errors create false positive variants?
! DISCUSSION: Why can a model mismatch create both false positives and false negatives?
```

### **9. Comparing Callers**

Two callers can disagree because they use different:

- candidate variant detection methods
- neural network models
- assumptions about ploidy
- filters
- indel representations
- handling of low-coverage or repetitive regions

Before comparing VCF files, variants are often normalized so that indels are represented consistently relative to the reference.

Useful tools:

- `bcftools view`: inspect and subset VCF files
- `bcftools stats`: summarize variant calls
- `bcftools query`: extract selected fields
- `bcftools norm`: normalize variant representation
- `bcftools isec`: compare intersections and differences between VCF files

```diff
+ Same BAM + same reference does not guarantee identical variant calls
+ Comparing callers teaches us where the evidence is strong and where it is ambiguous
```

### **10. Read-Based vs Assembly-Based Structural Variant Calling**

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

### **11. Practical Workflow Overview**

- Check that the sorted BAM and BAM index from lesson 5 are present
- Check that the reference genome is indexed
- Inspect coverage before variant calling
- Run Clair3 on the ONT BAM
- Run DeepVariant on the same BAM
- Inspect each VCF header and variant rows
- Count variants from each caller
- Filter variants
- Normalize VCF files
- Compare overlap and caller-specific calls
- Discuss structural variant signals in long-read alignments

```diff
+ The practical commands are in: lessons/6-longreads-variantcalling-exercise.md
```

### **12. Key Points**

- Variant calling identifies differences between a sample and a reference genome
- Long reads are useful for both small variants and structural variants
- Clair3 and DeepVariant are long-read-aware small variant callers
- The sequencing platform and model must match the data as closely as possible
- Bacterial datasets require attention to haploid calling assumptions
- VCF is the common output format for variant calls
- Filtering and comparison are necessary because not every call is reliable
- The quality of variant calls depends on read quality, alignment quality, model choice, and coverage

### **13. Reference**

- Zheng, Z., Li, S., Su, J. et al. [Symphonizing pileup and full-alignment for deep learning-based long-read variant calling](https://www.nature.com/articles/s43588-022-00387-x). *Nature Computational Science* 2, 797-803 (2022). https://doi.org/10.1038/s43588-022-00387-x
