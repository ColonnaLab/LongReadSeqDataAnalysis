[back to course home page ](../README.md)

## Advanced 5: Pangenomics Genotyping with PanGenie

This lesson is adapted from the [PanGenie workshop](https://pangenie-workshop.readthedocs.io/en/latest/index.html) by Jana Ebler, adapted to the course server and workspace.

The practical challenge is in [the pangenomics genotyping hackathon brief](11-pangenomics-genotyping-hackathon.md).

PanGenie is documented at [pangenie.readthedocs.io](https://pangenie.readthedocs.io/en/latest/index.html). The workshop material used here is by Jana Ebler and is available from the [PanGenie workshop repository](https://github.com/eblerjana/pangenie-workshop).

Dataset assumed for the practical:

- source tutorial: `https://pangenie-workshop.readthedocs.io/en/latest/index.html`
- server data folder: `/scratch/user1/data-pangenie`
- working directory: `/scratch/user1/pangenie-working`
- pangenome panel VCF: `panel_multi_chr5:50200000-50400000.vcf`
- reference FASTA: `reference_chr5:50200000-50400000.fasta`
- sample FASTA files: `NA19189`, `NA19190`, and `NA19191` chromosome 5 region FASTA files

```diff
+ We will follow the PanGenie workshop workflow
+ We will build a PanGenie index from a pangenome panel and reference
+ We will genotype three sample FASTA files
+ We will convert PanGenie bubble genotypes to biallelic VCF records
```

### **1. What to Read**

Use the original tutorial as the main lesson text:

- [PanGenie workshop](https://pangenie-workshop.readthedocs.io/en/latest/index.html)
- [PanGenie documentation](https://pangenie.readthedocs.io/en/latest/index.html)

This course page only records where the practical is located and how it has been adapted for the course machine.

### **2. Key Point**

PanGenie is a genotyping tool for variants already represented in a pangenome panel. It is not used here to discover new variants. In this practical, we run the workshop workflow on the prepared server dataset and inspect the generated genotype VCF files.

### **3. Bubble VCF and Biallelic Callset VCF**

PanGenie uses a graph-aware representation during genotyping, then the results can be converted into a more conventional per-variant representation.

|                         | **Bubble VCF** | **Biallelic callset VCF** |
| ----------------------- | -------------- | ------------------------- |
| What one row represents | One complete top-level graph bubble | One individual variant allele |
| Variant structure | Usually multiallelic | Biallelic |
| Nested variants | Multiple SNPs, indels, or structural variants may be contained inside one record | Each nested variant receives its own record |
| Genotype meaning | Which complete paths through the bubble the sample carries | Presence (`1`) or absence (`0`) of one particular variant allele |
| Role in PanGenie | Used as the main graph/variant input for indexing and genotyping | Used during post-processing to produce conventional per-variant genotypes |

In a pangenome graph, a bubble is a region where haplotypes can follow different paths. In the bubble VCF, one VCF record describes the whole bubble: `REF` is the reference path through that bubble, and each `ALT` is another complete haplotype path. One `ALT` path may contain several nested variants, such as two SNPs plus a deletion. PanGenie directly genotypes these complete bubble paths.

The biallelic callset VCF decomposes those complex bubble records into one record per individual variant. Each record asks a simpler question: is this particular allele absent (`0`) or present (`1`) on each haplotype?

Example:

```diff
Bubble paths:
Allele 0: reference sequence
Allele 1: SNP
Allele 2: deletion

Simplified bubble:

                   Allele 0 reference path
                        /--- A----\
                       /           \
       left flank --+--------T--------+-- right flank
                       \           /
                        \---------/
                     Allele 2: deletion


+ Individual haplotypes:
haplotype 1 carries Allele 1: SNP 
haplotype 2 carries Allele 2: deletion

+ PanGenie bubble genotype:
        1|2

+ Biallelic callset genotypes:
        SNP_A       1|0    haplotype 1 has SNP_A, haplotype 2 does not
        deletion_C  0|1    haplotype 1 does not have deletion_C, haplotype 2 has it
```

The biological result is the same, but the level of representation is different:

```diff
+ Bubble VCF: which complete graph paths are present?
+ Biallelic callset VCF: which individual variants are present?
```

### **4. Credit**

