# File Formats Used in This Course

This file lists the main file formats used in the lessons, with links to official specifications, documentation, or widely used reference pages.

## Sequence and Read Files

| Format | Used for | Specification / reference |
|---|---|---|
| FASTA (`.fa`, `.fasta`, `.fna`) | Reference genomes, assembled sequences, transcript or protein sequences | https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ |
| FASTQ (`.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`) | Sequencing reads with per-base quality scores | https://maq.sourceforge.net/fastq.shtml |
| SRA (`.sra`) | NCBI Sequence Read Archive data before conversion to FASTQ | https://github.com/ncbi/sra-tools/wiki |

## Alignment Files

| Format | Used for | Specification / reference |
|---|---|---|
| SAM (`.sam`) | Text format for read alignments against a reference genome | https://samtools.github.io/hts-specs/SAMv1.pdf |
| BAM (`.bam`) | Compressed binary version of SAM; common input for downstream analysis | https://samtools.github.io/hts-specs/SAMv1.pdf |
| CRAM (`.cram`) | Reference-based compressed alignment format | https://samtools.github.io/hts-specs/CRAMv3.pdf |
| PAF (`.paf`) | Pairwise mapping format often produced by minimap2 | https://github.com/lh3/miniasm/blob/master/PAF.md |

## Alignment Index Files

| Format | Used for | Specification / reference |
|---|---|---|
| FAI (`.fai`) | FASTA index created by `samtools faidx` | https://www.htslib.org/doc/faidx.html |
| BAI (`.bam.bai`, `.bai`) | BAM index for random access to genomic regions | https://samtools.github.io/hts-specs/SAMv1.pdf |
| CSI (`.csi`) | Coordinate-sorted index format for BAM/CRAM/VCF-like files | https://samtools.github.io/hts-specs/CSIv1.pdf |

## Variant Files

| Format | Used for | Specification / reference |
|---|---|---|
| VCF (`.vcf`, `.vcf.gz`) | Text format for small variants and many structural variant calls | https://samtools.github.io/hts-specs/VCFv4.5.pdf |
| BCF (`.bcf`) | Binary version of VCF used by BCFtools and other tools | https://samtools.github.io/hts-specs/BCFv2_qref.pdf |
| TBI (`.tbi`) | Tabix index for compressed VCF and other tab-delimited genomic files | https://www.htslib.org/doc/tabix.html |

## Annotation and Interval Files

| Format | Used for | Specification / reference |
|---|---|---|
| BED (`.bed`) | Genomic intervals such as regions, targets, and annotations | https://genome.ucsc.edu/FAQ/FAQformat.html#format1 |
| GFF3 (`.gff`, `.gff3`) | Genome feature annotations | https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md |
| GTF (`.gtf`) | Gene and transcript annotations, commonly used in RNA-seq workflows | https://useast.ensembl.org/info/website/upload/gff.html |

## Long-Read Signal and Run Files

| Format | Used for | Specification / reference |
|---|---|---|
| POD5 (`.pod5`) | Oxford Nanopore signal-level data | https://pod5-file-format.readthedocs.io/ |
| FAST5 (`.fast5`) | Older Oxford Nanopore signal-level data format | https://github.com/nanoporetech/ont_fast5_api |
| sequencing summary (`sequencing_summary.txt`) | Oxford Nanopore run-level read metadata and QC input for tools such as NanoPlot and pycoQC | https://github.com/wdecoster/NanoPlot |

## Common Compression

| Format | Used for | Specification / reference |
|---|---|---|
| gzip (`.gz`) | Compression of FASTQ, VCF, and other text files | https://www.gnu.org/software/gzip/ |
| BGZF (`.gz`, often block-compressed VCF/BED/GFF) | Blocked gzip compression that allows indexing and random access | https://samtools.github.io/hts-specs/SAMv1.pdf |
