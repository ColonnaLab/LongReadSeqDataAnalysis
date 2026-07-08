# Software Used in This Course

This file lists the main software used in the lessons, with links to the official source, repository, or installation page.

## Command Line Basics

| Software | Used for | Source / installation |
|---|---|---|
| Bash / Unix shell | Command-line work, navigation, pipes, redirects, scripts | https://www.gnu.org/software/bash/ |
| GNU coreutils | Common commands such as `ls`, `pwd`, `cat`, `head`, `wc`, `mkdir` | https://www.gnu.org/software/coreutils/ |
| gzip | Compression and decompression with `gzip`, `gunzip`, `zcat` | https://www.gnu.org/software/gzip/ |
| curl | Downloading files from URLs | https://curl.se/ |
| OpenSSH | Remote login and file transfer with `ssh` and `scp` | https://www.openssh.com/ |
| less | Viewing long text files in the terminal | https://www.greenwoodsoftware.com/less/ |

## Short-Read Analysis

| Software | Used for | Source / installation |
|---|---|---|
| FastQC | Quality control reports for FASTQ files | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Trimmomatic | Adapter and quality trimming of Illumina reads | http://www.usadellab.org/cms/?page=trimmomatic |
| Java | Running Trimmomatic `.jar` files | https://openjdk.org/ |
| BWA | Short-read alignment to a reference genome | https://github.com/lh3/bwa |
| SAMtools | SAM/BAM conversion, sorting, indexing, viewing, FASTA indexing |https://www.htslib.org/doc/samtools.html  https://github.com/samtools/samtools |
| FreeBayes | Small variant calling from aligned reads | https://github.com/freebayes/freebayes |
| BCFtools | Filtering and manipulating VCF/BCF files | https://github.com/samtools/bcftools |

## Long-Read Analysis

| Software | Used for | Source / installation |
|---|---|---|
| SRA Toolkit | Downloading SRA data with `prefetch` and `fasterq-dump` | https://github.com/ncbi/sra-tools |
| seqkit | FASTA/FASTQ summary statistics, read counts, lengths, N50 | https://github.com/shenwei356/seqkit |
| NanoPlot | Oxford Nanopore read QC plots and HTML reports | https://github.com/wdecoster/NanoPlot |
| pycoQC | Oxford Nanopore run-level QC from sequencing summary files | https://github.com/a-slide/pycoQC |
| minimap2 | Long-read alignment to a reference genome |https://lh3.github.io/minimap2/minimap2.html  https://github.com/lh3/minimap2 |
| SAMtools | BAM inspection, indexing, coverage and alignment summaries | https://github.com/samtools/samtools |
| BCFtools | VCF/BCF filtering and manipulation | https://github.com/samtools/bcftools |

## R and Reproducible Analysis

| Software | Used for | Source / installation |
|---|---|---|
| R | Statistical computing and reproducible data analysis | https://www.r-project.org/ |
| RStudio Desktop | Graphical IDE for R | https://posit.co/download/rstudio-desktop/ |

## Common Conda Installation Sources

Many bioinformatics tools above are available through Bioconda. When using conda/mamba, the usual channels are:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Useful package pages:

| Package | Bioconda page |
|---|---|
| FastQC | https://bioconda.github.io/recipes/fastqc/README.html |
| Trimmomatic | https://bioconda.github.io/recipes/trimmomatic/README.html |
| BWA | https://bioconda.github.io/recipes/bwa/README.html |
| SAMtools | https://bioconda.github.io/recipes/samtools/README.html |
| FreeBayes | https://bioconda.github.io/recipes/freebayes/README.html |
| BCFtools | https://bioconda.github.io/recipes/bcftools/README.html |
| SRA Toolkit | https://bioconda.github.io/recipes/sra-tools/README.html |
| seqkit | https://bioconda.github.io/recipes/seqkit/README.html |
| NanoPlot | https://bioconda.github.io/recipes/nanoplot/README.html |
| pycoQC | https://bioconda.github.io/recipes/pycoqc/README.html |
| minimap2 | https://bioconda.github.io/recipes/minimap2/README.html |

